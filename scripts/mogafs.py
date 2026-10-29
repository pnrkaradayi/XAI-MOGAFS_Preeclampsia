"""
XAI-MOGAFS — SHAP yönlendirmeli çok amaçlı genetik özellik seçimi (NSGA-II), sendrom özgüllüğü amacıyla.

Temsil      : K adaylık havuzdan (eğitim katlamasında ön-filtre) seçilen 1..max_k özellik (ikili vektör).
Sınıflandırıcı: L2 lojistik regresyon (n=~125, ≤5 özellik için uygun).
Amaçlar (hepsi en aza indirilir):
  f1 = 1 − AUC(PE vs kontrol)            eğitim katlamasında iç çapraz doğrulama (5 kat x 2 tekrar)
  f2 = max(0, AUC_off − 0.5)            AUC_off: geliştirme kohortunda (GSE149440) komplikasyon
                                         (sPTD + PPROM) vs kontrol; model eğitim katlamasının tamamına fit
                                         edilir ve geliştirme kohortuna uygulanır ("alarm başka sendromda
                                         ne kadar çalıyor")
  f3 = |S| / max_k
  f4 = SD(AUC) gebelik haftası dilimleri arasında (iç CV tahminlerinden)
  f5 = 1 − ortalama bootstrap seçim kararlılığı (isteğe bağlı; `stab` verilirse). Kararlılık: eğitim
       katlamasının B bootstrap örneğinde özelliğin ilk-N içinde kalma sıklığı (bootstrap_stability).
Tekrarlı ölçüm (isteğe bağlı; `groups`): aynı kadının birden çok örneği varsa iç CV kadın düzeyinde
  (StratifiedGroupKFold) bölünür; aynı kadın hem iç-eğitimde hem iç-testte olmaz.
Merkez değişmezliği (isteğe bağlı; `sites`) — XAI-MOGAFS-SI (Val-2 açıldıktan SONRA eklendi; keşifsel):
  f_si = 1 − min_merkez AUC (iç CV tahminlerinin her eğitim merkezindeki AUC'sinin en kötüsü); SHAP yönlendirmesi
  her merkezde ayrı hesaplanır ve özellik skoru = merkezler arası EN DÜŞÜK normalize |SHAP| (yalnız her merkezde
  önemli olan özellikler öne çıkar).
Eşlik eden değişken (isteğe bağlı; `covar`): örn. gebelik haftası dilimi; tüm lojistik fitlere sabit
  sütun olarak eklenir, panel boyutuna (≤ max_k) sayılmaz.
SHAP yönlendirmesi (popülasyon tohumlama + mutasyon):
  imp_pe  = ortalama |SHAP|, PE vs kontrol modeli (eğitim katlaması)
  imp_off = ortalama |SHAP|, komplikasyon vs kontrol modeli (geliştirme kohortu)
  özgüllük kolu: skor = imp_pe/max − λ·imp_off/max   (sınıf-çiftine özgü SHAP farkı)
  düz kol      : skor = imp_pe/max
  Doğrusal model için girişimsel (interventional) SHAP kesin olarak φ_ij = β_j (x_ij − E[x_j])'dir;
  shap.LinearExplainer ile aynı sonucu verir.
"""
import numpy as np
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.core.crossover import Crossover
from pymoo.core.mutation import Mutation
from pymoo.core.problem import Problem
from pymoo.core.repair import Repair
from pymoo.core.sampling import Sampling
from pymoo.optimize import minimize
from pymoo.util.nds.non_dominated_sorting import NonDominatedSorting
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedGroupKFold


def fit_lr(X, y, C=1.0):
    return LogisticRegression(C=C, solver="liblinear", max_iter=500).fit(X, y)


def linear_shap_importance(X, y, C=0.1):
    """Mean |SHAP| per feature for an L2 logistic model (exact interventional SHAP, log-odds scale)."""
    m = LogisticRegression(C=C, solver="lbfgs", max_iter=5000).fit(X, y)
    phi = m.coef_[0] * (X - X.mean(axis=0))
    return np.abs(phi).mean(axis=0)


def safe_auc(y, s):
    return roc_auc_score(y, s) if len(np.unique(y)) == 2 else np.nan


class SubsetEvaluator:
    """Evaluates objectives for a feature subset; caches by subset."""

    def __init__(self, X, y, ga_bin, Xdev=None, dev_label=None, use_spec=True,
                 max_k=5, n_splits=5, n_repeats=2, C=1.0, seed=0, covar=None, stab=None, groups=None, sites=None):
        self.X, self.y, self.ga_bin = X, y, ga_bin
        self.covar = None if covar is None else np.asarray(covar, float).reshape(len(y), -1)
        self.stab = stab
        self.sites = None if sites is None else np.asarray(sites)
        self.Xdev, self.dev_label, self.use_spec = Xdev, dev_label, use_spec
        self.max_k, self.C, self.n_rep, self.n_splits = max_k, C, n_repeats, n_splits
        if groups is None:
            self.splits = list(RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats,
                                                       random_state=seed).split(X, y))
        else:
            self.splits = [sp for r in range(n_repeats) for sp in
                           StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed + r).split(X, y, groups)]
        self.bins = [np.where(ga_bin == b)[0] for b in np.unique(ga_bin)]
        if Xdev is not None:
            self.off_rows = np.where(np.isin(dev_label, ["Control", "sPTD", "PPROM"]))[0]
            self.off_y = np.isin(dev_label[self.off_rows], ["sPTD", "PPROM"]).astype(int)
            self.pe_rows = np.where(np.isin(dev_label, ["PE", "sPTD", "PPROM"]))[0]
            self.pe_y = (dev_label[self.pe_rows] == "PE").astype(int)
        self.cache = {}

    def evaluate(self, idx):
        key = tuple(sorted(int(i) for i in idx))
        if key in self.cache:
            return self.cache[key]
        Xs = self.X[:, key]
        if self.covar is not None:
            Xs = np.column_stack([Xs, self.covar])
        oof = np.zeros((self.n_rep, len(self.y)))
        for r, (tr, te) in enumerate(self.splits):
            m = fit_lr(Xs[tr], self.y[tr], self.C)
            oof[r // self.n_splits, te] = m.decision_function(Xs[te])
        p = oof.mean(axis=0)
        auc = roc_auc_score(self.y, p)
        bin_aucs = [safe_auc(self.y[b], p[b]) for b in self.bins]
        rec = {"auc_cv": auc, "f1": 1 - auc, "f3": len(key) / self.max_k,
               "f4": float(np.nanstd(bin_aucs)), "bin_aucs": bin_aucs}
        if self.stab is not None:
            rec["f5"] = 1 - float(np.mean(self.stab[list(key)]))
        if self.sites is not None:
            site_aucs = [safe_auc(self.y[self.sites == st], p[self.sites == st]) for st in np.unique(self.sites)]
            rec["worst_site_auc"] = float(np.nanmin(site_aucs))
            rec["f_si"] = 1 - rec["worst_site_auc"]
        if self.Xdev is not None:
            full = fit_lr(Xs, self.y, self.C)
            s = full.decision_function(self.Xdev[:, key])
            rec["auc_off"] = roc_auc_score(self.off_y, s[self.off_rows])
            rec["auc_pe_vs_comp"] = roc_auc_score(self.pe_y, s[self.pe_rows])
            rec["f2"] = max(0.0, rec["auc_off"] - 0.5)
        self.cache[key] = rec
        return rec

    def F(self, rec):
        f = [rec["f1"], rec["f2"], rec["f3"], rec["f4"]] if self.use_spec else [rec["f1"], rec["f3"], rec["f4"]]
        if self.stab is not None:
            f = f + [rec["f5"]]
        if self.sites is not None:
            f = f + [rec["f_si"]]
        return f


# ---------------------------------------------------------------- pymoo bileşenleri
class SubsetProblem(Problem):
    def __init__(self, evaluator, n_var):
        self.ev = evaluator
        n_obj = (4 if evaluator.use_spec else 3) + (evaluator.stab is not None) + (evaluator.sites is not None)
        super().__init__(n_var=n_var, n_obj=n_obj, xl=0, xu=1, vtype=bool)

    def _evaluate(self, X, out, *args, **kwargs):
        out["F"] = np.array([self.ev.F(self.ev.evaluate(np.flatnonzero(x))) for x in X], dtype=float)


class GuidedSampling(Sampling):
    """İlk popülasyon: en yüksek skorlu 1..max_k özellik (deterministik tohumlar) + skor ağırlıklı rastgele kümeler."""

    def __init__(self, prob, max_k):
        super().__init__()
        self.p, self.max_k = prob, max_k

    def _do(self, problem, n_samples, *args, random_state=None, **kwargs):
        X = np.zeros((n_samples, problem.n_var), dtype=bool)
        order = np.argsort(-self.p)
        for i in range(n_samples):
            if i < self.max_k:
                X[i, order[: i + 1]] = True
            else:
                k = random_state.integers(1, self.max_k + 1)
                X[i, random_state.choice(problem.n_var, size=k, replace=False, p=self.p)] = True
        return X


class SubsetCrossover(Crossover):
    """Küme tabanlı çaprazlama: ebeveynlerin birleşiminden, ikisinde de olanlara iki kat ağırlıkla seçim."""

    def __init__(self, prob, max_k):
        super().__init__(2, 2, prob=0.9)
        self.p, self.max_k = prob, max_k

    def _do(self, problem, X, *args, random_state=None, **kwargs):
        _, n_mat, n_var = X.shape
        Y = np.zeros((2, n_mat, n_var), dtype=bool)
        for m in range(n_mat):
            a, b = X[0, m], X[1, m]
            union = np.flatnonzero(a | b)
            w = self.p[union] * (1 + (a & b)[union])
            w = w / w.sum()
            lo, hi = sorted([max(a.sum(), 1), max(b.sum(), 1)])
            for c in range(2):
                k = min(random_state.integers(lo, hi + 1), len(union), self.max_k)
                Y[c, m, random_state.choice(union, size=k, replace=False, p=w)] = True
        return Y


class GuidedMutation(Mutation):
    """SHAP yönlendirmeli mutasyon: ekle (skorla orantılı), çıkar (düşük skor öncelikli) veya değiştir."""

    def __init__(self, prob, max_k):
        super().__init__(prob=0.9)
        self.p, self.max_k = prob, max_k

    def _add(self, x, rs):
        free = np.flatnonzero(~x)
        w = self.p[free] / self.p[free].sum()
        x[rs.choice(free, p=w)] = True

    def _drop(self, x, rs):
        on = np.flatnonzero(x)
        w = (1 - self.p[on] / self.p.max()) + 1e-3
        x[rs.choice(on, p=w / w.sum())] = False

    def _do(self, problem, X, *args, random_state=None, **kwargs):
        X = X.copy()
        for x in X:
            k, u = x.sum(), random_state.random()
            if k < self.max_k and (u < 1 / 3 or k <= 1):
                self._add(x, random_state)
            elif k > 1 and u < 2 / 3:
                self._drop(x, random_state)
            else:
                self._drop(x, random_state)
                self._add(x, random_state)
        return X


class CardinalityRepair(Repair):
    def __init__(self, prob, max_k):
        super().__init__()
        self.p, self.max_k = prob, max_k

    def _do(self, problem, X, **kwargs):
        X = X.astype(bool).copy()
        for x in X:
            on = np.flatnonzero(x)
            if len(on) == 0:
                x[np.argmax(self.p)] = True
            elif len(on) > self.max_k:
                x[on[np.argsort(self.p[on])[: len(on) - self.max_k]]] = False
        return X


def guidance(X, y, Xdev=None, dev_label=None, lam=1.0, use_spec=True, floor=0.01, sites=None):
    """Özellik seçim olasılıkları (tohumlama + mutasyon) ve bileşenleri."""
    if sites is not None:
        per_site = []
        for st in np.unique(sites):
            m = sites == st
            imp = linear_shap_importance(X[m], y[m])
            per_site.append(imp / imp.max())
        a = np.min(per_site, axis=0)          # merkez-değişmez önem: her merkezde önemli olmalı
    else:
        imp_pe = linear_shap_importance(X, y)
        a = imp_pe / imp_pe.max()
    if use_spec:
        rows = np.isin(dev_label, ["Control", "sPTD", "PPROM"])
        yo = np.isin(dev_label[rows], ["sPTD", "PPROM"]).astype(int)
        imp_off = linear_shap_importance(Xdev[rows], yo)
        b = imp_off / imp_off.max()
        score = a - lam * b
    else:
        b = np.zeros_like(a)
        score = a
    prob = np.clip(score, 0, None) + floor
    return prob / prob.sum(), a, b


def bootstrap_stability(X, y, B=50, top=20, seed=0):
    """Her özelliğin B tabakalı bootstrap örneğinde |Welch t| ile ilk `top` içinde kalma sıklığı (0–1)."""
    rng = np.random.default_rng(seed)
    pos, neg = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    freq = np.zeros(X.shape[1])
    for _ in range(B):
        idx = np.r_[rng.choice(pos, len(pos)), rng.choice(neg, len(neg))]
        Xb, yb = X[idx], y[idx]
        a, b = Xb[yb == 1], Xb[yb == 0]
        t = (a.mean(0) - b.mean(0)) / np.sqrt(a.var(0, ddof=1) / len(a) + b.var(0, ddof=1) / len(b) + 1e-12)
        freq[np.argsort(-np.abs(t))[:top]] += 1
    return freq / B


def run_mogafs(X, y, ga_bin, Xdev=None, dev_label=None, use_spec=True, max_k=5,
               pop_size=100, n_gen=60, seed=0, lam=1.0, covar=None, stab=None, groups=None, sites=None):
    ev = SubsetEvaluator(X, y, ga_bin, Xdev, dev_label, use_spec=use_spec, max_k=max_k, seed=seed,
                         covar=covar, stab=stab, groups=groups, sites=sites)
    prob, a, b = guidance(X, y, Xdev, dev_label, lam=lam, use_spec=use_spec, sites=sites)
    algo = NSGA2(pop_size=pop_size, sampling=GuidedSampling(prob, max_k),
                 crossover=SubsetCrossover(prob, max_k), mutation=GuidedMutation(prob, max_k),
                 repair=CardinalityRepair(prob, max_k), eliminate_duplicates=True)
    minimize(SubsetProblem(ev, X.shape[1]), algo, ("n_gen", n_gen), seed=seed, verbose=False)
    # Arşiv: değerlendirilen TÜM kümeler üzerinden baskın olmayan cephe (son popülasyondan daha eksiksiz)
    keys = list(ev.cache)
    F = np.array([ev.F(ev.cache[k]) for k in keys])
    front = NonDominatedSorting().do(F, only_non_dominated_front=True)
    archive = [{"idx": keys[i], **{k: v for k, v in ev.cache[keys[i]].items() if k != "bin_aucs"},
                "pareto": i in set(front)} for i in range(len(keys))]
    return archive, {"prob": prob, "imp_pe": a, "imp_off": b, "n_evaluated": len(keys)}


def choose(archive, use_spec, tol=0.03, use_stab=False, use_si=False):
    """ÖN-KAYITLI seçim kuralı (Pareto cephesinden tek panel):
    iç CV AUC'si en iyinin `tol` kadar yakınındaki çözümler; özgüllük kolunda en düşük f2, sonra en az
    özellik, sonra en düşük f4; düz kolda en az özellik, sonra en düşük f4, sonra en yüksek AUC.
    Kararlılık amacı varsa (use_stab) düz kolda önce en düşük f5 (en kararlı), sonra aynı sıra."""
    front = [r for r in archive if r["pareto"]]
    best = max(r["auc_cv"] for r in front)
    cand = [r for r in front if r["auc_cv"] >= best - tol]
    if use_spec:
        cand.sort(key=lambda r: (round(r["f2"], 4), len(r["idx"]), r["f4"], -r["auc_cv"]))
    elif use_si:
        cand.sort(key=lambda r: (round(r["f_si"], 3), len(r["idx"]), r["f4"], -r["auc_cv"]))
    elif use_stab:
        cand.sort(key=lambda r: (round(r["f5"], 3), len(r["idx"]), r["f4"], -r["auc_cv"]))
    else:
        cand.sort(key=lambda r: (len(r["idx"]), r["f4"], -r["auc_cv"]))
    return cand[0]
