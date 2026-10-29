"""
sPTB paneli — GSE149440'ta geliştir, GSE59491'de BİR KEZ test et. Protokol: results/ptb/PROTOCOL.md (önceden yazıldı).
"""
import hashlib
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geo_utils import RAW, RESULTS, ROOT, pd, read_expression  # noqa: E402
from mogafs import bootstrap_stability, choose, fit_lr, run_mogafs  # noqa: E402
from panel_utils import boot_auc, delong_test, lasso_k, welch_t  # noqa: E402
from sklearn.linear_model import LogisticRegressionCV  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402
from sklearn.model_selection import RepeatedStratifiedKFold  # noqa: E402

OUT = RESULTS / "ptb"
CFG = dict(K=200, max_k=5, pop_size=100, n_gen=60, tol=0.03, seed=2026, folds=5, repeats=3)
METHODS = ["XAI-MOGAFS", "LASSO5", "TTEST5"]


def load():
    uni = pd.read_csv(ROOT / "data/processed/gene_universe.tsv", sep="\t", dtype={"entrez": str}).set_index("entrez")
    wb = uni.index[uni.wholeblood_all]
    # GSE149440: 17–23 hf penceresinde kadının en erken örneği (PE hariç)
    s = pd.read_csv(RESULTS / "inventory/GSE149440_samples.tsv", sep="\t", dtype={"woman_id": str})
    s = s[(s.ga_wk >= 17) & (s.ga_wk < 24) & s.group.isin(["Control", "sPTD", "PPROM"])]
    s = s.sort_values(["woman_id", "ga_wk"]).drop_duplicates("woman_id")
    E = read_expression(RAW / "GSE149440/GSE149440_series_matrix.txt.gz")[s.gsm]
    ent = E.index.str.replace("_at", "", regex=False)
    E = E.groupby(ent.values).median()
    E = E.loc[E.index.intersection(wb)].T
    E.index = s.woman_id.values
    yd = s.group.isin(["sPTD", "PPROM"]).astype(int).values
    # GSE59491: zaman noktası 1 (Adım 3 gen düzeyi matrisi)
    T = pd.read_csv(ROOT / "data/processed/expr/GSE59491.tsv.gz", sep="\t", dtype={"entrez": str}).set_index("entrez").T
    w = pd.read_csv(ROOT / "data/processed/selection/women_all.tsv", sep="\t", dtype=str)
    w = w[w.cohort == "GSE59491"].set_index("woman_id").loc[T.index]
    genes = sorted(set(E.columns) & set(T.columns), key=int)
    return E[genes], yd, s, T[genes], w, genes, uni.symbol.to_dict()


def zs(M, mu=None, sd=None):
    mu = M.mean(0) if mu is None else mu
    sd = M.std(0, ddof=1) if sd is None else sd
    sd = np.where(sd == 0, 1, sd)
    return (M - mu) / sd, mu, sd


def run_task(task, X, y, cfg):
    name, tr, te = task
    Ztr, mu, sd = zs(X[tr])
    t = welch_t(Ztr, y[tr])
    top = np.argsort(-np.abs(t))[: cfg["K"]]
    stab = bootstrap_stability(Ztr[:, top], y[tr], seed=cfg["seed"])
    arc, info = run_mogafs(Ztr[:, top], y[tr], np.zeros(len(tr), int), use_spec=False, max_k=cfg["max_k"],
                           pop_size=cfg["pop_size"], n_gen=cfg["n_gen"], seed=cfg["seed"], stab=stab)
    best = choose(arc, False, cfg["tol"], use_stab=True)
    sel = {"XAI-MOGAFS": top[list(best["idx"])], "LASSO5": lasso_k(Ztr, y[tr], cfg["max_k"]), "TTEST5": top[: cfg["max_k"]]}
    out = {"name": name, "te": None if te is None else te.tolist(), "scores": {}, "models": {}}
    for m in METHODS:
        S = np.asarray(sel[m])
        lr = fit_lr(Ztr[:, S], y[tr])
        out["models"][m] = {"idx": S.tolist(), "coef": lr.coef_[0].tolist(), "intercept": float(lr.intercept_[0])}
        if te is not None:
            Zte, _, _ = zs(X[te], mu, sd)
            out["scores"][m] = lr.decision_function(Zte[:, S])
    return out


def main():
    lines = []

    def say(m=""):
        print(m, flush=True)
        lines.append(str(m))

    E, yd, sdev, T, wt, genes, sym = load()
    yt = wt.group.isin(["sPTD", "PPROM"]).astype(int).values
    say(f"Geliştirme GSE149440 (17–23 hf): {len(yd)} kadın, {yd.sum()} sPTB (sPTD {int((sdev.group == 'sPTD').sum())}, "
        f"PPROM {int((sdev.group == 'PPROM').sum())}) / {(1 - yd).sum()} kontrol; hafta medyan {sdev.ga_wk.median():.1f}")
    say(f"Dış test GSE59491 (TP1): {len(yt)} kadın, {yt.sum()} sPTB / {(1 - yt).sum()} term; ortak gen {len(genes)}")
    X = E.values
    rskf = RepeatedStratifiedKFold(n_splits=CFG["folds"], n_repeats=CFG["repeats"], random_state=CFG["seed"])
    tasks = [(f"cv{i}", tr, te) for i, (tr, te) in enumerate(rskf.split(X, yd))] + [("FINAL", np.arange(len(yd)), None)]
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=min(14, len(tasks))) as ex:
        res = list(ex.map(run_task, tasks, [X] * len(tasks), [yd] * len(tasks), [CFG] * len(tasks)))
    say(f"Süre {time.time() - t0:.0f} sn")
    say("\nA. GELİŞTİRME — iç içe CV (GSE149440)")
    for m in METHODS:
        aucs = []
        for r in range(CFG["repeats"]):
            s = np.zeros(len(yd))
            for x in res[r * CFG["folds"]:(r + 1) * CFG["folds"]]:
                s[x["te"]] = x["scores"][m]
            aucs.append(roc_auc_score(yd, s))
        say(f"  {m:11s} AUC {np.mean(aucs):.3f} ({', '.join(f'{a:.3f}' for a in aucs)})")
    ridge = []
    for r in range(CFG["repeats"]):
        from sklearn.model_selection import StratifiedKFold, cross_val_predict
        from sklearn.pipeline import make_pipeline
        from sklearn.preprocessing import StandardScaler
        p = cross_val_predict(make_pipeline(StandardScaler(), LogisticRegressionCV(Cs=np.logspace(-5, -1, 8), max_iter=4000, scoring="roc_auc")),
                              X, yd, cv=StratifiedKFold(5, shuffle=True, random_state=r), method="predict_proba")[:, 1]
        ridge.append(roc_auc_score(yd, p))
    say(f"  ridge (tüm genler, referans) AUC {np.mean(ridge):.3f}")

    final = [x for x in res if x["name"] == "FINAL"][0]
    locked = {m: {"genes": [genes[i] for i in md["idx"]], "symbols": [sym.get(genes[i], genes[i]) for i in md["idx"]],
                  "coef": md["coef"], "intercept": md["intercept"]} for m, md in final["models"].items()}
    path = OUT / "locked_models.json"
    path.write_text(json.dumps(locked, indent=2, ensure_ascii=False), encoding="utf-8")
    dig = hashlib.sha256(path.read_bytes()).hexdigest()
    (OUT / "locked_models.sha256").write_text(dig + "  locked_models.json\n", encoding="utf-8")
    say(f"\nB. KİLİT (GSE149440 tamamı) SHA-256 {dig[:16]}…")
    for m, md in locked.items():
        say(f"  {m:11s} " + ", ".join(f"{s} ({c:+.2f})" for s, c in zip(md["symbols"], md["coef"])))

    say("\nC. DIŞ TEST — GSE59491 (bir kez)")
    Zt, _, _ = zs(T.values)
    gpos = {g: i for i, g in enumerate(genes)}
    sc = {}
    for m, md in locked.items():
        s = md["intercept"] + Zt[:, [gpos[g] for g in md["genes"]]] @ np.array(md["coef"])
        sc[m] = s
        lo, hi = boot_auc(yt, s)
        sub = []
        for g in ["sPTD", "PPROM"]:
            mk = wt.group.isin([g, "Control"]).values
            sub.append(f"{g} {roc_auc_score((wt.group.values[mk] == g).astype(int), s[mk]):.3f}")
        say(f"  {m:11s} AUC {roc_auc_score(yt, s):.3f} [%95 GA {lo:.3f}–{hi:.3f}] | " + " | ".join(sub))
    for b in ["LASSO5", "TTEST5"]:
        a1, a2, d, z, p = delong_test(yt, sc["XAI-MOGAFS"], sc[b])
        say(f"  DeLong XAI-MOGAFS vs {b}: {a1:.3f} vs {a2:.3f}, fark {d:+.3f}, p = {p:.3f}")
    (OUT / "report.txt").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
