"""Panel betiklerinin ortak yardımcıları (Welch t ön-filtresi, ≤k LASSO, bootstrap AUC GA)."""
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score


def welch_t(Z, y):
    a, b = Z[y == 1], Z[y == 0]
    return (a.mean(0) - b.mean(0)) / np.sqrt(a.var(0, ddof=1) / len(a) + b.var(0, ddof=1) / len(b))


def lasso_k(Z, y, k):
    """İkili arama: ≤k sıfırdan farklı katsayı veren en büyük C (L1-lojistik)."""
    lo, hi, best = -4.0, 1.0, None
    for _ in range(22):
        mid = (lo + hi) / 2
        m = LogisticRegression(penalty="l1", C=10 ** mid, solver="liblinear", max_iter=3000).fit(Z, y)
        nz = np.flatnonzero(m.coef_[0])
        if len(nz) <= k:
            if len(nz) >= 1 and (best is None or len(nz) >= len(best[0])):
                best = (nz[np.argsort(-np.abs(m.coef_[0][nz]))], mid)
            lo = mid
        else:
            hi = mid
    return best[0]


def boot_auc(y, s, n=2000, seed=0):
    """Tabakalı bootstrap %95 GA (vakalar ve kontroller ayrı yeniden örneklenir)."""
    rng = np.random.default_rng(seed)
    pos, neg = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    aucs = [roc_auc_score(np.r_[np.ones(len(pos)), np.zeros(len(neg))],
                          np.r_[s[rng.choice(pos, len(pos))], s[rng.choice(neg, len(neg))]]) for _ in range(n)]
    return np.percentile(aucs, [2.5, 97.5])


# ---------------------------------------------------------------- plazma cfRNA yardımcıları
COMPOSITION = {"COMP:RN7SL": r"^RN7SL", "COMP:RN7SK": r"^RN7SK", "COMP:RNA5": r"^RNA5",
               "COMP:RNY": r"^RNY", "COMP:SNOR": r"^(SNOR|SCARNA)", "COMP:RNU": r"^RNU", "COMP:MT": r"^MT-"}


def composition_features(counts):
    """counts: gen_adı x örnek ham sayım -> örnek x bileşim özelliği (log2 kütüphane yüzdesi)."""
    import pandas as pd
    tot = counts.sum(0)
    return pd.DataFrame({k: np.log2(counts[counts.index.str.contains(rx, regex=True)].sum(0) / tot * 100 + 0.01)
                         for k, rx in COMPOSITION.items()})


def cfrna_gene_level(raw, gene_info_path):
    """raw: GSE192902 sayım tablosu (gene_name, gene_num, örnekler...). Adım 3 ile aynı eşleme:
    Ensembl (eşsiz) -> Entrez, yoksa eşsiz sembol; aynı Entrez'de sayım toplamı; log2(CPM+1).
    Dönüş: Entrez x örnek."""
    import pandas as pd
    gi = pd.read_csv(gene_info_path, sep="\t", usecols=["GeneID", "Symbol", "dbXrefs"], dtype=str)
    sc = gi.Symbol.value_counts()
    sym2id = {s: i for i, s in zip(gi.GeneID, gi.Symbol) if sc[s] == 1}
    ens = gi.dbXrefs.str.extractall(r"Ensembl:(ENSG\d+)")[0].reset_index()
    ens["GeneID"] = gi.GeneID.values[ens.level_0]
    ec = ens[0].value_counts()
    ens2id = {e: g for e, g in zip(ens[0], ens.GeneID) if ec[e] == 1}
    ent = pd.Series([ens2id.get(e) or sym2id.get(s) for e, s in zip(raw.gene_num, raw.gene_name)], index=raw.index)
    counts = raw.drop(columns=["gene_name", "gene_num"]).assign(entrez=ent).dropna(subset=["entrez"]).groupby("entrez").sum()
    return np.log2(counts / counts.sum(0) * 1e6 + 1)


def delong_test(y, s1, s2):
    """İki eşli skorun AUC farkı için DeLong testi (Sun & Xu 2014 hızlı algoritması).
    Dönüş: (auc1, auc2, fark, z, iki yönlü p)."""
    from scipy import stats

    def midrank(x):
        order = np.argsort(x)
        xs = x[order]
        n = len(x)
        t = np.zeros(n)
        i = 0
        while i < n:
            j = i
            while j < n and xs[j] == xs[i]:
                j += 1
            t[i:j] = 0.5 * (i + j - 1) + 1
            i = j
        out = np.empty(n)
        out[order] = t
        return out

    y = np.asarray(y).astype(bool)
    m, n = y.sum(), (~y).sum()
    aucs, v01, v10 = [], [], []
    for s in (np.asarray(s1, float), np.asarray(s2, float)):
        pos, neg = s[y], s[~y]
        tx, ty, tz = midrank(pos), midrank(neg), midrank(np.r_[pos, neg])
        auc = (tz[:m].sum() - m * (m + 1) / 2) / (m * n)
        aucs.append(auc)
        v01.append((tz[:m] - tx) / n)
        v10.append(1 - (tz[m:] - ty) / m)
    sx, sy = np.cov(np.vstack(v01)), np.cov(np.vstack(v10))
    var = (sx[0, 0] + sx[1, 1] - 2 * sx[0, 1]) / m + (sy[0, 0] + sy[1, 1] - 2 * sy[0, 1]) / n
    diff = aucs[0] - aucs[1]
    z = diff / np.sqrt(var) if var > 0 else np.nan
    return aucs[0], aucs[1], diff, z, 2 * stats.norm.sf(abs(z))
