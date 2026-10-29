"""
Adım 4c — Pozitif kontrol: aynı işleme hattı, sinyal BEKLENEN yerde (plazma cfRNA, Moufarrej 2022)
PE sinyali buluyor mu? Ve tam kanda PE sinyali kohort içinde herhangi bir yerde var mı?
Her kohortta kohort içi 5x3 CV: ridge (tüm genler) ve tek-gen FDR özeti. Betimleyici.
"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from geo_utils import RESULTS, ROOT, pd  # noqa: E402
from scipy import stats  # noqa: E402
from sklearn.linear_model import LogisticRegressionCV  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402
from sklearn.model_selection import StratifiedKFold, cross_val_predict  # noqa: E402
from sklearn.pipeline import make_pipeline  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

w = pd.read_csv(ROOT / "data/processed/selection/women_all.tsv", sep="\t", dtype=str)
lines = []
def say(m=""):
    print(m, flush=True); lines.append(m)

tests = [("GSE85307", "GSE85307", "PE", None, "tam kan, eğitim"),
         ("GSE149440", "GSE149440", "PE", None, "tam kan"),
         ("GSE108497", "GSE108497", "PE", "SLE", "tam kan, SLE içi"),
         ("GSE192902", "GSE192902_cfRNA_discovery", "PE", None, "plazma cfRNA (pozitif kontrol)"),
         ("GSE192902", "GSE192902_cfRNA_val1", "PE", None, "plazma cfRNA (pozitif kontrol)")]
say("KOHORT İÇİ PE vs KONTROL — ridge (tüm genler, 5x3 CV) ve tek-gen Welch t")
for coh, f, g1, bg, lab in tests:
    E = pd.read_csv(ROOT / f"data/processed/expr/{f}.tsv.gz", sep="\t", dtype={"entrez": str}).set_index("entrez").T
    ww = w[(w.cohort == coh)].set_index("woman_id").loc[E.index]
    m = ww.group.isin([g1, "Control"]).values.copy()
    if bg:
        m &= (ww.background == bg).values
    X, y = E.values[m], (ww.group.values[m] == g1).astype(int)
    if "cfRNA" in f:
        X = X[:, (X > 1).mean(0) >= 0.5]          # düşük sayımlı genler: yarıdan fazlasında log2CPM>1
    X = X[:, X.std(0) > 0]
    aucs = []
    for r in range(3):
        mdl = make_pipeline(StandardScaler(), LogisticRegressionCV(Cs=np.logspace(-5, -1, 8), max_iter=4000, scoring="roc_auc"))
        p = cross_val_predict(mdl, X, y, cv=StratifiedKFold(5, shuffle=True, random_state=r), method="predict_proba")[:, 1]
        aucs.append(roc_auc_score(y, p))
    a, b = X[y == 1], X[y == 0]
    t = (a.mean(0) - b.mean(0)) / np.sqrt(a.var(0, ddof=1) / len(a) + b.var(0, ddof=1) / len(b))
    pv = 2 * stats.norm.sf(np.abs(t)); q = stats.false_discovery_control(pv)
    say(f"  {f:28s} {lab:32s} n={y.sum():3d}/{(1-y).sum():3d}  genler {X.shape[1]:6d}  ridge AUC {np.mean(aucs):.3f} ± {np.std(aucs):.3f}  "
        f"FDR<0.05: {(q<.05).sum():4d}  FDR<0.20: {(q<.2).sum():4d}  p<0.001: {(pv<.001).sum():4d} (null ≈ {X.shape[1]*.001:.0f})")
(RESULTS / "panel" / "positive_control_report.txt").write_text("\n".join(lines), encoding="utf-8")
