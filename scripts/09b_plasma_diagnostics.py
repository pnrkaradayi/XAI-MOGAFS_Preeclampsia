"""
Adım 4 (plazma) tanılaması — PE sinyali teknik bir iz mi?
Kilitli panellerde RN7SL (7SL/SRP RNA) psödogenleri baskın. Bu betik sınar:
  1. Kütüphanedeki RN7SL* ve diğer yapısal küçük RNA psödogeni payı PE / kontrol / hafta dilimine göre
  2. Gen türüne göre kohort içi sinyal (ridge 5x3 CV): tüm genler / yalnız protein kodlayan / yalnız psödogen
  3. Hafta dilimi içinde AUC (13–20 hf dengesizliği)
Betimleyicidir. Ham sayımlar (Discovery, Val-1) kullanılır; Val-2 açılmaz.
"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from geo_utils import RAW, RESULTS, ROOT, pd  # noqa: E402
from scipy import stats  # noqa: E402
from sklearn.linear_model import LogisticRegressionCV  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402
from sklearn.model_selection import StratifiedKFold, cross_val_predict  # noqa: E402
from sklearn.pipeline import make_pipeline  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

lines = []
def say(m=""):
    print(m, flush=True); lines.append(str(m))

w = pd.read_csv(ROOT / "data/processed/selection/women_all.tsv", sep="\t", dtype=str)
u = pd.read_csv(ROOT / "data/processed/gene_universe.tsv", sep="\t", dtype={"entrez": str}).set_index("entrez")
STRUCT = r"^(RN7SL|RN7SK|RNA5-8S|RNA5S|RNU\d|RNY|RNA18S|RNA28S|RNA45S|SNOR|SCARNA|MT-)"

def ridge_auc(X, y):
    a = []
    for r in range(3):
        m = make_pipeline(StandardScaler(), LogisticRegressionCV(Cs=np.logspace(-5, -1, 8), max_iter=4000, scoring="roc_auc"))
        p = cross_val_predict(m, X, y, cv=StratifiedKFold(5, shuffle=True, random_state=r), method="predict_proba")[:, 1]
        a.append(roc_auc_score(y, p))
    return np.mean(a), np.std(a)

for set_, fname in [("cfRNA_discovery", "GSE192902_counts_Discovery_postQC.csv.gz"),
                    ("cfRNA_val1", "GSE192902_counts_Validation1_postQC.csv.gz")]:
    ww = w[(w.cohort == "GSE192902") & (w.set == set_)].set_index("sample_id")
    c = pd.read_csv(RAW / "GSE192902" / fname).set_index("gene_name")[list(ww.index)]
    tot = c.sum(0)
    struct = c.index.str.contains(STRUCT, regex=True)
    rn7sl = c.index.str.startswith("RN7SL")
    y = (ww.group == "PE").astype(int).values
    late = (ww.ga_bin == "13–20").values
    say(f"=== {set_}: {len(y)} kadın ({y.sum()} PE); {len(c)} gen satırı; RN7SL* satırı {rn7sl.sum()}, yapısal küçük RNA satırı {struct.sum()}")
    for lab, mask in [("RN7SL* payı (%)", rn7sl), ("yapısal küçük RNA payı (%)", struct)]:
        frac = (c[mask].sum(0) / tot * 100).values
        say(f"  {lab:28s} PE medyan {np.median(frac[y==1]):6.2f} | kontrol {np.median(frac[y==0]):6.2f} | "
            f"AUC(PE vs kontrol) {roc_auc_score(y, frac):.2f} p={stats.mannwhitneyu(frac[y==1], frac[y==0]).pvalue:.1e} | "
            f"kontrollerde ≤12 vs 13–20 hf: {np.median(frac[(y==0)&~late]):.2f} / {np.median(frac[(y==0)&late]):.2f}")
    say(f"  toplam sayım (kütüphane): PE medyan {np.median(tot[y==1]):.0f} | kontrol {np.median(tot[y==0]):.0f} | AUC {roc_auc_score(y, tot.values):.2f}")
    E = pd.read_csv(ROOT / f"data/processed/expr/GSE192902_{set_}.tsv.gz", sep="\t", dtype={"entrez": str}).set_index("entrez")
    E = E[ww.woman_id.tolist()].T
    sym = u.symbol.reindex(E.columns).fillna("")
    typ = u.type_of_gene.reindex(E.columns).fillna("")
    pseudo = sym.str.contains(STRUCT, regex=True) | (typ == "pseudo")
    pc = (typ == "protein-coding") & ~pseudo
    if set_ == "cfRNA_discovery":
        for lab, cols in [("tüm genler", np.ones(len(sym), bool)), ("yalnız protein kodlayan", pc.values),
                          ("yalnız psödogen + yapısal küçük RNA", pseudo.values)]:
            a, s = ridge_auc(E.values[:, cols], y)
            say(f"  ridge AUC ({lab}, {cols.sum()} gen): {a:.3f} ± {s:.3f}")
        Xpc = E.values[:, pc.values]
        for b, lab in [(False, "≤12 hf"), (True, "13–20 hf")]:
            m = late == b
            if y[m].sum() >= 5:
                a, s = ridge_auc(Xpc[m], y[m])
                say(f"  protein kodlayan, yalnız {lab}: n={y[m].sum()}/{(1-y[m]).sum()}  ridge AUC {a:.3f} ± {s:.3f}")
        t = stats.ttest_ind(Xpc[y == 1], Xpc[y == 0], equal_var=False)
        q = stats.false_discovery_control(np.nan_to_num(t.pvalue, nan=1))
        top = np.argsort(t.pvalue)[:12]
        say(f"  protein kodlayan tek-gen: FDR<0.05 {int((q<.05).sum())}, FDR<0.20 {int((q<.2).sum())}; en güçlü: " +
            ", ".join(f"{sym[pc].iloc[i]}({'+' if t.statistic[i]>0 else '−'})" for i in top))
    say()
(RESULTS / "plasma" / "plasma_diagnostics_report.txt").write_text("\n".join(lines), encoding="utf-8")

# ---------------------------------------------------------------- 4. trombosit hipotezi
# 7SL RNA trombositlerde çok bol. Trombosite özgü mRNA'lar da PE'de yüksek mi ve RN7SL payıyla birlikte mi değişiyor?
PLT = ["PPBP", "PF4", "PF4V1", "ITGA2B", "GP9", "GP1BB", "GP6", "TUBB1", "CLU", "NRGN", "SPARC", "TREML1", "MYL9",
       "TPM4", "CAVIN2", "GNG11", "CLEC1B", "SELP", "TSPAN9", "MMRN1"]
say("4. TROMBOSİT HİPOTEZİ (trombosit mRNA skoru = genlerin log2CPM z ortalaması)")
for set_, fname in [("cfRNA_discovery", "GSE192902_counts_Discovery_postQC.csv.gz"),
                    ("cfRNA_val1", "GSE192902_counts_Validation1_postQC.csv.gz")]:
    ww = w[(w.cohort == "GSE192902") & (w.set == set_)].set_index("sample_id")
    c = pd.read_csv(RAW / "GSE192902" / fname).set_index("gene_name")[list(ww.index)]
    cpm = np.log2(c / c.sum(0) * 1e6 + 1)
    have = [g for g in PLT if g in cpm.index]
    z = cpm.loc[have].sub(cpm.loc[have].mean(1), axis=0).div(cpm.loc[have].std(1).replace(0, 1), axis=0)
    plt_score = z.mean(0).values
    rn7 = (c[c.index.str.startswith("RN7SL")].sum(0) / c.sum(0) * 100).values
    y = (ww.group == "PE").astype(int).values
    r = stats.spearmanr(plt_score, rn7)
    say(f"  {set_}: {len(have)}/{len(PLT)} trombosit geni mevcut; skor AUC(PE vs kontrol) {roc_auc_score(y, plt_score):.2f} "
        f"(p={stats.mannwhitneyu(plt_score[y==1], plt_score[y==0]).pvalue:.1e}); trombosit skoru ~ RN7SL payı Spearman ρ={r.statistic:.2f} (p={r.pvalue:.1e})")
    # RN7SL, trombosit skoruna göre düzeltildiğinde PE'yi hâlâ ayırıyor mu?
    import statsmodels.api as sm
    Xd = sm.add_constant(np.column_stack([(rn7 - rn7.mean()) / rn7.std(), (plt_score - plt_score.mean()) / plt_score.std()]))
    f = sm.Logit(y, Xd).fit(disp=0)
    say(f"    lojistik PE ~ RN7SL_z + trombosit_z: OR(RN7SL) {np.exp(f.params[1]):.2f} p={f.pvalues[1]:.3f}; "
        f"OR(trombosit) {np.exp(f.params[2]):.2f} p={f.pvalues[2]:.3f}")
(RESULTS / "plasma" / "plasma_diagnostics_report.txt").write_text("\n".join(lines), encoding="utf-8")
