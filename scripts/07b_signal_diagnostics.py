"""
Adım 4b — Sinyal tanılaması: GSE85307'de PE vs kontrol sinyali gerçekten yok mu, yoksa gizleniyor mu?

Adım 4'te tüm yöntemler dış CV'de AUC ≈ 0.4–0.5 verdi ve FDR<0.20 gen yok. Bu betik dört alternatif
açıklamayı sınar:
  A. Parti (hibridizasyon tarihi) / kovaryat gürültüsü: PCA–değişken ilişkisi; kovaryat düzeltmeli
     tek-gen analizi (ifade ~ PE + parti + hafta + etnik köken + astım) ve p-değeri dağılımı (λ).
  B. Etiketlerin geçerliliği: klinik risk faktörleri (BMI, yaş, etnik köken) PE'yi öngörüyor mu?
  C. Çok değişkenli sinyal: tüm genlerle ridge, ilk 20 PC, random forest (iç içe olmayan, 5x3 tekrarlı
     CV; seçim yok, yalnız ham sinyal var mı sorusu); parti merkezlemeli ve merkezlemesiz.
  D. miRNA katmanı ve literatür imzaları (nötrofil granül / olgunlaşmamış nötrofil — Hong 2019;
     interferon; plazma hücresi) — kohortlar arası.
Bu betik betimleyicidir; panel seçimini değiştirmez.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geo_utils import RESULTS, ROOT, pd  # noqa: E402
from scipy import stats  # noqa: E402
from sklearn.decomposition import PCA  # noqa: E402
from sklearn.ensemble import RandomForestClassifier  # noqa: E402
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV  # noqa: E402
from sklearn.model_selection import RepeatedStratifiedKFold, cross_val_predict  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402
from sklearn.pipeline import make_pipeline  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

OUT = RESULTS / "panel"
lines = []


def say(msg=""):
    print(msg, flush=True)
    lines.append(str(msg))


def cv_auc(model, X, y, seed=0, reps=3):
    aucs = []
    for r in range(reps):
        cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=1, random_state=seed + r)
        p = cross_val_predict(model, X, y, cv=cv, method="predict_proba")[:, 1]
        aucs.append(roc_auc_score(y, p))
    return np.mean(aucs), np.std(aucs)


def ols_t(Y, design, col=1):
    """Her gen için OLS; `col` sütunundaki katsayının t ve p değeri."""
    Xd = np.asarray(design, float)
    beta, *_ = np.linalg.lstsq(Xd, Y, rcond=None)
    resid = Y - Xd @ beta
    dof = Xd.shape[0] - np.linalg.matrix_rank(Xd)
    s2 = (resid ** 2).sum(0) / dof
    cov = np.linalg.pinv(Xd.T @ Xd)[col, col]
    t = beta[col] / np.sqrt(s2 * cov)
    return t, 2 * stats.t.sf(np.abs(t), dof)


def summarize_p(p, label):
    q = stats.false_discovery_control(p)
    lam = np.median(stats.chi2.isf(p, 1)) / stats.chi2.ppf(0.5, 1)
    say(f"  {label:44s} λ={lam:.2f}  p<0.001: {(p < 0.001).sum():4d} (null ≈ {len(p) * .001:.0f})  "
        f"FDR<0.05: {(q < .05).sum():4d}  FDR<0.20: {(q < .2).sum():4d}")


def main():
    w = pd.read_csv(ROOT / "data" / "processed" / "selection" / "women_all.tsv", sep="\t", dtype=str)
    u = pd.read_csv(ROOT / "data" / "processed" / "gene_universe.tsv", sep="\t", dtype={"entrez": str}).set_index("entrez")
    meta = pd.read_csv(RESULTS / "inventory" / "GSE85307_samples.tsv", sep="\t", dtype=str)
    meta = meta[(meta.layer == "mRNA")].set_index("woman_id")
    E = pd.read_csv(ROOT / "data" / "processed" / "expr" / "GSE85307.tsv.gz", sep="\t",
                    dtype={"entrez": str}).set_index("entrez").T
    wm = w[w.cohort == "GSE85307"].set_index("woman_id").loc[E.index]
    meta = meta.loc[E.index]
    y = (wm.group == "PE").astype(int).values
    ga = wm.ga_sample_wk.astype(float).values
    bmi = wm.bmi.astype(float).values
    age = wm.age_y.astype(float).values
    batch = meta.batch.values
    eth = meta.ethnicity.values
    asthma = meta.asthma.values

    # ------------------------------------------------------------------ A. parti ve kovaryatlar
    say("A. PARTİ VE KOVARYATLAR (GSE85307)")
    bt = pd.crosstab(pd.Series(batch, name="parti"), pd.Series(y, name="grup"))
    chi = stats.chi2_contingency(bt)
    say(f"  hibridizasyon partisi: {bt.shape[0]} parti; parti x PE ki-kare p = {chi.pvalue:.3f}")
    say("  " + bt.rename(columns={0: "kontrol", 1: "PE"}).T.to_string().replace("\n", "\n  "))
    X = E.values
    Xs = (X - X.mean(0)) / X.std(0, ddof=1)
    pca = PCA(n_components=10).fit(Xs)
    pcs = pca.transform(Xs)
    say("  İlk 10 PC: açıklanan varyans % " + " ".join(f"{v * 100:.1f}" for v in pca.explained_variance_ratio_))
    covs = {"parti": pd.get_dummies(batch).values.astype(float), "PE": y[:, None].astype(float),
            "hafta": ga[:, None], "etnik": pd.get_dummies(eth).values.astype(float),
            "astım": pd.get_dummies(asthma).values.astype(float), "BMI": bmi[:, None], "yaş": age[:, None]}
    say("  PC başına R² (değişken açıklama gücü):")
    header = "      " + "".join(f"{k:>8s}" for k in covs)
    say(header)
    for i in range(10):
        r2 = []
        for k, c in covs.items():
            ok = ~np.isnan(c).any(1)
            D = np.column_stack([np.ones(ok.sum()), c[ok]])
            b, *_ = np.linalg.lstsq(D, pcs[ok, i], rcond=None)
            res = pcs[ok, i] - D @ b
            r2.append(1 - res.var() / pcs[ok, i].var())
        say(f"  PC{i + 1:<3d}" + "".join(f"{v:8.2f}" for v in r2))
    say("  Tek-gen PE etkisi (p-değeri özeti):")
    one = np.ones((len(y), 1))
    _, p0 = ols_t(X, np.column_stack([one, y]))
    summarize_p(p0, "düzeltmesiz")
    Db = np.column_stack([one, y, pd.get_dummies(batch, drop_first=True).values.astype(float)])
    _, p1 = ols_t(X, Db)
    summarize_p(p1, "parti düzeltmeli")
    ok = ~np.isnan(bmi)
    Dfull = np.column_stack([one, y, pd.get_dummies(batch, drop_first=True).values.astype(float), ga,
                             pd.get_dummies(eth, drop_first=True).values.astype(float),
                             pd.get_dummies(asthma, drop_first=True).values.astype(float)])
    _, p2 = ols_t(X, Dfull)
    summarize_p(p2, "parti + hafta + etnik + astım düzeltmeli")
    Dpc = np.column_stack([one, y, pcs[:, :5]])
    _, p3 = ols_t(X, Dpc)
    summarize_p(p3, "ilk 5 PC düzeltmeli")
    # Permütasyon: gerçek etiketle p<0.001 sayısı, karıştırılmış etiket dağılımının neresinde?
    rng = np.random.default_rng(1)
    perm = []
    for _ in range(200):
        _, pp = ols_t(X, np.column_stack([one, rng.permutation(y)]))
        perm.append((pp < 0.001).sum())
    say(f"  Permütasyon (200): karıştırılmış etiketle p<0.001 gen sayısı medyan {np.median(perm):.0f} "
        f"[%95: {np.percentile(perm, 2.5):.0f}–{np.percentile(perm, 97.5):.0f}]; gerçek {(p0 < 0.001).sum()}")
    say()

    # ------------------------------------------------------------------ B. klinik etiket geçerliliği
    say("B. ETİKET GEÇERLİLİĞİ — klinik risk faktörleriyle PE öngörüsü (5x3 CV)")
    clin = pd.DataFrame({"bmi": bmi, "age": age}).assign(**{f"eth_{k}": (eth == k).astype(float) for k in np.unique(eth)})
    ok = clin.notna().all(1).values
    m = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000))
    a, s = cv_auc(m, clin[ok].values, y[ok])
    say(f"  BMI + yaş + etnik köken: AUC {a:.3f} ± {s:.3f} (n={ok.sum()})")
    a, s = cv_auc(m, clin[["bmi"]][ok].values, y[ok])
    say(f"  yalnız BMI: AUC {a:.3f} ± {s:.3f}")
    say(f"  PE grubunda doğum haftası: medyan {np.median(wm.ga_delivery_wk.astype(float)[y == 1]):.1f}; "
        f"<37 hf doğum {(wm.ga_delivery_wk.astype(float)[y == 1] < 37).sum()}/{y.sum()}; "
        f"<34 hf {(wm.ga_delivery_wk.astype(float)[y == 1] < 34).sum()}/{y.sum()}")
    say()

    # ------------------------------------------------------------------ C. çok değişkenli ham sinyal
    say("C. ÇOK DEĞİŞKENLİ HAM SİNYAL (5x3 CV; ölçekleme CV içinde)")
    Xc = X.copy()
    for b in np.unique(batch):
        Xc[batch == b] -= Xc[batch == b].mean(0)
    hv = np.argsort(-X.var(0))[:2000]
    models = {
        "ridge, tüm genler": (make_pipeline(StandardScaler(), LogisticRegressionCV(Cs=np.logspace(-5, -1, 8), max_iter=3000, scoring="roc_auc")), None),
        "ilk 20 PC + lojistik": (make_pipeline(StandardScaler(), PCA(20), LogisticRegression(max_iter=2000)), None),
        "random forest, en değişken 2000 gen": (RandomForestClassifier(500, min_samples_leaf=3, random_state=0, n_jobs=-1), hv),
    }
    for name, (mdl, cols) in models.items():
        for lab, M in [("ham", X), ("parti merkezlemeli", Xc)]:
            Z = M if cols is None else M[:, cols]
            a, s = cv_auc(mdl, Z, y)
            say(f"  {name:38s} {lab:20s} AUC {a:.3f} ± {s:.3f}")
    say()

    # ------------------------------------------------------------------ D. miRNA ve literatür imzaları
    say("D. miRNA KATMANI (GSE85307)")
    ct = pd.read_csv(ROOT / "data" / "processed" / "expr" / "GSE85307_miRNA_Ct.tsv.gz", sep="\t").set_index("mirna_assay").T
    ct = ct.loc[E.index]
    keep = ct.notna().mean() >= 0.8
    c = ct.loc[:, keep]
    c = c.fillna(c.max() + 1)
    dct = c.sub(c.mean(axis=1), axis=0)
    expr_m = -dct.values
    _, pm = ols_t(expr_m, np.column_stack([one, y]))
    summarize_p(pm, f"miRNA ({keep.sum()} test, global ortalama normalize)")
    a, s = cv_auc(make_pipeline(StandardScaler(), LogisticRegressionCV(Cs=np.logspace(-4, 0, 8), max_iter=3000, scoring="roc_auc")), expr_m, y)
    say(f"  ridge, tüm miRNA: AUC {a:.3f} ± {s:.3f}")
    top = np.argsort(pm)[:5]
    say("  en güçlü 5 miRNA: " + ", ".join(f"{c.columns[i]} (p={pm[i]:.1e})" for i in top))
    say()

    say("D2. LİTERATÜR İMZALARI — skor = genlerin kohort içi z ortalaması; AUC(grup vs kontrol)")
    sigs = {
        "nötrofil granül (Hong 2019, olgunlaşmamış nötrofil)": ["ELANE", "MPO", "DEFA4", "CEACAM8", "LTF", "CAMP", "MMP8", "LCN2", "ARG1", "OLFM4", "BPI", "CTSG", "AZU1"],
        "interferon": ["IFI27", "IFI44L", "IFI44", "RSAD2", "ISG15", "IFIT1", "IFIT3", "MX1", "OAS1", "OAS3", "SIGLEC1", "HERC5"],
        "plazma hücresi": ["JCHAIN", "TNFRSF17", "MZB1", "TXNDC5", "IGLL5", "POU2AF1", "CD38"],
    }
    sym2id = {s: i for i, s in u.symbol.dropna().items()}
    tests = [("GSE85307", "PE", "Control", None), ("GSE149440", "PE", "Control", None),
             ("GSE149440", "sPTD", "Control", None), ("GSE149440", "PPROM", "Control", None),
             ("GSE59491", "sPTD", "Control", None), ("GSE59491", "PPROM", "Control", None),
             ("GSE142974", "sPTD", "Control", None),
             ("GSE108497", "PE", "Control", "SLE"), ("GSE108497", "FGR_SGA_FD", "Control", "SLE")]
    cache = {}
    rows = []
    for sname, gl in sigs.items():
        ids = [sym2id[g] for g in gl if g in sym2id and u.loc[sym2id[g], "wholeblood_all"]]
        say(f"  {sname} ({len(ids)}/{len(gl)} gen ortak evrende)")
        for coh, g1, g0, bg in tests:
            if coh not in cache:
                cache[coh] = pd.read_csv(ROOT / "data" / "processed" / "expr" / f"{coh}.tsv.gz", sep="\t",
                                         dtype={"entrez": str}).set_index("entrez").T
            M = cache[coh]
            ww = w[w.cohort == coh].set_index("woman_id").loc[M.index]
            Z = (M[ids] - M[ids].mean()) / M[ids].std(ddof=1)
            sc = Z.mean(axis=1).values
            m1 = (ww.group == g1).values.copy()
            m0 = (ww.group == g0).values.copy()
            if bg:
                m0 &= (ww.background == bg).values
                m1 &= (ww.background == bg).values
            sel = m1 | m0
            auc = roc_auc_score(m1[sel].astype(int), sc[sel])
            pmw = stats.mannwhitneyu(sc[m1], sc[m0]).pvalue
            rows.append({"signature": sname, "cohort": coh, "group": g1, "n": int(m1.sum()), "n_ctrl": int(m0.sum()),
                         "auc": auc, "p_mw": pmw})
            say(f"    {coh:9s} {g1:10s} vs {g0}{' (' + bg + ')' if bg else '':6s}  n={m1.sum():3d}/{m0.sum():3d}  AUC {auc:.2f}  p={pmw:.3f}")
    pd.DataFrame(rows).to_csv(OUT / "diagnostics_signatures.tsv", sep="\t", index=False)
    (OUT / "diagnostics_report.txt").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
