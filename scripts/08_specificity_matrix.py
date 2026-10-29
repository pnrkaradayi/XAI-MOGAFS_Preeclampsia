"""
Adım 5 — Özgüllük matrisi: son paneller bağımsız kohortlarda "alarm başka sendromlarda da çalıyor mu?"

Uygulama: her kohortta panel genleri kohort içinde (tüm seçili kadınlar, etiketsiz) z-ölçeklenir;
skor = kesişim + Σ β_j z_j (katsayılar GSE85307 son modelinden, DEĞİŞTİRİLMEDEN).

Kohort rolleri:
  GSE142974, GSE59491, GSE108497   BAĞIMSIZ özgüllük testleri
  GSE149440                        GELİŞTİRME (f2 burada hesaplandı -> iyimser; ayrı gösterilir)
  GSE85307                         eğitim; satır dış-CV AUC'sinden (Adım 4) gelir

Her hücre: AUC(grup vs o kohortun kontrolleri) + stratified bootstrap %95 GA (2000).
GSE108497'de kontrol = SLE komplikasyonsuz (arka plan eşli, ana); sağlıklı gebe kontrol ikincil.
Ek: PE vs diğer sendrom (aynı kohort içinde), yalnız eğitim penceresi (10–18.9 hf) AUC'si ve
hafta düzeltmeli OR (vaka ~ skor_z + örnek haftası; skor kohort içinde z).

Çıktılar: results/specificity/matrix.tsv, pairwise.tsv, report.txt
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geo_utils import RESULTS, ROOT, pd  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402
import statsmodels.api as sm  # noqa: E402

OUT = RESULTS / "specificity"
OUT.mkdir(parents=True, exist_ok=True)
COHORTS = [("GSE142974", "bağımsız"), ("GSE59491", "bağımsız"), ("GSE108497", "bağımsız"), ("GSE149440", "geliştirme")]
GROUPS = ["PE", "sPTD", "PPROM", "FGR_SGA_FD"]
lines = []


def say(msg=""):
    print(msg, flush=True)
    lines.append(str(msg))


def boot_ci(y, s, n=2000, seed=0):
    rng = np.random.default_rng(seed)
    pos, neg = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    a = [roc_auc_score(np.r_[np.ones(len(pos)), np.zeros(len(neg))],
                       np.r_[s[rng.choice(pos, len(pos))], s[rng.choice(neg, len(neg))]]) for _ in range(n)]
    return np.percentile(a, [2.5, 97.5])


def adj_or(case, score, ga):
    """Hafta düzeltmeli OR (skorun 1 SD artışı başına)."""
    z = (score - score.mean()) / score.std(ddof=1)
    Xd = sm.add_constant(np.column_stack([z, ga]))
    try:
        f = sm.Logit(case, Xd).fit(disp=0)
        b, se = f.params[1], f.bse[1]
        return np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se), f.pvalues[1]
    except Exception:
        return (np.nan,) * 4


def main():
    panels = json.loads((RESULTS / "panel" / "final_panels.json").read_text(encoding="utf-8"))
    women = pd.read_csv(ROOT / "data" / "processed" / "selection" / "women_all.tsv", sep="\t", dtype=str)
    uni = pd.read_csv(ROOT / "data" / "processed" / "gene_universe.tsv", sep="\t", dtype={"entrez": str}).set_index("entrez")
    cv = pd.read_csv(RESULTS / "panel" / "outer_cv_summary.tsv", sep="\t").set_index("method")
    rows, pair = [], []
    for method, p in panels.items():
        genes, coef = p["genes"], np.array(p["coef"])
        low = [uni.loc[g, "symbol"] for g in genes if uni.loc[g, "GSE108497_detect_rate"] < 0.5]
        say(f"=== {method}: {', '.join(p['symbols'])}" + (f"  (GSE108497'de düşük saptama: {low})" if low else ""))
        r = cv.loc[method]
        rows.append({"method": method, "cohort": "GSE85307", "role": "eğitim (dış CV)", "control": "Control",
                     "group": "PE", "n_case": 47, "n_ctrl": 110, "auc": r.pooled_outer_auc,
                     "ci_lo": r.ci_lo, "ci_hi": r.ci_hi})
        for cohort, role in COHORTS:
            E = pd.read_csv(ROOT / "data" / "processed" / "expr" / f"{cohort}.tsv.gz", sep="\t",
                            dtype={"entrez": str}).set_index("entrez").loc[genes].T
            w = women[women.cohort == cohort].set_index("woman_id").loc[E.index]
            Z = (E - E.mean()) / E.std(ddof=1)
            score = p["intercept"] + Z.values @ coef
            ga = w.ga_sample_wk.astype(float).values
            tw = (w.in_train_window == "True").values
            ctrl_sets = {"Control": (w.group == "Control").values}
            if cohort == "GSE108497":
                ctrl_sets = {"Control (SLE)": ((w.group == "Control") & (w.background == "SLE")).values,
                             "Control (sağlıklı)": ((w.group == "Control") & (w.background == "healthy")).values}
            for cname, cmask in ctrl_sets.items():
                for g in GROUPS:
                    gmask = (w.group == g).values
                    if gmask.sum() < 5:
                        continue
                    sel = gmask | cmask
                    yy, ss = gmask[sel].astype(int), score[sel]
                    auc = roc_auc_score(yy, ss)
                    lo, hi = boot_ci(yy, ss)
                    rec = {"method": method, "cohort": cohort, "role": role, "control": cname, "group": g,
                           "n_case": int(gmask.sum()), "n_ctrl": int(cmask.sum()), "auc": auc, "ci_lo": lo, "ci_hi": hi}
                    tsel = sel & tw
                    if gmask[tsel].sum() >= 5 and cmask[tsel].sum() >= 5:
                        rec["auc_train_window"] = roc_auc_score(gmask[tsel].astype(int), score[tsel])
                        rec["n_case_tw"] = int(gmask[tsel].sum())
                    o, ol, oh, pv = adj_or(yy, ss, ga[sel])
                    rec.update({"or_ga_adj": o, "or_lo": ol, "or_hi": oh, "or_p": pv})
                    rows.append(rec)
            # PE vs diğer sendromlar (aynı kohort)
            if (w.group == "PE").sum() >= 5:
                for g in ["sPTD", "PPROM", "FGR_SGA_FD"]:
                    gmask = (w.group == g).values
                    if gmask.sum() < 5:
                        continue
                    pm = (w.group == "PE").values
                    sel = pm | gmask
                    yy, ss = pm[sel].astype(int), score[sel]
                    lo, hi = boot_ci(yy, ss)
                    pair.append({"method": method, "cohort": cohort, "role": role, "comparison": f"PE vs {g}",
                                 "n_pe": int(pm.sum()), "n_other": int(gmask.sum()),
                                 "auc": roc_auc_score(yy, ss), "ci_lo": lo, "ci_hi": hi})
    M = pd.DataFrame(rows)
    M.to_csv(OUT / "matrix.tsv", sep="\t", index=False)
    P = pd.DataFrame(pair)
    P.to_csv(OUT / "pairwise.tsv", sep="\t", index=False)

    say("\nÖZGÜLLÜK MATRİSİ — AUC(grup vs kontrol) [%95 GA]   (0.5 = alarm çalmıyor)")
    for method in panels:
        say(f"\n  {method}")
        sub = M[M.method == method]
        for _, r in sub.iterrows():
            extra = ""
            if pd.notna(r.get("auc_train_window", np.nan)):
                extra += f" | 10–18.9 hf: {r.auc_train_window:.2f} (n={int(r.n_case_tw)})"
            if pd.notna(r.get("or_ga_adj", np.nan)):
                extra += f" | hafta-düzeltmeli OR {r.or_ga_adj:.2f} [{r.or_lo:.2f}–{r.or_hi:.2f}]"
            say(f"    {r.cohort:9s} {r.role:15s} {r.group:10s} vs {r['control']:18s} n={r.n_case:3d}/{r.n_ctrl:3d}  "
                f"AUC {r.auc:.2f} [{r.ci_lo:.2f}–{r.ci_hi:.2f}]{extra}")
        ps = P[P.method == method]
        for _, r in ps.iterrows():
            say(f"    {r.cohort:9s} {r.role:15s} {r.comparison:22s} n={r.n_pe}/{r.n_other}  AUC {r.auc:.2f} [{r.ci_lo:.2f}–{r.ci_hi:.2f}]")
    (OUT / "report.txt").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
