"""
Adım 4 — Paneli kur (tam kan mRNA, GSE85307) ve iç içe çapraz doğrulama ile değerlendir.

Tasarım (kullanıcı kararları, 2026-10-22):
  * Eğitim: GSE85307 (47 PE / 110 kontrol), 5 tam kan kohortunun ortak 17.621 geni.
  * Özgüllük geliştirme kohortu: GSE149440 (f2 burada hesaplanır; bu yüzden bağımsız test DEĞİLDİR).
    Bağımsız özgüllük testleri Adım 5'te: GSE59491, GSE108497, GSE142974.
  * miRNA katmanı burada yok: özgüllük kohortlarında miRNA ölçülmediği için f2 hesaplanamaz.
    Katman başına temel performans (Şekil 2) ayrı betikte.

Yöntemler (hepsi aynı dış katmanlarda, aynı ön-işleme ile):
  MOGAFS_spec  XAI-MOGAFS, 4 amaç (f1 f2 f3 f4), sınıf-çiftine özgü SHAP yönlendirmesi
  MOGAFS       XAI-MOGAFS, 3 amaç (f1 f3 f4), düz SHAP yönlendirmesi  (özgüllük katkısının ablasyonu)
  LASSO5       L1-lojistik, ≤5 sıfırdan farklı katsayı veren en büyük C; seçilenlerle L2 yeniden fit
  TTEST5       Welch t ile en güçlü 5 gen + L2 lojistik

İç içe CV: dış 5 kat (yalnız değerlendirme). Her dış eğitim katlamasında: z-ölçekleme (katlama
istatistikleri) -> Welch t ön-filtresi (K=200) -> seçim (iç CV 5x2) -> L2 lojistik fit. Dış test
katlamasına aynı ölçekleme ve model uygulanır. GSE149440 kendi içinde (etiketsiz) z-ölçeklenir.
Son panel: aynı prosedür 157 kadının tamamında.

Çıktılar: results/panel/
"""
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geo_utils import RESULTS, ROOT, pd  # noqa: E402
from mogafs import choose, fit_lr, run_mogafs  # noqa: E402
from panel_utils import boot_auc, lasso_k, welch_t  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402
from sklearn.model_selection import StratifiedKFold  # noqa: E402

CFG = dict(K=200, max_k=5, pop_size=100, n_gen=60, tol=0.03, lam=1.0, outer_folds=5, seed=2026,
           ga_edges=[13, 16])
OUT = RESULTS / "panel"
OUT.mkdir(parents=True, exist_ok=True)
METHODS = ["MOGAFS_spec", "MOGAFS", "LASSO5", "TTEST5"]


def load():
    w = pd.read_csv(ROOT / "data" / "processed" / "selection" / "women_all.tsv", sep="\t", dtype=str)
    u = pd.read_csv(ROOT / "data" / "processed" / "gene_universe.tsv", sep="\t", dtype={"entrez": str}).set_index("entrez")
    wb = u.index[u.wholeblood_all]

    def mat(name):
        return pd.read_csv(ROOT / "data" / "processed" / "expr" / f"{name}.tsv.gz", sep="\t",
                           dtype={"entrez": str}).set_index("entrez").loc[wb].T

    X = mat("GSE85307")
    w85 = w[w.cohort == "GSE85307"].set_index("woman_id").loc[X.index]
    D = mat("GSE149440")
    w49 = w[w.cohort == "GSE149440"].set_index("woman_id").loc[D.index]
    ZD = ((D - D.mean()) / D.std(ddof=1).replace(0, 1)).values
    return dict(X=X.values, y=(w85.group == "PE").astype(int).values, ga=w85.ga_sample_wk.astype(float).values,
                women=X.index.tolist(), genes=wb.tolist(), symbols=u.loc[wb, "symbol"].tolist(),
                ZD=ZD, dev_label=w49.group.values)


def run_task(task, d, cfg):
    fold, tr, te = task
    t0 = time.time()
    X, y, ga = d["X"], d["y"], d["ga"]
    mu, sd = X[tr].mean(0), X[tr].std(0, ddof=1)
    sd[sd == 0] = 1
    Ztr = (X[tr] - mu) / sd
    Zte = (X[te] - mu) / sd if te is not None else None
    t = welch_t(Ztr, y[tr])
    top = np.argsort(-np.abs(t))[: cfg["K"]]
    gbin = np.digitize(ga[tr], cfg["ga_edges"])
    sel, extra = {}, {}
    for arm, spec in [("MOGAFS_spec", True), ("MOGAFS", False)]:
        archive, info = run_mogafs(Ztr[:, top], y[tr], gbin, d["ZD"][:, top], d["dev_label"], use_spec=spec,
                                   max_k=cfg["max_k"], pop_size=cfg["pop_size"], n_gen=cfg["n_gen"],
                                   seed=cfg["seed"] + (0 if fold == "final" else fold), lam=cfg["lam"])
        best = choose(archive, spec, cfg["tol"])
        sel[arm] = top[list(best["idx"])]
        for r in archive:
            r["genes"] = [d["genes"][top[i]] for i in r["idx"]]
        extra[arm] = {"archive": archive, "info": info, "chosen": best}
    sel["LASSO5"] = lasso_k(Ztr, y[tr], cfg["max_k"])
    sel["TTEST5"] = top[: cfg["max_k"]]
    off_rows = np.isin(d["dev_label"], ["Control", "sPTD", "PPROM"])
    off_y = np.isin(d["dev_label"][off_rows], ["sPTD", "PPROM"]).astype(int)
    pe_rows = np.isin(d["dev_label"], ["PE", "sPTD", "PPROM"])
    pe_y = (d["dev_label"][pe_rows] == "PE").astype(int)
    rows, scores, models = [], {}, {}
    for m in METHODS:
        S = np.asarray(sel[m])
        lr = fit_lr(Ztr[:, S], y[tr])
        sdev = lr.decision_function(d["ZD"][:, S])
        rec = {"fold": fold, "method": m, "n_features": len(S),
               "genes": ";".join(d["genes"][i] for i in S), "symbols": ";".join(str(d["symbols"][i]) for i in S),
               "dev_auc_off": roc_auc_score(off_y, sdev[off_rows]),
               "dev_auc_pe_vs_comp": roc_auc_score(pe_y, sdev[pe_rows]),
               "dev_auc_pe_vs_ctrl": roc_auc_score((d["dev_label"][np.isin(d["dev_label"], ["PE", "Control"])] == "PE").astype(int),
                                                   sdev[np.isin(d["dev_label"], ["PE", "Control"])])}
        if m.startswith("MOGAFS"):
            rec["inner_auc_cv"] = extra[m]["chosen"]["auc_cv"]
        if te is not None:
            s = lr.decision_function(Zte[:, S])
            rec["outer_auc"] = roc_auc_score(y[te], s)
            scores[m] = s
        models[m] = {"genes": [d["genes"][i] for i in S], "symbols": [str(d["symbols"][i]) for i in S],
                     "coef": lr.coef_[0].tolist(), "intercept": float(lr.intercept_[0]),
                     "train_mean": mu[S].tolist(), "train_sd": sd[S].tolist()}
        rows.append(rec)
    return {"fold": fold, "te": None if te is None else te.tolist(), "rows": rows, "scores": scores,
            "models": models, "extra": extra, "top": top.tolist(), "seconds": time.time() - t0}


def main():
    lines = []

    def say(msg=""):
        print(msg, flush=True)
        lines.append(str(msg))

    d = load()
    y = d["y"]
    say(f"Eğitim: {len(y)} kadın ({y.sum()} PE / {(1 - y).sum()} kontrol), {len(d['genes'])} gen; "
        f"geliştirme kohortu GSE149440: {pd.Series(d['dev_label']).value_counts().to_dict()}")
    say(f"Ayarlar: {CFG}")
    # Betimleyici (seçimde KULLANILMAZ): tüm veride tek gen farkları
    Z = (d["X"] - d["X"].mean(0)) / d["X"].std(0, ddof=1)
    t = welch_t(Z, y)
    from scipy import stats
    p = 2 * stats.t.sf(np.abs(t), df=len(y) - 2)
    q = stats.false_discovery_control(p)
    say(f"Betimleyici: max |t| = {np.abs(t).max():.2f}; FDR<0.05 gen sayısı = {(q < 0.05).sum()}; "
        f"FDR<0.20 = {(q < 0.20).sum()}; p<0.001 = {(p < 0.001).sum()} (null beklentisi ≈ {len(p) * 0.001:.0f})")

    outer = StratifiedKFold(CFG["outer_folds"], shuffle=True, random_state=CFG["seed"])
    tasks = [(i, tr, te) for i, (tr, te) in enumerate(outer.split(d["X"], y))]
    tasks.append(("final", np.arange(len(y)), None))
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=len(tasks)) as ex:
        res = list(ex.map(run_task, tasks, [d] * len(tasks), [CFG] * len(tasks)))
    say(f"Süre: {time.time() - t0:.0f} sn (görev başına {[round(r['seconds']) for r in res]})")

    folds = [r for r in res if r["fold"] != "final"]
    final = [r for r in res if r["fold"] == "final"][0]
    df = pd.DataFrame([row for r in res for row in r["rows"]])
    df.to_csv(OUT / "outer_cv_folds.tsv", sep="\t", index=False)

    # Havuzlanmış dış-katman skorları -> AUC + bootstrap %95 GA
    say("\nDIŞ ÇAPRAZ DOĞRULAMA (GSE85307, PE vs kontrol) ve GELİŞTİRME KOHORTU (GSE149440)")
    summ = []
    for m in METHODS:
        s = np.zeros(len(y))
        for r in folds:
            s[r["te"]] = r["scores"][m]
        auc = roc_auc_score(y, s)
        lo, hi = boot_auc(y, s)
        fm = df[(df.method == m) & (df.fold != "final")]
        fin = df[(df.method == m) & (df.fold == "final")].iloc[0]
        summ.append({"method": m, "pooled_outer_auc": auc, "ci_lo": lo, "ci_hi": hi,
                     "fold_auc_mean": fm.outer_auc.mean(), "fold_auc_sd": fm.outer_auc.std(),
                     "n_features_mean": fm.n_features.mean(),
                     "dev_auc_off_mean_folds": fm.dev_auc_off.mean(),
                     "final_n": fin.n_features, "final_symbols": fin.symbols,
                     "final_dev_auc_off": fin.dev_auc_off, "final_dev_auc_pe_vs_comp": fin.dev_auc_pe_vs_comp,
                     "final_dev_auc_pe_vs_ctrl": fin.dev_auc_pe_vs_ctrl})
    summ = pd.DataFrame(summ)
    summ.to_csv(OUT / "outer_cv_summary.tsv", sep="\t", index=False)
    for _, r in summ.iterrows():
        say(f"  {r.method:12s} dış AUC {r.pooled_outer_auc:.3f} [%95 GA {r.ci_lo:.3f}–{r.ci_hi:.3f}]  "
            f"(katman ort. {r.fold_auc_mean:.3f} ± {r.fold_auc_sd:.3f}; ort. {r.n_features_mean:.1f} özellik) | "
            f"GSE149440 alarm AUC(komplikasyon vs kontrol) katman ort. {r.dev_auc_off_mean_folds:.3f}")
    say("\nSON PANELLER (157 kadının tamamında)")
    for _, r in summ.iterrows():
        say(f"  {r.method:12s} {int(r.final_n)} gen: {r.final_symbols}")
        say(f"  {'':12s} GSE149440 (geliştirme): komplikasyon vs kontrol {r.final_dev_auc_off:.3f} | "
            f"PE vs komplikasyon {r.final_dev_auc_pe_vs_comp:.3f} | PE vs kontrol {r.final_dev_auc_pe_vs_ctrl:.3f}")

    say("\nSEÇİM KARARLILIĞI (dış katmanlarda kaç kez seçildi)")
    for m in METHODS:
        c = pd.Series([s for g in df[(df.method == m) & (df.fold != "final")].symbols for s in g.split(";")]).value_counts()
        say(f"  {m:12s} " + ", ".join(f"{k}×{v}" for k, v in c.head(8).items()))

    # Son panel modelleri ve Pareto arşivleri
    (OUT / "final_panels.json").write_text(json.dumps(final["models"], indent=2), encoding="utf-8")
    for arm in ["MOGAFS_spec", "MOGAFS"]:
        arc = pd.DataFrame(final["extra"][arm]["archive"])
        arc["idx"] = arc["idx"].map(lambda t: ";".join(map(str, t)))
        arc["genes"] = arc["genes"].map(";".join)
        arc.to_csv(OUT / f"archive_final_{arm}.tsv", sep="\t", index=False)
        info = final["extra"][arm]["info"]
        g = pd.DataFrame({"entrez": [d["genes"][i] for i in final["top"]],
                          "symbol": [d["symbols"][i] for i in final["top"]],
                          "shap_pe_vs_ctrl": info["imp_pe"], "shap_comp_vs_ctrl": info["imp_off"],
                          "selection_prob": info["prob"]})
        g.to_csv(OUT / f"guidance_final_{arm}.tsv", sep="\t", index=False)
        say(f"  {arm}: son çalışmada {info['n_evaluated']} farklı küme değerlendirildi; Pareto cephesi "
            f"{int(arc.pareto.sum())} çözüm")
    (OUT / "build_panel_report.txt").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
