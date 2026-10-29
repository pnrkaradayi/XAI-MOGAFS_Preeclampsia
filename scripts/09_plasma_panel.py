"""
Adım 4 (plazma kolu) — Plazma cfRNA'da XAI-MOGAFS paneli: geliştir, iç doğrula, KİLİTLE.

ÖN-KAYITLI PROTOKOL (Val-2 açılmadan önce sabitlendi, 2026-10-22):
  Veri: GSE192902 ≤20 hf, kadın başına en erken örnek (Adım 2). log2(CPM+1).
  Yöntemler: MOGAFS (3 amaç: f1 = 1−AUC, f3 = |S|/5, f4 = SD(AUC) ≤12 / 13–20 hf dilimleri arasında;
             plazmada başka sendrom olmadığı için f2 yok), LASSO5, TTEST5.
  Katlama içi ön-işleme: düşük ifade filtresi (eğitim kadınlarının ≥%50'sinde log2CPM > 1),
             z-ölçekleme (eğitim istatistikleri), Welch t ön-filtresi (K=200).
  A. Yöntem performansı: Discovery içinde iç içe CV, dış 5 kat x 3 tekrar.
  B. İç doğrulama: Discovery'nin tamamında kurulan model -> Validation-1 (kohort içi z-ölçekleme).
  C. KİLİT: Discovery + Validation-1 (her kohort kendi içinde z-ölçeklenip birleştirilir) üzerinde aynı
     prosedür -> locked_models.json + SHA-256.
  Gen evreni (--universe): all = tüm 7.025 gen; protein_coding = yalnız protein kodlayan, yapısal küçük
             RNA / psödogen hariç. Kronoloji: results/plasma/LOCK_REGISTRY.md.
  Varyant seçenekleri (results/plasma/VARIANTS_LOG.md):
    --composition  ham sayımlardan kütüphane payları (RN7SL, RN7SK, RNA5, RNY, SNOR, RNU, MT) aday özellik
    --stability    MOGAFS'a f5 (bootstrap seçim kararlılığı) amacı + kararlılık öncelikli seçim
    --ga-covar     hafta dilimi (13–20 = 1) tüm lojistik fitlere eşlik eden değişken (panel boyutuna sayılmaz)
    --tag          çıktı klasörü results/plasma/variants/<tag>/
    --tol          Pareto cephesinden seçim toleransı (varsayılan 0.03; 0 = en yüksek iç-CV AUC)
  Son test (Adım 7, ayrı betik): Validation-2 kendi içinde z-ölçeklenir, kilitli modeller DEĞİŞTİRİLMEDEN
  uygulanır; bir kez. Kilitli bir gen Val-2'de yoksa z = 0 (kohort ortalaması) alınır.
"""
import argparse
import hashlib
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geo_utils import RAW, RESULTS, ROOT, pd  # noqa: E402
from mogafs import bootstrap_stability, choose, fit_lr, run_mogafs  # noqa: E402
from panel_utils import boot_auc, lasso_k, welch_t  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402
from sklearn.model_selection import RepeatedStratifiedKFold  # noqa: E402

CFG = dict(K=200, max_k=5, pop_size=100, n_gen=60, tol=0.03, outer_folds=5, outer_repeats=3, seed=2026,
           min_expr=1.0, min_frac=0.5)
STRUCT = r"^(RN7SL|RN7SK|RNA5-8S|RNA5S|RNU\d|RNY|RNA18S|RNA28S|RNA45S|SNOR|SCARNA|MT-)"
COMPOSITION = {"COMP:RN7SL": r"^RN7SL", "COMP:RN7SK": r"^RN7SK", "COMP:RNA5": r"^RNA5",
               "COMP:RNY": r"^RNY", "COMP:SNOR": r"^(SNOR|SCARNA)", "COMP:RNU": r"^RNU", "COMP:MT": r"^MT-"}
COUNT_FILES = {"cfRNA_discovery": RAW / "GSE192902" / "GSE192902_counts_Discovery_postQC.csv.gz",
               "cfRNA_val1": RAW / "GSE192902" / "GSE192902_counts_Validation1_postQC.csv.gz"}
METHODS = ["MOGAFS", "LASSO5", "TTEST5"]


def composition_features(counts):
    """counts: gen_adı x örnek ham sayım -> örnek x bileşim özelliği (log2 kütüphane yüzdesi)."""
    tot = counts.sum(0)
    out = {k: np.log2(counts[counts.index.str.contains(rx, regex=True)].sum(0) / tot * 100 + 0.01)
           for k, rx in COMPOSITION.items()}
    return pd.DataFrame(out)


def load_set(set_, with_composition):
    w = pd.read_csv(ROOT / "data" / "processed" / "selection" / "women_all.tsv", sep="\t", dtype=str)
    E = pd.read_csv(ROOT / "data" / "processed" / "expr" / f"GSE192902_{set_}.tsv.gz", sep="\t",
                    dtype={"entrez": str}).set_index("entrez").T
    ww = w[(w.cohort == "GSE192902") & (w.set == set_)].set_index("woman_id").loc[E.index]
    comp = None
    if with_composition:
        counts = pd.read_csv(COUNT_FILES[set_]).set_index("gene_name")[ww.sample_id.tolist()]
        comp = composition_features(counts)
        comp.index = ww.index
    return E, comp, (ww.group == "PE").astype(int).values, (ww.ga_bin == "13–20").astype(int).values


def zscore(M, mu=None, sd=None):
    mu = M.mean(0) if mu is None else mu
    sd = M.std(0, ddof=1) if sd is None else sd
    sd = np.where(sd == 0, 1, sd)
    return (M - mu) / sd, mu, sd


def run_task(task, d, cfg):
    """task = (name, train_idx, test_idx or None, external) ; d = birleşik veri sözlüğü."""
    name, tr, te, external = task
    t0 = time.time()
    X, y, gb, coh = d["X"], d["y"], d["gb"], d["cohort"]
    # düşük ifade filtresi (eğitim kadınlarında); bileşim özellikleri her zaman tutulur
    keep = ((X[tr] > cfg["min_expr"]).mean(0) >= cfg["min_frac"]) | d["is_comp"]
    Xk = X[:, keep]
    genes = np.asarray(d["genes"])[keep]
    # z-ölçekleme: eğitim kümesinde her kohort kendi içinde (birleşik kilit eğitimi için)
    Ztr = np.empty((len(tr), Xk.shape[1]))
    for c in np.unique(coh[tr]):
        m = coh[tr] == c
        Ztr[m], _, _ = zscore(Xk[tr][m])
    if te is not None:
        if external:                                   # başka kohort: kendi içinde, etiketsiz
            Zte, _, _ = zscore(Xk[te])
        else:                                          # aynı kohort dış katmanı: eğitim istatistikleri
            _, mu, sd = zscore(Xk[tr])
            Zte, _, _ = zscore(Xk[te], mu, sd)
    t = welch_t(Ztr, y[tr])
    top = np.argsort(-np.abs(t))[: cfg["K"]]
    cov_tr = gb[tr].astype(float) if cfg["ga_covar"] else None
    stab = bootstrap_stability(Ztr[:, top], y[tr], seed=cfg["seed"]) if cfg["stability"] else None
    archive, info = run_mogafs(Ztr[:, top], y[tr], gb[tr], use_spec=False, max_k=cfg["max_k"],
                               pop_size=cfg["pop_size"], n_gen=cfg["n_gen"], seed=cfg["seed"],
                               covar=cov_tr, stab=stab)
    best = choose(archive, False, cfg["tol"], use_stab=cfg["stability"])
    sel = {"MOGAFS": top[list(best["idx"])], "LASSO5": lasso_k(Ztr, y[tr], cfg["max_k"]),
           "TTEST5": top[: cfg["max_k"]]}
    out = {"name": name, "te": None if te is None else te.tolist(), "rows": [], "scores": {}, "models": {},
           "n_genes_after_filter": int(keep.sum()), "seconds": None}

    def design(Z, rows, S):
        return np.column_stack([Z[:, S], gb[rows]]) if cfg["ga_covar"] else Z[:, S]

    for m in METHODS:
        S = np.asarray(sel[m])
        lr = fit_lr(design(Ztr, tr, S), y[tr])
        rec = {"task": name, "method": m, "n_features": len(S), "genes": ";".join(genes[S]),
               "symbols": ";".join(d["sym"].get(g, g) for g in genes[S])}
        if m == "MOGAFS":
            rec["inner_auc_cv"] = best["auc_cv"]
        if te is not None:
            s = lr.decision_function(design(Zte, te, S))
            out["scores"][m] = s
            rec["test_auc"] = roc_auc_score(y[te], s)
        out["rows"].append(rec)
        coef = lr.coef_[0]
        out["models"][m] = {"genes": genes[S].tolist(), "symbols": [d["sym"].get(g, g) for g in genes[S]],
                            "coef": coef[: len(S)].tolist(), "intercept": float(lr.intercept_[0])}
        if cfg["ga_covar"]:
            out["models"][m]["covariate"] = {"name": "ga_bin_13_20", "coef": float(coef[-1])}
    if name == "LOCK":
        arc = pd.DataFrame(archive)
        arc["genes"] = arc["idx"].map(lambda ix: ";".join(genes[top[list(ix)]]))
        arc["idx"] = arc["idx"].map(lambda ix: ";".join(map(str, ix)))
        out["archive"] = arc
        out["n_evaluated"] = info["n_evaluated"]
    out["seconds"] = time.time() - t0
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--universe", choices=["all", "protein_coding"], default="protein_coding")
    ap.add_argument("--composition", action="store_true")
    ap.add_argument("--stability", action="store_true")
    ap.add_argument("--ga-covar", action="store_true")
    ap.add_argument("--tag", default=None)
    ap.add_argument("--tol", type=float, default=CFG["tol"], help="Pareto seçim toleransı (0 = en yüksek iç-CV AUC)")
    args = ap.parse_args()
    universe = args.universe
    cfg = dict(CFG, composition=args.composition, stability=args.stability, ga_covar=args.ga_covar, tol=args.tol)
    if args.tag:
        OUT = RESULTS / "plasma" / "variants" / args.tag
    else:
        OUT = RESULTS / "plasma" / {"all": "all_genes", "protein_coding": "protein_coding"}[universe]
    OUT.mkdir(parents=True, exist_ok=True)
    lines = []

    def say(msg=""):
        print(msg, flush=True)
        lines.append(str(msg))

    Ed, Cd, yd, gd = load_set("cfRNA_discovery", args.composition)
    Ev, Cv, yv, gv = load_set("cfRNA_val1", args.composition)
    assert list(Ed.columns) == list(Ev.columns)
    u = pd.read_csv(ROOT / "data" / "processed" / "gene_universe.tsv", sep="\t", dtype={"entrez": str}).set_index("entrez")
    if universe == "protein_coding":
        sym = u.symbol.reindex(Ed.columns).fillna("")
        keep = ((u.type_of_gene.reindex(Ed.columns) == "protein-coding") & ~sym.str.contains(STRUCT, regex=True)).values
        Ed, Ev = Ed.loc[:, keep], Ev.loc[:, keep]
    if args.composition:
        Ed, Ev = Ed.join(Cd), Ev.join(Cv)
    genes = list(Ed.columns)
    d = {"X": np.vstack([Ed.values, Ev.values]), "y": np.r_[yd, yv], "gb": np.r_[gd, gv],
         "cohort": np.r_[np.zeros(len(yd), int), np.ones(len(yv), int)], "genes": genes,
         "is_comp": np.array([g.startswith("COMP:") for g in genes]), "sym": u.symbol.dropna().to_dict()}
    nd, nv = len(yd), len(yv)
    disc, val = np.arange(nd), np.arange(nd, nd + nv)
    say(f"Discovery: {nd} kadın ({yd.sum()} PE), Validation-1: {nv} ({yv.sum()} PE); {Ed.shape[1]} özellik "
        f"({int(d['is_comp'].sum())} bileşim)")
    say(f"Gen evreni: {universe}; varyant: bileşim={args.composition}, kararlılık={args.stability}, "
        f"hafta eşlik={args.ga_covar}, tol={args.tol}; Ayarlar: {CFG}")
    say(f"Hafta dilimi 13–20 oranı: Discovery PE {gd[yd == 1].mean():.2f} / kontrol {gd[yd == 0].mean():.2f}; "
        f"Val-1 PE {gv[yv == 1].mean():.2f} / kontrol {gv[yv == 0].mean():.2f}")

    rskf = RepeatedStratifiedKFold(n_splits=CFG["outer_folds"], n_repeats=CFG["outer_repeats"], random_state=CFG["seed"])
    tasks = [(f"cv{i}", disc[tr], disc[te], False) for i, (tr, te) in enumerate(rskf.split(Ed.values, yd))]
    tasks += [("DISC_FULL", disc, val, True), ("LOCK", np.r_[disc, val], None, False)]
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=min(14, len(tasks))) as ex:
        res = list(ex.map(run_task, tasks, [d] * len(tasks), [cfg] * len(tasks)))
    say(f"Süre: {time.time() - t0:.0f} sn")
    df = pd.DataFrame([r for x in res for r in x["rows"]])
    df.to_csv(OUT / "plasma_tasks.tsv", sep="\t", index=False)

    say("\nA. DISCOVERY İÇ İÇE CV (dış 5 kat x 3 tekrar) — PE vs kontrol")
    cv = [x for x in res if x["name"].startswith("cv")]
    summary = []
    for m in METHODS:
        per_rep, ci = [], []
        for r in range(CFG["outer_repeats"]):
            s = np.zeros(nd)
            for x in cv[r * CFG["outer_folds"]:(r + 1) * CFG["outer_folds"]]:
                s[x["te"]] = x["scores"][m]
            per_rep.append(roc_auc_score(yd, s))
            if r == 0:
                ci = boot_auc(yd, s)
        fa = df[(df.method == m) & df.task.str.startswith("cv")]
        say(f"  {m:8s} havuzlanmış dış AUC {np.mean(per_rep):.3f} (tekrarlar {', '.join(f'{a:.3f}' for a in per_rep)}; "
            f"1. tekrar %95 GA {ci[0]:.3f}–{ci[1]:.3f}); ort. {fa.n_features.mean():.1f} özellik")
        c = pd.Series([g for s_ in fa.symbols for g in s_.split(";")]).value_counts()
        say(f"  {'':8s} en sık seçilenler: " + ", ".join(f"{k}×{v}" for k, v in c.head(8).items()))
        summary.append({"method": m, "disc_cv_auc": float(np.mean(per_rep))})

    say("\nB. İÇ DOĞRULAMA — Discovery'de kurulan model, Validation-1'de")
    full = [x for x in res if x["name"] == "DISC_FULL"][0]
    for i, m in enumerate(METHODS):
        s = full["scores"][m]
        lo, hi = boot_auc(yv, s)
        rec = df[(df.task == "DISC_FULL") & (df.method == m)].iloc[0]
        summary[i]["val1_auc"] = float(roc_auc_score(yv, s))
        say(f"  {m:8s} {rec.symbols:45s} Val-1 AUC {roc_auc_score(yv, s):.3f} [%95 GA {lo:.3f}–{hi:.3f}] (n={yv.sum()}/{(1 - yv).sum()})")

    lock = [x for x in res if x["name"] == "LOCK"][0]
    say(f"\nC. KİLİTLİ MODELLER (Discovery + Validation-1, {len(d['y'])} kadın; filtre sonrası {lock['n_genes_after_filter']} özellik)")
    for m in METHODS:
        md = lock["models"][m]
        cv_txt = f"  [hafta eşlik katsayısı {md['covariate']['coef']:+.2f}]" if "covariate" in md else ""
        say(f"  {m:8s} " + ", ".join(f"{s} ({c:+.2f})" for s, c in zip(md["symbols"], md["coef"])) + cv_txt)
    locked = {"created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
              "protocol": "scripts/09_plasma_panel.py header; Val-2: within-cohort z-score, apply unchanged, open once",
              "universe": universe, "variant": {k: cfg[k] for k in ["composition", "stability", "ga_covar", "tol"]},
              "tag": args.tag, "config": CFG, "models": lock["models"]}
    path = OUT / "locked_models.json"
    path.write_text(json.dumps(locked, indent=2, ensure_ascii=False), encoding="utf-8")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    (OUT / "locked_models.sha256").write_text(f"{digest}  locked_models.json\n", encoding="utf-8")
    say(f"  kilit dosyası SHA-256: {digest[:16]}…  ({lock['n_evaluated']} küme değerlendirildi; Pareto "
        f"{int(lock['archive'].pareto.sum())} çözüm)")
    lock["archive"].to_csv(OUT / "archive_LOCK_MOGAFS.tsv", sep="\t", index=False)
    pd.DataFrame(summary).to_csv(OUT / "summary.tsv", sep="\t", index=False)
    (OUT / "plasma_report.txt").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
