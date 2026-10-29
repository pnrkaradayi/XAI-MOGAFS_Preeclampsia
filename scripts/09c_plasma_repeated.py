"""
Adım 4 (plazma, V7/V8) — Tekrarlı-örnek eğitimi, kadın düzeyinde sızıntısız bölme.

Moufarrej 2022 eğitimde kadın başına birden çok erken örnek kullanmış ve örnek düzeyinde raporlamış. Burada:
  * Eğitim: eğitim kadınlarının TÜM ≤20 hf örnekleri (≤12 ve 13–20 dilimleri, '.1' tekrarları dahil).
  * Bölmeler: dış CV StratifiedGroupKFold (kadın), 5 kat x 3 tekrar; iç CV (MOGAFS) de kadın gruplu.
  * Test: hiç görülmemiş kadınlar; (a) kadın düzeyi = kadının en erken örneği (Adım 2 seçimiyle aynı),
          (b) örnek düzeyi = test kadınlarının tüm ≤20 hf örnekleri (Moufarrej'in raporlama birimi).
  * Ön-işleme kuralları 09_plasma_panel.py ile aynı (düşük ifade filtresi, kohort içi z, Welch K=200, tol 0.03).
  * Val-1 karşılaştırması: Discovery'de kurulan modeller; aynı Val-1 örneklerinde V0 paneli (Discovery'de yeniden
    fit) ve Moufarrej 18 gen (yayımlanmış OR'lerle) — örnek ve kadın düzeyinde.
  * Kilit: Discovery + Val-1 tüm erken örnekleri; eğitim (apparent) AUC'si "eğitim" etiketiyle raporlanır.
  Validation-2 açılmaz.

Kullanım: python scripts/09c_plasma_repeated.py --tag V7 --max-k 5
          python scripts/09c_plasma_repeated.py --tag V8 --max-k 18
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
from mogafs import choose, fit_lr, run_mogafs  # noqa: E402
from panel_utils import boot_auc, cfrna_gene_level, composition_features, lasso_k, welch_t  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402
from sklearn.model_selection import StratifiedGroupKFold  # noqa: E402

CFG = dict(K=200, pop_size=100, n_gen=60, tol=0.03, outer_folds=5, outer_repeats=3, seed=2026,
           min_expr=1.0, min_frac=0.5)
FILES = {"Discovery": "GSE192902_counts_Discovery_postQC.csv.gz", "Validation 1": "GSE192902_counts_Validation1_postQC.csv.gz"}
GENE_INFO = ROOT / "data" / "reference" / "Homo_sapiens.gene_info.gz"
METHODS = ["MOGAFS", "LASSO", "TTEST"]


def build():
    meta = pd.read_csv(RESULTS / "inventory" / "GSE192902_samples.tsv", sep="\t", dtype=str)
    women = pd.read_csv(ROOT / "data" / "processed" / "selection" / "women_all.tsv", sep="\t", dtype=str)
    earliest = set(women.loc[women.cohort == "GSE192902", "sample_id"])
    blocks = []
    for ci, (coh, f) in enumerate(FILES.items()):
        m = meta[(meta.cohort == coh) & meta.tp_label.isin(["≤12", "13–20"])].set_index("sample_name")
        raw = pd.read_csv(RAW / "GSE192902" / f)
        g = cfrna_gene_level(raw, GENE_INFO)[m.index].T
        comp = composition_features(raw.set_index("gene_name").drop(columns=["gene_num"])[m.index])
        X = g.join(comp)
        blocks.append(pd.DataFrame({"sample": m.index, "woman": m.woman_id.values,
                                    "y": (m.disease != "control").astype(int).values,
                                    "gb": (m.tp_label == "13–20").astype(int).values, "cohort": ci,
                                    "earliest": m.index.isin(earliest)}).assign(_X=list(X.values)))
        cols = list(X.columns)
    info = pd.concat(blocks, ignore_index=True)
    X = np.vstack(info.pop("_X").values)
    return info, X, cols


def zscore(M, mu=None, sd=None):
    mu = M.mean(0) if mu is None else mu
    sd = M.std(0, ddof=1) if sd is None else sd
    sd = np.where(sd == 0, 1, sd)
    return (M - mu) / sd, mu, sd


def run_task(task, d, cfg):
    name, tr, te, external = task
    t0 = time.time()
    X, y, gb, coh, woman = d["X"], d["y"], d["gb"], d["cohort"], d["woman"]
    keep = ((X[tr] > cfg["min_expr"]).mean(0) >= cfg["min_frac"]) | d["is_comp"]
    Xk = X[:, keep]
    feats = np.asarray(d["cols"])[keep]
    Ztr = np.empty((len(tr), Xk.shape[1]))
    for c in np.unique(coh[tr]):
        m = coh[tr] == c
        Ztr[m], _, _ = zscore(Xk[tr][m])
    if te is not None:
        if external:
            Zte, _, _ = zscore(Xk[te])
        else:
            _, mu, sd = zscore(Xk[tr])
            Zte, _, _ = zscore(Xk[te], mu, sd)
    t = welch_t(Ztr, y[tr])
    top = np.argsort(-np.abs(t))[: cfg["K"]]
    archive, info = run_mogafs(Ztr[:, top], y[tr], gb[tr], use_spec=False, max_k=cfg["max_k"],
                               pop_size=cfg["pop_size"], n_gen=cfg["n_gen"], seed=cfg["seed"], groups=woman[tr])
    best = choose(archive, False, cfg["tol"])
    sel = {"MOGAFS": top[list(best["idx"])], "LASSO": lasso_k(Ztr, y[tr], cfg["max_k"]), "TTEST": top[: cfg["max_k"]]}
    out = {"name": name, "te": None if te is None else te.tolist(), "scores": {}, "models": {}, "rows": [],
           "seconds": None}
    for m in METHODS:
        S = np.asarray(sel[m])
        lr = fit_lr(Ztr[:, S], y[tr])
        rec = {"task": name, "method": m, "n_features": len(S), "symbols": ";".join(d["sym"].get(f, f) for f in feats[S]),
               "train_apparent_auc": roc_auc_score(y[tr], lr.decision_function(Ztr[:, S]))}
        if te is not None:
            out["scores"][m] = lr.decision_function(Zte[:, S])
        out["rows"].append(rec)
        out["models"][m] = {"features": feats[S].tolist(), "symbols": [d["sym"].get(f, f) for f in feats[S]],
                            "coef": lr.coef_[0].tolist(), "intercept": float(lr.intercept_[0])}
    if name == "LOCK":
        out["n_evaluated"] = info["n_evaluated"]
        out["pareto"] = int(sum(r["pareto"] for r in archive))
    out["seconds"] = time.time() - t0
    return out


def auc_ci(y, s):
    lo, hi = boot_auc(y, s, n=1000)
    return roc_auc_score(y, s), lo, hi


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--max-k", type=int, default=5)
    args = ap.parse_args()
    cfg = dict(CFG, max_k=args.max_k)
    OUT = RESULTS / "plasma" / "variants" / args.tag
    OUT.mkdir(parents=True, exist_ok=True)
    lines = []

    def say(msg=""):
        print(msg, flush=True)
        lines.append(str(msg))

    info, X, cols = build()
    u = pd.read_csv(ROOT / "data" / "processed" / "gene_universe.tsv", sep="\t", dtype={"entrez": str}).set_index("entrez")
    d = {"X": X, "y": info.y.values, "gb": info.gb.values, "cohort": info.cohort.values, "woman": info.woman.values,
         "cols": cols, "is_comp": np.array([c.startswith("COMP:") for c in cols]), "sym": u.symbol.dropna().to_dict()}
    disc = np.flatnonzero(info.cohort == 0)
    val = np.flatnonzero(info.cohort == 1)
    for nm, ix in [("Discovery", disc), ("Validation-1", val)]:
        sub = info.iloc[ix]
        say(f"{nm}: {len(ix)} örnek / {sub.woman.nunique()} kadın; PE örnek {sub.y.sum()} / kadın "
            f"{sub[sub.y == 1].woman.nunique()}; kadın başına örnek ort. {len(ix) / sub.woman.nunique():.2f}")
    say(f"Varyant {args.tag}: tüm genler + bileşim, ≤{args.max_k} özellik, tekrarlı-örnek eğitimi, kadın gruplu CV")

    tasks = []
    yw = info.iloc[disc].groupby("woman").y.first()
    for r in range(cfg["outer_repeats"]):
        sgk = StratifiedGroupKFold(n_splits=cfg["outer_folds"], shuffle=True, random_state=cfg["seed"] + r)
        for k, (a, b) in enumerate(sgk.split(disc, info.y.values[disc], info.woman.values[disc])):
            tasks.append((f"cv{r}_{k}", disc[a], disc[b], False))
    tasks += [("DISC_FULL", disc, val, True), ("LOCK", np.r_[disc, val], None, False)]
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=min(14, len(tasks))) as ex:
        res = list(ex.map(run_task, tasks, [d] * len(tasks), [cfg] * len(tasks)))
    say(f"Süre: {time.time() - t0:.0f} sn")
    rows = pd.DataFrame([r for x in res for r in x["rows"]])
    rows.to_csv(OUT / "tasks.tsv", sep="\t", index=False)

    say("\nA. DISCOVERY İÇ İÇE CV (kadın gruplu, 5 kat x 3 tekrar)")
    summary = []
    for m in METHODS:
        aw, asmp = [], []
        for r in range(cfg["outer_repeats"]):
            sc = np.full(len(info), np.nan)
            for x in res:
                if x["name"].startswith(f"cv{r}_"):
                    sc[x["te"]] = x["scores"][m]
            ew = disc[info.earliest.values[disc]]
            aw.append(roc_auc_score(info.y.values[ew], sc[ew]))
            asmp.append(roc_auc_score(info.y.values[disc], sc[disc]))
        tr_app = rows[(rows.method == m) & rows.task.str.startswith("cv")].train_apparent_auc.mean()
        say(f"  {m:7s} kadın düzeyi AUC {np.mean(aw):.3f} ({', '.join(f'{a:.3f}' for a in aw)}) | örnek düzeyi "
            f"{np.mean(asmp):.3f} | eğitim (apparent) {tr_app:.3f}")
        summary.append({"method": m, "disc_cv_auc_woman": np.mean(aw), "disc_cv_auc_sample": np.mean(asmp),
                        "train_apparent_auc": tr_app})

    say("\nB. VAL-1 (Discovery'de kurulan modeller; aynı örneklerde rakipler)")
    full = [x for x in res if x["name"] == "DISC_FULL"][0]
    vy = info.y.values[val]
    vearly = info.earliest.values[val]
    for i, m in enumerate(METHODS):
        s = full["scores"][m]
        a1, l1, h1 = auc_ci(vy[vearly], s[vearly])
        a2, l2, h2 = auc_ci(vy, s)
        summary[i].update({"val1_auc_woman": a1, "val1_auc_sample": a2})
        say(f"  {m:7s} {';'.join(full['models'][m]['symbols'])[:60]:60s} kadın {a1:.3f} [{l1:.2f}–{h1:.2f}] | "
            f"örnek {a2:.3f} [{l2:.2f}–{h2:.2f}]")
    # Rakipler aynı Val-1 örneklerinde: V0 paneli (Discovery'de yeniden fit) ve Moufarrej 18 gen
    sym2col = {d["sym"].get(c, c): j for j, c in enumerate(cols)}
    Zv, _, _ = zscore(X[val])
    Zd, _, _ = zscore(X[disc][info.earliest.values[disc]])
    v0 = [sym2col[s] for s in ["RN7SL230P", "HIP1R", "KIF13B"]]
    lr0 = fit_lr(Zd[:, v0], info.y.values[disc][info.earliest.values[disc]])
    s0 = lr0.decision_function(Zv[:, v0])
    cat = pd.read_csv(RESULTS / "literature" / "published_signatures.tsv", sep="\t", dtype=str)
    mo = cat[cat.signature_id == "Moufarrej2022_cfRNA_LR18"]
    gi = pd.read_csv(GENE_INFO, sep="\t", usecols=["GeneID", "dbXrefs"], dtype=str)
    ens = gi.dbXrefs.str.extractall(r"Ensembl:(ENSG\d+)")[0].reset_index()
    ens["GeneID"] = gi.GeneID.values[ens.level_0]
    e2i = ens.drop_duplicates(0, keep=False).set_index(0).GeneID.to_dict()
    colpos = {c: j for j, c in enumerate(cols)}
    w_, j_ = [], []
    for _, r in mo.iterrows():
        g = e2i.get(r.feature_id_if_given)
        if g in colpos:
            j_.append(colpos[g])
            w_.append(np.log(float(r.weight.replace("OR=", ""))))
    smo = Zv[:, j_] @ np.array(w_) / np.abs(w_).sum()
    for lab, s in [("V0 paneli (RN7SL230P, HIP1R, KIF13B)", s0), (f"Moufarrej 18 gen ({len(j_)}/18 mevcut, log OR)", smo)]:
        a1, l1, h1 = auc_ci(vy[vearly], s[vearly])
        a2, l2, h2 = auc_ci(vy, s)
        say(f"  {lab:68s} kadın {a1:.3f} [{l1:.2f}–{h1:.2f}] | örnek {a2:.3f} [{l2:.2f}–{h2:.2f}]")
    say(f"  (Val-1: {vearly.sum()} kadın ({vy[vearly].sum()} PE) / {len(val)} örnek ({vy.sum()} PE örneği); "
        f"Moufarrej yayımlanmış Val-1 AUC 0.71, ≤16 hf örnekleri)")

    lock = [x for x in res if x["name"] == "LOCK"][0]
    say(f"\nC. KİLİT (Discovery + Val-1, {len(info)} örnek / {info.woman.nunique()} kadın)")
    lr_app = rows[(rows.task == "LOCK")].set_index("method").train_apparent_auc
    for m in METHODS:
        md = lock["models"][m]
        say(f"  {m:7s} eğitim (apparent) AUC {lr_app[m]:.3f} | " + ", ".join(f"{s} ({c:+.2f})" for s, c in zip(md["symbols"], md["coef"])))
    locked = {"created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
              "protocol": "scripts/09c_plasma_repeated.py header; Val-2: within-cohort z-score over all ≤20 wk samples, apply unchanged, open once",
              "variant": {"tag": args.tag, "max_k": args.max_k, "repeated_samples": True, "composition": True},
              "config": cfg, "models": lock["models"]}
    path = OUT / "locked_models.json"
    path.write_text(json.dumps(locked, indent=2, ensure_ascii=False), encoding="utf-8")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    (OUT / "locked_models.sha256").write_text(f"{digest}  locked_models.json\n", encoding="utf-8")
    say(f"  kilit SHA-256: {digest[:16]}… ({lock['n_evaluated']} küme; Pareto {lock['pareto']})")
    pd.DataFrame(summary).to_csv(OUT / "summary.tsv", sep="\t", index=False)
    (OUT / "report.txt").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
