"""
Adım 8 (SONRADAN / KEŞİFSEL) — XAI-MOGAFS-SI: merkezler arası değişmez panel, merkez-dışarıda-bırakma (LOSO).

ÖNEMLİ: Bu analiz, önceden kaydedilmiş son test (Val-2, birincil V7 AUC 0.418) AÇILDIKTAN SONRA, başarısızlığın
nedenine (Stanford'a özgü RN7SL izi) yanıt olarak tasarlanmıştır. Makalede "post-hoc / keşifsel" olarak etiketlenir;
birincil sonucun yerine geçmez.

Merkezler: S1 Stanford (GSE192902 Discovery + Validation-1), S2 Validation-2 kurumları, S3 Del Vecchio (GSE154377).
Ön-işleme (her merkezde ayrı, etiketsiz): Entrez düzeyi ham sayım -> ORTAK protein kodlayan gen evreni (yapısal küçük
RNA/psödogen yok) -> bu evren üzerinden CPM -> log2(CPM+1) -> merkez içi z. Eğitimde kadın başına tüm ≤20 hf örnekleri
(S3: 1. + 2. trimester); test birimi kadının en erken örneği (S3: 1. trimester).

Yöntemler (her LOSO katında aynı eğitim/test):
  SI-MOGAFS   merkez-tutarlı ön-filtre (her eğitim merkezinde aynı yönde; skor = merkezler arası en küçük |t|, K=200)
              + NSGA-II amaçları f1 (havuzlanmış iç CV AUC), f3, f4, f_si (en kötü merkez AUC'si)
              + merkez-değişmez SHAP yönlendirmesi; seçim: tol içinde en iyi en-kötü-merkez AUC'si
  MOGAFS      havuzlanmış Welch ön-filtre + f1, f3, f4 (merkez bilgisi yok)
  META-TTEST5 merkez-tutarlı sıralamanın ilk 5 geni + lojistik (ön-filtrenin tek başına etkisi)
  LASSO5, TTEST5  havuzlanmış
  Moufarrej   yayımlanmış 18 gen, log OR ağırlıkları (S1'de bağımsız değil: eğitim kohortu)
Kadın gruplu iç CV; ≤5 özellik.
"""
import gzip
import io
import sys
import tarfile
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geo_utils import RAW, RESULTS, ROOT, pd, read_sample_annotation  # noqa: E402
from mogafs import choose, fit_lr, run_mogafs  # noqa: E402
from panel_utils import boot_auc, lasso_k, welch_t  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

OUT = RESULTS / "site_invariant"
OUT.mkdir(parents=True, exist_ok=True)
GENE_INFO = ROOT / "data" / "reference" / "Homo_sapiens.gene_info.gz"
STRUCT = r"^(RN7SL|RN7SK|RNA5-8S|RNA5S|RNU\d|RNY|RNA18S|RNA28S|RNA45S|SNOR|SCARNA|MT-)"
CFG = dict(K=200, max_k=5, pop_size=100, n_gen=60, tol=0.03, seed=2026, min_expr=1.0, min_frac=0.5)
METHODS = ["SI-MOGAFS", "MOGAFS", "META-TTEST5", "LASSO5", "TTEST5"]


def maps():
    gi = pd.read_csv(GENE_INFO, sep="\t", usecols=["GeneID", "Symbol", "dbXrefs", "type_of_gene"], dtype=str)
    ens = gi.dbXrefs.str.extractall(r"Ensembl:(ENSG\d+)")[0].reset_index()
    ens["GeneID"] = gi.GeneID.values[ens.level_0]
    e2i = ens.drop_duplicates(0, keep=False).set_index(0).GeneID.to_dict()
    sc = gi.Symbol.value_counts()
    s2i = {s: i for i, s in zip(gi.GeneID, gi.Symbol) if sc[s] == 1}
    pc = set(gi.GeneID[(gi.type_of_gene == "protein-coding") & ~gi.Symbol.str.contains(STRUCT, regex=True)])
    return e2i, s2i, pc, gi.set_index("GeneID").Symbol.to_dict()


def counts_192902(path, samples, e2i, s2i):
    raw = pd.read_csv(path)
    ent = [e2i.get(e) or s2i.get(s) for e, s in zip(raw.gene_num, raw.gene_name)]
    c = raw[list(samples)].assign(entrez=ent).dropna(subset=["entrez"]).groupby("entrez").sum()
    return c


def counts_delvecchio(gsms, e2i):
    cols = {}
    with tarfile.open(ROOT / "data" / "sealed" / "GSE154377_RAW.tar") as tf:
        for mem in tf.getmembers():
            g = mem.name.split("_")[0]
            if g in gsms:
                df = pd.read_csv(io.StringIO(gzip.decompress(tf.extractfile(mem).read()).decode()), sep="\t",
                                 header=None, index_col=0)
                cols[g] = df[~df.index.str.startswith("N_")].iloc[:, 0]
    c = pd.DataFrame(cols)
    ent = [e2i.get(e) for e in c.index.str.replace(r"\.\d+$", "", regex=True)]
    return c.assign(entrez=ent).dropna(subset=["entrez"]).groupby("entrez").sum()


def build():
    e2i, s2i, pc, sym = maps()
    meta = pd.read_csv(RESULTS / "inventory" / "GSE192902_samples.tsv", sep="\t", dtype=str)
    women = pd.read_csv(ROOT / "data" / "processed" / "selection" / "women_all.tsv", sep="\t", dtype=str)
    earliest = set(women.loc[women.cohort == "GSE192902", "sample_id"])
    sites = []
    for name, cohorts, files in [
        ("S1_Stanford", ["Discovery", "Validation 1"], [RAW / "GSE192902" / "GSE192902_counts_Discovery_postQC.csv.gz",
                                                        RAW / "GSE192902" / "GSE192902_counts_Validation1_postQC.csv.gz"]),
        ("S2_Val2", ["Validation 2"], [ROOT / "data" / "sealed" / "GSE192902_counts_Validation2_postQC.csv.gz"])]:
        parts, info = [], []
        for coh, f in zip(cohorts, files):
            m = meta[(meta.cohort == coh) & meta.tp_label.isin(["≤12", "13–20"])]
            parts.append(counts_192902(f, m.sample_name, e2i, s2i))
            info.append(pd.DataFrame({"sample": m.sample_name.values, "woman": name + ":" + m.woman_id.values,
                                      "group": np.where(m.disease == "control", "Control", "PE"),
                                      "gb": (m.tp_label == "13–20").astype(int).values,
                                      "earliest": m.sample_name.isin(earliest).values}))
        sites.append((name, pd.concat(parts, axis=1).fillna(0), pd.concat(info, ignore_index=True)))
    a = read_sample_annotation(RAW / "GSE154377" / "GSE154377_series_matrix.txt.gz")
    t = a.title.str.extract(r"^(?P<grp>.+?) Patient (?P<pid>\d+), (?P<tp>.+?) \(RNA-seq\)")
    dv = pd.DataFrame({"gsm": a.index}).join(t.reset_index(drop=True))
    dv = dv[dv.tp.isin(["1st Trimester", "2nd Trimester"])]
    dv["group"] = dv.grp.map({"Normal Pregnancy": "Control", "Preeclampsia/gestational hypertension": "PE",
                              "Gestational diabetes": "GDM", "Chronic  hypertension": "CHTN"})
    dv = dv[dv.group.isin(["Control", "PE", "GDM"])]
    c3 = counts_delvecchio(set(dv.gsm), e2i)
    sites.append(("S3_DelVecchio", c3[dv.gsm.tolist()],
                  pd.DataFrame({"sample": dv.gsm.values, "woman": "S3:" + dv.grp.values + dv.pid.values,
                                "group": dv.group.values, "gb": (dv.tp == "2nd Trimester").astype(int).values,
                                "earliest": (dv.tp == "1st Trimester").values})))
    # ortak protein kodlayan evren ve bu evren üzerinden CPM
    U = sorted(set.intersection(*[set(c.index) for _, c, _ in sites]) & pc, key=int)
    Xs, infos = [], []
    for k, (name, c, info) in enumerate(sites):
        cu = c.loc[U, info["sample"]]
        lg = np.log2(cu / cu.sum(0) * 1e6 + 1).T
        Xs.append(lg.values)
        infos.append(info.assign(site=k, site_name=name))
    info = pd.concat(infos, ignore_index=True)
    return np.vstack(Xs), info, U, sym


def youden(y, s):
    """Eğitim verisinde Youden J'yi en büyükleyen eşik."""
    from sklearn.metrics import roc_curve
    fpr, tpr, thr = roc_curve(y, s)
    return float(thr[np.argmax(tpr - fpr)])


def run_fold(task, d, cfg):
    held, = task
    t0 = time.time()
    X, info = d["X"], d["info"]
    trm = (info.site != held) & info.group.isin(["Control", "PE"])
    tr = np.flatnonzero(trm)
    y = (info.group == "PE").astype(int).values
    site, woman, gb = info.site.values, info.woman.values, info.gb.values
    keep = (X[tr] > cfg["min_expr"]).mean(0) >= cfg["min_frac"]
    # merkez içi z (etiketsiz; her merkezin tüm örnekleri)
    Z = np.empty_like(X[:, keep])
    for s_ in np.unique(site):
        m = site == s_
        mu, sd = X[m][:, keep].mean(0), X[m][:, keep].std(0, ddof=1)
        Z[m] = (X[m][:, keep] - mu) / np.where(sd == 0, 1, sd)
    genes = np.asarray(d["U"])[keep]
    Ztr, ytr = Z[tr], y[tr]
    # ön-filtreler
    t_pool = welch_t(Ztr, ytr)
    top_pool = np.argsort(-np.abs(t_pool))[: cfg["K"]]
    ts = np.array([welch_t(Ztr[site[tr] == s_], ytr[site[tr] == s_]) for s_ in np.unique(site[tr])])
    consistent = np.all(np.sign(ts) == np.sign(ts[0]), axis=0)
    si_score = np.where(consistent, np.min(np.abs(ts), axis=0), 0)
    top_si = np.argsort(-si_score)[: cfg["K"]]
    sel, extra = {}, {}
    for meth, top, sites_arg in [("SI-MOGAFS", top_si, site[tr]), ("MOGAFS", top_pool, None)]:
        arc, inf = run_mogafs(Ztr[:, top], ytr, gb[tr], use_spec=False, max_k=cfg["max_k"], pop_size=cfg["pop_size"],
                              n_gen=cfg["n_gen"], seed=cfg["seed"], groups=woman[tr], sites=sites_arg)
        best = choose(arc, False, cfg["tol"], use_si=sites_arg is not None)
        sel[meth] = top[list(best["idx"])]
        extra[meth] = {"inner_auc": best["auc_cv"], "worst_site_inner": best.get("worst_site_auc")}
    sel["META-TTEST5"] = top_si[: cfg["max_k"]]
    sel["LASSO5"] = lasso_k(Ztr, ytr, cfg["max_k"])
    sel["TTEST5"] = top_pool[: cfg["max_k"]]
    te_all = np.flatnonzero(info.site == held)
    rows, scores, thresholds = [], {}, {}
    for meth in METHODS:
        S = np.asarray(sel[meth])
        lr = fit_lr(Ztr[:, S], ytr)
        s = lr.decision_function(Z[te_all][:, S])
        scores[meth] = s
        thresholds[meth] = youden(ytr, lr.decision_function(Ztr[:, S]))  # eşik YALNIZ eğitim verisinden
        rows.append({"held_out": d["names"][held], "method": meth, "n_features": len(S),
                     "symbols": ";".join(d["sym"].get(g, g) for g in genes[S]), **extra.get(meth, {})})
    # Moufarrej (yayımlanmış, sabit)
    j, w = [], []
    gpos = {g: i for i, g in enumerate(genes)}
    for g, lw in d["mou"]:
        if g in gpos:
            j.append(gpos[g])
            w.append(lw)
    scores["Moufarrej18"] = Z[te_all][:, j] @ np.array(w) / np.abs(w).sum()
    rows.append({"held_out": d["names"][held], "method": "Moufarrej18", "n_features": len(j), "symbols": "(yayımlanmış)"})
    trs = Z[tr][:, j] @ np.array(w) / np.abs(w).sum()
    thresholds["Moufarrej18"] = youden(ytr, trs)
    return {"held": held, "te": te_all.tolist(), "scores": scores, "thresholds": thresholds, "rows": rows,
            "seconds": time.time() - t0}


def main():
    lines = []

    def say(m=""):
        print(m, flush=True)
        lines.append(str(m))

    X, info, U, sym = build()
    names = sorted(info.site_name.unique(), key=lambda n: info.loc[info.site_name == n, "site"].iloc[0])
    say("XAI-MOGAFS-SI — merkez-dışarıda-bırakma (SONRADAN / KEŞİFSEL; Val-2 açıldıktan sonra tasarlandı)")
    say(f"Ortak protein kodlayan evren: {len(U)} gen")
    for k, n in enumerate(names):
        sub = info[info.site == k]
        say(f"  {n}: {len(sub)} örnek / {sub.woman.nunique()} kadın; grup (kadın): "
            f"{sub.drop_duplicates('woman').group.value_counts().to_dict()}")
    cat = pd.read_csv(RESULTS / "literature" / "published_signatures.tsv", sep="\t", dtype=str)
    e2i, _, _, _ = maps()
    mo = cat[cat.signature_id == "Moufarrej2022_cfRNA_LR18"]
    mou = [(e2i.get(r.feature_id_if_given), np.log(float(r.weight.replace("OR=", "")))) for _, r in mo.iterrows()]
    d = {"X": X, "info": info, "U": U, "sym": sym, "names": names, "mou": [m for m in mou if m[0]]}
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=3) as ex:
        res = list(ex.map(run_fold, [(k,) for k in range(3)], [d] * 3, [CFG] * 3))
    say(f"Süre: {time.time() - t0:.0f} sn\n")
    rows = []
    for r in res:
        sub = info.iloc[r["te"]]
        for meth, s in r["scores"].items():
            for lab, mask, case in [("kadın", sub.earliest.values, "PE"), ("örnek", np.ones(len(sub), bool), "PE"),
                                    ("kadın-GDM", sub.earliest.values, "GDM")]:
                m = mask & sub.group.isin(["Control", case]).values
                yy = (sub.group.values[m] == case).astype(int)
                if yy.sum() < 3:
                    continue
                lo, hi = boot_auc(yy, s[m], n=1000)
                pred = (s[m] >= r["thresholds"][meth]).astype(int)
                sens, spec = pred[yy == 1].mean(), 1 - pred[yy == 0].mean()
                rows.append({"held_out": names[r["held"]], "method": meth, "level": lab, "n_case": int(yy.sum()),
                             "n_ctrl": int((1 - yy).sum()), "auc": roc_auc_score(yy, s[m]), "ci_lo": lo, "ci_hi": hi,
                             "sens": sens, "spec": spec, "acc": (pred == yy).mean(), "bacc": (sens + spec) / 2})
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "loso_results.tsv", sep="\t", index=False)
    # ROC eğrileri ve ek metrikler için test skorları (eşikler yalnız eğitimden)
    sc_rows = []
    for r in res:
        sub = info.iloc[r["te"]]
        for meth, s in r["scores"].items():
            sc_rows += [{"held_out": names[r["held"]], "method": meth, "sample": a, "woman": b, "group": c,
                         "earliest": e, "score": v, "threshold": r["thresholds"][meth]}
                        for a, b, c, e, v in zip(sub["sample"], sub.woman, sub.group, sub.earliest, s)]
    pd.DataFrame(sc_rows).to_csv(OUT / "loso_scores.tsv", sep="\t", index=False)
    P = pd.DataFrame([x for r in res for x in r["rows"]])
    P.to_csv(OUT / "loso_panels.tsv", sep="\t", index=False)
    for n in names:
        say(f"=== Dışarıda bırakılan merkez: {n}")
        for _, r in P[P.held_out == n].iterrows():
            say(f"  panel {r.method:12s} ({r.n_features}) {r.symbols}")
        for lab in ["kadın", "örnek", "kadın-GDM"]:
            sub = R[(R.held_out == n) & (R.level == lab)]
            if len(sub):
                say(f"  [{lab}] " + " | ".join(f"{x.method} {x.auc:.2f} [{x.ci_lo:.2f}–{x.ci_hi:.2f}]" for _, x in sub.iterrows())
                    + f"  (n={sub.n_case.iloc[0]}/{sub.n_ctrl.iloc[0]})")
        say()
    (OUT / "report.txt").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
