"""
Adım 7 — SON TEST: mühürlü GSE192902 Validation-2 ve GSE154377 (Del Vecchio) BİR KEZ açılır.

Protokol: results/plasma/LOCK_REGISTRY.md (açılmadan önce yazıldı).
  1. Mühürlü dosyaların SHA-256'sı manifest ile doğrulanır.
  2. İfade verisi ETİKETSİZ işlenir (özellik matrisleri + kalite kontrol); etiketler bu aşamada yüklenmez.
  3. Etiketler birleştirilir; kayıtlı her model DEĞİŞTİRİLMEDEN uygulanır; AUC + tabakalı bootstrap %95 GA,
     hafta dilimine göre AUC, hafta-düzeltmeli OR, birincil model vs Moufarrej için DeLong testi.
  4. Açılış zamanı results/final_test/OPENED.txt'e yazılır; yeniden çalıştırma --rerun-reason ister.

--dry-run: kodu Validation-1 üzerinde (Val-2 yerine) çalıştırır, mühürlü hiçbir dosyaya dokunmaz. Val-1 kilit
           eğitiminde kullanıldığı için bu rakamlar ANLAMSIZDIR; yalnız kod denemesidir.
"""
import argparse
import gzip
import hashlib
import io
import json
import re
import sys
import tarfile
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geo_utils import RAW, RESULTS, ROOT, pd, read_sample_annotation  # noqa: E402
from panel_utils import (boot_auc, cfrna_gene_level, composition_features,  # noqa: E402
                         delong_test)
from sklearn.metrics import roc_auc_score  # noqa: E402
import statsmodels.api as sm  # noqa: E402

warnings.filterwarnings("ignore")
PL = RESULTS / "plasma"
GENE_INFO = ROOT / "data" / "reference" / "Homo_sapiens.gene_info.gz"
SEALED_V2 = ROOT / "data" / "sealed" / "GSE192902_counts_Validation2_postQC.csv.gz"
SEALED_DV = ROOT / "data" / "sealed" / "GSE154377_RAW.tar"
LIT_SIGS = ["Moufarrej2022_cfRNA_LR18", "Weiner2022_EOP_NAMPT_APOA1", "Weiner2026_EOP_NAMPT",
            "Elovitz2025_PE_vs_all_top2", "Tarca2021_LDA4_H19_FN1_TUBB6_FPR3", "Myhrer2024_RPPH1_leukocytes",
            "Lekva2021_PBMC_pilotarray_top20"]


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def registry():
    """Kayıtlı kilitler: (etiket, json yolu, z-ölçekleme birimi)."""
    out = [("V7", PL / "variants/V7/locked_models.json", "all_samples"),
           ("V8", PL / "variants/V8/locked_models.json", "all_samples"),
           ("V0", PL / "all_genes/locked_models.json", "women"),
           ("V0pc", PL / "protein_coding/locked_models.json", "women")]
    out += [(f"V{i}", PL / f"variants/V{i}/locked_models.json", "women") for i in range(1, 7)]
    return out


# ---------------------------------------------------------------- etiketsiz ön-işleme
def gene_maps():
    gi = pd.read_csv(GENE_INFO, sep="\t", usecols=["GeneID", "Symbol", "dbXrefs"], dtype=str)
    ens = gi.dbXrefs.str.extractall(r"Ensembl:(ENSG\d+)")[0].reset_index()
    ens["GeneID"] = gi.GeneID.values[ens.level_0]
    e2i = ens.drop_duplicates(0, keep=False).set_index(0).GeneID.to_dict()
    i2s = gi.set_index("GeneID").Symbol.to_dict()
    return e2i, i2s


def prep_gse192902(path, sample_names):
    raw = pd.read_csv(path)
    missing = set(sample_names) - set(raw.columns)
    assert not missing, f"sayım dosyasında olmayan örnekler: {sorted(missing)[:5]}"
    g = cfrna_gene_level(raw[["gene_name", "gene_num"] + list(sample_names)], GENE_INFO).T
    comp = composition_features(raw.set_index("gene_name")[list(sample_names)])
    rn7 = raw.set_index("gene_name")[list(sample_names)]
    frac = rn7[rn7.index.str.startswith("RN7SL")].sum(0) / rn7.sum(0)
    return g.join(comp), frac, {"genes_rows": len(raw), "samples": len(sample_names)}


def prep_delvecchio(gsms):
    e2i, i2s = gene_maps()
    cols = {}
    with tarfile.open(SEALED_DV) as tf:
        for mem in tf.getmembers():
            gsm = mem.name.split("_")[0]
            if gsm not in gsms:
                continue
            txt = gzip.decompress(tf.extractfile(mem).read()).decode()
            df = pd.read_csv(io.StringIO(txt), sep="\t", header=None, index_col=0)
            cols[gsm] = df[~df.index.str.startswith("N_")]
    assert set(cols) == set(gsms), f"tar'da eksik örnek: {set(gsms) - set(cols)}"
    count_cols = [c for c in next(iter(cols.values())).columns]  # STAR: 3 iplik sütunu; bu veri setinde tek sütun
    tot = {c: sum(cols[g][c].sum() for g in cols) for c in count_cols}
    strand = max(tot, key=tot.get)
    counts = pd.DataFrame({g: cols[g][strand] for g in gsms})
    ens = counts.index.str.replace(r"\.\d+$", "", regex=True)
    ent = pd.Series([e2i.get(e) for e in ens], index=counts.index)
    sym = ent.map(i2s).fillna("")
    comp = composition_features(counts.set_axis(sym.values))
    gl = counts.assign(entrez=ent.values).dropna(subset=["entrez"]).groupby("entrez").sum()
    g = np.log2(gl / gl.sum(0) * 1e6 + 1).T
    named = sym.values
    frac = counts[pd.Series(named, index=counts.index).str.startswith("RN7SL").values].sum(0) / counts.sum(0)
    qc = {"strand_column": int(strand), "assigned_by_column": tot, "genes_rows": len(counts),
          "mapped_entrez": int(ent.notna().sum()), "samples": len(gsms)}
    return g.join(comp), frac, qc


# ---------------------------------------------------------------- değerlendirme
def zscore(M):
    sd = M.std(ddof=1).replace(0, 1)
    return (M - M.mean()) / sd


def apply_model(md, Z, ga_bin):
    feats = md.get("features") or md["genes"]
    coef = np.array(md["coef"], float)
    X = np.column_stack([Z[f].values if f in Z.columns else np.zeros(len(Z)) for f in feats])
    s = md["intercept"] + X @ coef
    if "covariate" in md:
        s = s + md["covariate"]["coef"] * ga_bin
    return s, sum(f in Z.columns for f in feats), len(feats)


def score_signature(sig, Z, e2i):
    w = []
    for _, r in sig.iterrows():
        g = e2i.get(r.feature_id_if_given) if isinstance(r.feature_id_if_given, str) and r.feature_id_if_given.startswith("ENSG") else None
        w.append((g, r))
    return w


def metrics(y, s, ga_bin=None, n_boot=2000):
    y = np.asarray(y).astype(int)
    auc = roc_auc_score(y, s)
    lo, hi = boot_auc(y, s, n=n_boot)
    out = {"auc": auc, "ci_lo": lo, "ci_hi": hi, "n_case": int(y.sum()), "n_ctrl": int((1 - y).sum())}
    if ga_bin is not None:
        for b, lab in [(0, "≤12"), (1, "13–20")]:
            m = ga_bin == b
            if y[m].sum() >= 3 and (1 - y[m]).sum() >= 3:
                out[f"auc_{lab}"] = roc_auc_score(y[m], s[m])
        try:
            z = (s - s.mean()) / s.std(ddof=1)
            f = sm.Logit(y, sm.add_constant(np.column_stack([z, ga_bin]))).fit(disp=0)
            out.update({"or_ga_adj": float(np.exp(f.params[1])), "or_p": float(f.pvalues[1])})
        except Exception:
            pass
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--rerun-reason", default=None)
    args = ap.parse_args()
    OUT = RESULTS / ("final_test_dryrun" if args.dry_run else "final_test")
    OUT.mkdir(parents=True, exist_ok=True)
    opened = OUT / "OPENED.txt"
    if not args.dry_run and opened.exists() and not args.rerun_reason:
        sys.exit("Son test zaten açıldı. Yalnız teknik hata düzeltmesi için --rerun-reason ile yeniden çalıştırılabilir.")
    lines = []

    def say(m=""):
        print(m, flush=True)
        lines.append(str(m))

    # 1. mühür doğrulaması
    man = pd.read_csv(ROOT / "logs" / "download_manifest.tsv", sep="\t")
    if not args.dry_run:
        for role, p in [("SEALED_final_test", SEALED_V2), ("SEALED2_final_test", SEALED_DV)]:
            rec = man.loc[man.role == role, "sha256"].iloc[0]
            ok = sha(p) == rec
            say(f"Mühür {p.name}: SHA-256 {'AYNI' if ok else 'FARKLI'} ({rec[:16]}…)")
            assert ok, "mühür bozulmuş"
        for tag, path, _ in registry():
            ref = (path.parent / "locked_models.sha256").read_text().split()[0]
            assert sha(path) == ref, f"kilit dosyası değişmiş: {tag}"
        say("Tüm kilit dosyalarının SHA-256'sı kayıtla aynı.")

    # 2. etiketsiz ön-işleme
    meta = pd.read_csv(RESULTS / "inventory" / "GSE192902_samples.tsv", sep="\t", dtype=str)
    coh = "Validation 1" if args.dry_run else "Validation 2"
    set_ = "cfRNA_val1" if args.dry_run else "FINAL_TEST"
    path = RAW / "GSE192902" / "GSE192902_counts_Validation1_postQC.csv.gz" if args.dry_run else SEALED_V2
    m2 = meta[(meta.cohort == coh) & meta.tp_label.isin(["≤12", "13–20"])]
    women = pd.read_csv(ROOT / "data" / "processed" / "selection" / "women_all.tsv", sep="\t", dtype=str)
    w2 = women[(women.cohort == "GSE192902") & (women.set == set_)]
    F2, rn7_2, qc2 = prep_gse192902(path, m2.sample_name.tolist())
    say(f"Val-2 ön-işleme (etiketsiz): {qc2}; özellik {F2.shape[1]}")
    cohorts = {"Val2": (F2, rn7_2, m2, w2)}
    if not args.dry_run:
        a = read_sample_annotation(RAW / "GSE154377" / "GSE154377_series_matrix.txt.gz")
        t = a.title.str.extract(r"^(?P<grp>.+?) Patient (?P<pid>\d+), (?P<tp>.+?) \(RNA-seq\)")
        dv = pd.DataFrame({"gsm": a.index}).join(t.reset_index(drop=True))
        dv = dv[dv.tp == "1st Trimester"]
        F3, rn7_3, qc3 = prep_delvecchio(set(dv.gsm))
        say(f"Del Vecchio ön-işleme (etiketsiz): {qc3}; özellik {F3.shape[1]}")
    (OUT / "preprocessing_qc.json").write_text(json.dumps({"val2": qc2, "delvecchio": qc3 if not args.dry_run else None},
                                                          indent=2, default=str), encoding="utf-8")
    if not args.dry_run:
        opened.write_text(f"Etiketler birleştirildi: {datetime.now(timezone.utc).isoformat()}\n"
                          f"rerun_reason: {args.rerun_reason}\n", encoding="utf-8")

    # 3. etiketleri birleştir ve değerlendir
    say("\n=== ETİKETLER BİRLEŞTİRİLDİ ===")
    e2i, _ = gene_maps()
    cat = pd.read_csv(RESULTS / "literature" / "published_signatures.tsv", sep="\t", dtype=str)
    gi = pd.read_csv(GENE_INFO, sep="\t", usecols=["GeneID", "Symbol"], dtype=str).drop_duplicates("Symbol", keep=False)
    s2i = gi.set_index("Symbol").GeneID.to_dict()
    rows, scores_keep = [], {}

    def evaluate(cname, F, rn7, y_sample, ga_sample, woman_mask):
        Zw = zscore(F[woman_mask])
        Za = zscore(F)
        yw, gw = y_sample[woman_mask], ga_sample[woman_mask]
        for tag, p, unit in registry():
            L = json.loads(p.read_text(encoding="utf-8"))
            for meth, md in L["models"].items():
                if unit == "all_samples":
                    sa, nh, nf = apply_model(md, Za, ga_sample)
                    sw = sa[woman_mask]
                    rows.append({"cohort": cname, "model": f"{tag}_{meth}", "level": "örnek", "n_feat_present": f"{nh}/{nf}",
                                 **metrics(y_sample, sa, ga_sample)})
                else:
                    sw, nh, nf = apply_model(md, Zw, gw)
                rows.append({"cohort": cname, "model": f"{tag}_{meth}", "level": "kadın", "n_feat_present": f"{nh}/{nf}",
                             **metrics(yw, sw, gw)})
                scores_keep[(cname, f"{tag}_{meth}")] = sw
        rf = rn7.values[woman_mask]
        rows.append({"cohort": cname, "model": "RN7SL_fraction", "level": "kadın", "n_feat_present": "1/1", **metrics(yw, rf, gw)})
        scores_keep[(cname, "RN7SL_fraction")] = rf
        for sid in LIT_SIGS + ["Moufarrej2022_cfRNA_LR18_logOR"]:
            base = sid.replace("_logOR", "")
            sg = cat[cat.signature_id == base]
            feats, wts = [], []
            for _, r in sg.iterrows():
                fid = r.feature_id_if_given
                g = e2i.get(fid) if isinstance(fid, str) and fid.startswith("ENSG") else s2i.get(str(r.feature).strip())
                if g is None or g not in Zw.columns or r.direction_in_pe not in ("up", "down"):
                    continue
                feats.append(g)
                wts.append(np.log(float(r.weight.replace("OR=", ""))) if sid.endswith("_logOR") else (1.0 if r.direction_in_pe == "up" else -1.0))
            if len(feats) / max(len(sg), 1) < 0.5:
                rows.append({"cohort": cname, "model": f"LIT_{sid}", "level": "kadın", "n_feat_present": f"{len(feats)}/{len(sg)}", "auc": np.nan})
                continue
            s = Zw[feats].values @ np.array(wts) / np.abs(wts).sum()
            rows.append({"cohort": cname, "model": f"LIT_{sid}", "level": "kadın", "n_feat_present": f"{len(feats)}/{len(sg)}",
                         **metrics(yw, s, gw)})
            scores_keep[(cname, f"LIT_{sid}")] = s
        return yw

    lab2 = m2.set_index("sample_name").loc[F2.index]
    y2 = (lab2.disease != "control").astype(int).values
    g2 = (lab2.tp_label == "13–20").astype(int).values
    wm2 = F2.index.isin(w2.sample_id).astype(bool)
    assert wm2.sum() == len(w2), "kadın düzeyi örnek eşleşmesi tutmadı"
    yw2 = evaluate("Val2", F2, rn7_2, y2, g2, wm2)
    if not args.dry_run:
        lab3 = dv.set_index("gsm").loc[F3.index]
        grp = lab3.grp.map({"Normal Pregnancy": "Control", "Preeclampsia/gestational hypertension": "PE_GH",
                            "Gestational diabetes": "GDM", "Chronic  hypertension": "CHTN"})
        for case in ["PE_GH", "GDM"]:
            msk = grp.isin(["Control", case]).values
            evaluate(f"DelVecchio_{case}", F3[msk], rn7_3[msk], (grp[msk] == case).astype(int).values,
                     np.zeros(msk.sum(), int), np.ones(msk.sum(), bool))

    R = pd.DataFrame(rows)
    R.to_csv(OUT / "final_results.tsv", sep="\t", index=False)
    say("\nBİRİNCİL SONUÇ")
    prim = R[(R.cohort == "Val2") & (R.model == "V7_MOGAFS") & (R.level == "kadın")].iloc[0]
    say(f"  V7 MOGAFS, Val-2 kadın düzeyi: AUC {prim.auc:.3f} [%95 GA {prim.ci_lo:.3f}–{prim.ci_hi:.3f}] "
        f"(n={int(prim.n_case)}/{int(prim.n_ctrl)}; özellik {prim.n_feat_present})")
    for k in ["Moufarrej2022_cfRNA_LR18_logOR", "Moufarrej2022_cfRNA_LR18"]:
        if ("Val2", f"LIT_{k}") in scores_keep:
            a1, a2, dlt, z, p = delong_test(yw2, scores_keep[("Val2", "V7_MOGAFS")], scores_keep[("Val2", f"LIT_{k}")])
            say(f"  DeLong V7 vs {k}: {a1:.3f} vs {a2:.3f}, fark {dlt:+.3f}, p = {p:.3f}")
    say("\nTÜM SONUÇLAR (kadın düzeyi; örnek düzeyi V7/V8 için ayrıca)")
    for _, r in R.sort_values(["cohort", "level", "auc"], ascending=[True, True, False]).iterrows():
        if pd.isna(r.auc):
            say(f"  {r.cohort:18s} {r.level:5s} {r.model:42s} kapsam yetersiz ({r.n_feat_present})")
            continue
        extra = "".join(f" | {k.replace('auc_', '')}: {r[k]:.2f}" for k in ["auc_≤12", "auc_13–20"] if k in r and pd.notna(r[k]))
        say(f"  {r.cohort:18s} {r.level:5s} {r.model:42s} AUC {r.auc:.3f} [{r.ci_lo:.2f}–{r.ci_hi:.2f}] "
            f"n={int(r.n_case)}/{int(r.n_ctrl)} özellik {r.n_feat_present}{extra}")
    (OUT / "final_report.txt").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
