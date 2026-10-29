"""
Makale için ek değerlendirme metrikleri ve şekil verileri (yeni model EĞİTİLMEZ; kilitli modeller yeniden uygulanır).

Çıktılar (results/paper/):
  wb_panels_by_syndrome.tsv  kilitli tam kan panellerimiz (PE: 4 yöntem; sPTB: 3 yöntem) x her kohort/sendrom hücresi
  val2_full_metrics.tsv      Val-2 kadın düzeyi: AUC, %95 GA, AUPRC, accuracy, dengeli accuracy, duyarlılık, özgüllük,
                             PPV, NPV, F1, MCC, LR+, LR− (eşikler YALNIZ eğitim verisinden)
  val2_scores.tsv            Val-2 kadın düzeyi skorlar (ROC eğrileri için)
  loso_full_metrics.tsv      LOSO (keşifsel) tüm yöntemler x merkez, aynı metrikler
  loso_delong.tsv            LOSO katlarında yöntem çiftleri için DeLong testi
  rn7sl_fraction.tsv         örnek başına RN7SL* kütüphane payı (kohort x grup)
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geo_utils import RAW, RESULTS, ROOT, pd  # noqa: E402
from panel_utils import boot_auc, delong_test  # noqa: E402
from scipy.stats import mannwhitneyu  # noqa: E402
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score, roc_curve  # noqa: E402

OUT = RESULTS / "paper"
OUT.mkdir(parents=True, exist_ok=True)
WB_CELLS = [("GSE85307", "PE", "eğitim (apparent)"),
            ("GSE149440", "PE", "geliştirme"), ("GSE149440", "sPTD", "geliştirme"), ("GSE149440", "PPROM", "geliştirme"),
            ("GSE59491", "sPTD", "bağımsız"), ("GSE59491", "PPROM", "bağımsız"),
            ("GSE142974", "sPTD", "bağımsız"),
            ("GSE108497", "PE", "bağımsız"), ("GSE108497", "FGR_SGA_FD", "bağımsız")]


def load_module(name, file):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).resolve().parent / file)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def youden(y, s):
    fpr, tpr, thr = roc_curve(y, s)
    return float(thr[np.argmax(tpr - fpr)])


def full_metrics(y, s, thr=None, n_boot=2000):
    """Eşikten bağımsız (AUC, GA, AUPRC) + eşik verildiyse sınıflandırma metrikleri."""
    y = np.asarray(y).astype(int)
    s = np.asarray(s, float)
    lo, hi = boot_auc(y, s, n=n_boot)
    out = {"n_case": int(y.sum()), "n_ctrl": int((1 - y).sum()), "prevalence": y.mean(),
           "auc": roc_auc_score(y, s), "ci_lo": lo, "ci_hi": hi, "auprc": average_precision_score(y, s)}
    if thr is None:
        return out
    p = (s >= thr).astype(int)
    tp, fp = int(((p == 1) & (y == 1)).sum()), int(((p == 1) & (y == 0)).sum())
    tn, fn = int(((p == 0) & (y == 0)).sum()), int(((p == 0) & (y == 1)).sum())
    sens = tp / max(tp + fn, 1)
    spec = tn / max(tn + fp, 1)
    ppv = tp / (tp + fp) if tp + fp else np.nan
    npv = tn / (tn + fn) if tn + fn else np.nan
    out.update({"threshold": thr, "tp": tp, "fp": fp, "tn": tn, "fn": fn, "acc": (tp + tn) / len(y),
                "bacc": (sens + spec) / 2, "sens": sens, "spec": spec, "ppv": ppv, "npv": npv,
                "f1": 2 * tp / (2 * tp + fp + fn) if tp else 0.0, "mcc": matthews_corrcoef(y, p),
                "lr_pos": sens / (1 - spec) if spec < 1 else np.inf,
                "lr_neg": (1 - sens) / spec if spec > 0 else np.inf})
    return out


def zscore(M):
    return (M - M.mean()) / M.std(ddof=1).replace(0, 1)


# ------------------------------------------------------------------ 1. tam kan panelleri x sendrom
def wb_panels_by_syndrome():
    w = pd.read_csv(ROOT / "data/processed/selection/women_all.tsv", sep="\t", dtype=str)
    panels = {}
    for m, md in json.loads((RESULTS / "panel/final_panels.json").read_text(encoding="utf-8")).items():
        panels[("PE paneli", m)] = md
    for m, md in json.loads((RESULTS / "ptb/locked_models.json").read_text(encoding="utf-8")).items():
        panels[("sPTB paneli", m)] = md
    rows = []
    for coh in dict.fromkeys(c for c, _, _ in WB_CELLS):
        E = pd.read_csv(ROOT / f"data/processed/expr/{coh}.tsv.gz", sep="\t", dtype={"entrez": str}).set_index("entrez").T
        meta = w[w.cohort == coh].set_index("woman_id").loc[E.index]
        Z = zscore(E)
        for (fam, m), md in panels.items():
            have = [g in Z.columns for g in md["genes"]]
            G = [g for g in md["genes"] if g in Z.columns]
            c = np.array([cf for cf, h in zip(md["coef"], have) if h])
            s = Z[G].values @ c
            for c2, grp, role in WB_CELLS:
                if c2 != coh:
                    continue
                case = (meta.group == grp).values.copy()
                ctrl = (meta.group == "Control").values.copy()
                if coh == "GSE108497":
                    case &= (meta.background == "SLE").values
                    ctrl &= (meta.background == "SLE").values
                if fam == "sPTB paneli":
                    role = {"GSE149440": "eğitim (apparent)", "GSE59491": "kilitli dış test"}.get(coh, "bağımsız")
                mk = case | ctrl
                yy = case[mk].astype(int)
                p = mannwhitneyu(s[mk][yy == 1], s[mk][yy == 0]).pvalue
                rows.append({"panel": fam, "method": m, "cohort": coh, "group": grp, "role": role,
                             "coverage": f"{len(G)}/{len(md['genes'])}", "p_mw": p, **full_metrics(yy, s[mk], n_boot=1000)})
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "wb_panels_by_syndrome.tsv", sep="\t", index=False)
    print("wb_panels_by_syndrome:", len(R), "satır")


# ------------------------------------------------------------------ 2. Val-2 metrikleri + skorlar
def val2():
    m09c = load_module("p09c", "09c_plasma_repeated.py")
    m11 = load_module("p11", "11_final_test.py")
    meta = pd.read_csv(RESULTS / "inventory/GSE192902_samples.tsv", sep="\t", dtype=str)
    women = pd.read_csv(ROOT / "data/processed/selection/women_all.tsv", sep="\t", dtype=str)
    m2 = meta[(meta.cohort == "Validation 2") & meta.tp_label.isin(["≤12", "13–20"])]
    F2, rn7, _ = m11.prep_gse192902(m11.SEALED_V2, m2.sample_name.tolist())
    lab2 = m2.set_index("sample_name").loc[F2.index]
    w2 = women[(women.cohort == "GSE192902") & (women.set == "FINAL_TEST")]
    wm = F2.index.isin(w2.sample_id)
    y = (lab2.disease != "control").astype(int).values[wm]
    gb = (lab2.tp_label == "13–20").astype(int).values
    Za, Zw = zscore(F2), zscore(F2[wm])

    # V7 eşiği: Stanford eğitim örnekleri (13_comparison_table.py ile aynı)
    info, X, cols = m09c.build()
    F = pd.DataFrame(X, columns=cols)
    Ztr = pd.concat([zscore(F[info.cohort.values == c]) for c in [0, 1]]).loc[F.index]
    scores, rows = {}, []
    L7 = json.loads((RESULTS / "plasma/variants/V7/locked_models.json").read_text(encoding="utf-8"))["models"]
    for meth, md in L7.items():
        s_tr = md["intercept"] + Ztr[md["features"]].values @ np.array(md["coef"])
        thr = youden(info.y.values, s_tr)
        s, _, _ = m11.apply_model(md, Za, gb)
        scores[f"V7_{meth}"] = s[wm]
        rows.append({"model": f"V7_{meth}", "design": "önceden kayıtlı, mühürlü dış test",
                     **full_metrics(y, s[wm], thr)})
    # diğer kilitler: eşikten bağımsız metrikler
    for tag, p, unit in m11.registry():
        if tag == "V7":
            continue
        for meth, md in json.loads(p.read_text(encoding="utf-8"))["models"].items():
            if unit == "all_samples":
                s, _, _ = m11.apply_model(md, Za, gb)
                s = s[wm]
            else:
                s, _, _ = m11.apply_model(md, Zw, gb[wm])
            scores[f"{tag}_{meth}"] = s
            rows.append({"model": f"{tag}_{meth}", "design": "ikincil kilit, mühürlü dış test", **full_metrics(y, s)})
    # Moufarrej 18 gen: eşik kendi eğitim kohortu Discovery'den
    cat = pd.read_csv(RESULTS / "literature/published_signatures.tsv", sep="\t", dtype=str)
    e2i, _ = m11.gene_maps()
    mo = cat[cat.signature_id == "Moufarrej2022_cfRNA_LR18"]
    feats = [(e2i.get(r.feature_id_if_given), np.log(float(r.weight.replace("OR=", "")))) for _, r in mo.iterrows()]
    feats = [(g, wt) for g, wt in feats if g]
    E = pd.read_csv(ROOT / "data/processed/expr/GSE192902_cfRNA_discovery.tsv.gz", sep="\t",
                    dtype={"entrez": str}).set_index("entrez").T
    wd = women[(women.cohort == "GSE192902") & (women.set == "cfRNA_discovery")].set_index("woman_id").loc[E.index]
    gd = [g for g, _ in feats if g in E.columns]
    wts = np.array([wt for g, wt in feats if g in E.columns])
    thr_mo = youden((wd.group == "PE").astype(int).values, zscore(E[gd]).values @ wts / np.abs(wts).sum())
    gv = [g for g, _ in feats if g in Zw.columns]
    wv = np.array([wt for g, wt in feats if g in Zw.columns])
    smo = Zw[gv].values @ wv / np.abs(wv).sum()
    scores["Moufarrej18_logOR"] = smo
    rows.append({"model": "Moufarrej18_logOR", "design": "yayımlanmış imza, aynı kadınlar", **full_metrics(y, smo, thr_mo)})
    scores["RN7SL_fraction"] = rn7.values[wm]
    rows.append({"model": "RN7SL_fraction", "design": "tek özellik", **full_metrics(y, rn7.values[wm])})
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "val2_full_metrics.tsv", sep="\t", index=False)
    S = pd.DataFrame(scores, index=F2.index[wm]).assign(y=y, ga_bin=gb[wm])
    S.to_csv(OUT / "val2_scores.tsv", sep="\t")
    print("val2:", len(R), "model; V7 AUC", round(R.loc[R.model == "V7_MOGAFS", "auc"].iloc[0], 3))


# ------------------------------------------------------------------ 3. LOSO tüm metrikler + DeLong
def loso():
    S = pd.read_csv(RESULTS / "site_invariant/loso_scores.tsv", sep="\t")
    rows, dl = [], []
    for (held, meth), d in S.groupby(["held_out", "method"], sort=False):
        for lab, case in [("PE", "PE"), ("GDM", "GDM")]:
            m = d.earliest.astype(bool) & d.group.isin(["Control", case])
            yy = (d.group[m] == case).astype(int).values
            if yy.sum() < 3:
                continue
            rows.append({"held_out": held, "method": meth, "contrast": f"{lab} vs kontrol",
                         **full_metrics(yy, d.score[m].values, d.threshold.iloc[0])})
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "loso_full_metrics.tsv", sep="\t", index=False)
    for held, d in S[S.earliest.astype(bool) & S.group.isin(["Control", "PE"])].groupby("held_out", sort=False):
        P = d.pivot_table(index="sample", columns="method", values="score")
        yy = (d.drop_duplicates("sample").set_index("sample").loc[P.index].group == "PE").astype(int).values
        for b in ["MOGAFS", "META-TTEST5", "LASSO5", "TTEST5", "Moufarrej18"]:
            a1, a2, diff, z, p = delong_test(yy, P["SI-MOGAFS"].values, P[b].values)
            dl.append({"held_out": held, "comparison": f"SI-MOGAFS vs {b}", "auc_1": a1, "auc_2": a2, "diff": diff,
                       "z": z, "p": p})
    pd.DataFrame(dl).to_csv(OUT / "loso_delong.tsv", sep="\t", index=False)
    print("loso:", len(R), "satır")


# ------------------------------------------------------------------ 4. RN7SL payı
def rn7sl():
    m11 = load_module("p11", "11_final_test.py")
    meta = pd.read_csv(RESULTS / "inventory/GSE192902_samples.tsv", sep="\t", dtype=str)
    women = pd.read_csv(ROOT / "data/processed/selection/women_all.tsv", sep="\t", dtype=str)
    rows = []
    for coh, f, set_ in [("Discovery", RAW / "GSE192902/GSE192902_counts_Discovery_postQC.csv.gz", "cfRNA_discovery"),
                         ("Validation 1", RAW / "GSE192902/GSE192902_counts_Validation1_postQC.csv.gz", "cfRNA_val1"),
                         ("Validation 2", m11.SEALED_V2, "FINAL_TEST")]:
        m = meta[(meta.cohort == coh) & meta.tp_label.isin(["≤12", "13–20"])]
        ws = set(women.loc[(women.cohort == "GSE192902") & (women.set == set_), "sample_id"])
        m = m[m.sample_name.isin(ws)]
        raw = pd.read_csv(f).set_index("gene_name")[m.sample_name.tolist()]
        frac = raw[raw.index.str.startswith("RN7SL")].sum(0) / raw.sum(0)
        g = m.set_index("sample_name").loc[frac.index].disease.ne("control").map({True: "PE", False: "Control"})
        rows += [{"cohort": coh, "group": gg, "frac": v} for v, gg in zip(frac.values, g.values)]
    a = m11.read_sample_annotation(RAW / "GSE154377/GSE154377_series_matrix.txt.gz")
    t = a.title.str.extract(r"^(?P<grp>.+?) Patient (?P<pid>\d+), (?P<tp>.+?) \(RNA-seq\)")
    dv = pd.DataFrame({"gsm": a.index}).join(t.reset_index(drop=True))
    dv = dv[dv.tp == "1st Trimester"]
    _, fr3, _ = m11.prep_delvecchio(set(dv.gsm))
    grp = dv.set_index("gsm").loc[fr3.index].grp.map({"Normal Pregnancy": "Control", "Preeclampsia/gestational hypertension": "PE",
                                                      "Gestational diabetes": "GDM", "Chronic  hypertension": "CHTN"})
    rows += [{"cohort": "Del Vecchio", "group": g, "frac": v} for v, g in zip(fr3.values, grp.values)]
    pd.DataFrame(rows).to_csv(OUT / "rn7sl_fraction.tsv", sep="\t", index=False)
    print("rn7sl:", len(rows), "örnek")


if __name__ == "__main__":
    which = sys.argv[1:] or ["wb", "val2", "loso", "rn7sl"]
    for k in which:
        {"wb": wb_panels_by_syndrome, "val2": val2, "loso": loso, "rn7sl": rn7sl}[k]()
