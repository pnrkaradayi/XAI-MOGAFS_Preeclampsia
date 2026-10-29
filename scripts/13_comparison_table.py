"""
Literatür karşılaştırma tablosu — AUC ve doğruluk (accuracy, dengeli doğruluk, duyarlılık, özgüllük).

Kurallar:
  * Literatür satırları: makalede BİLDİRİLEN değerler (results/literature_specificity/reported_performance.tsv'deki
    kaynaklardan); dengeli doğruluk = (duyarlılık + özgüllük) / 2 yalnız makale ikisini de bildirdiyse.
  * Bizim yeniden testlerimiz: aynı kohortlarda AUC + %95 GA.
  * Eşikler YALNIZ eğitim verisinden (Youden): V7 -> Stanford (Discovery + Val-1) eğitim örnekleri;
    Moufarrej -> kendi eğitim kohortu Discovery; SI-MOGAFS -> her LOSO katının eğitim merkezleri.
  * Accuracy sınıf oranına bağlıdır (Val-2'de herkese "kontrol" demek %70 verir); karşılaştırma dengeli doğrulukla yapılır.
Çıktı: results/comparison/literature_comparison.tsv ve .md
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geo_utils import RESULTS, ROOT, pd  # noqa: E402
from panel_utils import boot_auc  # noqa: E402
from sklearn.metrics import roc_auc_score, roc_curve  # noqa: E402

OUT = RESULTS / "comparison"
OUT.mkdir(parents=True, exist_ok=True)


def load_module(name, file):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).resolve().parent / file)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def youden(y, s):
    fpr, tpr, thr = roc_curve(y, s)
    return float(thr[np.argmax(tpr - fpr)])


def cls_metrics(y, s, thr):
    y = np.asarray(y).astype(int)
    pred = (s >= thr).astype(int)
    sens, spec = pred[y == 1].mean(), 1 - pred[y == 0].mean()
    lo, hi = boot_auc(y, s, n=2000)
    return {"auc": roc_auc_score(y, s), "ci_lo": lo, "ci_hi": hi, "acc": (pred == y).mean(),
            "bacc": (sens + spec) / 2, "sens": sens, "spec": spec, "n": f"{y.sum()}/{(1 - y).sum()}"}


def zscore(M):
    return (M - M.mean()) / M.std(ddof=1).replace(0, 1)


def main():
    m09c = load_module("p09c", "09c_plasma_repeated.py")
    m11 = load_module("p11", "11_final_test.py")
    rows = []

    # ---- V7 (birincil kilit) — eşik Stanford eğitim örneklerinden
    info, X, cols = m09c.build()
    F = pd.DataFrame(X, columns=cols)
    Ztr = pd.concat([zscore(F[info.cohort.values == c]) for c in [0, 1]]).loc[F.index]
    L = json.loads((RESULTS / "plasma" / "variants" / "V7" / "locked_models.json").read_text(encoding="utf-8"))
    md = L["models"]["MOGAFS"]
    s_tr = md["intercept"] + Ztr[md["features"]].values @ np.array(md["coef"])
    thr_v7 = youden(info.y.values, s_tr)
    meta = pd.read_csv(RESULTS / "inventory" / "GSE192902_samples.tsv", sep="\t", dtype=str)
    women = pd.read_csv(ROOT / "data" / "processed" / "selection" / "women_all.tsv", sep="\t", dtype=str)
    m2 = meta[(meta.cohort == "Validation 2") & meta.tp_label.isin(["≤12", "13–20"])]
    F2, _, _ = m11.prep_gse192902(m11.SEALED_V2, m2.sample_name.tolist())
    Z2 = zscore(F2)
    lab2 = m2.set_index("sample_name").loc[F2.index]
    w2 = women[(women.cohort == "GSE192902") & (women.set == "FINAL_TEST")]
    wm = F2.index.isin(w2.sample_id)
    y2 = (lab2.disease != "control").astype(int).values[wm]
    s2 = md["intercept"] + Z2[md["features"]].values[wm] @ np.array(md["coef"])
    rows.append({"grup": "Bu çalışma", "model": "XAI-MOGAFS V7 (önceden kayıtlı BİRİNCİL)", "bölme": "plazma cfRNA",
                 "değerlendirme": "Val-2, mühürlü dış test (bir kez)", **cls_metrics(y2, s2, thr_v7)})

    # ---- Moufarrej 18 gen — eşik Discovery'den (kendi eğitim kohortu)
    cat = pd.read_csv(RESULTS / "literature" / "published_signatures.tsv", sep="\t", dtype=str)
    e2i, _ = m11.gene_maps()
    mo = cat[cat.signature_id == "Moufarrej2022_cfRNA_LR18"]
    feats = [(e2i.get(r.feature_id_if_given), np.log(float(r.weight.replace("OR=", "")))) for _, r in mo.iterrows()]
    feats = [(g, w) for g, w in feats if g]
    E = pd.read_csv(ROOT / "data" / "processed" / "expr" / "GSE192902_cfRNA_discovery.tsv.gz", sep="\t",
                    dtype={"entrez": str}).set_index("entrez").T
    wd = women[(women.cohort == "GSE192902") & (women.set == "cfRNA_discovery")].set_index("woman_id").loc[E.index]
    gd = [g for g, _ in feats if g in E.columns]
    wts = np.array([w for g, w in feats if g in E.columns])
    sd = zscore(E[gd]).values @ wts / np.abs(wts).sum()
    thr_mo = youden((wd.group == "PE").astype(int).values, sd)
    Zw2 = zscore(F2[wm])
    gv = [g for g, _ in feats if g in Zw2.columns]
    wv = np.array([w for g, w in feats if g in Zw2.columns])
    smo = Zw2[gv].values @ wv / np.abs(wv).sum()
    rows.append({"grup": "Literatür (bizim yeniden testimiz)", "model": "Moufarrej 2022, 18 gen (log OR)",
                 "bölme": "plazma cfRNA", "değerlendirme": "Val-2, aynı kadınlar", **cls_metrics(y2, smo, thr_mo)})

    # ---- SI-MOGAFS LOSO (keşifsel)
    R = pd.read_csv(RESULTS / "site_invariant" / "loso_results.tsv", sep="\t")
    for _, r in R[(R.level == "kadın") & R.method.isin(["SI-MOGAFS", "LASSO5", "Moufarrej18"])].iterrows():
        rows.append({"grup": "Bu çalışma — keşifsel (Val-2 sonrası)" if r.method != "Moufarrej18" else "Literatür (bizim yeniden testimiz)",
                     "model": {"SI-MOGAFS": "XAI-MOGAFS-SI", "LASSO5": "LASSO (≤5)", "Moufarrej18": "Moufarrej 2022, 18 gen"}[r.method]
                     + f" — merkez-dışarıda: {r.held_out}", "bölme": "plazma cfRNA",
                     "değerlendirme": "merkez-dışarıda-bırakma" + (" (eğitim kohortu, bağımsız değil)" if r.method == "Moufarrej18" and r.held_out == "S1_Stanford" else ""),
                     "auc": r.auc, "ci_lo": r.ci_lo, "ci_hi": r.ci_hi, "acc": r.acc, "bacc": r.bacc,
                     "sens": r.sens, "spec": r.spec, "n": f"{r.n_case}/{r.n_ctrl}"})

    T = pd.DataFrame(rows)
    T.to_csv(OUT / "computed_metrics.tsv", sep="\t", index=False)
    print(T.round(3).to_string(index=False))
    print(f"\nEşikler: V7 {thr_v7:.3f} (Stanford eğitim), Moufarrej {thr_mo:.3f} (Discovery)")


if __name__ == "__main__":
    main()
