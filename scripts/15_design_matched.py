"""
Literatürdeki değerlendirme tasarımlarıyla EŞLEŞTİRİLMİŞ karşılaştırma (keşifsel; her satır tasarım adıyla etiketlenir).

"Castillo-Marco tarzı dış değerlendirme": imzanın GEN KÜMESİ sabit tutulur, katsayılar dış kohortun içinde
birini-dışarıda-bırak (LOOCV) lojistik regresyonla yeniden öğrenilir; AUC LOO tahminlerinden hesaplanır.
Bu bir bağımsız doğrulama DEĞİLDİR (dış kohortun etiketleri ağırlıkları belirler); gen kümesinin taşınabilirliğini
ölçer. Önceden kaydedilmiş birincil sonuç (V7, Val-2 0.418) değişmez ve ayrıca raporlanır.
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geo_utils import RESULTS, ROOT, pd  # noqa: E402
from panel_utils import boot_auc  # noqa: E402
from sklearn.linear_model import LogisticRegression  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402
from sklearn.model_selection import LeaveOneOut, cross_val_predict  # noqa: E402

OUT = RESULTS / "comparison"


def load_module(name, file):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).resolve().parent / file)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def loocv_auc(X, y):
    p = cross_val_predict(LogisticRegression(C=1.0, solver="liblinear"), X, y, cv=LeaveOneOut(),
                          method="decision_function")
    lo, hi = boot_auc(y, p, n=2000)
    return roc_auc_score(y, p), lo, hi


def main():
    m11 = load_module("p11", "11_final_test.py")
    meta = pd.read_csv(RESULTS / "inventory" / "GSE192902_samples.tsv", sep="\t", dtype=str)
    women = pd.read_csv(ROOT / "data" / "processed" / "selection" / "women_all.tsv", sep="\t", dtype=str)
    m2 = meta[(meta.cohort == "Validation 2") & meta.tp_label.isin(["≤12", "13–20"])]
    F2, _, _ = m11.prep_gse192902(m11.SEALED_V2, m2.sample_name.tolist())
    w2 = women[(women.cohort == "GSE192902") & (women.set == "FINAL_TEST")]
    wm = F2.index.isin(w2.sample_id)
    lab = m2.set_index("sample_name").loc[F2.index]
    y = (lab.disease != "control").astype(int).values[wm]
    Z = F2[wm]
    Z = (Z - Z.mean()) / Z.std(ddof=1).replace(0, 1)

    sets = {}
    for tag, p in [("V7 (birincil)", RESULTS / "plasma/variants/V7/locked_models.json"),
                   ("V0", RESULTS / "plasma/all_genes/locked_models.json"),
                   ("V0pc (protein kodlayan)", RESULTS / "plasma/protein_coding/locked_models.json")]:
        md = json.loads(p.read_text(encoding="utf-8"))["models"]["MOGAFS"]
        sets["XAI-MOGAFS " + tag] = md.get("features") or md["genes"]
    P = pd.read_csv(RESULTS / "site_invariant/loso_panels.tsv", sep="\t")
    si = P[(P.held_out == "S2_Val2") & (P.method == "SI-MOGAFS")].symbols.iloc[0].split(";")
    uni = pd.read_csv(ROOT / "data/processed/gene_universe.tsv", sep="\t", dtype={"entrez": str}).set_index("entrez")
    s2e = {s: e for e, s in uni.symbol.dropna().items()}
    sets["XAI-MOGAFS-SI (Val-2'siz eğitilmiş genler)"] = [s2e[s] for s in si if s in s2e]
    cat = pd.read_csv(RESULTS / "literature/published_signatures.tsv", sep="\t", dtype=str)
    e2i, _ = m11.gene_maps()
    mo = cat[cat.signature_id == "Moufarrej2022_cfRNA_LR18"]
    sets["Moufarrej 2022 (18 gen)"] = [e2i.get(r.feature_id_if_given) for _, r in mo.iterrows() if e2i.get(r.feature_id_if_given)]

    rows = []
    for name, feats in sets.items():
        have = [f for f in feats if f in Z.columns]
        a, lo, hi = loocv_auc(Z[have].values, y)
        rows.append({"gen_kumesi": name, "n_gen": f"{len(have)}/{len(feats)}", "tasarım": "Castillo-Marco tarzı: Val-2 içinde LOOCV yeniden ayar",
                     "auc": a, "ci_lo": lo, "ci_hi": hi})
        print(f"  {name:48s} {len(have)}/{len(feats)} gen | Val-2 LOOCV yeniden ayar AUC {a:.3f} [{lo:.2f}–{hi:.2f}]")
    pd.DataFrame(rows).to_csv(OUT / "design_matched_castillo_style.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
