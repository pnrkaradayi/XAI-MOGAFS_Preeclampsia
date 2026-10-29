"""
Adım 1 — GSE85307 envanteri ve doğrulama.

Beklenen (devir dosyası §3, §10):
  * 157 kadın: 47 preeklampsi / 110 kontrol
  * kan alımı haftası 10,0–18,9 (medyan 13,4)
  * mRNA (GPL6244) ve miRNA (GPL19162) aynı kadınlardan

Çıktılar results/inventory/ altına yazılır:
  GSE85307_samples.tsv   örnek düzeyi tablo (her iki katman)
  GSE85307_women.tsv     kadın düzeyi tablo (katman eşleşmesiyle)
  GSE85307_report.txt    doğrulama raporu
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geo_utils import RAW, RESULTS, read_expression, read_sample_annotation, to_num  # noqa: E402

SRC = RAW / "GSE85307"
OUT = RESULTS / "inventory"
OUT.mkdir(parents=True, exist_ok=True)

LAYERS = {
    "mRNA": SRC / "GSE85307-GPL6244_series_matrix.txt.gz",
    "miRNA": SRC / "GSE85307-GPL19162_series_matrix.txt.gz",
}
EXPECTED = {"n_women": 157, "PE": 47, "control": 110, "ga_min": 10.0, "ga_max": 18.9}

lines = []


def say(msg=""):
    print(msg)
    lines.append(msg)


def check(label, ok):
    say(f"  [{'TAMAM' if ok else 'UYUMSUZ'}] {label}")
    return ok


import pandas as pd  # noqa: E402

tables = []
for layer, path in LAYERS.items():
    ann = read_sample_annotation(path)
    expr = read_expression(path)
    t = pd.DataFrame({
        "gsm": ann.index,
        "layer": layer,
        "woman_id": ann["ch:recruitment_id"].values,
        "condition": ann["ch:pregnancy_condition"].values,
        "ga_draw_wk": to_num(ann["ch:blood_draw_at_gestation_week"]).values,
        "ga_delivery_wk": to_num(ann["ch:gestation_weeks"]).values,
        "age_y": to_num(ann["ch:age_years"]).values,
        "bmi": to_num(ann["ch:bmi_kg/m2"]).values,
        "ethnicity": ann["ch:ethnicity"].values,
        "asthma": ann["ch:asthma_0=healthy_1=asthmatic"].values,
        "batch": (ann.get("ch:microarray_hybridization_date", ann.get("ch:batch_id"))).values,
    })
    t["in_matrix"] = t["gsm"].isin(expr.columns)
    # Teknik tekrarlar: başlıkta "Replicate k" geçer (aynı kan örneği, farklı plaka).
    rep = ann["title"].str.extract(r"Replicate (\d+)", expand=False)
    t["replicate"] = pd.to_numeric(rep, errors="coerce").values
    t["primary"] = t["replicate"].isna() | (t["replicate"] == 1)
    tables.append(t)

    say(f"=== {layer}  ({path.name})")
    say(f"  matris boyutu: {expr.shape[0]} özellik x {expr.shape[1]} örnek")
    say(f"  eksik değer oranı: {expr.isna().mean().mean():.3f}")
    say(f"  değer aralığı: {expr.min().min():.2f} – {expr.max().max():.2f}")
    say(f"  anotasyondaki örnek: {len(t)}, matristeki: {t['in_matrix'].sum()}")
    say(f"  benzersiz kadın: {t['woman_id'].nunique()}  (kadın başına örnek maks: "
        f"{t['woman_id'].value_counts().max()})")
    say(f"  grup (örnek düzeyi): {t['condition'].value_counts().to_dict()}")
    reps = t[t["replicate"].notna()]
    if len(reps):
        say(f"  teknik tekrar: {reps['woman_id'].nunique()} kadında {len(reps)} örnek "
            f"({', '.join(sorted(reps['woman_id'].unique()))}); yalnız 'Replicate 1' tutulur")
    say()

samples = pd.concat(tables, ignore_index=True)
samples.to_csv(OUT / "GSE85307_samples.tsv", sep="\t", index=False)

# ---- Kadın düzeyi: katmanlar recruitment_id ile eşleştirilir (yalnız birincil örnekler)
primary = samples[samples["primary"]]
assert not primary.duplicated(["layer", "woman_id"]).any(), "kadın başına birden çok birincil örnek"
per_layer = {l: g.set_index("woman_id") for l, g in primary.groupby("layer")}
women = pd.DataFrame(index=sorted(set(samples["woman_id"])))
women.index.name = "woman_id"
for l, g in per_layer.items():
    women[f"gsm_{l}"] = g["gsm"]
    women[f"condition_{l}"] = g["condition"]
    women[f"ga_draw_{l}"] = g["ga_draw_wk"]
women["has_mRNA"] = women["gsm_mRNA"].notna()
women["has_miRNA"] = women["gsm_miRNA"].notna()

say("=== Katmanlar arası tutarlılık (aynı kadın)")
both = women[women.has_mRNA & women.has_miRNA]
say(f"  yalnız mRNA: {(women.has_mRNA & ~women.has_miRNA).sum()}, "
    f"yalnız miRNA: {(~women.has_mRNA & women.has_miRNA).sum()}, ikisi birden: {len(both)}")
cond_mismatch = (both["condition_mRNA"] != both["condition_miRNA"]).sum()
ga_diff = (both["ga_draw_mRNA"] - both["ga_draw_miRNA"]).abs()
say(f"  grup etiketi uyuşmazlığı: {cond_mismatch}")
say(f"  kan alım haftası farkı (maks): {ga_diff.max():.3f} hf  -> aynı kan örneği: {bool((ga_diff < 1e-6).all())}")
say()

women["condition"] = women["condition_mRNA"].fillna(women["condition_miRNA"])
women["ga_draw_wk"] = women["ga_draw_mRNA"].fillna(women["ga_draw_miRNA"])
meta = primary.drop_duplicates("woman_id").set_index("woman_id")
for c in ["ga_delivery_wk", "age_y", "bmi", "ethnicity", "asthma"]:
    women[c] = meta[c]
women.to_csv(OUT / "GSE85307_women.tsv", sep="\t")

# ---- Devir dosyasındaki rakamlarla karşılaştırma
say("=== Devir dosyasıyla karşılaştırma (kadın düzeyi)")
counts = women["condition"].value_counts()
n_pe = int(counts.get("Preeclampsia", 0))
n_ctrl = int(counts.get("Normal", 0))
ga = women["ga_draw_wk"]
ok = all([
    check(f"kadın sayısı = {len(women)} (beklenen {EXPECTED['n_women']})", len(women) == EXPECTED["n_women"]),
    check(f"preeklampsi = {n_pe} (beklenen {EXPECTED['PE']})", n_pe == EXPECTED["PE"]),
    check(f"kontrol = {n_ctrl} (beklenen {EXPECTED['control']})", n_ctrl == EXPECTED["control"]),
    check(f"kan alımı en erken = {ga.min():.2f} hf (beklenen ≥ {EXPECTED['ga_min']})", ga.min() >= EXPECTED["ga_min"] - 1e-9),
    check(f"kan alımı en geç = {ga.max():.2f} hf (beklenen ≤ {EXPECTED['ga_max']})", ga.max() <= EXPECTED["ga_max"] + 0.05),
])
say(f"  kan alımı medyanı = {ga.median():.2f} hf (devir dosyası: 13,4)")
say()

say("=== Grup x kan alım haftası (kadın sayısı)")
bins = pd.cut(ga, bins=[10, 12, 14, 16, 18, 20], right=False,
              labels=["10–11.9", "12–13.9", "14–15.9", "16–17.9", "18–19.9"])
xt = pd.crosstab(women["condition"], bins, margins=True, margins_name="Toplam")
say(xt.to_string())
say()

say("=== Grup bazında tanımlayıcılar (medyan [min–maks])")
for c in ["ga_draw_wk", "ga_delivery_wk", "age_y", "bmi"]:
    parts = []
    for grp, g in women.groupby("condition"):
        v = g[c].dropna()
        parts.append(f"{grp}: {v.median():.1f} [{v.min():.1f}–{v.max():.1f}] (n={len(v)})")
    say(f"  {c:15s} " + " | ".join(parts))
say()

say("SONUÇ: " + ("Tüm beklenen rakamlar doğrulandı." if ok else "Bazı rakamlar devir dosyasıyla UYUŞMUYOR — yukarıya bakın."))

(OUT / "GSE85307_report.txt").write_text("\n".join(lines), encoding="utf-8")
