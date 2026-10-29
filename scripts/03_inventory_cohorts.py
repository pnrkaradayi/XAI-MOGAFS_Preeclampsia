"""
Adım 1 (devam) — Özgüllük ve son test kohortlarının envanteri.

Protokol §5.1: her seride kadın sayısı, grup dağılımı, örnekleme haftası tablosu;
GSM kesişimleriyle mükerrer kontrolü. Devir dosyasındaki rakamlar EXPECTED olarak
yazıldı; her biri TAMAM / UYUMSUZ olarak raporlanır.

Mühürlü dosya (data/sealed/) AÇILMAZ; yalnızca SHA-256'sı manifest ile karşılaştırılır.

Çıktılar (results/inventory/):
  <GSE>_samples.tsv        örnek düzeyi temiz anotasyon
  cohorts_report.txt       doğrulama raporu
"""
import hashlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geo_utils import RAW, RESULTS, ROOT, pd, read_sample_annotation, to_num  # noqa: E402

OUT = RESULTS / "inventory"
OUT.mkdir(parents=True, exist_ok=True)
lines = []


def say(msg=""):
    print(msg)
    lines.append(str(msg))


def check(label, got, expected):
    ok = got == expected
    say(f"  [{'TAMAM' if ok else 'UYUMSUZ'}] {label}: {got} (beklenen {expected})")
    return ok


def ann(gse, fname=None):
    return read_sample_annotation(RAW / gse / (fname or f"{gse}_series_matrix.txt.gz"))


results = []

# =====================================================================
# GSE149440 — Özgüllük 1 (Tarca 2021)
# =====================================================================
say("=" * 70)
say("GSE149440 — ÖZGÜLLÜK 1 (tam kan, GPL28460)")
say("=" * 70)
a = ann("GSE149440")
s149440 = pd.DataFrame({
    "gsm": a.index,
    "woman_id": a["ch:individual"].values,
    "group": a["ch:group"].values,
    "ga_wk": to_num(a["ch:gestational age"]).values,
    "ga_delivery_wk": to_num(a["ch:gadelivery"]).values,
    "subgroup": a["ch:subgroup"].values,
    "batch": a["ch:batch"].values,
    "orig_train_flag": a["ch:train"].values,
})
s149440.to_csv(OUT / "GSE149440_samples.tsv", sep="\t", index=False)
w = s149440.groupby("woman_id").agg(group=("group", "first"), n_groups=("group", "nunique"),
                                    n=("gsm", "size"), ga_first=("ga_wk", "min"))
say(f"  örnek: {len(s149440)}, kadın başına örnek medyan {w.n.median():.0f} (maks {w.n.max()})")
say(f"  grup etiketi kadın içinde tutarsız: {(w.n_groups > 1).sum()}")
results += [
    check("kadın sayısı", len(w), 149),
    check("örnek sayısı", len(s149440), 735),
]
gc = w.group.value_counts()
for g, e in [("Early_Preeclampsia", 13), ("PPROM", 37), ("sPTD", 34), ("Control", 65)]:
    results.append(check(f"{g} (kadın)", int(gc.get(g, 0)), e))
results.append(check("<20 hf örneği olan kadın", int((w.ga_first < 20).sum()), 123))
say()
say("  Pencereye göre kadın sayısı (kadının o penceredeki EN ERKEN örneği):")
rows = {}
for lab, lo, hi in [("<20 hf", 0, 20), ("10–18.9 hf (eğitim penceresi)", 10, 19), ("17–23 hf (GSE59491 penceresi)", 17, 24)]:
    sub = s149440[(s149440.ga_wk >= lo) & (s149440.ga_wk < hi)]
    rows[lab] = sub.sort_values("ga_wk").drop_duplicates("woman_id")["group"].value_counts()
say(pd.DataFrame(rows).fillna(0).astype(int).T.to_string())
say()

# =====================================================================
# GSE59491 — Özgüllük 2 (All Our Babies)
# =====================================================================
say("=" * 70)
say("GSE59491 — ÖZGÜLLÜK 2 (tam kan, GPL18964)")
say("=" * 70)
a = ann("GSE59491")
outcome_map = {"1": "sPTL", "2": "PPROM", "3": "Term"}
s59491 = pd.DataFrame({
    "gsm": a.index,
    "woman_id": a["ch:individual id"].values,
    "timepoint": to_num(a["ch:sample collected at time pt (17-23weeks 1, 27-33weeks 2)"]).values,
    "ga_wk": to_num(a["ch:gest age at sampling (weeks)"]).values,
    "ga_delivery_wk": to_num(a["ch:gest age at delivery (weeks)"]).values,
    "group": a["ch:birthing outcome (sptl 1, pprom 2, term 3)"].map(outcome_map).values,
})
s59491.to_csv(OUT / "GSE59491_samples.tsv", sep="\t", index=False)
w = s59491.groupby("woman_id").agg(group=("group", "first"), n_groups=("group", "nunique"), n=("gsm", "size"))
say(f"  grup etiketi kadın içinde tutarsız: {(w.n_groups > 1).sum()}")
results += [check("örnek sayısı", len(s59491), 326), check("kadın sayısı", len(w), 165)]
gc = w.group.value_counts()
for g, e in [("sPTL", 15), ("PPROM", 36), ("Term", 114)]:
    results.append(check(f"{g} (kadın)", int(gc.get(g, 0)), e))
tp1 = s59491[s59491.timepoint == 1]
say(f"  zaman noktası 1: {len(tp1)} örnek / {tp1.woman_id.nunique()} kadın, "
    f"hafta {tp1.ga_wk.min():.0f}–{tp1.ga_wk.max():.0f}")
say(f"  zaman noktası 1 kadın başına tek örnek mi: {not tp1.woman_id.duplicated().any()}")
say("  Zaman noktası 1, grup x örnekleme haftası (kadın):")
say(pd.crosstab(tp1.group, tp1.ga_wk.astype(int), margins=True, margins_name="Top").to_string())
ov = tp1[tp1.ga_wk <= 18]
say(f"  [§9 açık soru] eğitim penceresiyle (≤18.9 hf) kesişen TP1 örnekleri: {len(ov)} kadın "
    f"{ov.group.value_counts().to_dict()}")
say("  NOT: GSE59491'de preeklampsi grubu yok; özgüllük testi yalnız sPTL/PPROM vs term.")
say()

# =====================================================================
# GSE192902 — Bölme + son test (Moufarrej 2022, plazma cfRNA)
# =====================================================================
say("=" * 70)
say("GSE192902 — BÖLME + SON TEST (plazma cfRNA, GPL24676)")
say("=" * 70)
a = ann("GSE192902")
tp_map = {"1": "≤12", "2": "13–20", "3": "≥23", "4": "postpartum"}
title = a["title"].astype(str)
s192902 = pd.DataFrame({
    "gsm": a.index,
    "sample_name": title.values,
    "woman_id": title.str.split("_").str[0].values,
    "timepoint": title.str.split("_").str[1].str.replace(r"\.\d+$", "", regex=True).values,
    "within_window_repeat": title.str.contains(r"_\d+\.\d+$").values,
    "window": a["ch:sampling time group"].values,
    "disease": a["ch:disease"].values,
    "cohort": a["ch:cohort"].values,
})
s192902["tp_label"] = s192902["timepoint"].map(tp_map)
s192902.to_csv(OUT / "GSE192902_samples.tsv", sep="\t", index=False)
results.append(check("örnek sayısı", len(s192902), 404))
mism = (s192902["tp_label"].str.replace("–", "-") != s192902["window"]
        .str.extract(r"^(≤12|13-20|≥23|Post)")[0].str.replace("Post", "postpartum")).sum()
say(f"  başlıktaki zaman kodu ile anotasyon penceresi uyuşmazlığı: {mism}")
say(f"  aynı pencerede tekrar örnek ('.1' eki): {int(s192902.within_window_repeat.sum())} "
    f"({', '.join(s192902.loc[s192902.within_window_repeat, 'sample_name'])})")
say(f"  birden çok kohortta yer alan kadın: "
    f"{int((s192902.groupby('woman_id').cohort.nunique() > 1).sum())}")
say()
say("  Örnek düzeyi: pencere x tanı")
say(pd.crosstab(s192902.tp_label, s192902.disease).loc[["≤12", "13–20", "≥23", "postpartum"]].to_string())
xt = pd.crosstab(s192902.tp_label, s192902.disease)
for win, exp in [("≤12", (96, 18, 14)), ("13–20", (84, 19, 18))]:
    got = tuple(int(xt.loc[win, d]) for d in ["control", "pre-eclampsia", "severe pre-eclampsia"])
    results.append(check(f"{win} hf örnek (kontrol, PE, şiddetli PE)", got, exp))
say()
early = s192902[s192902.tp_label.isin(["≤12", "13–20"])]
say("  ≤20 hf, kohort x tanı — ÖRNEK düzeyi (devir dosyasındaki tablo bu düzeyde):")
xs = pd.crosstab(early.cohort, early.disease)
say(xs.to_string())
for coh, exp in [("Discovery", (73, 14, 16)), ("Validation 1", (46, 8, 3)), ("Validation 2", (61, 15, 13))]:
    got = tuple(int(xs.loc[coh, d]) for d in ["control", "pre-eclampsia", "severe pre-eclampsia"])
    results.append(check(f"{coh} ≤20 hf örnek", got, exp))
say()
say("  ≤20 hf, kohort x tanı — KADIN düzeyi (kadının en erken örneği; analizde kullanılacak):")
first = early.sort_values(["woman_id", "timepoint", "within_window_repeat"]).drop_duplicates("woman_id")
say(pd.crosstab(first.cohort, first.disease, margins=True, margins_name="Top").to_string())
say()
say("  Kadının en erken örneğinin penceresi (≤20 hf, kohort x pencere):")
say(pd.crosstab(first.cohort, first.tp_label).to_string())
say()

# Sayım dosyaları: sütunlar anotasyonla eşleşiyor mu (mühürlü dosya hariç)
for coh, fname in [("Discovery", "GSE192902_counts_Discovery_postQC.csv.gz"),
                   ("Validation 1", "GSE192902_counts_Validation1_postQC.csv.gz")]:
    cols = pd.read_csv(RAW / "GSE192902" / fname, nrows=0).columns
    samp = pd.Index(cols[2:])
    meta = s192902.loc[s192902.cohort == coh, "sample_name"]
    n_genes = len(pd.read_csv(RAW / "GSE192902" / fname, usecols=[0]))
    say(f"  {fname}: {n_genes} gen x {len(samp)} örnek; anotasyonda bu kohort: {len(meta)}; "
        f"eşleşen: {samp.isin(meta).sum()}; sütunlar {list(cols[:2])}")

# Mühür kontrolü
sealed = ROOT / "data" / "sealed" / "GSE192902_counts_Validation2_postQC.csv.gz"
man = pd.read_csv(ROOT / "logs" / "download_manifest.tsv", sep="\t")
rec = man.loc[man.role == "SEALED_final_test", "sha256"].iloc[0]
h = hashlib.sha256(sealed.read_bytes()).hexdigest()
say(f"  MÜHÜR: Validation-2 dosyası açılmadı; SHA-256 manifest ile "
    f"{'AYNI' if h == rec else 'FARKLI !!!'} ({h[:16]}…)")
say()

# =====================================================================
# GSE149436 — Hafta referansı
# =====================================================================
say("=" * 70)
say("GSE149436 — HAFTA REFERANSI (normal gebelik, GPL28460)")
say("=" * 70)
a = ann("GSE149436")
s149436 = pd.DataFrame({
    "gsm": a.index,
    "woman_id": a["ch:individual"].values,
    "group": a["ch:group"].values,
    "ga_wk": to_num(a["ch:gestational age"]).values,
    "batch": a["ch:batch"].values,
})
s149436.to_csv(OUT / "GSE149436_samples.tsv", sep="\t", index=False)
results.append(check("örnek sayısı", len(s149436), 261))
say(f"  kadın: {s149436.woman_id.nunique()}, grup: {s149436.group.value_counts().to_dict()}, "
    f"hafta {s149436.ga_wk.min():.1f}–{s149436.ga_wk.max():.1f}")
shared_w = set(s149436.woman_id) & set(s149440.woman_id)
say(f"  GSE149440 ile ortak kadın kimliği: {len(shared_w)}")
say()

# =====================================================================
# Seriler arası mükerrer kontrolü (GSM kesişimi)
# =====================================================================
say("=" * 70)
say("SERİLER ARASI MÜKERRER KONTROLÜ (GSM kesişimi)")
say("=" * 70)
s85307 = pd.read_csv(OUT / "GSE85307_samples.tsv", sep="\t")
gsms = {
    "GSE85307_mRNA": set(s85307.loc[s85307.layer == "mRNA", "gsm"]),
    "GSE85307_miRNA": set(s85307.loc[s85307.layer == "miRNA", "gsm"]),
    "GSE149440": set(s149440.gsm), "GSE59491": set(s59491.gsm),
    "GSE192902": set(s192902.gsm), "GSE149436": set(s149436.gsm),
}
names = list(gsms)
mat = pd.DataFrame([[len(gsms[i] & gsms[j]) for j in names] for i in names], index=names, columns=names)
say(mat.to_string())
say()

n_ok, n_all = sum(results), len(results)
say(f"SONUÇ: {n_ok}/{n_all} devir dosyası rakamı doğrulandı.")
(OUT / "cohorts_report.txt").write_text("\n".join(lines), encoding="utf-8")
