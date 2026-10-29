"""
Adım 1 (ek) — Sistematik taramayla eklenen iki serinin envanteri.

  GSE142974  VDAART, GPL6244, 10–18 hf: sPTB (<32 hf) vs normal  -> özgüllük 0
  GSE108497  Hong 2019, GPL10558, SLE + sağlıklı gebelik, boylamsal -> özgüllük 3 (FGR/SGA dahil)

EXPECTED değerleri, tarama sırasında GEO'nun GSM kayıtlarından (acc.cgi) okunan
sayılardır; burada indirilen dosyalardan bağımsız olarak yeniden doğrulanır.

Çıktılar (results/inventory/):
  GSE142974_samples.tsv, GSE108497_samples.tsv, added_report.txt
"""
import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geo_utils import RAW, RESULTS, pd, read_expression, read_sample_annotation, to_num  # noqa: E402

OUT = RESULTS / "inventory"
lines = []


def say(msg=""):
    print(msg)
    lines.append(str(msg))


def check(label, got, expected):
    ok = got == expected
    say(f"  [{'TAMAM' if ok else 'UYUMSUZ'}] {label}: {got} (beklenen {expected})")
    return ok


results = []
w85307 = pd.read_csv(OUT / "GSE85307_women.tsv", sep="\t")

# =====================================================================
# GSE142974 — Özgüllük 0 (aynı kohort, aynı platform, aynı pencere)
# =====================================================================
say("=" * 70)
say("GSE142974 — ÖZGÜLLÜK 0 (VDAART tam kan, GPL6244)")
say("=" * 70)
path = RAW / "GSE142974" / "GSE142974_series_matrix.txt.gz"
a = read_sample_annotation(path)
expr = read_expression(path)
grp = a["ch:pregnancy condition"].map({"Normal Pregnancy": "Normal", "Spontaneous Preterm Births": "sPTB"})
s142974 = pd.DataFrame({
    "gsm": a.index,
    "woman_id": a["ch:recruitment id"].values,
    "group": grp.values,
    "ga_draw_wk": to_num(a["ch:blood draw at gestation week"]).values,
    "ga_delivery_wk": to_num(a["ch:delivery gestation age week"]).values,
    "age_y": to_num(a["ch:age years"]).values,
    "ethnicity": a["ch:ethnicity"].values,
    "site": a["ch:study site"].values,
})
s142974.to_csv(OUT / "GSE142974_samples.tsv", sep="\t", index=False)
say(f"  matris: {expr.shape[0]} prob x {expr.shape[1]} örnek; değer aralığı "
    f"{expr.min().min():.2f}–{expr.max().max():.2f}")
results += [
    check("kadın sayısı", s142974.woman_id.nunique(), 24),
    check("kadın başına tek örnek", bool(not s142974.woman_id.duplicated().any()), True),
    check("sPTB", int((s142974.group == "sPTB").sum()), 8),
    check("normal", int((s142974.group == "Normal").sum()), 16),
    check("GSE85307 ile ortak kadın", int(s142974.woman_id.isin(w85307.woman_id).sum()), 0),
]
say(f"  kan alımı: {s142974.ga_draw_wk.min():.2f}–{s142974.ga_draw_wk.max():.2f} hf "
    f"(GSE85307: 10.00–18.86)")
for g, sub in s142974.groupby("group"):
    say(f"  {g}: doğum haftası {sub.ga_delivery_wk.min():.1f}–{sub.ga_delivery_wk.max():.1f}, "
        f"kan alımı medyan {sub.ga_draw_wk.median():.1f}")
# GSE85307 mRNA matrisiyle aynı prob kümesi ve benzer ölçek mi (aynı işleme hattı)?
ref = read_expression(RAW / "GSE85307" / "GSE85307-GPL6244_series_matrix.txt.gz")
same_probes = set(ref.index) == set(expr.index)
say(f"  GSE85307 ile prob kümesi aynı: {same_probes}; medyan ifade GSE85307 {ref.median().median():.2f} "
    f"/ GSE142974 {expr.median().median():.2f}")
say()

# =====================================================================
# GSE108497 — Özgüllük 3 (SLE + sağlıklı gebelik)
# =====================================================================
say("=" * 70)
say("GSE108497 — ÖZGÜLLÜK 3 (tam kan, GPL10558 Illumina HT-12)")
say("=" * 70)
a = read_sample_annotation(RAW / "GSE108497" / "GSE108497_series_matrix.txt.gz")
grp_tp = a["ch:grp_p_tp"]
s108497 = pd.DataFrame({
    "gsm": a.index,
    "chip_id": a["description"].values,
    "manifest_note": a["description.1"].values,
    "woman_id": a["ch:donor_id"].values,
    "group": grp_tp.str.rsplit("_", n=1).str[0].values,
    "tp": to_num(a["ch:tp"]).values,
    "ga_wk": to_num(a["ch:ga_at_collection"]).values,
    "ga_end_wk": to_num(a["ch:ga_at_end_of_pregnancy"]).values,
    "sle": a["ch:sle"].values,
    "pe": a["ch:pe"].values, "pe_onset_36": a["ch:if_pe_before_or_after_36_weeks"].values,
    "iugr": a["ch:iugr"].values, "sga": a["ch:sga"].values, "fd": a["ch:fd"].values,
    "nnd": a["ch:nnd"].values, "pl_insuff": a["ch:pl_insuff"].values,
    "apl": a["ch:apl"].values, "lac": a["ch:lac"].values,
    "age_y": to_num(a["ch:age"]).values, "race": a["ch:race"].values,
    "batch": a["ch:batch"].values,
})
s108497.to_csv(OUT / "GSE108497_samples.tsv", sep="\t", index=False)

w = s108497.groupby("woman_id").agg(group=("group", "first"), n_groups=("group", "nunique"), n=("gsm", "size"))
results += [check("örnek sayısı", len(s108497), 512), check("kadın sayısı", len(w), 158)]
say(f"  grup etiketi kadın içinde tutarsız: {(w.n_groups > 1).sum()}; kadın başına örnek maks {w.n.max()}")
dup_tp = s108497[s108497.tp.notna()].duplicated(["woman_id", "tp"], keep=False)
say(f"  aynı kadın + aynı zaman noktasında birden çok örnek: {int(dup_tp.sum())}")
gc = w.group.value_counts()
for g, e in [("HC_P", 43), ("SLE_P_NC", 46), ("SLE_P_PE", 24), ("SLE_P_OC", 22), ("HC_NP", 23)]:
    results.append(check(f"{g} (kadın)", int(gc.get(g, 0)), e))
say()

say("  Parti (batch) x grup — örnek sayısı:")
say(pd.crosstab(s108497.group, s108497.batch, margins=True, margins_name="Top").to_string())
say()

early = s108497[s108497.tp.isin([1, 2])].sort_values(["woman_id", "ga_wk"]).drop_duplicates("woman_id")
say("  <24 hf (zaman noktası 1–2), kadının en erken örneği:")
say(pd.crosstab(early.group, early.tp.astype(int), margins=True, margins_name="Top").to_string())
for g, e in [("HC_P", 42), ("SLE_P_NC", 46), ("SLE_P_OC", 21), ("SLE_P_PE", 23)]:
    results.append(check(f"{g} <24 hf örneği olan kadın", int((early.group == g).sum()), e))
say(f"  zaman noktası 1 hafta aralığı: {early[early.tp == 1].ga_wk.min():.1f}–"
    f"{early[early.tp == 1].ga_wk.max():.1f} (medyan {early[early.tp == 1].ga_wk.median():.1f})")
win = early[(early.ga_wk >= 10) & (early.ga_wk < 19)]
say(f"  eğitim penceresi (10–18.9 hf) içinde: {win.group.value_counts().to_dict()}")
say()
say("  Sonlanım bayrakları (<24 hf örneği olan kadınlar, değer=1 sayısı):")
flags = ["pe", "iugr", "sga", "fd", "nnd", "pl_insuff"]
say(early.groupby("group")[flags].apply(lambda d: (d == "1").sum()).to_string())
say(f"  PE başlangıcı (36 hf öncesi/sonrası): "
    f"{early.loc[early.group == 'SLE_P_PE', 'pe_onset_36'].value_counts().to_dict()}")
oc = early[early.group == "SLE_P_OC"]
say(f"  SLE_P_OC içinde FGR veya SGA: {int(((oc.iugr == '1') | (oc.sga == '1')).sum())} / {len(oc)}; "
    f"hiçbir bayrağı olmayan: {int((oc[flags] != '1').all(axis=1).sum())}")
say()

# İfade dosyaları: çip barkodu -> GSM eşlemesi. Ana dosyada 510 çip var; 9269325029_J ve _K
# GEO'da ayrı bir dosyada yayımlanmış.
def read_illumina(path):
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
    idx = [i for i, h in enumerate(header) if i > 0 and h != "Detection Pval"]
    cols = [header[i] for i in idx]
    df = pd.read_csv(path, sep="\t", usecols=[0] + idx, index_col=0, dtype={c: "float32" for c in cols})
    return df


parts = {}
for fname in ["GSE108497_normalized_data.txt.gz", "GSE108497_9269325029_J_and_K_normalized.txt.gz"]:
    parts[fname] = read_illumina(RAW / "GSE108497" / fname)
    say(f"  {fname}: {parts[fname].shape[1]} çip x {parts[fname].shape[0]} prob")
main, jk = parts.values()
say(f"  ortak prob: {len(main.index.intersection(jk.index))} (yalnız ana dosyada "
    f"{len(main.index.difference(jk.index))}, yalnız J/K dosyasında {len(jk.index.difference(main.index))})")
vals = main.join(jk, how="inner")
matched = s108497.chip_id.isin(vals.columns)
results.append(check("anotasyondaki örneklerin ifade verisinde bulunan", int(matched.sum()), 512))
if (~matched).any():
    say(f"  eşleşmeyen GSM: {s108497.loc[~matched, ['gsm', 'chip_id', 'group']].to_dict('records')[:10]}")
say(f"  ifade verisinde olup anotasyonda olmayan çip: {len(set(vals.columns) - set(s108497.chip_id))}")
say(f"  birleşik matris: {vals.shape[0]} prob x {vals.shape[1]} çip; değer aralığı "
    f"{vals.min().min():.1f}–{vals.max().max():.1f}; negatif değer oranı {(vals < 0).mean().mean():.3f}")
med = vals.median()
say(f"  J/K çiplerinin medyan yoğunluğu {med[jk.columns].round(1).tolist()} / diğer çipler medyanı "
    f"{med.drop(jk.columns).median():.1f} (aynı ölçekte mi kontrolü)")
say("  NOT: değerler log ölçeğinde değil (arka plan çıkarılmış ham yoğunluk). Kohort içinde "
    "offset + log2 + normalizasyon gerekecek; dış test seti olduğu için eğitimden bağımsız yapılabilir.")
say(f"  erken manifest notu taşıyan örnek: {int(s108497.manifest_note.str.contains('Early', na=False).sum())}")
say()

# =====================================================================
# Seriler arası mükerrer kontrolü
# =====================================================================
say("=" * 70)
say("MÜKERRER KONTROLÜ (yeni seriler x mevcut seriler, GSM kesişimi)")
say("=" * 70)
others = {
    "GSE85307": set(pd.read_csv(OUT / "GSE85307_samples.tsv", sep="\t").gsm),
    "GSE149440": set(pd.read_csv(OUT / "GSE149440_samples.tsv", sep="\t").gsm),
    "GSE59491": set(pd.read_csv(OUT / "GSE59491_samples.tsv", sep="\t").gsm),
    "GSE192902": set(pd.read_csv(OUT / "GSE192902_samples.tsv", sep="\t").gsm),
}
for new, gs in [("GSE142974", set(s142974.gsm)), ("GSE108497", set(s108497.gsm))]:
    say(f"  {new}: " + ", ".join(f"{k} {len(gs & v)}" for k, v in others.items()))
say()

say(f"SONUÇ: {sum(results)}/{len(results)} kontrol tuttu.")
(OUT / "added_report.txt").write_text("\n".join(lines), encoding="utf-8")
