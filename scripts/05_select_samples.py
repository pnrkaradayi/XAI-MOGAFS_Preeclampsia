"""
Adım 2 — Ayıklama: kadın başına tek örnek, uyumlaştırılmış gruplar, akış şeması.

ÖN-KAYITLI SEÇİM KURALI (analizden önce sabitlendi):
  Ana seçim   : her kadının, kohortun erken gebelik penceresindeki EN ERKEN örneği.
                  GSE85307 / GSE142974  kayıt örneği (10–18.9 hf, kadın başına tek)
                  GSE149440             gebelik haftası < 20
                  GSE59491              zaman noktası 1 (17–23 hf)
                  GSE108497             zaman noktası 1–2 (< 24 hf); gebe olmayan kontroller (HC_NP) dışlanır
                  GSE192902             ≤12 ve 13–20 hf dilimleri; aynı dilimde '.1' tekrarı ikinci tercihtir
  Duyarlılık  : in_train_window = örnek haftası 10.0 ≤ hf < 19.0 (eğitim penceresi; GSE192902'de
                  kesin hafta yok -> NA).

Uyumlaştırılmış gruplar: PE, Control, sPTD, PPROM, FGR_SGA_FD (GSE108497 'diğer komplikasyon').
GSE108497'de tüm komplikasyonlar SLE hastalarındadır -> background sütunu (SLE / healthy).

GSE192902 Validation-2 için seçim YALNIZCA anotasyondan yapılır; ifade dosyası açılmaz.
Seçimin şimdi yapılması, son testte hangi örneğin kullanılacağını önceden sabitler.

Çıktılar:
  data/processed/selection/women_all.tsv   tüm kohortlar, kadın başına bir satır
  results/flow/flow.tsv                   akış şeması sayıları
  results/flow/table1.tsv                 grup bazında tanımlayıcılar
  results/flow/flow_report.txt
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geo_utils import RESULTS, ROOT, pd  # noqa: E402

INV = RESULTS / "inventory"
SEL = ROOT / "data" / "processed" / "selection"
FLOW = RESULTS / "flow"
SEL.mkdir(parents=True, exist_ok=True)
FLOW.mkdir(parents=True, exist_ok=True)

COLS = ["cohort", "set", "woman_id", "sample_id", "gsm", "group", "subtype", "background",
        "ga_sample_wk", "ga_bin", "in_train_window", "ga_delivery_wk", "age_y", "bmi", "sealed"]
lines, flow = [], []


def say(msg=""):
    print(msg)
    lines.append(str(msg))


def step(cohort, label, n, note=""):
    flow.append({"cohort": cohort, "step": label, "n": n, "note": note})


def train_window(ga):
    return ((ga >= 10) & (ga < 19)).where(ga.notna())


def finish(df, cohort, **const):
    for k, v in const.items():
        df[k] = v
    for c in COLS:
        if c not in df:
            df[c] = pd.NA
    df["cohort"] = cohort
    return df[COLS]


def inv(name):
    return pd.read_csv(INV / name, sep="\t", dtype={"woman_id": str})


out = []

# ---------------------------------------------------------------- GSE85307 (eğitim)
w = inv("GSE85307_women.tsv")
step("GSE85307", "örnek (mRNA / miRNA)", 157 + 159, "miRNA'da 2 kadının teknik tekrarı")
step("GSE85307", "kadın", len(w))
d = pd.DataFrame({
    "woman_id": w.woman_id, "sample_id": w.gsm_mRNA, "gsm": w.gsm_mRNA,
    "sample_id_miRNA": w.gsm_miRNA,
    "group": w.condition.map({"Preeclampsia": "PE", "Normal": "Control"}),
    "ga_sample_wk": w.ga_draw_wk, "ga_delivery_wk": w.ga_delivery_wk, "age_y": w.age_y, "bmi": w.bmi,
})
d["in_train_window"] = train_window(d.ga_sample_wk)
step("GSE85307", "seçilen (kadın başına 1)", len(d))
d85307_mirna = d[["woman_id", "sample_id_miRNA"]]
out.append(finish(d, "GSE85307", set="train", background="none", sealed=False))

# ---------------------------------------------------------------- GSE142974 (özgüllük 0)
s = inv("GSE142974_samples.tsv")
step("GSE142974", "örnek", len(s))
step("GSE142974", "kadın", s.woman_id.nunique())
d = pd.DataFrame({
    "woman_id": s.woman_id, "sample_id": s.gsm, "gsm": s.gsm,
    "group": s.group.map({"sPTB": "sPTD", "Normal": "Control"}), "subtype": s.group.map({"sPTB": "sPTB<32", "Normal": pd.NA}),
    "ga_sample_wk": s.ga_draw_wk, "ga_delivery_wk": s.ga_delivery_wk, "age_y": s.age_y,
})
d["in_train_window"] = train_window(d.ga_sample_wk)
step("GSE142974", "seçilen (kadın başına 1)", len(d))
out.append(finish(d, "GSE142974", set="spec0", background="none", sealed=False))

# ---------------------------------------------------------------- GSE149440 (özgüllük 1)
s = inv("GSE149440_samples.tsv")
step("GSE149440", "örnek", len(s))
step("GSE149440", "kadın", s.woman_id.nunique())
e = s[s.ga_wk < 20].sort_values(["woman_id", "ga_wk"])
n_no = s.woman_id.nunique() - e.woman_id.nunique()
step("GSE149440", "dışlanan: < 20 hf örneği yok", n_no)
e = e.drop_duplicates("woman_id")
step("GSE149440", "seçilen (en erken < 20 hf)", len(e))
d = pd.DataFrame({
    "woman_id": e.woman_id, "sample_id": e.gsm, "gsm": e.gsm,
    "group": e.group.map({"Early_Preeclampsia": "PE", "PPROM": "PPROM", "sPTD": "sPTD", "Control": "Control"}),
    "subtype": e.group.map({"Early_Preeclampsia": "early-onset"}),
    "ga_sample_wk": e.ga_wk, "ga_delivery_wk": e.ga_delivery_wk,
})
d["in_train_window"] = train_window(d.ga_sample_wk)
out.append(finish(d, "GSE149440", set="spec1", background="none", sealed=False))

# ---------------------------------------------------------------- GSE59491 (özgüllük 2)
s = inv("GSE59491_samples.tsv")
step("GSE59491", "örnek", len(s))
step("GSE59491", "kadın", s.woman_id.nunique())
e = s[s.timepoint == 1]
step("GSE59491", "dışlanan: zaman noktası 2 örnekleri (27–33 hf)", int((s.timepoint == 2).sum()), "örnek")
step("GSE59491", "seçilen (zaman noktası 1)", e.woman_id.nunique())
d = pd.DataFrame({
    "woman_id": e.woman_id, "sample_id": e.gsm, "gsm": e.gsm,
    "group": e.group.map({"sPTL": "sPTD", "PPROM": "PPROM", "Term": "Control"}),
    "ga_sample_wk": e.ga_wk, "ga_delivery_wk": e.ga_delivery_wk,
})
d["in_train_window"] = train_window(d.ga_sample_wk)
out.append(finish(d, "GSE59491", set="spec2", background="none", sealed=False))

# ---------------------------------------------------------------- GSE108497 (özgüllük 3)
s = inv("GSE108497_samples.tsv")
step("GSE108497", "örnek", len(s))
step("GSE108497", "kadın", s.woman_id.nunique())
np_w = s.loc[s.group == "HC_NP", "woman_id"].nunique()
step("GSE108497", "dışlanan: gebe olmayan kontrol (HC_NP)", np_w)
s = s[s.group != "HC_NP"]
e = s[s.tp.isin([1, 2])].sort_values(["woman_id", "ga_wk"])
step("GSE108497", "dışlanan: < 24 hf örneği yok", s.woman_id.nunique() - e.woman_id.nunique())
e = e.drop_duplicates("woman_id")
step("GSE108497", "seçilen (en erken < 24 hf)", len(e))
sub = []
for _, r in e.iterrows():
    if r.group == "SLE_P_PE":
        sub.append("PE<36" if r.pe_onset_36 == "Before" else "PE≥36")
    elif r.group == "SLE_P_OC":
        tags = [t for t, c in [("FGR", "iugr"), ("SGA", "sga"), ("FD", "fd"), ("NND", "nnd"), ("PlInsuf", "pl_insuff")]
                if str(r[c]) == "1"]
        sub.append("+".join(tags))
    else:
        sub.append(pd.NA)
d = pd.DataFrame({
    "woman_id": e.woman_id, "sample_id": e.chip_id, "gsm": e.gsm,
    "group": e.group.map({"SLE_P_PE": "PE", "SLE_P_OC": "FGR_SGA_FD", "SLE_P_NC": "Control", "HC_P": "Control"}),
    "subtype": sub, "background": e.group.map(lambda g: "healthy" if g.startswith("HC") else "SLE"),
    "ga_sample_wk": e.ga_wk, "ga_delivery_wk": e.ga_end_wk, "age_y": e.age_y,
})
d["in_train_window"] = train_window(d.ga_sample_wk)
out.append(finish(d, "GSE108497", set="spec3", sealed=False))

# ---------------------------------------------------------------- GSE192902 (bölme + son test)
s = inv("GSE192902_samples.tsv")
set_map = {"Discovery": "cfRNA_discovery", "Validation 1": "cfRNA_val1", "Validation 2": "FINAL_TEST"}
for coh, g in s.groupby("cohort"):
    tag = f"GSE192902 {coh}"
    step(tag, "örnek", len(g))
    step(tag, "kadın", g.woman_id.nunique())
    e = g[g.tp_label.isin(["≤12", "13–20"])].sort_values(["woman_id", "timepoint", "within_window_repeat"])
    step(tag, "dışlanan: ≤20 hf örneği yok", g.woman_id.nunique() - e.woman_id.nunique())
    rep = int(e.duplicated("woman_id").sum())
    e = e.drop_duplicates("woman_id")
    step(tag, "seçilen (en erken ≤20 hf)", len(e), f"{rep} ek erken örnek (sonraki dilim/tekrar) kullanılmadı")
    d = pd.DataFrame({
        "woman_id": e.woman_id, "sample_id": e.sample_name, "gsm": e.gsm,
        "group": e.disease.map({"control": "Control", "pre-eclampsia": "PE", "severe pre-eclampsia": "PE"}),
        "subtype": e.disease.map({"severe pre-eclampsia": "severe", "pre-eclampsia": "non-severe"}),
        "ga_bin": e.tp_label,
    })
    out.append(finish(d, "GSE192902", set=set_map[coh], background="none", sealed=(coh == "Validation 2")))

women = pd.concat(out, ignore_index=True)
assert not women.duplicated(["cohort", "woman_id"]).any()
women.to_csv(SEL / "women_all.tsv", sep="\t", index=False)
d85307_mirna.to_csv(SEL / "GSE85307_miRNA_samples.tsv", sep="\t", index=False)
fl = pd.DataFrame(flow)
fl.to_csv(FLOW / "flow.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- rapor
say("AKIŞ ŞEMASI")
for c, g in fl.groupby("cohort", sort=False):
    say(f"  {c}")
    for _, r in g.iterrows():
        say(f"    {r.step:48s} {r.n:5d}  {r.note}")
say()
say("SEÇİLEN KADINLAR — kohort x grup (ana seçim)")
order = ["PE", "Control", "sPTD", "PPROM", "FGR_SGA_FD"]
lab = (women.cohort + " | " + women.set).rename("kohort | set")
xt = pd.crosstab(lab, women.group).reindex(columns=order, fill_value=0)
xt["Toplam"] = xt.sum(axis=1)
say(xt.to_string())
say()
say("DUYARLILIK — yalnız eğitim penceresi (10–18.9 hf)")
tw = women[women.in_train_window == True]  # noqa: E712
xt2 = pd.crosstab(tw.cohort, tw.group).reindex(columns=order, fill_value=0)
say(xt2.to_string())
say()
say("GSE108497 arka plan x grup")
g8 = women[women.cohort == "GSE108497"]
say(pd.crosstab(g8.background, g8.group).to_string())
say("GSE108497 'diğer komplikasyon' alt tipleri: " + str(g8.loc[g8.group == "FGR_SGA_FD", "subtype"].value_counts().to_dict()))
say()

# Tablo 1: grup bazında tanımlayıcılar
rows = []
for (c, st, g), sub in women.groupby(["cohort", "set", "group"], sort=False):
    row = {"cohort": c, "set": st, "group": g, "n": len(sub)}
    for v in ["ga_sample_wk", "ga_delivery_wk", "age_y", "bmi"]:
        x = pd.to_numeric(sub[v], errors="coerce").dropna()
        row[v] = f"{x.median():.1f} [{x.quantile(.25):.1f}–{x.quantile(.75):.1f}]" if len(x) else ""
    if c == "GSE192902":
        row["ga_sample_wk"] = " / ".join(f"{k}: {v}" for k, v in sub.ga_bin.value_counts().sort_index().items())
    rows.append(row)
t1 = pd.DataFrame(rows)
t1.to_csv(FLOW / "table1.tsv", sep="\t", index=False)
say("TABLO 1 — medyan [ÇAA]")
say(t1.to_string(index=False))
(FLOW / "flow_report.txt").write_text("\n".join(lines), encoding="utf-8")
