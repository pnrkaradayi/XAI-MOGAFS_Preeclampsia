"""
Adım 3 — Ortak dil: prob -> Entrez gen, kohort başına gen düzeyi matris, platformlar arası kesişim.

Kurallar (devir dosyası §4.7):
  * Anahtar Entrez Gene ID (NCBI gene_info ile güncel sembol eklenir). Emekli Entrez kimlikleri
    yalnızca güncel sembol eşsizse kurtarılır.
  * Birden çok gene eşlenen problar atılır; aynı gene düşen problarda MEDYAN alınır.
  * Platformlar arası KESİŞİM (birleşim değil).
  * Yalnız etiketsiz, kohort içi dönüşümler yapılır (log, CPM, dış kohortta quantile). Eğitim verisine
    fit edilen ölçekleme / filtre / seçim Adım 4'te çapraz doğrulamanın İÇİNDE yapılır.
  * Yalnız Adım 2'de seçilen örnekler (kadın başına bir) işlenir; sütun adı = woman_id.
  * GSE192902 Validation-2 (mühürlü) AÇILMAZ.

Kohort başına dönüşüm:
  GSE85307 / GSE142974 / GSE149440 / GSE59491  GEO'daki RMA log2 değerleri olduğu gibi
  GSE108497   arka planı çıkarılmış yoğunluk -> max(x,0)+16 -> log2 -> kohort içi quantile
              (+16 ofseti limma::neqc ile aynı); çip başına saptama p-değerinden gen saptama oranı
  GSE192902   ham sayım -> aynı Entrez'de toplam -> log2(CPM+1)
  GSE85307 miRNA  ham Ct olduğu gibi (normalizasyon çapraz doğrulama içinde); yalnız hsa- problar

Çıktılar:
  data/processed/expr/<kohort>.tsv.gz             gen (Entrez) x kadın
  data/processed/expr/GSE85307_miRNA_Ct.tsv.gz    miRNA x kadın
  data/processed/gene_universe.tsv                gen x kohort varlık tablosu + sembol + tür
  results/genes/gene_level_report.txt
"""
import gzip
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geo_utils import RAW, RESULTS, ROOT, pd, read_expression  # noqa: E402

PLAT = ROOT / "data" / "platforms"
EXPR = ROOT / "data" / "processed" / "expr"
OUTR = RESULTS / "genes"
EXPR.mkdir(parents=True, exist_ok=True)
OUTR.mkdir(parents=True, exist_ok=True)
lines = []


def say(msg=""):
    print(msg, flush=True)
    lines.append(str(msg))


# ------------------------------------------------------------------ gen referansı
gi = pd.read_csv(ROOT / "data" / "reference" / "Homo_sapiens.gene_info.gz", sep="\t",
                 usecols=["GeneID", "Symbol", "dbXrefs", "type_of_gene", "chromosome"], dtype=str)
gi = gi.set_index("GeneID")
current = set(gi.index)
sym_counts = gi.Symbol.value_counts()
sym2id = {s: i for i, s in gi.Symbol.items() if sym_counts[s] == 1}
ens = gi.dbXrefs.str.extractall(r"Ensembl:(ENSG\d+)")[0].reset_index()
ens_counts = ens[0].value_counts()
ens2id = {e: g for g, e in zip(ens.GeneID, ens[0]) if ens_counts[e] == 1}


def resolve(entrez, symbol=None):
    """Güncel Entrez kimliğini döndürür; emekliyse ve sembol eşsizse sembolden kurtarır."""
    if entrez in current:
        return entrez, "current"
    if symbol is not None and symbol in sym2id:
        return sym2id[symbol], "rescued_by_symbol"
    return None, "dropped_retired"


def build_map(ids, entrez, symbols=None, name=""):
    """prob -> Entrez eşlemesi; çok-genli (///) ve boş problar atılır."""
    df = pd.DataFrame({"probe": ids, "entrez": entrez, "symbol": symbols if symbols is not None else None})
    n0 = len(df)
    empty = df.entrez.isna() | df.entrez.astype(str).str.strip().isin(["", "nan", "None", "<NA>"])
    df["entrez"] = df.entrez.astype(str).str.strip()
    multi = df.entrez.str.contains("///")
    df = df[~empty & ~multi].copy()
    df["symbol"] = df.symbol.astype(str).str.split("///").str[0].str.strip() if symbols is not None else None
    res = [resolve(e, s if symbols is not None else None) for e, s in zip(df.entrez, df.symbol)]
    df["entrez"] = [r[0] for r in res]
    df["status"] = [r[1] for r in res]
    st = df.status.value_counts().to_dict()
    df = df[df.entrez.notna()]
    say(f"  {name}: {n0} prob -> boş {int(empty.sum())}, çok-genli {int((multi & ~empty).sum())} atıldı; "
        f"{st}; kalan {len(df)} prob / {df.entrez.nunique()} gen")
    return df.set_index("probe")["entrez"]


def read_annot(path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("ID\t"):
                cols = line.rstrip("\n").split("\t")
                break
        df = pd.read_csv(f, sep="\t", names=cols, usecols=["ID", "Gene symbol", "Gene ID"], dtype=str,
                         comment="!", on_bad_lines="skip")
    return df


def read_table(path):
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("!platform_table_begin"):
                break
        df = pd.read_csv(f, sep="\t", dtype=str, comment="!")
    return df


def collapse(expr, pmap):
    """prob x örnek -> gen x örnek (aynı gene düşen problarda medyan)."""
    sub = expr.loc[expr.index.intersection(pmap.index)]
    return sub.groupby(pmap.loc[sub.index].values).median()


def quantile_normalize(df):
    arr = df.to_numpy(dtype=np.float64)
    order = np.argsort(arr, axis=0)
    ranks = np.empty_like(order)
    np.put_along_axis(ranks, order, np.arange(arr.shape[0])[:, None], axis=0)
    mean_sorted = np.sort(arr, axis=0).mean(axis=1)
    return pd.DataFrame(mean_sorted[ranks], index=df.index, columns=df.columns)


def save(df, name):
    df.index.name = "entrez"
    df.round(5).to_csv(EXPR / f"{name}.tsv.gz", sep="\t", compression="gzip")
    say(f"  -> {name}: {df.shape[0]} gen x {df.shape[1]} kadın; değer {df.min().min():.2f}–{df.max().max():.2f}")


women = pd.read_csv(ROOT / "data" / "processed" / "selection" / "women_all.tsv", sep="\t", dtype=str)


def selected(cohort, set_=None):
    w = women[women.cohort == cohort]
    if set_:
        w = w[w.set == set_]
    return w.set_index("sample_id")["woman_id"]


# ------------------------------------------------------------------ prob eşlemeleri
say("PROB -> ENTREZ EŞLEMELERİ")
a = read_annot(PLAT / "GPL6244.annot.gz")
map6244 = build_map(a["ID"], a["Gene ID"], a["Gene symbol"], "GPL6244 (HuGene 1.0; GSE85307, GSE142974)")
a = read_annot(PLAT / "GPL10558.annot.gz")
map10558 = build_map(a["ID"], a["Gene ID"], a["Gene symbol"], "GPL10558 (HT-12; GSE108497)")
t = read_table(PLAT / "GPL28460_table.txt")
map28460 = build_map(t["ID"], t["ENTREZ_GENE_ID"], None, "GPL28460 (HTA2 ENTREZG; GSE149440)")
t = read_table(PLAT / "GPL18964_table.txt")
map18964 = build_map(t["ID"], t["ENTREZ_GENE_ID"], None, "GPL18964 (HuGene2.1 ENTREZG; GSE59491)")
say()

genes = {}
say("GEN DÜZEYİ MATRİSLER")
# ------------------------------------------------------------------ RMA dizileri
for cohort, path, pmap in [
    ("GSE85307", RAW / "GSE85307" / "GSE85307-GPL6244_series_matrix.txt.gz", map6244),
    ("GSE142974", RAW / "GSE142974" / "GSE142974_series_matrix.txt.gz", map6244),
    ("GSE149440", RAW / "GSE149440" / "GSE149440_series_matrix.txt.gz", map28460),
    ("GSE59491", RAW / "GSE59491" / "GSE59491_series_matrix.txt.gz", map18964),
]:
    sel = selected(cohort)
    expr = read_expression(path)
    missing = set(sel.index) - set(expr.columns)
    assert not missing, f"{cohort}: matriste olmayan seçili örnek {list(missing)[:5]}"
    expr = expr[sel.index].rename(columns=sel)
    say(f"  {cohort}: {expr.shape[0]} prob, {expr.isna().sum().sum()} eksik değer")
    g = collapse(expr, pmap)
    save(g, cohort)
    genes[cohort] = g

# ------------------------------------------------------------------ GSE108497 (Illumina)
sel = selected("GSE108497")
vals, pvals = [], []
for fname in ["GSE108497_normalized_data.txt.gz", "GSE108497_9269325029_J_and_K_normalized.txt.gz"]:
    p = RAW / "GSE108497" / fname
    with gzip.open(p, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
    raw = pd.read_csv(p, sep="\t", index_col=0, header=None, skiprows=1, dtype={0: str})
    raw.columns = range(1, raw.shape[1] + 1)
    vcols = [i for i, h in enumerate(header) if i > 0 and h != "Detection Pval"]
    names = [header[i] for i in vcols]
    v = raw[vcols].astype("float32")
    v.columns = names
    pv = raw[[i + 1 for i in vcols]].astype("float32")
    pv.columns = names
    vals.append(v)
    pvals.append(pv)
v = vals[0].join(vals[1], how="inner")[sel.index]
pv = pvals[0].join(pvals[1], how="inner")[sel.index]
say(f"  GSE108497: {v.shape[0]} ortak prob; seçili {v.shape[1]} çip")
lv = np.log2(v.clip(lower=0) + 16)
qn = quantile_normalize(lv).rename(columns=sel)
g = collapse(qn, map10558)
save(g, "GSE108497")
genes["GSE108497"] = g
det = (pv < 0.01).rename(columns=sel)
det_gene = det.loc[det.index.intersection(map10558.index)].groupby(map10558.loc[det.index.intersection(map10558.index)].values).max()
det_rate = det_gene.mean(axis=1).rename("GSE108497_detect_rate")
say(f"  GSE108497 saptama: örneklerin ≥%50'sinde saptanan (p<0.01) gen {int((det_rate >= .5).sum())} / {len(det_rate)}")

# ------------------------------------------------------------------ GSE192902 (cfRNA sayım)
for set_, fname in [("cfRNA_discovery", "GSE192902_counts_Discovery_postQC.csv.gz"),
                    ("cfRNA_val1", "GSE192902_counts_Validation1_postQC.csv.gz")]:
    sel = selected("GSE192902", set_)
    c = pd.read_csv(RAW / "GSE192902" / fname)
    ent = pd.Series([ens2id.get(e) or sym2id.get(s) for e, s in zip(c.gene_num, c.gene_name)],
                    index=c.index, name="entrez")
    c = pd.concat([c, ent], axis=1)
    n_un = int(pd.isna(c.entrez).sum())
    c = c.dropna(subset=["entrez"])
    counts = c.groupby("entrez")[list(sel.index)].sum()
    cpm = counts / counts.sum(axis=0) * 1e6
    g = np.log2(cpm + 1).rename(columns=sel)
    say(f"  GSE192902 {set_}: {len(c) + n_un} gen satırı, eşlenemeyen {n_un}")
    save(g, f"GSE192902_{set_}")
    genes[f"GSE192902_{set_}"] = g

# ------------------------------------------------------------------ GSE85307 miRNA (ham Ct)
m = read_expression(RAW / "GSE85307" / "GSE85307-GPL19162_series_matrix.txt.gz")
msel = pd.read_csv(ROOT / "data" / "processed" / "selection" / "GSE85307_miRNA_samples.tsv", sep="\t", dtype=str)
msel = msel.set_index("sample_id_miRNA")["woman_id"]
m = m[msel.index].rename(columns=msel)
m = m[m.index.str.startswith("hsa-")]
m.index.name = "mirna_assay"
m.round(4).to_csv(EXPR / "GSE85307_miRNA_Ct.tsv.gz", sep="\t", compression="gzip")
det_m = m.notna().mean(axis=1)
say(f"  -> GSE85307_miRNA_Ct: {m.shape[0]} insan miRNA testi x {m.shape[1]} kadın; "
    f"kadınların ≥%50'sinde saptanan {int((det_m >= .5).sum())}, ≥%80'inde {int((det_m >= .8).sum())}")
say()

# ------------------------------------------------------------------ kesişim
WB = ["GSE85307", "GSE142974", "GSE149440", "GSE59491", "GSE108497"]
uni = pd.DataFrame(index=sorted(set().union(*[set(g.index) for g in genes.values()]), key=int))
uni.index.name = "entrez"
for k, g in genes.items():
    uni[k] = uni.index.isin(g.index)
uni["symbol"] = gi.Symbol.reindex(uni.index).values
uni["type_of_gene"] = gi.type_of_gene.reindex(uni.index).values
uni["chromosome"] = gi.chromosome.reindex(uni.index).values
uni["wholeblood_all"] = uni[WB].all(axis=1)
uni["plasma_discovery_and_wb"] = uni.wholeblood_all & uni["GSE192902_cfRNA_discovery"]
uni = uni.join(det_rate)
uni.to_csv(ROOT / "data" / "processed" / "gene_universe.tsv", sep="\t")

say("KESİŞİM (Entrez gen sayısı)")
for k in genes:
    say(f"  {k:28s} {int(uni[k].sum()):6d}")
pc = uni.type_of_gene == "protein-coding"
say(f"  Tam kan 5 kohort kesişimi: {int(uni.wholeblood_all.sum())} gen "
    f"({int((uni.wholeblood_all & pc).sum())} protein kodlayan)")
say(f"  + plazma cfRNA (Discovery) ile kesişim: {int(uni.plasma_discovery_and_wb.sum())} gen")
say("  Ardışık kesişim (sırayla eklenince kalan):")
acc = None
for k in WB + ["GSE192902_cfRNA_discovery"]:
    acc = uni[k] if acc is None else acc & uni[k]
    say(f"    + {k:26s} {int(acc.sum()):6d}")
say()

# ------------------------------------------------------------------ akıl sağlığı kontrolleri
say("KONTROLLER")
core = uni.index[uni.wholeblood_all]
means = pd.DataFrame({k: genes[k].loc[core].mean(axis=1) for k in WB})
say("  Gen ortalamalarının kohortlar arası Spearman korelasyonu (tam kan kesişimi):")
say(means.corr(method="spearman").round(2).to_string())
pl = genes["GSE192902_cfRNA_discovery"]
both = core.intersection(pl.index)
r = pd.Series(pl.loc[both].mean(axis=1)).corr(means.loc[both, "GSE85307"], method="spearman")
say(f"  Plazma cfRNA ile GSE85307 tam kan gen ortalaması Spearman: {r:.2f} (bölmeler farklı; düşük beklenir)")
for gname in ["HBB", "XIST", "RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY"]:
    gid = sym2id.get(gname)
    parts = []
    for k, g in list(genes.items()):
        if gid in g.index:
            pct = (g.mean(axis=1).rank(pct=True)[gid] * 100)
            parts.append(f"{k} %{pct:.0f}")
    say(f"  {gname:7s} ortalama ifade yüzdeliği: " + ", ".join(parts))
y_ids = [sym2id[s] for s in ["DDX3Y", "KDM5D", "UTY", "EIF1AY"]]
yd = uni.loc[uni.index.intersection(y_ids), "GSE108497_detect_rate"]
say(f"  GSE108497 Y geni saptama oranı (RPS4Y1 hariç): {yd.round(2).to_dict()}")
say("  Yorum: tüm örnekler kadın; HBB en üstte, XIST yüksek, Y genleri tabanda. GSE108497'de RPS4Y1'in"
    " zayıf sinyali Illumina probunun RPS4X ile çapraz hibridizasyonudur (diğer Y genleri saptanmıyor).")

(OUTR / "gene_level_report.txt").write_text("\n".join(lines), encoding="utf-8")
