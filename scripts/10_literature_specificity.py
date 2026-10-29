"""
Adım 5 (ana şekil) — Yayımlanmış PE imzalarının sendrom özgüllüğü.

Katalog: results/literature/published_signatures.tsv (alt ajan; Moufarrej 2022 listesi Extended Data
Tablo 4 görüntüsüyle elle doğrulandı). Yalnız mRNA/lncRNA/kodlamayan gen imzaları (miRNA imzaları ayrı).

ÖN-YAZILI PUANLAMA KURALI:
  * Özellik -> Entrez: Ensembl kimliği verilmişse onunla, değilse güncel sembol, değilse eşsiz eş-ad.
  * Her kohortta imza genleri kohort içinde (seçili tüm kadınlar, etiketsiz) z-ölçeklenir.
  * Skor = Σ yön_j · z_j / n (yayımlanmış yön; up = +1, down = −1). Ağırlık kullanılmaz (platformlar arası
    ölçek farkı nedeniyle tüm imzalara aynı kural).
  * Yön yayımlanmamışsa: cfRNA/plazma imzaları için GSE192902 Discovery'de, tam kan/PBMC imzaları için
    GSE85307'de Welch t işareti. O "yön kalibrasyon" kohortundaki PE AUC'si iyimserdir ve işaretlenir.
  * Yön kalibrasyonunda gen birincil kalibrasyon kohortunda yoksa diğerinden (plazma <-> GSE85307) alınır;
    kullanılan her kalibrasyon kohortu o imza için "bağımsız değil" sayılır.
  * VDAART imzaları (Mirzakhani 2016, Yadama 2020): kaydedilen yönler kendi keşif verileri GSE85307 ile
    uyumsuz (işaret uyumu %52 ve %22; ajan da başlık çelişkisi bildirdi) -> yön GSE85307'de yeniden
    kalibre edilir (kontrol: scripts/10 çıktısı report.txt).
  * Moufarrej 2022 ayrıca yayımlanmış OR'lerle ağırlıklı puanlanır (Σ log(OR)_j z_j / Σ|log(OR)_j|).
  * Kapsam: kohortta bulunan imza özelliği oranı ≥ 0.5 değilse hücre NA.
  * Bağımsızlık: imzanın türetildiği kohort (Hong 2019 -> GSE108497; Tarca 2021 LDA -> GSE149440;
    Mirzakhani/Yadama -> GSE85307; Moufarrej -> GSE192902 Discovery) ve yön kalibrasyon kohortları özetten çıkar.
  * Her hücre: Mann–Whitney iki yönlü p; bağımsız hücreler arasında BH-FDR (q).
  * Karşılaştırmalar: her grup vs aynı kohortun kontrolleri (GSE108497'de SLE kontrolleri);
    AUC + tabakalı bootstrap %95 GA (1000). AUC > 0.5 = imza o grupta PE yönünde "alarm" veriyor.
  * Özet: PE AUC (bağımsız kohortlar ortalaması) ve alarm AUC'si (sPTD, PPROM, FGR/SGA/FD; bağımsız kohortlar
    ortalaması); özgüllük indeksi = PE AUC − alarm AUC.

Çıktılar: results/literature_specificity/{matrix.tsv, summary.tsv, coverage.tsv, report.txt}
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geo_utils import RESULTS, ROOT, pd  # noqa: E402
from panel_utils import boot_auc  # noqa: E402
from scipy import stats  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

OUT = RESULTS / "literature_specificity"
OUT.mkdir(parents=True, exist_ok=True)
lines = []

# kohort, ifade dosyası, bölme, karşılaştırmalar (grup, kontrol arka planı)
COHORTS = [
    ("GSE85307", "GSE85307", "whole_blood", [("PE", None)]),
    ("GSE149440", "GSE149440", "whole_blood", [("PE", None), ("sPTD", None), ("PPROM", None)]),
    ("GSE59491", "GSE59491", "whole_blood", [("sPTD", None), ("PPROM", None)]),
    ("GSE142974", "GSE142974", "whole_blood", [("sPTD", None)]),
    ("GSE108497", "GSE108497", "whole_blood", [("PE", "SLE"), ("FGR_SGA_FD", "SLE")]),
    ("GSE192902_disc", "GSE192902_cfRNA_discovery", "plasma", [("PE", None)]),
    ("GSE192902_val1", "GSE192902_cfRNA_val1", "plasma", [("PE", None)]),
]
SET_OF = {"GSE192902_disc": "cfRNA_discovery", "GSE192902_val1": "cfRNA_val1"}


def say(m=""):
    print(m, flush=True)
    lines.append(str(m))


def main():
    cat = pd.read_csv(RESULTS / "literature" / "published_signatures.tsv", sep="\t", dtype=str)
    cat = cat[~cat.molecule.str.contains("miRNA", case=False, na=False) | cat.molecule.str.contains("mRNA", na=False)]
    cat = cat[~cat.compartment.str.contains("miRNA", case=False, na=False)]
    gi = pd.read_csv(ROOT / "data" / "reference" / "Homo_sapiens.gene_info.gz", sep="\t",
                     usecols=["GeneID", "Symbol", "Synonyms", "dbXrefs"], dtype=str)
    sym2id = gi.drop_duplicates("Symbol", keep=False).set_index("Symbol").GeneID.to_dict()
    syn = gi.assign(s=gi.Synonyms.str.split("|")).explode("s")
    syn = syn[syn.s.notna() & (syn.s != "-")].drop_duplicates("s", keep=False)
    syn2id = syn.set_index("s").GeneID.to_dict()
    ens = gi.dbXrefs.str.extractall(r"Ensembl:(ENSG\d+)")[0].reset_index()
    ens = ens.assign(GeneID=gi.GeneID.values[ens.level_0]).drop_duplicates(0, keep=False)
    ens2id = ens.set_index(0).GeneID.to_dict()

    def to_entrez(feat, fid):
        if isinstance(fid, str) and fid.startswith("ENSG"):
            return ens2id.get(fid.split(".")[0])
        f = str(feat).strip()
        return sym2id.get(f) or syn2id.get(f)

    cat["entrez"] = [to_entrez(f, i) for f, i in zip(cat.feature, cat.feature_id_if_given)]
    cat["sign"] = cat.direction_in_pe.map({"up": 1.0, "down": -1.0})

    w = pd.read_csv(ROOT / "data" / "processed" / "selection" / "women_all.tsv", sep="\t", dtype=str)
    expr, meta = {}, {}
    for key, f, comp, _ in COHORTS:
        E = pd.read_csv(ROOT / "data" / "processed" / "expr" / f"{f}.tsv.gz", sep="\t", dtype={"entrez": str}).set_index("entrez").T
        coh = "GSE192902" if key.startswith("GSE192902") else key
        ww = w[w.cohort == coh]
        if key in SET_OF:
            ww = ww[ww.set == SET_OF[key]]
        meta[key] = ww.set_index("woman_id").loc[E.index]
        Z = (E - E.mean()) / E.std(ddof=1).replace(0, np.nan)
        expr[key] = Z

    def direction_from(key, ids):
        Z, m = expr[key], meta[key]
        pe, ct = (m.group == "PE").values, (m.group == "Control").values
        if key == "GSE108497":
            pe &= (m.background == "SLE").values
        out = {}
        for g in ids:
            if g in Z.columns:
                a, b = Z.loc[pe, g].dropna(), Z.loc[ct, g].dropna()
                out[g] = np.sign(a.mean() - b.mean()) or 1.0
        return out

    RECALIBRATE = {"Mirzakhani2016_VDAART_348", "Yadama2020_VDAART_overlap238"}
    mou = cat[cat.signature_id == "Moufarrej2022_cfRNA_LR18"].copy()
    mou["signature_id"] = "Moufarrej2022_cfRNA_LR18_logOR"
    mou["weight_num"] = np.log(mou.weight.str.replace("OR=", "", regex=False).astype(float))
    cat = pd.concat([cat, mou], ignore_index=True)

    rows, cov_rows = [], []
    for sid, s in cat.groupby("signature_id", sort=False):
        s = s.dropna(subset=["entrez"]).drop_duplicates("entrez")
        n_total = cat[cat.signature_id == sid].shape[0]
        comp = s.compartment.iloc[0] if len(s) else ""
        plasma_sig = "plasma" in comp or "serum" in comp
        calib = []
        if s.sign.isna().all() or sid in RECALIBRATE:
            order = ["GSE192902_disc", "GSE85307"] if plasma_sig else ["GSE85307", "GSE192902_disc"]
            d = {}
            for ck in order:
                got = {g: v for g, v in direction_from(ck, s.entrez.tolist()).items() if g not in d}
                if got:
                    calib.append(ck)
                    d.update(got)
            s = s.assign(sign=s.entrez.map(d))
        if "weight_num" in s and s.weight_num.notna().all():
            s = s.assign(sign=s.weight_num / s.weight_num.abs().mean())
        s = s.dropna(subset=["sign"])
        for key, f, comp_c, contrasts in COHORTS:
            Z, m = expr[key], meta[key]
            have = [g for g in s.entrez if g in Z.columns and Z[g].notna().all()]
            cov = len(have) / n_total if n_total else 0
            cov_rows.append({"signature_id": sid, "cohort": key, "n_features": n_total, "n_mapped": len(s), "n_present": len(have), "coverage": cov})
            if cov < 0.5 or len(have) == 0:
                continue
            sg = s.set_index("entrez").loc[have, "sign"].astype(float).values
            score = (Z[have].values * sg).mean(axis=1)
            for grp, bg in contrasts:
                g1 = (m.group == grp).values
                g0 = (m.group == "Control").values
                if bg:
                    g1, g0 = g1 & (m.background == bg).values, g0 & (m.background == bg).values
                sel = g1 | g0
                yy, ss = g1[sel].astype(int), score[sel]
                lo, hi = boot_auc(yy, ss, n=500)
                pmw = stats.mannwhitneyu(ss[yy == 1], ss[yy == 0]).pvalue
                rows.append({"signature_id": sid, "study": s.study.iloc[0], "sig_compartment": comp,
                             "at_diagnosis": s.at_diagnosis.iloc[0], "cohort": key, "cohort_compartment": comp_c,
                             "group": grp, "n_case": int(g1.sum()), "n_ctrl": int(g0.sum()), "coverage": round(cov, 2),
                             "auc": roc_auc_score(yy, ss), "ci_lo": lo, "ci_hi": hi, "p_mw": pmw,
                             "direction_calibrated_here": key in calib,
                             "direction_source": ("published" if not calib else "calibrated:" + "+".join(calib))
                             + ("+logOR" if sid.endswith("_logOR") else "")})
    M = pd.DataFrame(rows)
    M.to_csv(OUT / "matrix.tsv", sep="\t", index=False)
    pd.DataFrame(cov_rows).to_csv(OUT / "coverage.tsv", sep="\t", index=False)

    # Özet: bağımsız hücreler (yön kalibrasyon kohortu hariç; Mirzakhani/Yadama GSE85307'den türetildiği
    # için GSE85307 onlar için bağımsız değildir)
    def derived(sid):
        if sid.startswith(("Mirzakhani2016", "Yadama2020")):
            return "GSE85307"
        if sid.startswith("Moufarrej2022"):
            return "GSE192902_disc"
        if sid.startswith("Hong2019"):
            return "GSE108497"
        if sid.startswith("Tarca2021_LDA"):
            return "GSE149440"
        return None
    M["independent"] = ~M.direction_calibrated_here & (M.signature_id.map(derived) != M.cohort)
    q = np.full(len(M), np.nan)
    q[M.independent.values] = stats.false_discovery_control(M.loc[M.independent, "p_mw"].values)
    M["q_bh_independent"] = q
    M.to_csv(OUT / "matrix.tsv", sep="	", index=False)
    ind = M[M.independent]
    summ = []
    for sid, g in ind.groupby("signature_id", sort=False):
        pe = g[g.group == "PE"]
        off = g[g.group.isin(["sPTD", "PPROM", "FGR_SGA_FD"])]
        summ.append({"signature_id": sid, "study": g.study.iloc[0], "sig_compartment": g.sig_compartment.iloc[0],
                     "at_diagnosis": g.at_diagnosis.iloc[0], "direction_source": g.direction_source.iloc[0],
                     "pe_auc_mean": pe.auc.mean(), "pe_cells": len(pe),
                     "pe_auc_wholeblood": pe[pe.cohort_compartment == "whole_blood"].auc.mean(),
                     "pe_auc_plasma": pe[pe.cohort_compartment == "plasma"].auc.mean(),
                     "alarm_auc_mean": off.auc.mean(), "alarm_cells": len(off),
                     "alarm_cells_ci_above_05": int((off.ci_lo > 0.5).sum()),
                     "pe_cells_ci_above_05": int((pe.ci_lo > 0.5).sum()),
                     "alarm_cells_q_lt_010": int((off.q_bh_independent < 0.10).sum())})
    S = pd.DataFrame(summ)
    S["specificity_index"] = S.pe_auc_mean - S.alarm_auc_mean
    S.to_csv(OUT / "summary.tsv", sep="\t", index=False)

    say(f"{M.signature_id.nunique()} imza değerlendirildi ({len(M)} hücre; bağımsız {int(M.independent.sum())}).")
    say("Kapsamı < 0.5 olduğu için hiç değerlendirilemeyen imzalar: " +
        ", ".join(sorted(set(cat.signature_id) - set(M.signature_id))))
    say("\nÖZET (bağımsız hücreler): PE AUC ortalaması, alarm AUC'si (sPTD/PPROM/FGR vs kontrol)")
    for _, r in S.sort_values("pe_auc_mean", ascending=False).iterrows():
        say(f"  {r.signature_id:40s} {r.sig_compartment[:14]:14s} PE {r.pe_auc_mean:.2f} (n={r.pe_cells}, "
            f"GA>0.5: {r.pe_cells_ci_above_05}) | alarm {r.alarm_auc_mean:.2f} (n={r.alarm_cells}, GA>0.5: "
            f"{r.alarm_cells_ci_above_05}) | özgüllük {r.specificity_index:+.2f}  [{r.direction_source}]")
    sig = M[(M.independent) & (M.ci_lo > 0.5)]
    say(f"\n%95 GA'sı 0.5'in üzerinde olan bağımsız hücreler ({len(sig)}):")
    for _, r in sig.sort_values(["group", "auc"], ascending=[True, False]).iterrows():
        say(f"  {r.signature_id:40s} {r.cohort:15s} {r.group:10s} AUC {r.auc:.2f} [{r.ci_lo:.2f}–{r.ci_hi:.2f}] n={r.n_case}/{r.n_ctrl}  p={r.p_mw:.4f} q={r.q_bh_independent:.3f}")
    say(f"\nBağımsız hücrelerde FDR q<0.10: {int((M.q_bh_independent < 0.10).sum())}; q<0.20: {int((M.q_bh_independent < 0.20).sum())} "
        f"(toplam {int(M.independent.sum())} test)")
    mo = M[M.signature_id.str.startswith("Moufarrej2022")]
    say("\nMoufarrej 2022 — aynı kohortlarda yeniden puanlama (yayımlanmış: Val-1 0.71, Val-2 0.72; örnek tanımı farklı):")
    for _, r in mo.iterrows():
        say(f"  {r.signature_id:32s} {r.cohort:15s} {r.group:10s} AUC {r.auc:.2f} [{r.ci_lo:.2f}–{r.ci_hi:.2f}] n={r.n_case}/{r.n_ctrl}"
            + ("  (eğitim kohortu)" if r.cohort == "GSE192902_disc" else ""))
    (OUT / "report.txt").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
