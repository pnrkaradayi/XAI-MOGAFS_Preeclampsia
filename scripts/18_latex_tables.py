"""
Makale tabloları (paper/tables/*.tex) — tüm sayılar kaydedilmiş sonuç dosyalarından okunur; elle sayı girilmez
(yalnız literatürde BİLDİRİLEN değerler results/literature_specificity/reported_performance.tsv kaynaklıdır).
"""
import json
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geo_utils import RESULTS, ROOT, pd  # noqa: E402

TAB = ROOT / "paper" / "tables"
TAB.mkdir(parents=True, exist_ok=True)
P = RESULTS / "paper"

SIG = {"Moufarrej2022_cfRNA_LR18": ("Moufarrej 2022, 18 genes (sign)", "moufarrej2022"),
       "Moufarrej2022_cfRNA_LR18_logOR": ("Moufarrej 2022, 18 genes (log-OR)", "moufarrej2022"),
       "Rasmussen2022_cfRNA_7gene": ("Rasmussen 2022, 7 genes", "rasmussen2022"),
       "CastilloMarco2025_EOPE_T1": ("Castillo-Marco 2025, EOPE-T1", "castillomarco2025"),
       "CastilloMarco2025_EOPE_T2": ("Castillo-Marco 2025, EOPE-T2", "castillomarco2025"),
       "Hong2019_SLE_PE_P1_58probe": ("Hong 2019, SLE-PE 58 probes", "hong2019"),
       "Hong2019_M5.15_neutrophil_module": ("Hong 2019, neutrophil M5.15", "hong2019"),
       "Mirzakhani2016_VDAART_348": ("Mirzakhani 2016, VDAART 348", "mirzakhani2016"),
       "Yadama2020_VDAART_overlap238": ("Yadama 2020, VDAART 238", "yadama2020scirep"),
       "Chaiworapongsa2013_EOPE_atDx": ("Chaiworapongsa 2013, EOPE", "chaiworapongsa2013"),
       "Chaiworapongsa2013_LOPE_atDx": ("Chaiworapongsa 2013, LOPE", "chaiworapongsa2013"),
       "Tarca2021_untargeted16_wholeblood": ("Tarca 2021, untargeted 16", "tarca2021jmfnm"),
       "Tarca2021_LDA4_H19_FN1_TUBB6_FPR3": ("Tarca 2021, LDA-4", "tarca2021jmfnm"),
       "Tarca2021_EVT20_wholeblood": ("Tarca 2021, EVT 20", "tarca2021jmfnm"),
       "Tarca2021_LDA3_MMP11_SLC6A2_AIF1L": ("Tarca 2021, LDA-3a", "tarca2021jmfnm"),
       "Tarca2021_LDA3_MMP11_SLC6A2_IL18BP": ("Tarca 2021, LDA-3b", "tarca2021jmfnm"),
       "Elovitz2025_placentalHDP_6gene": ("Elovitz 2025, placental HDP 6", "elovitz2025"),
       "Elovitz2025_placental_vs_immune_6gene": ("Elovitz 2025, placental vs immune 6", "elovitz2025"),
       "Elovitz2025_PE_vs_all_top2": ("Elovitz 2025, PE vs all 2", "elovitz2025"),
       "Weiner2022_EOP_NAMPT_APOA1": ("Weiner 2022, NAMPT + APOA1", "weiner2022"),
       "Weiner2026_EOP_NAMPT": ("Weiner 2026, NAMPT", "weiner2026"),
       "Carter2025_inflammatory_cfRNA8": ("Carter 2025, inflammatory 8", "carter2025"),
       "Lekva2021_TERC_PBMC": ("Lekva 2021, TERC", "lekva2021"),
       "Lekva2021_PBMC_pilotarray_top20": ("Lekva 2021, PBMC top-20", "lekva2021"),
       "Myhrer2024_RPPH1_leukocytes": ("Myhrer 2024, RPPH1", "myhrer2024"),
       "Froystad2026_PBMC_up10_22_24wk": ("Fr{\\o}ystad 2026, PBMC 10", "froystad2026")}
COMP = {"whole_blood": "WB", "plasma_cfRNA": "cfRNA", "PBMC/leukocyte": "PBMC", "PBMC/leukocytes": "PBMC"}
COH = {"GSE192902_disc": "GSE192902-D", "GSE192902_val1": "GSE192902-V1"}


def f(x, d=2):
    if x is None or (isinstance(x, float) and (np.isnan(x) or np.isinf(x))):
        return "--"
    return f"{x:.{d}f}"


def ci(a, lo, hi, d=2):
    return "--" if pd.isna(a) else f"{a:.{d}f} ({lo:.{d}f}--{hi:.{d}f})"


def pf(p):
    if pd.isna(p):
        return "--"
    return "$<$0.001" if p < 0.001 else f"{p:.3f}"


def sig(sid, cite=True):
    n, k = SIG.get(sid, (sid.replace("_", " "), None))
    return f"{n}~\\cite{{{k}}}" if (cite and k) else n


def write(name, caption, label, colspec, header, rows, notes=None, size="\\footnotesize", env="tableorg",
          pos="!htbp", full=True, midrules=(), tabcolsep="4pt"):
    colspec = colspec.replace("p{", ">{\\raggedright\\arraybackslash}p{")  # sola yaslı metin sütunları
    out = [f"\\begin{{{env}}}[{pos}]" if pos else f"\\begin{{{env}}}", "\\centering", f"\\caption{{{caption}}}\\label{{{label}}}", f"{size}",
           f"\\setlength{{\\tabcolsep}}{{{tabcolsep}}}"]
    if full:
        mw = "\\textheight" if env.startswith("sideways") else "\\linewidth"
        out.append(f"\\adjustbox{{max width={mw}}}{{\\begin{{tabular}}{{@{{}}{colspec}@{{}}}}")
    else:
        out.append(f"\\begin{{tabular}}{{@{{}}{colspec}@{{}}}}")
    out.append("\\toprule")
    out += [h if h.endswith("\\\\") or h.startswith("\\cmidrule") else h + " \\\\" for h in header]
    out.append("\\midrule")
    for i, r in enumerate(rows):
        if i in midrules:
            out.append("\\midrule")
        out.append(r if isinstance(r, str) else " & ".join(r) + " \\\\")
    out.append("\\botrule")
    out.append("\\end{tabular}}" if full else "\\end{tabular}")
    if notes:
        out.append(f"\\tabnote{{{notes}}}")
    out.append(f"\\end{{{env}}}")
    txt = "\n".join(out) + "\n"
    for a, b in [("≤", "$\\le$"), ("≥", "$\\ge$"), ("–", "--"), ("−", "$-$"), ("×", "$\\times$"), ("±", "$\\pm$")]:
        txt = txt.replace(a, b)
    (TAB / f"{name}.tex").write_text(txt, encoding="utf-8")
    print("  ", name, len(rows), "satır")


# ---------------------------------------------------------------------------------------------- kohort özellikleri
def tab_cohort_chars():
    T = pd.read_csv(RESULTS / "flow/table1.tsv", sep="\t").fillna("")
    role = {"train": "PE training", "spec0": "specificity", "spec1": "development ($f_2$)", "spec2": "specificity",
            "spec3": "specificity", "cfRNA_discovery": "plasma training", "cfRNA_val1": "same-site validation",
            "FINAL_TEST": "sealed external test"}
    rows, prev = [], None
    mids = []
    for _, r in T.iterrows():
        coh = r.cohort if r.cohort != prev else ""
        if r.cohort != prev and prev is not None:
            mids.append(len(rows))
        prev = r.cohort
        g = r.group.replace("FGR_SGA_FD", "FGR/SGA/FD")
        rows.append([coh, role.get(r.set, r.set) if coh else "", g, str(r.n), r.ga_sample_wk.replace("–", "--"),
                     r.ga_delivery_wk.replace("–", "--") or "--", r.age_y.replace("–", "--") or "--", r.bmi.replace("–", "--") or "--"])
    write("tab_cohort_chars", "Characteristics of the women analysed (one sample per woman). Values are median "
          "[interquartile range]. For GSE192902 only the gestational-age bin of sampling is public.",
          "tab:chars", "llllllll",
          ["Cohort & Role & Group & $n$ & GA at sampling (wk) & GA at delivery (wk) & Age (y) & BMI (kg/m$^2$)"],
          rows, notes="GA, gestational age; FGR/SGA/FD, fetal growth restriction, small for gestational age or fetal "
          "death; -- not reported in the public metadata.", size="\\scriptsize", midrules=mids, tabcolsep="3pt")


# ---------------------------------------------------------------------------------------------- tam kan CV
def tab_wb_cv():
    S = pd.read_csv(RESULTS / "panel/outer_cv_summary.tsv", sep="\t")
    name = {"MOGAFS_spec": "XAI-MOGAFS ($f_1$--$f_4$, with $f_2$)", "MOGAFS": "XAI-MOGAFS ($f_1,f_3,f_4$)",
            "LASSO5": "LASSO-5", "TTEST5": "$t$-test-5"}
    rows = []
    for _, r in S.iterrows():
        rows.append([name[r.method], ci(r.pooled_outer_auc, r.ci_lo, r.ci_hi), f"{r.fold_auc_mean:.3f} $\\pm$ {r.fold_auc_sd:.3f}",
                     f(r.n_features_mean, 1), f(r.dev_auc_off_mean_folds, 3), r.final_symbols.replace(";", ", "),
                     f(r.final_dev_auc_off, 3), f(r.final_dev_auc_pe_vs_comp, 3)])
    write("tab_wb_cv", "Whole blood, preeclampsia (PE) vs control in GSE85307: nested cross-validation (5 outer folds) "
          "and final panels. Off-target AUC is complication (sPTD + PPROM) vs control in the development cohort GSE149440.",
          "tab:wbcv", "lccccp{3.6cm}cc",
          ["Method & Pooled outer AUC (95\\% CI) & Fold AUC (mean $\\pm$ SD) & $|S|$ & Off-target AUC (folds) & "
           "Final panel (all 157 women) & Final off-target AUC & Final PE vs other complications"],
          rows, size="\\scriptsize", tabcolsep="3pt")


# ---------------------------------------------------------------------------------------------- literatür özeti
def tab_lit_summary():
    S = pd.read_csv(RESULTS / "literature_specificity/summary.tsv", sep="\t")
    S = S.sort_values("specificity_index", ascending=False)
    rows = []
    for _, r in S.iterrows():
        rows.append([sig(r.signature_id), COMP.get(r.sig_compartment, r.sig_compartment), r.direction_source.split(":")[0],
                     f(r.pe_auc_mean), str(int(r.pe_cells)), f(r.alarm_auc_mean), str(int(r.alarm_cells)),
                     f"{int(r.pe_cells_ci_above_05)}/{int(r.alarm_cells_ci_above_05)}", f"{r.specificity_index:+.3f}"])
    write("tab_lit_summary", "Independent re-test of 26 published early-pregnancy maternal-blood signatures: mean AUC in "
          "preeclampsia cells and in off-target (sPTD, PPROM, FGR) cells, and the specificity index "
          "(PE AUC $-$ off-target AUC). Only cells independent of each signature's development data are counted.",
          "tab:litsum", "lllcccccc",
          ["Signature & Comp. & Direction & PE AUC & $k_{\\mathrm{PE}}$ & Off-target AUC & $k_{\\mathrm{off}}$ & "
           "CI$>$0.5 (PE/off) & Specificity index"],
          rows, notes="Comp., compartment (WB whole blood; cfRNA plasma cell-free RNA; PBMC peripheral blood mononuclear "
          "cells). Direction: published (signs from the paper) or calibrated (signs estimated in a training cohort, which "
          "is then excluded from that signature's independent cells). Across all 228 independent cells no test reached "
          "FDR $q<0.20$ (Benjamini--Hochberg).", size="\\scriptsize", tabcolsep="3pt")


# ---------------------------------------------------------------------------------------------- PE vs her sendrom
def tab_pe_vs(syndrome, cells, label, title, fname):
    M = pd.read_csv(RESULTS / "literature_specificity/matrix.tsv", sep="\t")
    M["indep"] = M.independent.astype(str) == "True"
    pe_cells = [("GSE85307", "PE"), ("GSE149440", "PE"), ("GSE108497", "PE"), ("GSE192902_disc", "PE"), ("GSE192902_val1", "PE")]
    rows, deltas = [], []
    for sid in SIG:
        d = M[M.signature_id == sid]
        if d.empty:
            continue
        pe = d[d.group.eq("PE") & d.indep]
        sy = d[d.group.eq(syndrome) & d.indep]
        cols = []
        for c, g in cells:
            x = d[(d.cohort == c) & (d.group == g)]
            if x.empty or pd.isna(x.auc.iloc[0]):
                cols.append("--")
            else:
                x = x.iloc[0]
                cols.append(ci(x.auc, x.ci_lo, x.ci_hi) + ("" if x.indep else "$^\\dagger$"))
        pe_mean = pe.auc.mean() if len(pe) else np.nan
        sy_mean = sy.auc.mean() if len(sy) else np.nan
        delta = pe_mean - sy_mean
        deltas.append(delta)
        best_pe = pe.loc[pe.auc.idxmax()] if len(pe) else None
        rows.append([sig(sid), f(pe_mean) + (f" ({len(pe)})" if len(pe) else ""),
                     (f"{best_pe.auc:.2f} ({COH.get(best_pe.cohort, best_pe.cohort)})" if best_pe is not None else "--"),
                     *cols, f(sy_mean), "--" if np.isnan(delta) else f"{delta:+.2f}"])
    hdr = " & ".join([f"{COH.get(c, c)}" for c, _ in cells])
    short = syndrome.replace("FGR_SGA_FD", "FGR")
    write(fname, f"Preeclampsia (PE) vs {title}: per-signature AUC (95\\% CI) for {short} vs control in each cohort, "
          f"compared with the same signature's AUC for PE. $\\Delta$ = mean independent PE AUC $-$ mean independent "
          f"{short} AUC; a PE-specific signature should have PE AUC $>$ 0.5, {short} AUC $\\approx$ 0.5 and "
          f"$\\Delta>0$.", label, "l" + "c" * (len(cells) + 4),
          [f"Signature & PE AUC (k) & Best PE cell & {hdr} & Mean {syndrome.replace('FGR_SGA_FD', 'FGR')} & $\\Delta$"],
          rows, notes="$^\\dagger$ Not independent (cohort used to calibrate the signature's directions). "
          f"Median $\\Delta$ across signatures: {np.nanmedian(deltas):+.3f}; signatures with $\\Delta>0.05$: "
          f"{int(np.nansum(np.array(deltas) > 0.05))} of {int(np.sum(~np.isnan(deltas)))}.",
          size="\\scriptsize", tabcolsep="2.5pt")


def tab_pe_vs_gdm():
    R = pd.read_csv(RESULTS / "final_test/final_results.tsv", sep="\t")
    R = R[R.level == "kadın"]
    pe = R[R.cohort == "DelVecchio_PE_GH"].set_index("model")
    gd = R[R.cohort == "DelVecchio_GDM"].set_index("model")
    models = [m for m in pe.index if m.startswith("LIT_") or m in ("RN7SL_fraction",) or m.endswith("_MOGAFS")]
    rows = []
    for m in models:
        a, b = pe.loc[m], gd.loc[m]
        if pd.isna(a.auc):
            continue
        lab = sig(m.replace("LIT_", "")) if m.startswith("LIT_") else (
            "RN7SL* library fraction" if m == "RN7SL_fraction" else f"XAI-MOGAFS lock {m.split('_')[0]}")
        rows.append([lab, a.n_feat_present, ci(a.auc, a.ci_lo, a.ci_hi), ci(b.auc, b.ci_lo, b.ci_hi), f"{a.auc - b.auc:+.2f}"])
    L = pd.read_csv(P / "loso_full_metrics.tsv", sep="\t")
    L = L[L.held_out == "S3_DelVecchio"]
    mids = [len(rows)]
    for m, lab in [("SI-MOGAFS", "XAI-MOGAFS-SI (LOSO)"), ("MOGAFS", "XAI-MOGAFS (LOSO)"), ("LASSO5", "LASSO-5 (LOSO)"),
                   ("TTEST5", "$t$-test-5 (LOSO)"), ("Moufarrej18", "Moufarrej 2022 (LOSO frame)")]:
        a = L[(L.method == m) & L.contrast.str.startswith("PE")].iloc[0]
        b = L[(L.method == m) & L.contrast.str.startswith("GDM")].iloc[0]
        rows.append([lab, "--", ci(a.auc, a.ci_lo, a.ci_hi), ci(b.auc, b.ci_lo, b.ci_hi), f"{a.auc - b.auc:+.2f}"])
    write("tab_pe_vs_gdm", "Preeclampsia/gestational hypertension vs gestational diabetes (GDM) in the sealed Del Vecchio "
          "plasma cohort (GSE154377, first trimester; 8 PE/GH, 7 GDM, 8 controls). A PE-specific score should separate "
          "PE/GH from controls but not GDM from controls.", "tab:gdm", "lcccc",
          ["Score & Features present & PE/GH vs control AUC (95\\% CI) & GDM vs control AUC (95\\% CI) & $\\Delta$"],
          rows, notes="Upper block: locked models and published signatures applied once (pre-registered sealed test). "
          "Lower block: post hoc leave-site-out (LOSO) fold in which Del Vecchio was held out (both trimesters used for "
          "training sites only). With $n\\le 8$ per group, CIs are wide.", size="\\scriptsize", midrules=mids)


# ---------------------------------------------------------------------------------------------- bizim paneller x sendrom
def tab_ours_syndrome():
    R = pd.read_csv(P / "wb_panels_by_syndrome.tsv", sep="\t")
    role = {"eğitim (apparent)": "training", "geliştirme": "development", "bağımsız": "independent",
            "kilitli dış test": "locked test"}
    mname = {"MOGAFS_spec": "XAI-MOGAFS ($f_2$)", "MOGAFS": "XAI-MOGAFS", "LASSO5": "LASSO-5", "TTEST5": "$t$-test-5",
             "XAI-MOGAFS": "XAI-MOGAFS"}
    for panel, fname, title in [("PE paneli", "tab_ours_pe", "PE panels"), ("sPTB paneli", "tab_ours_ptb", "sPTB panels")]:
        d = R[R.panel == panel]
        rows, mids, prev = [], [], None
        for _, r in d.iterrows():
            if prev is not None and r.method != prev:
                mids.append(len(rows))
            rows.append([mname[r.method] if r.method != prev else "", r.cohort, r.group.replace("FGR_SGA_FD", "FGR"),
                         role[r.role], f"{r.n_case}/{r.n_ctrl}", ci(r.auc, r.ci_lo, r.ci_hi), f(r.auprc), pf(r.p_mw)])
            prev = r.method
        what = ("whole-blood PE panels (locked on GSE85307)" if panel == "PE paneli"
                else "whole-blood spontaneous preterm birth (sPTB) panels (locked on GSE149440)")
        write(fname, f"Syndrome-by-syndrome evaluation of our {what}: each panel scored in every cohort and contrasted "
              "case vs control for preeclampsia (PE), spontaneous preterm delivery (sPTD), preterm premature rupture of "
              "membranes (PPROM) and fetal growth restriction (FGR).", f"tab:{fname[4:]}", "lllllccc",
              ["Method & Cohort & Syndrome & Role & Cases/controls & AUC (95\\% CI) & AUPRC & $p$ (Mann--Whitney)"],
              rows, size="\\scriptsize", midrules=mids,
              notes="Role: training = apparent performance on the panel's own training data; development = cohort used "
              "for the $f_2$ objective; locked test = pre-registered external test; independent = never used.")


# ---------------------------------------------------------------------------------------------- plazma varyantları
def tab_variants():
    V = pd.read_csv(P / "val2_full_metrics.tsv", sep="\t").set_index("model")
    F = pd.read_csv(RESULTS / "final_test/final_results.tsv", sep="\t")
    F = F[(F.cohort == "Val2") & (F.level == "kadın")].set_index("model")
    spec = {"V0": ("all", "--", "--", "--", "--"), "V0pc": ("protein-coding", "--", "--", "--", "--"),
            "V1": ("all", "\\checkmark", "--", "--", "--"), "V2": ("all", "\\checkmark", "\\checkmark", "--", "--"),
            "V3": ("all", "\\checkmark", "\\checkmark", "\\checkmark", "--"),
            "V4": ("protein-coding", "\\checkmark", "\\checkmark", "\\checkmark", "--"),
            "V5": ("all", "\\checkmark", "\\checkmark", "--", "--"), "V6": ("protein-coding", "\\checkmark", "\\checkmark", "\\checkmark", "--"),
            "V7": ("all", "\\checkmark", "--", "--", "\\checkmark"), "V8": ("all", "\\checkmark", "--", "--", "\\checkmark")}
    rows = []
    for code in ["V0", "V0pc", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8"]:
        if code in ("V0", "V0pc"):
            key = "all_genes" if code == "V0" else "protein_coding"
            t = (RESULTS / f"plasma/{key}/plasma_report.txt").read_text(encoding="utf-8")
            cv = float(re.search(r"MOGAFS\s+havuzlanmış dış AUC ([\d.]+)", t).group(1))
            v1 = float(re.search(r"^\s+MOGAFS\s+.*?Val-1 AUC ([\d.]+)", t, re.M).group(1))
            lcv = float(re.search(r"LASSO5\s+havuzlanmış dış AUC ([\d.]+)", t).group(1))
        else:
            s = pd.read_csv(RESULTS / f"plasma/variants/{code}/summary.tsv", sep="\t").set_index("method")
            m = s.loc["MOGAFS"]
            cv = m.get("disc_cv_auc_woman", m.get("disc_cv_auc"))
            v1 = m.get("val1_auc_woman", m.get("val1_auc"))
            lk = "LASSO" if "LASSO" in s.index else "LASSO5"
            lcv = s.loc[lk].get("disc_cv_auc_woman", s.loc[lk].get("disc_cv_auc"))
        v2 = V.loc[f"{code}_MOGAFS"]
        lk2 = f"{code}_LASSO" if f"{code}_LASSO" in V.index else f"{code}_LASSO5"
        tol = "0" if code in ("V5", "V6") else "0.03"
        kmax = "18" if code == "V8" else "5"
        u, comp, stab, ga, rep = spec[code]
        name = code + (" (primary)" if code == "V7" else "")
        rows.append([name, u, comp, stab, ga, rep, kmax, tol, f(cv, 3), f(v1, 3), ci(v2.auc, v2.ci_lo, v2.ci_hi), f(v2.auprc),
                     f(F.loc[f"{code}_MOGAFS", "auc_≤12"]), f(F.loc[f"{code}_MOGAFS", "auc_13–20"]), f(lcv, 3), f(V.loc[lk2, "auc"])])
    write("tab_variants", "Plasma cfRNA: every XAI-MOGAFS variant tried before the sealed test, with its development "
          "estimates and its result on the sealed Validation-2 cohort (26 PE / 61 control). All locks were registered "
          "(SHA-256) before Validation-2 was opened; the pre-written rule selected V7 as primary.", "tab:variants",
          "lllllllccccccccc",
          ["Lock & Genes & Comp. & $f_5$ & GA cov. & Rep. & $k_{\\max}$ & tol & Disc. CV & Val-1 & Val-2 AUC (95\\% CI) & "
           "Val-2 AUPRC & Val-2 $\\le$12 wk & Val-2 13--20 wk & LASSO Disc. CV & LASSO Val-2"],
          rows, notes="Comp., library-composition features; $f_5$, bootstrap-stability objective; GA cov., gestational-age "
          "bin covariate; Rep., all early samples per woman used in training (woman-grouped CV); tol, Pareto-selection "
          "tolerance. Disc. CV = Discovery nested cross-validation (woman level); Val-1 = same-site held-out validation.",
          size="\\scriptsize", tabcolsep="2pt", env="sidewaystableorg", pos="")


# ---------------------------------------------------------------------------------------------- Val-2 tüm metrikler
METR = [("auc", "AUC"), ("auprc", "AUPRC"), ("acc", "Acc."), ("bacc", "BAcc"), ("sens", "Sens."), ("spec", "Spec."),
        ("ppv", "PPV"), ("npv", "NPV"), ("f1", "F1"), ("mcc", "MCC"), ("lr_pos", "LR$+$"), ("lr_neg", "LR$-$")]


def tab_val2_metrics():
    V = pd.read_csv(P / "val2_full_metrics.tsv", sep="\t").set_index("model")
    L = pd.read_csv(P / "loso_full_metrics.tsv", sep="\t")
    L = L[(L.held_out == "S2_Val2") & L.contrast.str.startswith("PE")].set_index("method")
    rows = []
    for key, lab in [("V7_MOGAFS", "XAI-MOGAFS V7 (primary)"), ("V7_LASSO", "LASSO-5 (V7 space)"),
                     ("V7_TTEST", "$t$-test-5 (V7 space)"), ("Moufarrej18_logOR", "Moufarrej 2022 (18)")]:
        r = V.loc[key]
        rows.append([lab, "sealed"] + [f(r[m]) for m, _ in METR] + [f"{int(r.tp)}/{int(r.fp)}/{int(r.tn)}/{int(r.fn)}"])
    mids = [len(rows)]
    for key, lab in [("SI-MOGAFS", "XAI-MOGAFS-SI"), ("MOGAFS", "XAI-MOGAFS (pooled)"), ("META-TTEST5", "meta-$t$-test-5"),
                     ("LASSO5", "LASSO-5"), ("TTEST5", "$t$-test-5"), ("Moufarrej18", "Moufarrej 2022, 18 genes")]:
        r = L.loc[key]
        rows.append([lab, "LOSO"] + [f(r[m]) for m, _ in METR] + [f"{int(r.tp)}/{int(r.fp)}/{int(r.tn)}/{int(r.fn)}"])
    write("tab_val2_metrics", "Complete evaluation-metric panel on the same 87 Validation-2 women (26 PE / 61 control). "
          "Classification thresholds were fixed on training data only (Youden index); no threshold was tuned on "
          "Validation-2.", "tab:val2metrics", "ll" + "c" * (len(METR) + 1),
          ["Model & Design & " + " & ".join(m[1] for m in METR) + " & TP/FP/TN/FN"],
          rows, notes="Upper block: pre-registered sealed test (V7 thresholds from Stanford training samples; Moufarrej "
          "threshold from its own Discovery cohort). Lower block: post hoc leave-site-out fold with the Validation-2 "
          "sites held out (thresholds from the two training sites). Acc., accuracy; BAcc, balanced accuracy; PPV/NPV, "
          "positive/negative predictive value; MCC, Matthews correlation coefficient; LR, likelihood ratio. Prevalence "
          "0.30: a constant `control' prediction attains accuracy 0.70 and BAcc 0.50.",
          size="\\scriptsize", tabcolsep="2.2pt", midrules=mids)


def tab_val2_all():
    V = pd.read_csv(P / "val2_full_metrics.tsv", sep="\t").set_index("model")
    F = pd.read_csv(RESULTS / "final_test/final_results.tsv", sep="\t")
    F = F[(F.cohort == "Val2") & (F.level == "kadın")].set_index("model")
    rows = []
    for m in V.index:
        fk = {"Moufarrej18_logOR": "LIT_Moufarrej2022_cfRNA_LR18_logOR"}.get(m, m)
        r, q = V.loc[m], F.loc[fk] if fk in F.index else None
        lab = m.replace("_", " ").replace("MOGAFS", "XAI-MOGAFS").replace("TTEST5", "$t$-test-5").replace("TTEST", "$t$-test") \
            .replace("LASSO5", "LASSO-5").replace("Moufarrej18 logOR", "Moufarrej 2022 (18, log-OR)").replace("RN7SL fraction", "RN7SL* fraction")
        rows.append([lab, q.n_feat_present if q is not None else "--", ci(r.auc, r.ci_lo, r.ci_hi), f(r.auprc),
                     f(q["auc_≤12"]) if q is not None else "--", f(q["auc_13–20"]) if q is not None else "--",
                     f(q.or_ga_adj) if q is not None else "--", pf(q.or_p) if q is not None else "--"])
    for k in ["LIT_Elovitz2025_PE_vs_all_top2", "LIT_Weiner2022_EOP_NAMPT_APOA1", "LIT_Weiner2026_EOP_NAMPT",
              "LIT_Tarca2021_LDA4_H19_FN1_TUBB6_FPR3", "LIT_Myhrer2024_RPPH1_leukocytes", "LIT_Moufarrej2022_cfRNA_LR18"]:
        q = F.loc[k]
        rows.append([sig(k[4:]), q.n_feat_present, ci(q.auc, q.ci_lo, q.ci_hi), "--", f(q["auc_≤12"]), f(q["auc_13–20"]),
                     f(q.or_ga_adj), pf(q.or_p)])
    write("tab_val2_all", "Sealed Validation-2 (woman level, 26 PE / 61 control): every registered lock and every "
          "published signature with sufficient gene coverage, applied once without re-training.", "tab:val2all",
          "llcccccc", ["Model & Features present & AUC (95\\% CI) & AUPRC & AUC $\\le$12 wk & AUC 13--20 wk & "
                       "GA-adjusted OR & $p$"], rows, size="\\scriptsize", tabcolsep="3pt",
          notes="GA-adjusted OR: odds ratio per SD of the score in a logistic model with the gestational-age bin. "
          "DeLong test, V7 vs Moufarrej (sign weights): $\\Delta$AUC $=-0.316$, $p=0.004$.")


# ---------------------------------------------------------------------------------------------- LOSO
def tab_loso():
    L = pd.read_csv(P / "loso_full_metrics.tsv", sep="\t")
    L = L[L.contrast.str.startswith("PE")]
    site = {"S1_Stanford": "Stanford (29/72)", "S2_Val2": "Val-2 sites (26/61)", "S3_DelVecchio": "Del Vecchio (8/8)"}
    rows, mids, prev = [], [], None
    for s in ["S2_Val2", "S1_Stanford", "S3_DelVecchio"]:
        for _, r in L[L.held_out == s].iterrows():
            if prev is not None and s != prev:
                mids.append(len(rows))
            rows.append([site[s] if s != prev else "", {"SI-MOGAFS": "XAI-MOGAFS-SI", "MOGAFS": "XAI-MOGAFS",
                                                        "META-TTEST5": "meta-$t$-test-5", "LASSO5": "LASSO-5",
                                                        "TTEST5": "$t$-test-5", "Moufarrej18": "Moufarrej 2022"}[r.method],
                         ci(r.auc, r.ci_lo, r.ci_hi), f(r.auprc), f(r.acc), f(r.bacc), f(r.sens), f(r.spec), f(r.ppv),
                         f(r.npv), f(r.f1), f(r.mcc)])
            prev = s
    write("tab_loso", "Post hoc leave-site-out (LOSO) evaluation of XAI-MOGAFS-SI and comparators (woman level, PE vs "
          "control). Each row is trained on the other two sites and applied once to the held-out site.", "tab:loso",
          "llcccccccccc", ["Held-out site (PE/control) & Method & AUC (95\\% CI) & AUPRC & Acc. & BAcc & Sens. & Spec. & "
                           "PPV & NPV & F1 & MCC"], rows, size="\\scriptsize", midrules=mids, tabcolsep="2.5pt",
          notes="Moufarrej 2022 uses the published genes and log-odds weights; at Stanford it is evaluated on its own "
          "training cohort and is therefore not independent there.")
    D = pd.read_csv(P / "loso_delong.tsv", sep="\t")
    rows, mids, prev = [], [], None
    for _, r in D.iterrows():
        if prev is not None and r.held_out != prev:
            mids.append(len(rows))
        rows.append([site[r.held_out] if r.held_out != prev else "", r.comparison.replace("SI-MOGAFS", "XAI-MOGAFS-SI")
                     .replace("META-TTEST5", "meta-$t$-test-5").replace("TTEST5", "$t$-test-5").replace("LASSO5", "LASSO-5")
                     .replace("Moufarrej18", "Moufarrej 2022").replace("vs MOGAFS", "vs XAI-MOGAFS"),
                     f(r.auc_1, 3), f(r.auc_2, 3), f"{r['diff']:+.3f}", f(r.z, 2), pf(r.p)])
        prev = r.held_out
    write("tab_loso_delong", "Paired DeLong tests of XAI-MOGAFS-SI against each comparator on the same held-out women.",
          "tab:delong", "llccccc", ["Held-out site & Comparison & AUC$_1$ & AUC$_2$ & $\\Delta$AUC & $z$ & $p$"], rows,
          size="\\scriptsize", midrules=mids)
    Pn = pd.read_csv(RESULTS / "site_invariant/loso_panels.tsv", sep="\t")
    rows, mids, prev = [], [], None
    for _, r in Pn.iterrows():
        if r.method == "Moufarrej18":
            continue
        if prev is not None and r.held_out != prev:
            mids.append(len(rows))
        inner = f(r.inner_auc, 3) if "inner_auc" in r and pd.notna(r.inner_auc) else "--"
        worst = f(r.worst_site_inner, 3) if "worst_site_inner" in r and pd.notna(r.worst_site_inner) else "--"
        rows.append([site[r.held_out].split(" (")[0] if r.held_out != prev else "", {"SI-MOGAFS": "XAI-MOGAFS-SI",
                     "MOGAFS": "XAI-MOGAFS", "META-TTEST5": "meta-$t$-test-5", "LASSO5": "LASSO-5", "TTEST5": "$t$-test-5"}[r.method],
                     str(r.n_features), r.symbols.replace(";", ", "), inner, worst])
        prev = r.held_out
    write("tab_loso_panels", "Panels selected in each leave-site-out fold (protein-coding universe; structural small RNAs "
          "excluded).", "tab:losopanels", "llclcc", ["Held-out site & Method & $|S|$ & Genes & Inner-CV AUC & Worst-site inner AUC"],
          rows, size="\\scriptsize", midrules=mids)


# ---------------------------------------------------------------------------------------------- sPTB
def tab_ptb():
    t = (RESULTS / "ptb/report.txt").read_text(encoding="utf-8")
    dev = {m: (a, r) for m, a, r in re.findall(r"^\s+(XAI-MOGAFS|LASSO5|TTEST5)\s+AUC ([\d.]+) \(([^)]+)\)", t, re.M)}
    ext = {m: (a, lo, hi, s1, s2) for m, a, lo, hi, s1, s2 in
           re.findall(r"^\s+(XAI-MOGAFS|LASSO5|TTEST5)\s+AUC ([\d.]+) \[%95 GA ([\d.]+)–([\d.]+)\] \| sPTD ([\d.]+) \| PPROM ([\d.]+)", t, re.M)}
    lock = {m: g for m, g in re.findall(r"^\s+(XAI-MOGAFS|LASSO5|TTEST5)\s+([A-Z0-9].*\(.*)$", t, re.M)}
    dl = dict(re.findall(r"DeLong XAI-MOGAFS vs (\w+): .*p = ([\d.]+)", t))
    rows = []
    for m, lab in [("XAI-MOGAFS", "XAI-MOGAFS"), ("LASSO5", "LASSO-5"), ("TTEST5", "$t$-test-5")]:
        a, lo, hi, s1, s2 = ext[m]
        rows.append([lab, f"{dev[m][0]} ({dev[m][1].replace(', ', '/')})", f"{a} ({lo}--{hi})", s1, s2,
                     dl.get(m, "--"), lock[m].replace("(", "(").strip()])
    ridge = re.search(r"ridge .*?AUC ([\d.]+)", t).group(1)
    write("tab_ptb", "Pre-registered spontaneous preterm birth (sPTB) experiment in whole blood (17--23 weeks): "
          "development by nested CV in GSE149440 (68 sPTB / 50 term) and a single locked test in GSE59491 "
          "(51 sPTB / 114 term).", "tab:ptb", "lcccccp{4.6cm}",
          ["Method & Dev. nested-CV AUC (repeats) & External AUC (95\\% CI) & sPTD & PPROM & DeLong $p$ vs XAI-MOGAFS & "
           "Locked panel (coefficient)"], rows, size="\\scriptsize", tabcolsep="3pt",
          notes=f"Ridge regression on all 17,621 genes reached a development CV AUC of {ridge}, i.e. no learnable signal.")


# ---------------------------------------------------------------------------------------------- paneller
def tab_panels():
    gi = pd.read_csv(ROOT / "data/processed/gene_universe.tsv", sep="\t", dtype={"entrez": str}).set_index("entrez").symbol.to_dict()
    gi2 = pd.read_csv(ROOT / "data/reference/Homo_sapiens.gene_info.gz", sep="\t", usecols=["GeneID", "Symbol"], dtype=str) \
        .set_index("GeneID").Symbol.to_dict()

    def symb(x):
        return gi.get(x) or gi2.get(x) or x.replace("_", "\\_")

    rows = []
    fp = json.loads((RESULTS / "panel/final_panels.json").read_text(encoding="utf-8"))
    for m, lab in [("MOGAFS_spec", "XAI-MOGAFS ($f_2$)"), ("MOGAFS", "XAI-MOGAFS"), ("LASSO5", "LASSO-5"), ("TTEST5", "$t$-test-5")]:
        md = fp[m]
        rows.append(["Whole blood, PE (GSE85307)", lab, ", ".join(f"{s} ({c:+.2f})" for s, c in zip(md["symbols"], md["coef"]))])
    mids = [len(rows)]
    for key, code in [("all_genes", "V0"), ("protein_coding", "V0pc"), ("variants/V4", "V4"), ("variants/V7", "V7")]:
        L = json.loads((RESULTS / f"plasma/{key}/locked_models.json").read_text(encoding="utf-8"))["models"]
        for m, md in L.items():
            feats = md.get("features") or md["genes"]
            lab = {"MOGAFS": "XAI-MOGAFS", "LASSO": "LASSO-5", "LASSO5": "LASSO-5", "TTEST": "$t$-test-5", "TTEST5": "$t$-test-5"}[m]
            rows.append([f"Plasma cfRNA lock {code}" + (" (primary)" if code == "V7" and m == "MOGAFS" else ""), lab,
                         ", ".join(f"{symb(s)} ({c:+.2f})" for s, c in zip(feats, md["coef"]))])
    mids.append(len(rows))
    L = json.loads((RESULTS / "ptb/locked_models.json").read_text(encoding="utf-8"))
    for m, md in L.items():
        rows.append(["Whole blood, sPTB (GSE149440)", {"XAI-MOGAFS": "XAI-MOGAFS", "LASSO5": "LASSO-5", "TTEST5": "$t$-test-5"}[m],
                     ", ".join(f"{s} ({c:+.2f})" for s, c in zip(md["symbols"], md["coef"]))])
    rows = [[a, b, c.replace("_", "\\_")] for a, b, c in rows]
    write("tab_panels", "Locked panels: genes (or composition features) and standardized logistic-regression coefficients.",
          "tab:panels", "llp{8.2cm}", ["Arm & Method & Features (coefficient)"], rows, size="\\scriptsize", midrules=mids,
          notes="RN7SL*P, RN7SL pseudogene loci (7SL/SRP RNA); comp\\_* features are library-composition fractions. Positive "
          "coefficients increase the predicted probability of the case class.")


# ---------------------------------------------------------------------------------------------- tasarım eşleştirmeli
def tab_design():
    V7 = pd.read_csv(RESULTS / "plasma/variants/V7/summary.tsv", sep="\t").set_index("method").loc["MOGAFS"]
    V0 = (RESULTS / "plasma/all_genes/plasma_report.txt").read_text(encoding="utf-8")
    v0_val1 = float(re.search(r"^\s+MOGAFS\s+.*?Val-1 AUC ([\d.]+)", V0, re.M).group(1))
    D = pd.read_csv(RESULTS / "comparison/design_matched_castillo_style.tsv", sep="\t")
    V = pd.read_csv(P / "val2_full_metrics.tsv", sep="\t").set_index("model")
    L = pd.read_csv(P / "loso_full_metrics.tsv", sep="\t")
    L = L[L.contrast.str.startswith("PE")].set_index(["held_out", "method"])
    dl = pd.read_csv(P / "loso_delong.tsv", sep="\t")
    p_si_mo = dl[(dl.held_out == "S2_Val2") & (dl.comparison == "SI-MOGAFS vs Moufarrej18")].p.iloc[0]
    rw = {r.gen_kumesi: r for _, r in D.iterrows()}
    si_rw = [r for k, r in rw.items() if k.startswith("XAI-MOGAFS-SI")][0]
    mo_rw = [r for k, r in rw.items() if k.startswith("Moufarrej")][0]
    v7_rw = [r for k, r in rw.items() if k.startswith("XAI-MOGAFS V7")][0]
    mo, si, v7 = V.loc["Moufarrej18_logOR"], L.loc[("S2_Val2", "SI-MOGAFS")], V.loc["V7_MOGAFS"]
    rows = [
        ["Training data (apparent)", "Moufarrej 2022, reported~\\cite{moufarrej2022}", "0.99", "XAI-MOGAFS V7", f(V7.train_apparent_auc),
         "4 vs 18 features", "--"],
        ["Internal cross-validation", "Rasmussen 2022, LOOCV, reported~\\cite{rasmussen2022}", "0.82",
         "XAI-MOGAFS V7, nested CV", f(V7.disc_cv_auc_woman), "nested (stricter) vs LOOCV", "--"],
        ["Same-site held-out", "Moufarrej 2022 Val-1, reported", "0.71", "XAI-MOGAFS V7 / V0, Val-1",
         f"{V7.val1_auc_woman:.2f} / {v0_val1:.2f}", "same Stanford Val-1 women", "--"],
        ["External, gene set re-fitted (LOOCV)", "Castillo-Marco 2025 design~\\cite{castillomarco2025}: Moufarrej genes",
         ci(mo_rw.auc, mo_rw.ci_lo, mo_rw.ci_hi), "XAI-MOGAFS-SI genes (5)", ci(si_rw.auc, si_rw.ci_lo, si_rw.ci_hi),
         "same 87 Val-2 women; V7 genes " + f(v7_rw.auc), "--"],
        ["External, locked, leave-site-out", "Moufarrej 2022, published weights", ci(mo.auc, mo.ci_lo, mo.ci_hi),
         "XAI-MOGAFS-SI (post hoc)", ci(si.auc, si.ci_lo, si.ci_hi), f"BAcc {si.bacc:.2f} vs {mo.bacc:.2f}; MCC {si.mcc:.2f} vs {mo.mcc:.2f}",
         pf(p_si_mo)],
        ["External, locked, pre-registered", "Moufarrej 2022, published weights", ci(mo.auc, mo.ci_lo, mo.ci_hi),
         "XAI-MOGAFS V7 (primary)", ci(v7.auc, v7.ci_lo, v7.ci_hi), f"BAcc {v7.bacc:.2f} vs {mo.bacc:.2f}", "0.004"],
        ["External, different laboratory", "Moufarrej 2022 at Del Vecchio, reported (90\\% CI)", "0.74",
         "XAI-MOGAFS-SI / V7", f"{L.loc[('S3_DelVecchio', 'SI-MOGAFS'), 'auc']:.2f} / 0.41",
         f"our re-test of Moufarrej: {L.loc[('S3_DelVecchio', 'Moufarrej18'), 'auc']:.2f}", "--"]]
    # portre düzen: AUC'yi yöntem hücresine kalın yaz
    rows = [[d, f"{lit}: \\textbf{{{la}}}", f"{ours}: \\textbf{{{oa}}}", note, p] for d, lit, la, ours, oa, note, p in rows]
    write("tab_design", "Design-matched comparison with the literature. Each row compares like with like: the same "
          "evaluation design, and where possible the same women. Published values are taken from the original papers; "
          "all other values were computed here. AUC values in bold.", "tab:design", "p{2.7cm}p{3.6cm}p{3.3cm}p{2.9cm}c",
          ["Evaluation design & Literature comparator: AUC & This study: AUC & Same data / note & $p$"], rows,
          size="\\scriptsize", tabcolsep="3pt",
          notes="Pre-registered result: V7 (row 6) is the study's primary endpoint and is reported unchanged. "
          "XAI-MOGAFS-SI was designed after Validation-2 had been opened and is exploratory (rows 4, 5, 7). $p$: paired "
          "DeLong test on the same women.")


# ============================================================================ ANA METİN ÖZET TABLOLARI (Results)
def bold_best(vals, fmt, higher=True):
    """Sütundaki en iyi değeri kalın yaz."""
    arr = np.array([np.nan if v is None else v for v in vals], float)
    best = np.nanmax(arr) if higher else np.nanmin(arr)
    return [("\\textbf{" + fmt(v) + "}") if (v is not None and not np.isnan(v) and np.isclose(v, best)) else fmt(v)
            for v in vals]


def tab_syndrome_summary():
    """Literatür imzaları: PE ve her diğer sendrom için tek satırlık özet (ayrıntı: Ek tablolar)."""
    M = pd.read_csv(RESULTS / "literature_specificity/matrix.tsv", sep="\t")
    M["indep"] = M.independent.astype(str) == "True"
    pe = M[M.group.eq("PE") & M.indep]
    pe_means = pe.groupby("signature_id").auc.mean()
    rows = [["Preeclampsia (target)", "GSE85307, GSE149440, GSE108497, GSE192902-D/V1",
             f"{pe.signature_id.nunique()}", f(pe_means.median()), f"{int((pe.ci_lo > 0.5).sum())} / {len(pe)}", "--", "--"]]
    for syn, lab, cohorts, tab in [("sPTD", "Spontaneous preterm delivery", "GSE149440, GSE59491, GSE142974", "tab:pe_sptd"),
                                   ("PPROM", "PPROM", "GSE149440, GSE59491", "tab:pe_pprom"),
                                   ("FGR_SGA_FD", "Fetal growth restriction", "GSE108497 (SLE background)", "tab:pe_fgr")]:
        sy = M[M.group.eq(syn) & M.indep]
        sy_means = sy.groupby("signature_id").auc.mean()
        common = pe_means.index.intersection(sy_means.index)
        d = (pe_means[common] - sy_means[common])
        rows.append([f"{lab} (Table~\\ref{{{tab}}})", cohorts, f"{sy.signature_id.nunique()}", f(sy_means.median()),
                     f"{int((sy.ci_lo > 0.5).sum())} / {len(sy)}", f"{d.median():+.3f}", f"{int((d > 0.05).sum())} / {len(d)}"])
    R = pd.read_csv(RESULTS / "final_test/final_results.tsv", sep="\t")
    R = R[(R.level == "kadın") & R.model.str.startswith("LIT_")]
    a = R[R.cohort == "DelVecchio_PE_GH"].set_index("model").auc
    b = R[R.cohort == "DelVecchio_GDM"].set_index("model")
    common = a.dropna().index.intersection(b.auc.dropna().index)
    d = a[common] - b.auc[common]
    rows.append(["Gestational diabetes (Table~\\ref{tab:gdm})", "GSE154377 (sealed)", f"{len(common)}",
                 f(b.auc[common].median()), f"{int((b.ci_lo[common] > 0.5).sum())} / {len(common)}", f"{d.median():+.3f}",
                 f"{int((d > 0.05).sum())} / {len(d)}"])
    write("tab_syndrome_summary", "Published early-pregnancy PE signatures re-tested in independent cohorts, summarised by "
          "syndrome. A PE-specific signature should discriminate PE (AUC $>$ 0.5) and stay silent in the other syndromes "
          "($\\Delta>0$). Per-signature results are given in Supplementary Tables~\\ref{tab:pe_sptd}--\\ref{tab:gdm}.", "tab:synsum", "p{3.6cm}p{3.6cm}ccccc",
          ["Syndrome (vs control) & Cohorts & Signatures & Median AUC & Cells with CI $>$ 0.5 & Median $\\Delta$ & "
           "Signatures with $\\Delta>0.05$"], rows, size="\\footnotesize", midrules=[1],
          notes="How to read: AUC 0.5 = no better than chance, 1.0 = perfect. Median AUC = typical performance of the "
          "signatures in that condition (median over signatures of each signature's mean AUC across independent "
          "cohorts). Cells with CI $>$ 0.5 = number of tests in which a signature was clearly better than chance / number "
          "of tests. $\\Delta$ = PE AUC $-$ AUC in the condition; a PE-specific signature has $\\Delta$ clearly above 0 "
          "(for GDM: PE/gestational hypertension vs GDM in the same cohort). No independent test reached FDR $q<0.20$ "
          "(228 tests).")


def tab_spec_methods():
    """Tam kan: XAI-MOGAFS vs LASSO / t-test — özgüllük karşılaştırması (bağımsız kohortlar)."""
    R = pd.read_csv(P / "wb_panels_by_syndrome.tsv", sep="\t")
    R = R[R.panel == "PE paneli"]
    S = pd.read_csv(RESULTS / "panel/outer_cv_summary.tsv", sep="\t").set_index("method")
    off = [("GSE59491", "sPTD"), ("GSE59491", "PPROM"), ("GSE142974", "sPTD"), ("GSE108497", "FGR_SGA_FD")]
    meths = [("MOGAFS_spec", "XAI-MOGAFS (with $f_2$)"), ("MOGAFS", "XAI-MOGAFS (no $f_2$)"), ("LASSO5", "LASSO-5"),
             ("TTEST5", "$t$-test-5")]
    cols = {k: [] for k in ["size", "pe", "fgr", "margin", "mean_off", "max_off", "sig"]}
    for m, _ in meths:
        d = R[R.method == m]
        g = lambda c, s: d[(d.cohort == c) & (d.group == s)].iloc[0]
        dev = [abs(g(c, s).auc - 0.5) for c, s in off]
        cols["size"].append(int(S.loc[m, "final_n"]))
        cols["pe"].append(g("GSE108497", "PE").auc)
        cols["fgr"].append(g("GSE108497", "FGR_SGA_FD").auc)
        cols["margin"].append(g("GSE108497", "PE").auc - g("GSE108497", "FGR_SGA_FD").auc)
        cols["mean_off"].append(float(np.mean(dev)))
        cols["max_off"].append(float(np.max(dev)))
        cols["sig"].append(int(sum(g(c, s).p_mw < 0.05 for c, s in off)))
    fm2 = lambda v: f"{v:.2f}"
    fm3 = lambda v: f"{v:.3f}"
    pe_b = bold_best(cols["pe"], fm2)
    mar_b = bold_best(cols["margin"], lambda v: f"{v:+.2f}")
    mean_b = bold_best(cols["mean_off"], fm3, higher=False)
    max_b = bold_best(cols["max_off"], fm3, higher=False)
    rows = []
    for i, (m, lab) in enumerate(meths):
        rows.append([lab, str(cols["size"][i]), pe_b[i], fm2(cols["fgr"][i]), mar_b[i], mean_b[i], max_b[i],
                     f"{cols['sig'][i]} / 4"])
    write("tab_spec_methods", "Syndrome specificity of locked whole-blood PE panels in independent cohorts. XAI-MOGAFS "
          "panels carry the least off-target signal; the panel trained with the specificity objective $f_2$ shows the "
          "largest separation between PE and FGR in the same disease background. Best value per column in bold.",
          "tab:specmethods", "lccccccc",
          ["Selector & Genes & PE AUC & FGR AUC & PE $-$ FGR & Mean $|\\AUC-0.5|$ & Max $|\\AUC-0.5|$ & Off-target cells"
           " \\\\ & & \\multicolumn{3}{c}{GSE108497 (SLE background)} & \\multicolumn{3}{c}{4 independent off-target cells}"],
          rows, size="\\footnotesize",
          notes="Off-target cells: sPTD and PPROM in GSE59491, sPTD in GSE142974 and FGR in GSE108497; "
          "``off-target cells'' counts cells with Mann--Whitney $p<0.05$ (a false alarm or an inverted signal). "
          "Full cell-by-cell results: Tables~\\ref{tab:ours_pe}--\\ref{tab:ours_ptb}.")


def tab_headtohead():
    """Val-2 kadınlarında (LOSO, Val-2 dışarıda) tüm metrikler + DeLong — tek ana karşılaştırma tablosu."""
    L = pd.read_csv(P / "loso_full_metrics.tsv", sep="\t")
    L = L[(L.held_out == "S2_Val2") & L.contrast.str.startswith("PE")].set_index("method")
    D = pd.read_csv(P / "loso_delong.tsv", sep="\t")
    D = D[D.held_out == "S2_Val2"].set_index("comparison").p
    meths = [("SI-MOGAFS", "XAI-MOGAFS-SI (5 genes)"), ("Moufarrej18", "Moufarrej 2022 (18 genes)"),
             ("MOGAFS", "XAI-MOGAFS, pooled"), ("TTEST5", "$t$-test-5"), ("LASSO5", "LASSO-5"),
             ("META-TTEST5", "meta-$t$-test-5")]
    keys = ["auc", "auprc", "acc", "bacc", "sens", "spec", "ppv", "npv", "f1", "mcc"]
    colvals = {k: [L.loc[m, k] for m, _ in meths] for k in keys}
    colfmt = {k: bold_best(colvals[k], lambda v: f"{v:.2f}") for k in keys}
    rows = []
    for i, (m, lab) in enumerate(meths):
        r = L.loc[m]
        p = "--" if m == "SI-MOGAFS" else pf(D.get(f"SI-MOGAFS vs {m}"))
        auc_cell = colfmt["auc"][i] + f" ({r.ci_lo:.2f}--{r.ci_hi:.2f})"
        rows.append([lab, auc_cell] + [colfmt[k][i] for k in keys[1:]] + [p])
    write("tab_headtohead", "Head-to-head comparison on the same 87 external women (Validation-2 institutions; 26 PE / 61 "
          "control). All models were trained without these women; thresholds were fixed on the training sites. Best value "
          "per column in bold.", "tab:h2h", "lccccccccccc",
          ["Model & AUC (95\\% CI) & AUPRC & Acc. & BAcc & Sens. & Spec. & PPV & NPV & F1 & MCC & $p$ vs SI"],
          rows, size="\\footnotesize", tabcolsep="3pt",
          notes="Post hoc leave-site-out analysis (designed after the sealed test). $p$: paired DeLong test against "
          "XAI-MOGAFS-SI. Acc., accuracy; BAcc, balanced accuracy; PPV/NPV, positive/negative predictive value; MCC, "
          "Matthews correlation coefficient. The other two leave-site-out folds are reported in Table~\\ref{tab:loso}.")


def tab_locks():
    items = [("V0", RESULTS / "plasma/all_genes"), ("V0pc", RESULTS / "plasma/protein_coding")] + \
            [(f"V{i}", RESULTS / f"plasma/variants/V{i}") for i in range(1, 9)] + [("sPTB panels", RESULTS / "ptb")]
    rows = []
    for name, d in items:
        sha = (d / "locked_models.sha256").read_text().split()[0]
        try:
            created = json.loads((d / "locked_models.json").read_text(encoding="utf-8")).get("created_utc", "")[:16].replace("T", " ")
        except Exception:
            created = ""
        rows.append([name + (" (primary)" if name == "V7" else ""), f"\\texttt{{{sha[:32]}}}", created or "--"])
    man = pd.read_csv(ROOT / "logs/download_manifest.tsv", sep="\t")
    for role, lab in [("SEALED_final_test", "Sealed data: GSE192902 Validation-2 counts"),
                      ("SEALED2_final_test", "Sealed data: GSE154377 raw archive")]:
        rows.append([lab, f"\\texttt{{{man.loc[man.role == role, 'sha256'].iloc[0][:32]}}}", "--"])
    opened = (RESULTS / "final_test/OPENED.txt").read_text(encoding="utf-8").splitlines()[0].split(": ", 1)[1]
    write("tab_locks", "Lock registry: SHA-256 digests (first 32 hexadecimal characters) of every locked model file and "
          "of the sealed datasets. The final-test script verified all digests before reading any label.", "tab:locks",
          "llc", ["Item & SHA-256 (prefix) & Created (UTC)"], rows, size="\\scriptsize",
          notes=f"Labels of the sealed data were merged once, at {opened.replace('T', ' ')[:19]} UTC, after all locks above.")


if __name__ == "__main__":
    tab_syndrome_summary()
    tab_spec_methods()
    tab_headtohead()
    tab_locks()
    tab_cohort_chars()
    tab_wb_cv()
    tab_lit_summary()
    tab_pe_vs("sPTD", [("GSE149440", "sPTD"), ("GSE59491", "sPTD"), ("GSE142974", "sPTD")], "tab:pe_sptd",
              "spontaneous preterm delivery (sPTD)", "tab_pe_vs_sptd")
    tab_pe_vs("PPROM", [("GSE149440", "PPROM"), ("GSE59491", "PPROM")], "tab:pe_pprom",
              "preterm premature rupture of membranes (PPROM)", "tab_pe_vs_pprom")
    tab_pe_vs("FGR_SGA_FD", [("GSE108497", "FGR_SGA_FD")], "tab:pe_fgr", "fetal growth restriction (FGR)", "tab_pe_vs_fgr")
    tab_pe_vs_gdm()
    tab_ours_syndrome()
    tab_variants()
    tab_val2_metrics()
    tab_val2_all()
    tab_loso()
    tab_ptb()
    tab_panels()
    tab_design()
