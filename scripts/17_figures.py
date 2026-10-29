"""
Makale şekilleri (paper/figures/*.pdf). Yalnız kaydedilmiş sonuç dosyalarını okur; yeni model eğitmez.
(Tek istisna: SHAP şekli için Val-2-dışarıda LOSO katının SON lojistik modeli, kaydedilen panel genleriyle yeniden
fit edilir — 12_site_invariant.py'deki fit_lr ile aynı, deterministik.)

Renkler: dataviz referans paleti; her varlık (yöntem/grup) tüm şekillerde aynı rengi taşır.
"""
import json
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geo_utils import RESULTS, ROOT, pd  # noqa: E402

import matplotlib as mpl  # noqa: E402
mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from sklearn.metrics import roc_auc_score, roc_curve  # noqa: E402

FIG = ROOT / "paper" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
P = RESULTS / "paper"

INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e6e4df", "#fcfcfb"
NEUTRAL = "#9a9893"
C = {"SI-MOGAFS": "#2a78d6", "Moufarrej": "#eb6834", "MOGAFS": "#1baf7a", "LASSO": "#eda100",
     "TTEST": "#e87ba4", "META": "#008300", "MOGAFS_spec": "#6250d6", "RN7SL": "#e34948"}
G = {"Control": NEUTRAL, "PE": "#2a78d6", "sPTD": "#eb6834", "PPROM": "#1baf7a", "FGR": "#eda100",
     "GDM": "#e87ba4", "CHTN": "#6250d6"}
DIV = LinearSegmentedColormap.from_list("div", ["#184f95", "#3987e5", "#9ec5f4", "#f0efec", "#f4b7b1", "#e34948", "#9e2a2a"])
SEQ = LinearSegmentedColormap.from_list("seq", ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])

mpl.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7.5, "axes.titlesize": 8.5, "axes.labelsize": 7.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6.8, "axes.spines.top": False,
    "axes.spines.right": False, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.color": INK2, "ytick.color": INK2, "axes.labelcolor": INK, "axes.edgecolor": INK2, "text.color": INK,
    "legend.frameon": False, "pdf.fonttype": 42, "savefig.bbox": "tight", "savefig.dpi": 300,
    "axes.titleweight": "bold", "axes.titlelocation": "left", "lines.linewidth": 1.6})
W = 6.6  # tam metin genişliği (inç)


def save(fig, name):
    fig.savefig(FIG / f"{name}.pdf")
    fig.savefig(FIG / f"{name}.png", dpi=200)
    plt.close(fig)
    print("  ", name)


def chance(ax, axis="x"):
    (ax.axvline if axis == "x" else ax.axhline)(0.5, color=INK2, lw=0.7, ls=(0, (3, 2)), zorder=0)


def hgrid(ax, axis="x"):
    ax.grid(axis=axis, color=GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)


# ------------------------------------------------------------------------------------------ 1. kohortlar
def fig_cohorts():
    w = pd.read_csv(ROOT / "data/processed/selection/women_all.tsv", sep="\t", dtype=str)
    lab = {("GSE85307", "train"): "GSE85307 — whole blood\n(PE training)",
           ("GSE149440", "spec1"): "GSE149440 — whole blood\n(development, $f_2$)",
           ("GSE59491", "spec2"): "GSE59491 — whole blood",
           ("GSE142974", "spec0"): "GSE142974 — whole blood",
           ("GSE108497", "spec3"): "GSE108497 — whole blood\n(SLE + healthy)",
           ("GSE192902", "cfRNA_discovery"): "GSE192902 Discovery — plasma",
           ("GSE192902", "cfRNA_val1"): "GSE192902 Val-1 — plasma",
           ("GSE192902", "FINAL_TEST"): "GSE192902 Val-2 — plasma\n(sealed)"}
    w["lab"] = [lab.get((a, b)) for a, b in zip(w.cohort, w.set)]
    w["grp"] = w.group.replace({"FGR_SGA_FD": "FGR"})
    T = pd.crosstab(w.lab, w.grp)
    rn = pd.read_csv(P / "rn7sl_fraction.tsv", sep="\t")
    dv = rn[rn.cohort == "Del Vecchio"].group.value_counts()
    T.loc["GSE154377 Del Vecchio — plasma\n(sealed)"] = dv.reindex(T.columns).fillna(0)
    for g in ["GDM", "CHTN"]:
        if g not in T.columns:
            T[g] = 0
            T.loc["GSE154377 Del Vecchio — plasma\n(sealed)", g] = dv.get(g, 0)
    order = list(lab.values()) + ["GSE154377 Del Vecchio — plasma\n(sealed)"]
    cols = ["Control", "PE", "sPTD", "PPROM", "FGR", "GDM", "CHTN"]
    T = T.loc[order[::-1], cols].fillna(0).astype(int)
    fig, ax = plt.subplots(figsize=(W, 3.6))
    left = np.zeros(len(T))
    for g in cols:
        v = T[g].values
        ax.barh(range(len(T)), v, left=left, color=G[g], edgecolor=SURF, lw=1.0, height=0.68, label=g)
        for i, (l, x) in enumerate(zip(left, v)):
            if x >= 6:
                ax.text(l + x / 2, i, str(x), ha="center", va="center", fontsize=6.3,
                        color="white" if g in ("PE", "PPROM", "CHTN") else INK)
        left += v
    for i, t in enumerate(left):
        ax.text(t + 2, i, f"n = {int(t)}", va="center", fontsize=6.5, color=INK2)
    ax.set_yticks(range(len(T)), T.index)
    ax.set_xlabel("Women (one sample per woman in the primary analysis)")
    ax.legend(ncol=7, loc="lower left", bbox_to_anchor=(0, 1.0), handlelength=1, columnspacing=1)
    hgrid(ax)
    ax.set_xlim(0, left.max() * 1.12)
    save(fig, "fig_cohorts")


# ------------------------------------------------------------------------------------------ 2. tam kan null
def fig_wholeblood():
    S = pd.read_csv(RESULTS / "panel/outer_cv_summary.tsv", sep="\t")
    name = {"MOGAFS_spec": "XAI-MOGAFS (with $f_2$)", "MOGAFS": "XAI-MOGAFS (no $f_2$)", "LASSO5": "LASSO-5",
            "TTEST5": "t-test-5"}
    col = {"MOGAFS_spec": C["MOGAFS_spec"], "MOGAFS": C["MOGAFS"], "LASSO5": C["LASSO"], "TTEST5": C["TTEST"]}
    txt = (RESULTS / "panel/positive_control_report.txt").read_text(encoding="utf-8")
    pc = re.findall(r"^\s+(\S+)\s+(.+?)\s+n=\s*(\d+)/\s*(\d+)\s+genler\s+(\d+)\s+ridge AUC ([\d.]+) ± ([\d.]+)", txt, re.M)
    fig, (a, b) = plt.subplots(1, 2, figsize=(W, 2.5), gridspec_kw={"width_ratios": [1, 1.1], "wspace": 0.55})
    for i, r in enumerate(S.itertuples()):
        a.errorbar(r.pooled_outer_auc, i, xerr=[[r.pooled_outer_auc - r.ci_lo], [r.ci_hi - r.pooled_outer_auc]],
                   fmt="o", ms=5, color=col[r.method], ecolor=col[r.method], elinewidth=1.4, capsize=0)
        a.text(r.ci_hi + 0.015, i, f"{r.pooled_outer_auc:.2f}", va="center", fontsize=6.5, color=INK2)
    ridge = [float(x[5]) for x in pc if x[0] == "GSE85307"][0]
    a.axvline(ridge, color=INK2, lw=0.8, ls=":", zorder=0)
    a.text(ridge + 0.008, len(S) - 0.45, f"ridge, all genes {ridge:.2f}", fontsize=6, color=INK2, ha="left", va="bottom")
    chance(a)
    a.set_yticks(range(len(S)), [name[m] for m in S.method])
    a.set_xlim(0.25, 0.75)
    a.set_ylim(-0.6, len(S) + 0.4)
    a.set_xlabel("Nested-CV AUC, PE vs control (95% CI)")
    a.set_title("a  GSE85307, 47 PE / 110 control")
    hgrid(a)
    labs, vals, sds, cols = [], [], [], []
    for coh, desc, nc, nt, ng, auc, sd in pc:
        plasma = "plazma" in desc
        labs.append(coh.replace("GSE192902_cfRNA_discovery", "GSE192902 Disc.").replace("GSE192902_cfRNA_val1", "GSE192902 Val-1")
                    + f"\n({'plasma' if plasma else 'whole blood'}; {nc}/{nt})")
        vals.append(float(auc))
        sds.append(float(sd))
        cols.append(C["SI-MOGAFS"] if plasma else NEUTRAL)
    y = np.arange(len(vals))[::-1]
    b.barh(y, vals, xerr=sds, color=cols, height=0.6, edgecolor=SURF, error_kw={"elinewidth": 0.8, "ecolor": INK2})
    for yy, v in zip(y, vals):
        b.text(0.305, yy, f"{v:.2f}", va="center", fontsize=6.5, color="white")
    chance(b)
    b.set_xlim(0.3, 0.85)
    b.set_yticks(y, labs)
    b.set_xlabel("Ridge (all genes) CV AUC ± SD, 5×3 CV")
    b.set_title("b  Same pipeline: whole blood vs plasma cfRNA")
    hgrid(b)
    save(fig, "fig_wholeblood")


# ------------------------------------------------------------------------------------------ 3. literatür ısı haritası
CELLS = [("GSE85307", "PE"), ("GSE149440", "PE"), ("GSE108497", "PE"), ("GSE192902_disc", "PE"),
         ("GSE192902_val1", "PE"), ("GSE149440", "sPTD"), ("GSE59491", "sPTD"), ("GSE142974", "sPTD"),
         ("GSE149440", "PPROM"), ("GSE59491", "PPROM"), ("GSE108497", "FGR_SGA_FD")]
SHORT = {"GSE192902_disc": "GSE192902 D", "GSE192902_val1": "GSE192902 V1"}


def sig_label(s):
    s = s.replace("_", " ").replace("2022 cfRNA LR18", " 2022 (18 genes)").replace(" logOR", ", log-OR")
    return s


def fig_lit_heatmap():
    M = pd.read_csv(RESULTS / "literature_specificity/matrix.tsv", sep="\t")
    M["indep"] = M.independent.astype(str) == "True"
    A = M.pivot_table(index="signature_id", columns=["cohort", "group"], values="auc")
    I = M.pivot_table(index="signature_id", columns=["cohort", "group"], values="indep", aggfunc="first")
    A = A.reindex(columns=pd.MultiIndex.from_tuples(CELLS))
    I = I.reindex(columns=pd.MultiIndex.from_tuples(CELLS))
    S = pd.read_csv(RESULTS / "literature_specificity/summary.tsv", sep="\t").set_index("signature_id")
    order = S.sort_values(["sig_compartment", "pe_auc_mean"], ascending=[False, False]).index
    order = [o for o in order if o in A.index]
    A, I = A.loc[order], I.loc[order]
    fig, ax = plt.subplots(figsize=(W, 6.6))
    im = ax.imshow(A.values, cmap=DIV, vmin=0.2, vmax=0.8, aspect="auto")
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            v = A.values[i, j]
            if np.isnan(v):
                ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, fill=False, hatch="////", ec=GRID, lw=0))
                continue
            dep = I.values[i, j] is False or I.values[i, j] == 0
            ax.text(j, i, f"{v:.2f}" + ("*" if dep else ""), ha="center", va="center", fontsize=5.6,
                    color="white" if abs(v - 0.5) > 0.2 else INK, style="italic" if dep else "normal")
    ax.set_xticks(range(len(CELLS)), [f"{SHORT.get(c, c)}\n{g.replace('FGR_SGA_FD', 'FGR')}" for c, g in CELLS],
                  rotation=90, fontsize=6.3)
    ax.set_yticks(range(len(order)), [f"{sig_label(o)}  [{S.loc[o, 'sig_compartment'].replace('_', ' ')}]" for o in order],
                  fontsize=6.0)
    for x in [4.5, 7.5, 9.5]:
        ax.axvline(x, color=SURF, lw=2.5)
    ax.text(2, -1.2, "Preeclampsia", ha="center", fontsize=7, weight="bold")
    ax.text(6, -1.2, "sPTD", ha="center", fontsize=7, weight="bold")
    ax.text(8.5, -1.2, "PPROM", ha="center", fontsize=7, weight="bold")
    ax.text(10, -1.2, "FGR", ha="center", fontsize=7, weight="bold")
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
    cb.set_label("AUC (case vs control); 0.5 = no signal", fontsize=6.5)
    cb.outline.set_visible(False)
    save(fig, "fig_lit_heatmap")


# ------------------------------------------------------------------------------------------ 4. iyimserlik farkı
def fig_optimism():
    S = pd.read_csv(RESULTS / "literature_specificity/summary.tsv", sep="\t").set_index("signature_id")
    # Bildirilen değerler: results/literature_specificity/reported_performance.tsv (PubMed / PMC tam metin)
    rep = [("Moufarrej 2022 (18 genes), Val-2", "Moufarrej2022_cfRNA_LR18", 0.72, "AUC"),
           ("Rasmussen 2022 (7 genes), LOOCV", "Rasmussen2022_cfRNA_7gene", 0.82, "AUC"),
           ("Castillo-Marco 2025 EOPE-T1", "CastilloMarco2025_EOPE_T1", 0.88, "AUC"),
           ("Castillo-Marco 2025 EOPE-T2", "CastilloMarco2025_EOPE_T2", 0.89, "AUC"),
           ("Elovitz 2025 placental HDP", "Elovitz2025_placentalHDP_6gene", 0.88, "AUC"),
           ("Weiner 2022 NAMPT + APOA1", "Weiner2022_EOP_NAMPT_APOA1", 0.89, "AUC"),
           ("Tarca 2021 3-gene LDA (11–17 wk)", "Tarca2021_LDA3_MMP11_SLC6A2_IL18BP", (0.83 + 0.94) / 2, "BAcc"),
           ("Tarca 2021 4-gene LDA (22–28 wk)", "Tarca2021_LDA4_H19_FN1_TUBB6_FPR3", (0.85 + 0.92) / 2, "BAcc")]
    fig, ax = plt.subplots(figsize=(W, 3.1))
    for i, (lab, sid, r, kind) in enumerate(rep[::-1]):
        ind = S.loc[sid, "pe_auc_mean"]
        ax.plot([ind, r], [i, i], color=GRID, lw=3, solid_capstyle="round", zorder=1)
        ax.scatter(r, i, s=36, color=C["Moufarrej"], marker="o" if kind == "AUC" else "D", zorder=3,
                   edgecolor=SURF, lw=0.8)
        ax.scatter(ind, i, s=36, color=C["SI-MOGAFS"], zorder=3, edgecolor=SURF, lw=0.8)
        ax.text(r + 0.012, i, f"{r:.2f}", va="center", fontsize=6.3, color=INK2)
        ax.text(ind - 0.012, i, f"{ind:.2f} (k={int(S.loc[sid, 'pe_cells'])})", va="center", ha="right", fontsize=6.3, color=INK2)
    chance(ax)
    ax.set_yticks(range(len(rep)), [x[0] for x in rep[::-1]])
    ax.set_xlim(0.2, 1.02)
    ax.set_xlabel("Discrimination")
    hgrid(ax)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=C["Moufarrej"], label="Reported in the original paper (AUC)"),
                       Line2D([], [], marker="D", ls="", color=C["Moufarrej"], label="Reported (balanced accuracy)"),
                       Line2D([], [], marker="o", ls="", color=C["SI-MOGAFS"], label="Independent re-test here, mean PE AUC over k cohorts")],
              loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3)
    save(fig, "fig_optimism")


# ------------------------------------------------------------------------------------------ 5. bizim panellerimiz x sendrom
WB = [("GSE85307", "PE"), ("GSE149440", "PE"), ("GSE108497", "PE"), ("GSE149440", "sPTD"), ("GSE59491", "sPTD"),
      ("GSE142974", "sPTD"), ("GSE149440", "PPROM"), ("GSE59491", "PPROM"), ("GSE108497", "FGR_SGA_FD")]


def fig_syndrome_ours():
    R = pd.read_csv(P / "wb_panels_by_syndrome.tsv", sep="\t")
    rows = [("PE paneli", "MOGAFS_spec", "PE panel · XAI-MOGAFS ($f_2$)"), ("PE paneli", "MOGAFS", "PE panel · XAI-MOGAFS"),
            ("PE paneli", "LASSO5", "PE panel · LASSO-5"), ("PE paneli", "TTEST5", "PE panel · t-test-5"),
            ("sPTB paneli", "XAI-MOGAFS", "sPTB panel · XAI-MOGAFS"), ("sPTB paneli", "LASSO5", "sPTB panel · LASSO-5"),
            ("sPTB paneli", "TTEST5", "sPTB panel · t-test-5")]
    A = np.full((len(rows), len(WB)), np.nan)
    role = np.empty(A.shape, dtype=object)
    for i, (p, m, _) in enumerate(rows):
        for j, (c, g) in enumerate(WB):
            r = R[(R.panel == p) & (R.method == m) & (R.cohort == c) & (R.group == g)]
            if len(r):
                A[i, j] = r.auc.iloc[0]
                role[i, j] = r.role.iloc[0]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(W, 3.0), gridspec_kw={"width_ratios": [2.3, 1], "wspace": 0.08})
    im = ax.imshow(A, cmap=DIV, vmin=0.2, vmax=0.8, aspect="auto")
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            rr = role[i, j] or ""
            mark = "ᵗ" if "apparent" in rr else ("ᵈ" if "geli" in rr else ("ᴸ" if "kilitli" in rr else ""))
            ax.text(j, i, f"{A[i, j]:.2f}{mark}", ha="center", va="center", fontsize=6,
                    color="white" if abs(A[i, j] - 0.5) > 0.2 else INK)
    ax.set_xticks(range(len(WB)), [f"{c}\n{g.replace('FGR_SGA_FD', 'FGR')}" for c, g in WB], rotation=90, fontsize=6.2)
    ax.set_yticks(range(len(rows)), [r[2] for r in rows], fontsize=6.4)
    ax.axhline(3.5, color=SURF, lw=3)
    for x in [2.5, 5.5, 7.5]:
        ax.axvline(x, color=SURF, lw=2.5)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title("a  AUC per syndrome (case vs control)", fontsize=7.3)
    # b: hedef-dışı sinyal büyüklüğü — bağımsız hedef-dışı hücreler
    off = [("GSE59491", "sPTD"), ("GSE59491", "PPROM"), ("GSE142974", "sPTD"), ("GSE108497", "FGR_SGA_FD")]
    pe_rows = rows[:4]
    cols = [C["MOGAFS_spec"], C["MOGAFS"], C["LASSO"], C["TTEST"]]
    for i, ((p, m, lab), col) in enumerate(zip(pe_rows, cols)):
        d = np.array([abs(R[(R.panel == p) & (R.method == m) & (R.cohort == c) & (R.group == g)].auc.iloc[0] - 0.5) for c, g in off])
        y = len(pe_rows) - 1 - i
        bx.barh(y, d.mean(), color=col, height=0.6, edgecolor=SURF)
        bx.scatter(d, np.full(len(d), y), s=9, color=INK, zorder=3, lw=0)
        bx.text(0.345, y + 0.34, f"mean {d.mean():.3f}", va="center", ha="right", fontsize=6.0, color=INK2)
    bx.set_yticks(range(4), [r[2].replace("PE panel · ", "") for r in pe_rows[::-1]], fontsize=6.4)
    bx.yaxis.tick_right()
    bx.spines["right"].set_visible(True)
    bx.spines["left"].set_visible(False)
    bx.set_xlim(0, 0.35)
    bx.set_xlabel("|AUC − 0.5| in independent\noff-target cells (bar = mean)")
    bx.set_title("b  Off-target signal", fontsize=7.3)
    hgrid(bx)
    save(fig, "fig_syndrome_ours")


# ------------------------------------------------------------------------------------------ 6. Pareto cepheleri
def fig_pareto():
    fig, axs = plt.subplots(1, 3, figsize=(W, 2.5), gridspec_kw={"wspace": 0.42})
    fp = json.loads((RESULTS / "panel/final_panels.json").read_text(encoding="utf-8"))
    A = pd.read_csv(RESULTS / "panel/archive_final_MOGAFS_spec.tsv", sep="\t", dtype={"genes": str})
    A["k"] = A.idx.astype(str).str.count(";") + 1
    sel = set(fp["MOGAFS_spec"]["genes"])
    A["sel"] = [set(str(g).split(";")) == sel for g in A.genes]
    ax = axs[0]
    nd = A[~A.pareto]
    ax.scatter(nd.auc_off, nd.auc_cv, s=3, color=GRID, lw=0, rasterized=True)
    pf = A[A.pareto]
    sc = ax.scatter(pf.auc_off, pf.auc_cv, s=14, c=pf.k, cmap=SEQ, vmin=0.5, vmax=5, edgecolor=INK2, lw=0.3, zorder=3)
    s = A[A.sel]
    ax.scatter(s.auc_off, s.auc_cv, s=90, marker="*", color=C["MOGAFS_spec"], edgecolor=INK, lw=0.5, zorder=4)
    ax.axvline(0.5, color=INK2, lw=0.6, ls=(0, (3, 2)))
    ax.set_xlabel("Off-target AUC on development cohort\n(sPTD+PPROM vs control; $f_2$ penalises > 0.5)")
    ax.set_ylabel("Inner-CV AUC (PE vs control)")
    ax.set_title("a  Whole blood, with $f_2$")
    for k, (path, lab, key) in enumerate([(RESULTS / "plasma/all_genes/archive_LOCK_MOGAFS.tsv", "b  Plasma, all genes (V0)", "all_genes"),
                                           (RESULTS / "plasma/protein_coding/archive_LOCK_MOGAFS.tsv", "c  Plasma, protein-coding", "protein_coding")]):
        ax = axs[k + 1]
        B = pd.read_csv(path, sep="\t", dtype={"genes": str})
        B["k"] = B.idx.astype(str).str.count(";") + 1
        L = json.loads((RESULTS / f"plasma/{key}/locked_models.json").read_text(encoding="utf-8"))["models"]["MOGAFS"]
        sel = set(L.get("features") or L["genes"])
        B["sel"] = [set(str(g).split(";")) == sel for g in B.genes]
        rng = np.random.default_rng(0)
        jit = rng.uniform(-0.18, 0.18, len(B))
        nd = B[~B.pareto]
        ax.scatter(nd.k + jit[~B.pareto.values], nd.auc_cv, s=3, color=GRID, lw=0, rasterized=True)
        pf = B[B.pareto]
        sc = ax.scatter(pf.k + jit[B.pareto.values], pf.auc_cv, s=14, c=pf.f4, cmap=SEQ, edgecolor=INK2, lw=0.3, zorder=3)
        s = B[B.sel]
        ax.scatter(s.k, s.auc_cv, s=90, marker="*", color=C["MOGAFS"], edgecolor=INK, lw=0.5, zorder=4)
        ax.set_xticks(range(1, 6))
        ax.set_xlabel("Panel size |S| ($f_3$)")
        ax.set_title(lab)
        if k == 1:
            cb = fig.colorbar(sc, ax=ax, fraction=0.06, pad=0.02)
            cb.set_label("$f_4$ (SD of AUC across GA bins)", fontsize=6)
            cb.outline.set_visible(False)
    axs[1].set_ylabel("Inner-CV AUC")
    fig.legend(handles=[Line2D([], [], marker="o", ls="", color=GRID, label="evaluated subset"),
                        Line2D([], [], marker="o", ls="", mfc="#3987e5", mec=INK2, label="non-dominated (archive front)"),
                        Line2D([], [], marker="*", ls="", ms=10, color=INK2, label="selected by pre-registered rule")],
               loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3)
    save(fig, "fig_pareto")


# ------------------------------------------------------------------------------------------ 7. varyant gidişatı
def fig_variants():
    V = pd.read_csv(P / "val2_full_metrics.tsv", sep="\t").set_index("model")
    rows = []
    for key, code in [("all_genes", "V0"), ("protein_coding", "V0pc")]:
        t = (RESULTS / f"plasma/{key}/plasma_report.txt").read_text(encoding="utf-8")
        cv = float(re.search(r"MOGAFS\s+havuzlanmış dış AUC ([\d.]+)", t).group(1))
        v1 = float(re.search(r"^\s+MOGAFS\s+.*?Val-1 AUC ([\d.]+)", t, re.M).group(1))
        rows.append((code, cv, v1))
    for v in ["V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8"]:
        s = pd.read_csv(RESULTS / f"plasma/variants/{v}/summary.tsv", sep="\t").set_index("method").loc["MOGAFS"]
        rows.append((v, s.get("disc_cv_auc_woman", s.get("disc_cv_auc")), s.get("val1_auc_woman", s.get("val1_auc"))))
    fig, ax = plt.subplots(figsize=(W * 0.62, 3.2))
    x = [0, 1, 2]
    labels = []
    for code, cv, v1 in rows:
        v2 = V.loc[f"{code}_MOGAFS", "auc"]
        prim = code == "V7"
        pc = code in ("V0pc", "V4", "V6")
        col = C["MOGAFS"] if prim else ("#1c5cab" if pc else NEUTRAL)
        ax.plot(x, [cv, v1, v2], color=col, lw=2.2 if prim else 0.9, marker="o", ms=4 if prim else 2.5, zorder=3 if prim else 2,
                alpha=1 if (prim or pc) else 0.8)
        labels.append([v2, code + (" (primary)" if prim else ""), col if (prim or pc) else INK2])
    # aynı Val-2 değerli varyantları birleştir, sonra dikey çakışmayı aç
    labels.sort(key=lambda t: t[0])
    merged = []
    for v2, lab, col in labels:
        if merged and abs(merged[-1][0] - v2) < 0.004 and merged[-1][2] == col:
            merged[-1][1] += " = " + lab
        else:
            merged.append([v2, lab, col])
    ypos = []
    for v2, lab, col in merged:
        y = max(v2, ypos[-1] + 0.016) if ypos else v2
        ypos.append(y)
        ax.text(2.06, y, f"{lab} {v2:.2f}", va="center", fontsize=5.8, color=col)
    ax.plot(x, [0.82, 0.64, 0.730], color=C["Moufarrej"], lw=1.6, ls="--", marker="s", ms=3.5)
    ax.text(2.06, 0.745, "Moufarrej 2022\n(published weights)", va="bottom", fontsize=6, color=C["Moufarrej"])
    chance(ax, "y")
    ax.set_xticks(x, ["Discovery\nnested CV", "Validation-1\n(same site)", "Validation-2\n(sealed, external)"])
    ax.set_xlim(-0.15, 2.65)
    ax.set_ylabel("AUC (woman level)")
    hgrid(ax, "y")
    ax.legend(handles=[Line2D([], [], color=C["MOGAFS"], lw=2.2, label="V7, pre-registered primary"),
                       Line2D([], [], color="#1c5cab", lw=0.9, label="protein-coding variants"),
                       Line2D([], [], color=NEUTRAL, lw=0.9, label="other locked variants"),
                       Line2D([], [], color=C["Moufarrej"], ls="--", label="Moufarrej (Discovery = its training set)")],
              loc="lower left", fontsize=6)
    save(fig, "fig_variants")


# ------------------------------------------------------------------------------------------ 8. RN7SL
def fig_rn7sl():
    d = pd.read_csv(P / "rn7sl_fraction.tsv", sep="\t")
    d = d[d.group.isin(["Control", "PE", "GDM"])]
    cohorts = ["Discovery", "Validation 1", "Validation 2", "Del Vecchio"]
    fig, ax = plt.subplots(figsize=(W, 2.6))
    pos, ticks, labs = 0, [], []
    rng = np.random.default_rng(1)
    for coh in cohorts:
        sub = d[d.cohort == coh]
        groups = [g for g in ["Control", "PE", "GDM"] if (sub.group == g).any()]
        start = pos
        for g in groups:
            v = sub[sub.group == g].frac.values * 100
            bp = ax.boxplot(v, positions=[pos], widths=0.55, patch_artist=True, showfliers=False,
                            medianprops={"color": INK, "lw": 1}, whiskerprops={"color": INK2, "lw": 0.7},
                            capprops={"color": INK2, "lw": 0.7}, boxprops={"lw": 0})
            bp["boxes"][0].set_facecolor(G[g])
            bp["boxes"][0].set_alpha(0.35)
            ax.scatter(pos + rng.uniform(-0.17, 0.17, len(v)), v, s=6, color=G[g], lw=0, zorder=3)
            ticks.append(pos)
            labs.append(f"{g}\n(n={len(v)})")
            pos += 1
        y = (sub.group == "PE").astype(int)[sub.group.isin(["Control", "PE"])]
        auc = roc_auc_score(y, sub.frac[sub.group.isin(["Control", "PE"])])
        txt = f"{coh}\nAUC PE vs control {auc:.2f}"
        if "GDM" in groups:
            m = sub.group.isin(["Control", "GDM"])
            txt += f"\nGDM vs control {roc_auc_score((sub.group[m] == 'GDM').astype(int), sub.frac[m]):.2f}"
        ax.text((start + pos - 1) / 2, 30, txt, ha="center", va="top", fontsize=6.3, color=INK)
        pos += 0.8
    ax.set_yscale("log")
    ax.set_ylim(0.005, 40)
    ax.set_xticks(ticks, labs, fontsize=6.3)
    ax.set_ylabel("RN7SL* share of library (%)")
    hgrid(ax, "y")
    save(fig, "fig_rn7sl")


# ------------------------------------------------------------------------------------------ 9. ROC eğrileri
def roc_line(ax, y, s, col, lab, ls="-", lw=1.6):
    f, t, _ = roc_curve(y, s)
    ax.plot(f, t, color=col, lw=lw, ls=ls, label=f"{lab} ({roc_auc_score(y, s):.2f})", drawstyle="steps-post")


def fig_roc():
    S = pd.read_csv(P / "val2_scores.tsv", sep="\t", index_col=0)
    L = pd.read_csv(RESULTS / "site_invariant/loso_scores.tsv", sep="\t")
    fig, (a, b) = plt.subplots(1, 2, figsize=(W, 4.1), gridspec_kw={"wspace": 0.28})
    for col, key, lab, ls in [("V7_MOGAFS", C["MOGAFS"], "V7 XAI-MOGAFS, pre-registered primary", "-"),
                              ("V0pc_MOGAFS", "#1c5cab", "V0pc XAI-MOGAFS (protein-coding lock)", "-"),
                              ("V4_LASSO5", C["LASSO"], "V4 LASSO-5 (protein-coding lock)", "-"),
                              ("RN7SL_fraction", C["RN7SL"], "RN7SL fraction", ":"),
                              ("Moufarrej18_logOR", C["Moufarrej"], "Moufarrej 2022, 18 genes", "--")]:
        roc_line(a, S.y.values, S[col].values, key, lab, ls)
    a.set_title("a  Sealed Validation-2 (26 PE / 61 control)")
    d = L[(L.held_out == "S2_Val2") & L.earliest.astype(bool) & L.group.isin(["Control", "PE"])]
    for meth, key, lab, ls in [("SI-MOGAFS", C["SI-MOGAFS"], "XAI-MOGAFS-SI", "-"), ("MOGAFS", C["MOGAFS"], "XAI-MOGAFS", "-"),
                               ("LASSO5", C["LASSO"], "LASSO-5", "-"), ("TTEST5", C["TTEST"], "t-test-5", "-"),
                               ("META-TTEST5", C["META"], "meta-t-test-5", "-"), ("Moufarrej18", C["Moufarrej"], "Moufarrej 2022", "--")]:
        x = d[d.method == meth]
        roc_line(b, (x.group == "PE").astype(int).values, x.score.values, key, lab, ls, 2.0 if meth == "SI-MOGAFS" else 1.2)
    b.set_title("b  Leave-site-out: Val-2 sites held out (post hoc)")
    for ax in (a, b):
        ax.plot([0, 1], [0, 1], color=INK2, lw=0.6, ls=(0, (3, 2)))
        ax.set_xlabel("1 − specificity")
        ax.set_ylabel("Sensitivity")
        ax.set_aspect("equal")
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), fontsize=6, title="AUC in parentheses", title_fontsize=6)
        hgrid(ax, "both")
    save(fig, "fig_roc")


# ------------------------------------------------------------------------------------------ 10. LOSO
MNAME = {"SI-MOGAFS": "XAI-MOGAFS-SI", "MOGAFS": "XAI-MOGAFS", "META-TTEST5": "meta-t-test-5", "LASSO5": "LASSO-5",
         "TTEST5": "t-test-5", "Moufarrej18": "Moufarrej 2022 (18)"}
MCOL = {"SI-MOGAFS": C["SI-MOGAFS"], "MOGAFS": C["MOGAFS"], "META-TTEST5": C["META"], "LASSO5": C["LASSO"],
        "TTEST5": C["TTEST"], "Moufarrej18": C["Moufarrej"]}


def fig_loso():
    R = pd.read_csv(P / "loso_full_metrics.tsv", sep="\t")
    R = R[R.contrast.str.startswith("PE")]
    sites = [("S2_Val2", "Val-2 sites held out\n(26 PE / 61 control)"), ("S1_Stanford", "Stanford held out\n(29 PE / 72 control)"),
             ("S3_DelVecchio", "Del Vecchio held out\n(8 PE / 8 control)")]
    fig, axs = plt.subplots(1, 3, figsize=(W, 2.5), sharey=True, gridspec_kw={"wspace": 0.08})
    meths = list(MNAME)
    for ax, (s, lab) in zip(axs, sites):
        for i, m in enumerate(meths):
            r = R[(R.held_out == s) & (R.method == m)].iloc[0]
            y = len(meths) - 1 - i
            ax.errorbar(r.auc, y, xerr=[[r.auc - r.ci_lo], [r.ci_hi - r.auc]], fmt="o", ms=4.5, color=MCOL[m],
                        elinewidth=1.3, capsize=0)
            ax.text(0.985, y, f"{r.auc:.2f}", va="center", ha="right", fontsize=6.2, color=INK2)
        chance(ax)
        ax.set_xlim(0, 1)
        ax.set_xticks([0.2, 0.4, 0.6, 0.8])
        ax.set_title(lab, fontsize=7.2)
        ax.set_xlabel("AUC (95% CI)")
        hgrid(ax)
    axs[0].set_yticks(range(len(meths)), [MNAME[m] for m in meths[::-1]])
    save(fig, "fig_loso")


# ------------------------------------------------------------------------------------------ 11. çoklu metrik
METRICS = [("auc", "AUC"), ("auprc", "AUPRC"), ("acc", "Accuracy"), ("bacc", "Balanced acc."), ("sens", "Sensitivity"),
           ("spec", "Specificity"), ("ppv", "PPV"), ("npv", "NPV"), ("f1", "F1"), ("mcc", "MCC")]


def fig_metrics():
    V = pd.read_csv(P / "val2_full_metrics.tsv", sep="\t").set_index("model")
    R = pd.read_csv(P / "loso_full_metrics.tsv", sep="\t")
    R = R[(R.held_out == "S2_Val2") & R.contrast.str.startswith("PE")].set_index("method")
    fig, (a, b) = plt.subplots(1, 2, figsize=(W, 3.2), sharey=True, gridspec_kw={"wspace": 0.06})
    ys = np.arange(len(METRICS))[::-1]
    for ax, sets, title in [(a, [("V7_MOGAFS", V, C["MOGAFS"], "V7 XAI-MOGAFS (pre-registered)"),
                                 ("Moufarrej18_logOR", V, C["Moufarrej"], "Moufarrej 2022")],
                             "a  Sealed test, pre-registered"),
                            (b, [(m, R, MCOL[m], MNAME[m]) for m in ["SI-MOGAFS", "Moufarrej18", "MOGAFS", "LASSO5", "TTEST5"]],
                             "b  Same women, leave-site-out (post hoc)")]:
        n = len(sets)
        for k, (key, T, col, lab) in enumerate(sets):
            off = (k - (n - 1) / 2) * 0.13
            vals = [T.loc[key, m] for m, _ in METRICS]
            ax.scatter(vals, ys + off, s=16 if key in ("SI-MOGAFS", "V7_MOGAFS") else 11, color=col, label=lab, zorder=3,
                       marker="o" if "Mou" not in key else "s", edgecolor=SURF, lw=0.4)
        ax.axvline(0, color=INK2, lw=0.5)
        ax.set_xlim(-0.3, 1.0)
        ax.set_title(title)
        ax.set_xlabel("Value (thresholds fixed on training data)")
        ax.legend(loc="lower left", bbox_to_anchor=(0, 1.06), ncol=2, fontsize=5.9)
        hgrid(ax)
        for y in ys:
            ax.axhline(y, color=GRID, lw=0.4, zorder=0)
    a.set_yticks(ys, [m[1] for m in METRICS])
    save(fig, "fig_metrics")


# ------------------------------------------------------------------------------------------ 12. tasarım eşleştirmeli
def fig_design():
    V7 = pd.read_csv(RESULTS / "plasma/variants/V7/summary.tsv", sep="\t").set_index("method").loc["MOGAFS"]
    D = pd.read_csv(RESULTS / "comparison/design_matched_castillo_style.tsv", sep="\t")
    si_rw = D[D.gen_kumesi.str.startswith("XAI-MOGAFS-SI")].auc.iloc[0]
    mo_rw = D[D.gen_kumesi.str.startswith("Moufarrej")].auc.iloc[0]
    V = pd.read_csv(P / "val2_full_metrics.tsv", sep="\t").set_index("model")
    R = pd.read_csv(P / "loso_full_metrics.tsv", sep="\t")
    R = R[R.contrast.str.startswith("PE")].set_index(["held_out", "method"])
    rows = [("Training data (apparent)", 0.99, "Moufarrej, reported", V7.train_apparent_auc, "V7", "ours"),
            ("Internal cross-validation", 0.82, "Rasmussen, LOOCV, reported", V7.disc_cv_auc_woman, "V7, nested CV", "ours"),
            ("Same-site held-out (Stanford Val-1)", 0.71, "Moufarrej, reported", V7.val1_auc_woman, "V7", "ours"),
            ("External, gene set re-fitted by LOOCV\n(Castillo-Marco design; Val-2)", mo_rw, "Moufarrej genes", si_rw,
             "XAI-MOGAFS-SI genes", "posthoc"),
            ("External, locked, leave-site-out (Val-2)", R.loc[("S2_Val2", "Moufarrej18"), "auc"], "Moufarrej, re-tested",
             R.loc[("S2_Val2", "SI-MOGAFS"), "auc"], "XAI-MOGAFS-SI", "posthoc"),
            ("External, locked, pre-registered (Val-2)", V.loc["Moufarrej18_logOR", "auc"], "Moufarrej, re-tested",
             V.loc["V7_MOGAFS", "auc"], "V7", "prereg"),
            ("External, other laboratory (Del Vecchio, 8/8)", R.loc[("S3_DelVecchio", "Moufarrej18"), "auc"], "Moufarrej, re-tested",
             R.loc[("S3_DelVecchio", "SI-MOGAFS"), "auc"], "XAI-MOGAFS-SI", "posthoc")]
    fig, ax = plt.subplots(figsize=(W, 3.3))
    for i, (lab, lit, litl, our, ourl, kind) in enumerate(rows[::-1]):
        col = {"ours": C["MOGAFS"], "posthoc": C["SI-MOGAFS"], "prereg": C["MOGAFS"]}[kind]
        ax.plot([lit, our], [i, i], color=GRID, lw=3, zorder=1, solid_capstyle="round")
        ax.scatter(lit, i, s=40, marker="s", color=C["Moufarrej"], zorder=3, edgecolor=SURF, lw=0.6)
        ax.scatter(our, i, s=46, color=col, zorder=3, edgecolor=SURF, lw=0.6)
        hi, lo = (our, lit) if our >= lit else (lit, our)
        ax.text(hi + 0.012, i + 0.17, f"{(our if our >= lit else lit):.2f} {ourl if our >= lit else litl}", fontsize=5.8,
                color=INK2, va="center")
        ax.text(lo - 0.012, i - 0.17, f"{(lit if our >= lit else our):.2f} {litl if our >= lit else ourl}", fontsize=5.8,
                color=INK2, va="center", ha="right")
    chance(ax)
    ax.set_yticks(range(len(rows)), [r[0] for r in rows[::-1]], fontsize=6.4)
    ax.set_xlim(0.0, 1.18)
    ax.set_xlabel("AUC")
    hgrid(ax)
    ax.legend(handles=[Line2D([], [], marker="s", ls="", color=C["Moufarrej"], label="Literature (best published comparator)"),
                       Line2D([], [], marker="o", ls="", color=C["MOGAFS"], label="This study, XAI-MOGAFS"),
                       Line2D([], [], marker="o", ls="", color=C["SI-MOGAFS"], label="This study, XAI-MOGAFS-SI (post hoc)")],
              loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3)
    save(fig, "fig_design")


# ------------------------------------------------------------------------------------------ 13. sPTB
def fig_ptb():
    t = (RESULTS / "ptb/report.txt").read_text(encoding="utf-8")
    dev = dict(re.findall(r"^\s+(XAI-MOGAFS|LASSO5|TTEST5)\s+AUC ([\d.]+) \(", t, re.M))
    ext = dict(re.findall(r"^\s+(XAI-MOGAFS|LASSO5|TTEST5)\s+AUC ([\d.]+) \[", t, re.M))
    ridge = float(re.search(r"ridge .*?AUC ([\d.]+)", t).group(1))
    L = json.loads((RESULTS / "ptb/locked_models.json").read_text(encoding="utf-8"))
    E = pd.read_csv(ROOT / "data/processed/expr/GSE59491.tsv.gz", sep="\t", dtype={"entrez": str}).set_index("entrez").T
    w = pd.read_csv(ROOT / "data/processed/selection/women_all.tsv", sep="\t", dtype=str)
    w = w[w.cohort == "GSE59491"].set_index("woman_id").loc[E.index]
    y = w.group.isin(["sPTD", "PPROM"]).astype(int).values
    fig, (a, b) = plt.subplots(1, 2, figsize=(W, 2.7), gridspec_kw={"wspace": 0.35, "width_ratios": [1.2, 1]})
    meths = [("XAI-MOGAFS", C["MOGAFS"], "XAI-MOGAFS"), ("LASSO5", C["LASSO"], "LASSO-5"), ("TTEST5", C["TTEST"], "t-test-5")]
    xs = np.arange(len(meths))
    for k, (m, col, lab) in enumerate(meths):
        a.bar(k - 0.18, float(dev[m]), 0.34, color=col, alpha=0.45, edgecolor=SURF)
        a.bar(k + 0.18, float(ext[m]), 0.34, color=col, edgecolor=SURF)
        a.text(k - 0.18, float(dev[m]) + 0.01, dev[m][:5], ha="center", fontsize=6.2)
        a.text(k + 0.18, float(ext[m]) + 0.01, ext[m][:5], ha="center", fontsize=6.2)
        Z = (E - E.mean()) / E.std(ddof=1)
        s = Z[L[m]["genes"]].values @ np.array(L[m]["coef"])
        roc_line(b, y, s, col, lab)
    a.axhline(ridge, color=INK2, ls=":", lw=0.8)
    a.text(2.45, ridge + 0.01, f"ridge (dev.) {ridge:.2f}", fontsize=6, color=INK2, ha="right")
    chance(a, "y")
    a.set_xticks(xs, [m[2] for m in meths])
    a.set_ylim(0.3, 0.7)
    a.set_ylabel("AUC")
    a.set_title("a  Development CV (light) vs locked external (dark)", fontsize=7.4)
    a.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=NEUTRAL, alpha=0.45, label="GSE149440 nested CV (118 women)"),
                      plt.Rectangle((0, 0), 1, 1, color=NEUTRAL, label="GSE59491 locked test (51 sPTB / 114 term)")],
             loc="upper left", fontsize=5.8)
    hgrid(a, "y")
    b.plot([0, 1], [0, 1], color=INK2, lw=0.6, ls=(0, (3, 2)))
    b.set_aspect("equal")
    b.set_xlabel("1 − specificity")
    b.set_ylabel("Sensitivity")
    b.set_title("b  GSE59491, sPTB vs term", fontsize=7.4)
    b.legend(loc="lower right", fontsize=6)
    hgrid(b, "both")
    save(fig, "fig_ptb")


# ------------------------------------------------------------------------------------------ 14. merkez-değişmez SHAP
def fig_shap():
    import importlib.util
    spec = importlib.util.spec_from_file_location("p12", Path(__file__).resolve().parent / "12_site_invariant.py")
    m12 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m12)
    from mogafs import fit_lr
    X, info, U, sym = m12.build()
    held, cfg = 1, m12.CFG
    trm = (info.site != held) & info.group.isin(["Control", "PE"])
    tr = np.flatnonzero(trm)
    y = (info.group == "PE").astype(int).values
    site = info.site.values
    keep = (X[tr] > cfg["min_expr"]).mean(0) >= cfg["min_frac"]
    Z = np.empty_like(X[:, keep])
    for s_ in np.unique(site):
        m = site == s_
        mu, sd = X[m][:, keep].mean(0), X[m][:, keep].std(0, ddof=1)
        Z[m] = (X[m][:, keep] - mu) / np.where(sd == 0, 1, sd)
    genes = np.asarray(U)[keep]
    gsym = np.array([sym.get(g, g) for g in genes])
    Pn = pd.read_csv(RESULTS / "site_invariant/loso_panels.tsv", sep="\t")
    rows = []
    fig, axs = plt.subplots(1, 3, figsize=(W, 2.7), gridspec_kw={"wspace": 0.45, "width_ratios": [1.25, 1, 1]})
    for meth in ["SI-MOGAFS", "MOGAFS"]:
        panel = Pn[(Pn.held_out == "S2_Val2") & (Pn.method == meth)].symbols.iloc[0].split(";")
        idx = [int(np.flatnonzero(gsym == s)[0]) for s in panel]
        Ztr = Z[tr][:, idx]
        lr = fit_lr(Ztr, y[tr])
        phi = lr.coef_[0] * (Ztr - Ztr.mean(0))
        st = site[tr]
        for j, s in enumerate(panel):
            imp = {k: np.abs(phi[st == k, j]).mean() for k in np.unique(st)}
            rows.append({"method": meth, "gene": s, "coef": lr.coef_[0][j], **{f"site{k}": v for k, v in imp.items()}})
        if meth == "SI-MOGAFS":
            ax = axs[0]
            rng = np.random.default_rng(2)
            for j, s in enumerate(panel):
                yj = len(panel) - 1 - j
                for k, col, mk in [(0, C["SI-MOGAFS"], "o"), (2, C["Moufarrej"], "^")]:
                    v = phi[st == k, j]
                    ax.scatter(v, yj + rng.uniform(-0.2, 0.2, len(v)) + (0.12 if k else -0.12), s=4, color=col, lw=0,
                               marker=mk, alpha=0.8)
            ax.axvline(0, color=INK2, lw=0.6)
            ax.set_yticks(range(len(panel)), panel[::-1])
            ax.set_xlabel("SHAP value φ (log-odds)")
            ax.set_title("a  XAI-MOGAFS-SI panel", fontsize=7.4)
            hgrid(ax)
    T = pd.DataFrame(rows)
    T.to_csv(P / "shap_site_importance.tsv", sep="\t", index=False)
    for ax, meth, lab in [(axs[1], "SI-MOGAFS", "b  SI panel: mean |φ| by site"), (axs[2], "MOGAFS", "c  Pooled panel: mean |φ| by site")]:
        t = T[T.method == meth]
        yy = np.arange(len(t))[::-1]
        ax.barh(yy + 0.18, t.site0, 0.34, color=C["SI-MOGAFS"], edgecolor=SURF, label="Stanford")
        ax.barh(yy - 0.18, t.site2, 0.34, color=C["Moufarrej"], edgecolor=SURF, label="Del Vecchio")
        ax.set_yticks(yy, t.gene)
        ax.set_xlabel("mean |SHAP|")
        ax.set_title(lab, fontsize=7.4)
        hgrid(ax)
        ratio = (np.minimum(t.site0, t.site2) / np.maximum(t.site0, t.site2)).mean()
        ax.text(0.98, 0.02, f"mean min/max ratio {ratio:.2f}", transform=ax.transAxes, ha="right", fontsize=5.8, color=INK2)
    axs[1].legend(loc="upper right", fontsize=5.8)
    save(fig, "fig_shap")


# ------------------------------------------------------------------------------------------ 15. özgüllük haritası
def fig_spec_landscape():
    """26 yayımlanmış imza + 4 panelimiz: bağımsız hedef-dışı hücrelerde ortalama ve en kötü |AUC-0.5|."""
    M = pd.read_csv(RESULTS / "literature_specificity/matrix.tsv", sep="\t")
    M = M[(M.independent.astype(str) == "True") & M.group.isin(["sPTD", "PPROM", "FGR_SGA_FD"])]
    lit = M.assign(dev=(M.auc - 0.5).abs()).groupby("signature_id").dev.agg(["mean", "max"])
    R = pd.read_csv(P / "wb_panels_by_syndrome.tsv", sep="\t")
    R = R[R.panel == "PE paneli"]
    cells = [("GSE59491", "sPTD"), ("GSE59491", "PPROM"), ("GSE142974", "sPTD"), ("GSE108497", "FGR_SGA_FD")]
    ours = {}
    for m in ["MOGAFS", "MOGAFS_spec", "LASSO5", "TTEST5"]:
        d = R[R.method == m]
        dv = [abs(d[(d.cohort == c) & (d.group == g)].auc.iloc[0] - 0.5) for c, g in cells]
        ours[m] = (float(np.mean(dv)), float(np.max(dv)))
    fig, (a, b) = plt.subplots(1, 2, figsize=(W, 3.2), gridspec_kw={"width_ratios": [1.3, 1], "wspace": 0.5})
    a.scatter(lit["mean"], lit["max"], s=22, color=NEUTRAL, alpha=0.75, lw=0, label="published signatures (26)", zorder=2)
    lab = {"MOGAFS": ("XAI-MOGAFS", C["MOGAFS"]), "MOGAFS_spec": ("XAI-MOGAFS ($f_2$)", C["MOGAFS_spec"]),
           "LASSO5": ("LASSO-5", C["LASSO"]), "TTEST5": ("$t$-test-5", C["TTEST"])}
    for m, (mx, mm) in ours.items():
        n, col = lab[m]
        a.scatter(mx, mm, s=70, color=col, edgecolor=INK, lw=0.6, zorder=4, marker="D" if "MOGAFS" in m else "s")
        right = mx > 0.1
        a.annotate(n, (mx, mm), xytext=(-7 if right else 7, 4), textcoords="offset points", fontsize=6.5, color=INK,
                   ha="right" if right else "left")
    a.set_xlabel("Mean off-target signal, mean $|\\mathrm{AUC}-0.5|$")
    a.set_ylabel("Worst-case off-target signal, max $|\\mathrm{AUC}-0.5|$")
    a.set_title("a  Off-target signal of 30 PE signatures and panels", fontsize=7.4)
    a.text(0.02, 0.97, "better", transform=a.transAxes, fontsize=6.5, color=INK2, va="top")
    a.annotate("", xy=(0.02, 0.86), xytext=(0.12, 0.96), xycoords="axes fraction",
               arrowprops={"arrowstyle": "-|>", "color": INK2, "lw": 0.7})
    a.legend(loc="upper center", fontsize=6)
    hgrid(a, "both")
    # b: en kötü durum sıralaması
    short = {"Tarca2021_LDA3_MMP11_SLC6A2_IL18BP": "Tarca 2021 LDA-3b", "Tarca2021_LDA3_MMP11_SLC6A2_AIF1L": "Tarca 2021 LDA-3a",
             "Tarca2021_LDA4_H19_FN1_TUBB6_FPR3": "Tarca 2021 LDA-4", "Tarca2021_EVT20_wholeblood": "Tarca 2021 EVT-20",
             "Tarca2021_untargeted16_wholeblood": "Tarca 2021 untargeted", "Lekva2021_TERC_PBMC": "Lekva 2021 TERC",
             "Lekva2021_PBMC_pilotarray_top20": "Lekva 2021 top-20", "Elovitz2025_placental_vs_immune_6gene": "Elovitz 2025 (6, b)",
             "Elovitz2025_placentalHDP_6gene": "Elovitz 2025 (6, a)", "Elovitz2025_PE_vs_all_top2": "Elovitz 2025 (2)",
             "Hong2019_SLE_PE_P1_58probe": "Hong 2019 SLE-PE", "Hong2019_M5.15_neutrophil_module": "Hong 2019 M5.15",
             "Moufarrej2022_cfRNA_LR18": "Moufarrej 2022", "Moufarrej2022_cfRNA_LR18_logOR": "Moufarrej 2022 (log-OR)",
             "Weiner2022_EOP_NAMPT_APOA1": "Weiner 2022", "Weiner2026_EOP_NAMPT": "Weiner 2026",
             "Froystad2026_PBMC_up10_22_24wk": "Frøystad 2026", "Myhrer2024_RPPH1_leukocytes": "Myhrer 2024"}
    allv = [(short.get(s, s.split("_")[0][:-4] + " " + s.split("_")[0][-4:]), v, NEUTRAL) for s, v in lit["max"].items()]
    allv += [(lab[m][0], v[1], lab[m][1]) for m, v in ours.items()]
    allv.sort(key=lambda t: t[1])
    top = allv[:12]
    y = np.arange(len(top))[::-1]
    b.barh(y, [t[1] for t in top], color=[t[2] for t in top], height=0.68, edgecolor=SURF)
    for yy, (n, v, col) in zip(y, top):
        b.text(v + 0.003, yy, f"{v:.3f}", va="center", fontsize=6, color=INK2)
    b.set_yticks(y, [t[0] for t in top], fontsize=6)
    b.set_xlabel("Worst-case off-target signal (lower is better)")
    b.set_title("b  The 12 most syndrome-neutral of 30", fontsize=7.4)
    rank = [t[0] for t in allv]
    b.text(0.98, 0.99, f"of 30: LASSO-5 ranks {rank.index('LASSO-5') + 1},\n$t$-test-5 ranks {rank.index('$t$-test-5') + 1}",
           transform=b.transAxes, ha="right", va="top", fontsize=6, color=INK2)
    b.set_xlim(0, 0.16)
    hgrid(b)
    save(fig, "fig_spec_landscape")


# ------------------------------------------------------------------------------------------ 16. dış kadınlarda üstünlük
def fig_superiority():
    L = pd.read_csv(RESULTS / "site_invariant/loso_scores.tsv", sep="\t")
    d = L[(L.held_out == "S2_Val2") & L.earliest.astype(bool) & L.group.isin(["Control", "PE"])]
    Pv = d.pivot_table(index="sample", columns="method", values="score")
    y = (d.drop_duplicates("sample").set_index("sample").loc[Pv.index].group == "PE").astype(int).values
    thr = d.groupby("method").threshold.first()
    rng = np.random.default_rng(2026)
    pos, neg = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    comps = [("Moufarrej18", "Moufarrej 2022 (18 genes)"), ("MOGAFS", "XAI-MOGAFS, pooled"), ("TTEST5", "$t$-test-5"),
             ("LASSO5", "LASSO-5"), ("META-TTEST5", "meta-$t$-test-5")]
    boots = {k: [] for k, _ in comps}
    for _ in range(2000):
        idx = np.r_[rng.choice(pos, len(pos)), rng.choice(neg, len(neg))]
        yy = y[idx]
        a_si = roc_auc_score(yy, Pv["SI-MOGAFS"].values[idx])
        for k, _ in comps:
            boots[k].append(a_si - roc_auc_score(yy, Pv[k].values[idx]))
    rows = []
    for k, labl in comps:
        bt = np.array(boots[k])
        rows.append({"comparator": labl, "delta": roc_auc_score(y, Pv["SI-MOGAFS"]) - roc_auc_score(y, Pv[k]),
                     "lo": np.percentile(bt, 2.5), "hi": np.percentile(bt, 97.5), "p_sup": float((bt > 0).mean())})
    T = pd.DataFrame(rows)
    T.to_csv(P / "si_paired_bootstrap.tsv", sep="\t", index=False)
    fig = plt.figure(figsize=(W, 4.6))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.05], hspace=0.62, wspace=0.85)
    a = fig.add_subplot(gs[0, :2])
    for i, r in enumerate(T.itertuples()):
        yy = len(T) - 1 - i
        col = C["Moufarrej"] if "Moufarrej" in r.comparator else C["SI-MOGAFS"]
        a.errorbar(r.delta, yy, xerr=[[r.delta - r.lo], [r.hi - r.delta]], fmt="o", ms=5, color=col, elinewidth=1.4, capsize=0)
        a.text(0.52, yy, f"$\\Delta$ {r.delta:+.2f}   P(SI better) {r.p_sup:.2f}", va="center", fontsize=6.2, color=INK2)
    a.axvline(0, color=INK2, lw=0.7, ls=(0, (3, 2)))
    a.set_yticks(range(len(T)), [f"vs {c}" for c in T.comparator[::-1]])
    a.set_xlim(-0.2, 0.9)
    a.set_xlabel("AUC difference, XAI-MOGAFS-SI $-$ comparator (paired bootstrap 95% CI)")
    a.set_title("a  XAI-MOGAFS-SI vs each comparator (same 87 women)", fontsize=7.4)
    hgrid(a)
    # b: gen sayısı vs dış AUC
    b = fig.add_subplot(gs[0, 2])
    V = pd.read_csv(P / "val2_full_metrics.tsv", sep="\t").set_index("model")
    LF = pd.read_csv(P / "loso_full_metrics.tsv", sep="\t")
    LF = LF[(LF.held_out == "S2_Val2") & LF.contrast.str.startswith("PE")].set_index("method")
    pts = [("XAI-MOGAFS-SI", 5, LF.loc["SI-MOGAFS", "auc"], C["SI-MOGAFS"], "o"),
           ("Moufarrej 2022", 16, V.loc["Moufarrej18_logOR", "auc"], C["Moufarrej"], "s"),
           ("XAI-MOGAFS pooled", 4, LF.loc["MOGAFS", "auc"], C["MOGAFS"], "o"),
           ("$t$-test-5", 5, LF.loc["TTEST5", "auc"], C["TTEST"], "o"), ("LASSO-5", 5, LF.loc["LASSO5", "auc"], C["LASSO"], "o"),
           ("V7 (pre-reg.)", 4, V.loc["V7_MOGAFS", "auc"], C["MOGAFS"], "^"),
           ("V0pc", 4, V.loc["V0pc_MOGAFS", "auc"], "#1c5cab", "^"), ("V8 LASSO", 18, V.loc["V8_LASSO", "auc"], C["LASSO"], "^")]
    for n, g, auc, col, mk in pts:
        b.scatter(g, auc, s=40 if n == "XAI-MOGAFS-SI" else 24, color=col, marker=mk, edgecolor=INK if n == "XAI-MOGAFS-SI" else SURF,
                  lw=0.6, zorder=3)
        off = {"XAI-MOGAFS pooled": (5, -6), "LASSO-5": (5, -4), "$t$-test-5": (5, 2), "V0pc": (-4, 3)}.get(n, (4, 3))
        b.annotate(n, (g, auc), xytext=off, textcoords="offset points", fontsize=5.4,
                   ha="right" if off[0] < 0 else "left")
    chance(b, "y")
    b.set_xlabel("Genes in the panel")
    b.set_ylabel("External AUC (87 women)")
    b.set_xlim(2, 21)
    b.set_ylim(0.35, 0.85)
    b.set_title("b  AUC vs panel size", fontsize=7.4)
    hgrid(b, "both")
    # c, d: karmaşıklık matrisleri
    for j, (k, name, col) in enumerate([("SI-MOGAFS", "XAI-MOGAFS-SI", "#1c5cab"), ("Moufarrej18", "Moufarrej 2022", "#b4471d")]):
        ax = fig.add_subplot(gs[1, j])
        pred = (Pv[k].values >= thr[k]).astype(int)
        cm = np.array([[((pred == 1) & (y == 1)).sum(), ((pred == 0) & (y == 1)).sum()],
                       [((pred == 1) & (y == 0)).sum(), ((pred == 0) & (y == 0)).sum()]])
        frac = cm / cm.sum(1, keepdims=True)
        ax.imshow(frac, cmap=LinearSegmentedColormap.from_list("c", ["#f7f6f3", col]), vmin=0, vmax=1)
        for (r, c_), v in np.ndenumerate(cm):
            ax.text(c_, r, f"{v}\n({frac[r, c_]:.0%})", ha="center", va="center", fontsize=6.6,
                    color="white" if frac[r, c_] > 0.55 else INK)
        ax.set_xticks([0, 1], ["test +", "test $-$"])
        ax.set_yticks([0, 1], ["PE", "control"])
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_title(f"{'c' if j == 0 else 'd'}  {name}", fontsize=7.4)
    ax = fig.add_subplot(gs[1, 2])
    ax.axis("off")
    si, mo = LF.loc["SI-MOGAFS"], LF.loc["Moufarrej18"]
    txt = ("Thresholds fixed on training sites\n\n"
           f"False alarms in controls:\n  XAI-MOGAFS-SI {int(si.fp)}/61, Moufarrej {int(mo.fp)}/61\n"
           f"PE cases detected:\n  XAI-MOGAFS-SI {int(si.tp)}/26, Moufarrej {int(mo.tp)}/26\n"
           f"Specificity {si.spec:.2f} vs {mo.spec:.2f}\nPPV {si.ppv:.2f} vs {mo.ppv:.2f}\nMCC {si.mcc:.2f} vs {mo.mcc:.2f}")
    ax.text(0.0, 0.95, txt, va="top", fontsize=6.4, color=INK, linespacing=1.35)
    save(fig, "fig_superiority")


if __name__ == "__main__":
    which = sys.argv[1:] or ["cohorts", "wholeblood", "lit", "optimism", "syndrome", "pareto", "variants", "rn7sl", "roc",
                             "loso", "metrics", "design", "ptb", "shap", "speclandscape", "superiority"]
    F = {"cohorts": fig_cohorts, "wholeblood": fig_wholeblood, "lit": fig_lit_heatmap, "optimism": fig_optimism,
         "syndrome": fig_syndrome_ours, "pareto": fig_pareto, "variants": fig_variants, "rn7sl": fig_rn7sl, "roc": fig_roc,
         "loso": fig_loso, "metrics": fig_metrics, "design": fig_design, "ptb": fig_ptb, "shap": fig_shap,
         "speclandscape": fig_spec_landscape, "superiority": fig_superiority}
    for k in which:
        F[k]()
