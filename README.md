# Syndrome specificity and cross-site transportability of early-pregnancy blood RNA signatures for preeclampsia

Analysis code and result files for the manuscript

> Karadayı Ataş P., Vanneschi L. *Syndrome specificity and cross-site transportability of early-pregnancy
> maternal-blood transcriptomic signatures for preeclampsia: a multi-cohort benchmark with explainable multi-objective
> genetic feature selection* (manuscript in preparation).

The repository contains **XAI-MOGAFS**, a SHAP-guided multi-objective genetic feature selector (NSGA-II) with explicit
syndrome-specificity, stability and site-invariance objectives. It also contains the full pipeline used in the study:

| Arm | Question |
|---|---|
| A | Do published early-pregnancy PE signatures recognise PE in independent cohorts, and do they stay silent in other great obstetrical syndromes (sPTD, PPROM, FGR, GDM)? |
| B | Do XAI-MOGAFS whole-blood panels carry less off-target signal than LASSO / *t*-test panels? |
| C | Do plasma cell-free RNA panels, locked before testing, transport to sealed external cohorts? |
| D | Negative control: spontaneous preterm birth, developed in one cohort and tested once in another. |

## Repository layout

```
scripts/        analysis code, run in numerical order (comments are in Turkish)
  mogafs.py       XAI-MOGAFS: objectives, SHAP-guided sampling / crossover / mutation / repair, selection rule
  panel_utils.py  shared helpers (Welch t pre-filter, LASSO with <=k features, bootstrap AUC CI)
  geo_utils.py    GEO series-matrix readers
results/        all result files used for the tables and figures (TSV / JSON / Markdown)
  plasma/LOCK_REGISTRY.md   lock registry: SHA-256 of every locked model and of the sealed data
  plasma/VARIANTS_LOG.md    every plasma variant tried before the sealed test
  final_test/               single-use sealed test (OPENED.txt = opening time stamp)
logs/download_manifest.tsv  source URL, size and SHA-256 of every downloaded file
requirements.txt            exact package versions (Python 3.11)
```

Raw data are **not** included. All of them are public in NCBI GEO and are downloaded by `scripts/01_download.py`.

## Installation

```bash
py -3.11 -m venv .venv                      # or: python3.11 -m venv .venv
.venv/Scripts/python -m pip install -r requirements.txt   # Linux/macOS: .venv/bin/python
```

## Data

```bash
python scripts/01_download.py
```

The script downloads the processed series matrices, count tables and platform annotations to `data/`:

- Whole blood: GSE85307, GSE149440, GSE59491, GSE142974, GSE108497.
- Plasma cfRNA: GSE192902, GSE154377.
- NCBI gene info.

Compare the SHA-256 digests with `logs/download_manifest.tsv`. The two external test sets go to `data/sealed/`:

- GSE192902 Validation-2;
- GSE154377 expression data.

## Pipeline

| Script | Purpose | Main output |
|---|---|---|
| `00_dataset_search.py`, `00b_inspect_candidates.py` | systematic GEO search and inspection of candidate series | `results/dataset_search/` |
| `01_download.py` | download and SHA-256 manifest | `data/`, `logs/download_manifest.tsv` |
| `02`–`04_inventory_*.py` | cohort inventories and checks | `results/inventory/` |
| `05_select_samples.py` | one sample per woman (pre-specified rule), harmonised groups, flow | `results/flow/` |
| `06_gene_level.py` | probe → Entrez gene, cross-platform gene universe | `data/processed/` |
| `07_build_panel.py` | Arm B: whole-blood panels (XAI-MOGAFS ± specificity objective, LASSO-5, *t*-test-5), nested CV | `results/panel/` |
| `07b_signal_diagnostics.py`, `07c_positive_control.py` | whole-blood signal diagnostics; plasma positive control | `results/panel/` |
| `08_specificity_matrix.py` | off-target behaviour of the locked whole-blood panels | `results/panel/` |
| `09_plasma_panel.py`, `09b_plasma_diagnostics.py`, `09c_plasma_repeated.py` | Arm C: plasma variants V0–V8, development and locking; 7SL-RNA diagnostics | `results/plasma/` |
| `10_literature_specificity.py` | Arm A: syndrome-by-syndrome re-test of 26 published signatures | `results/literature_specificity/` |
| `11_final_test.py` | **single-use sealed test** of all registered locks and published signatures | `results/final_test/` |
| `12_site_invariant.py` | *post hoc*: XAI-MOGAFS-SI, leave-one-site-out evaluation | `results/site_invariant/` |
| `13_comparison_table.py`, `15_design_matched.py` | literature comparison by evaluation design (*exploratory*) | `results/comparison/` |
| `14_ptb_panel.py` | Arm D: spontaneous preterm birth negative control | `results/ptb/` |
| `16_paper_metrics.py` | additional metrics on locked models (no re-training) | `results/paper/` |
| `17_figures.py`, `18_latex_tables.py` | manuscript figures and tables, read only from the result files | `paper/figures/`, `paper/tables/` |

### Sealed test

- `11_final_test.py` was run **once**, after every model had been locked.
- It verifies all digests and processes the sealed data without labels.
- It then writes `results/final_test/OPENED.txt` and only after that merges the labels.
- It refuses to run again unless a documented technical fix is given.
- Scripts `12` and `15` were designed **after** the sealed test was opened and are exploratory.

## Reproducibility

- All stochastic components use seed 2026.
- Python 3.11 with the package versions in `requirements.txt`.
- Running times on a 14-core workstation:
  - whole-blood nested CV: ~8 min;
  - plasma leave-site-out: ~3 min;
  - preterm-birth experiment: ~6 min.

## Citation

See `CITATION.cff`. The manuscript reference will be added on publication.
