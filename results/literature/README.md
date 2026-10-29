# Published preeclampsia (PE) signatures in maternal blood

Catalogue of **published** PE prediction/association signatures measured in maternal blood
(whole blood, PBMC/leukocytes, plasma or serum cell-free RNA, circulating/EV miRNA), focused on
sampling **before ~24 weeks** of gestation. Signatures sampled at diagnosis are included and flagged
`at_diagnosis = yes`.

Main table: **`published_signatures.tsv`** — one row per (signature, feature).
Downloaded supplementary files: **`supp/<FirstAuthorYear>/`**.

**Rule applied throughout:** only feature names actually read in the paper, its tables, figure
legends/panels or supplementary files were recorded. Nothing was completed from memory. Where a list
could not be fully retrieved, the signature is marked `verified = partial` and the note says exactly
what is missing.

Totals: **51 signatures from 31 studies, 1270 feature rows** (1149 `yes`, 121 `partial`).

---

## 1. Summary table (one line per signature)

| signature_id | study (year) | compartment | molecule | sampling window | at_dx | #features | verified |
|---|---|---|---|---|---|---|---|
| Moufarrej2022_cfRNA_LR18 | Moufarrej, Nature (2022) | plasma cfRNA | mRNA +2 Y_RNA +1 pseudogene | 5–16 wk | no | 18 | yes |
| Rasmussen2022_cfRNA_7gene | Rasmussen, Nature / Mirvie (2022) | plasma cfRNA | mRNA +1 snoRNA | 16–27 wk | no | 7 | yes |
| CastilloMarco2025_EOPE_T1 | Castillo-Marco, Nat Commun (2025) | plasma cfRNA | mixed cfRNA | T1 9–14 wk | no | 36 | yes |
| CastilloMarco2025_LOPE_T1 | Castillo-Marco, Nat Commun (2025) | plasma cfRNA | mixed cfRNA | T1 9–14 wk | no | 19 | yes |
| CastilloMarco2025_EOPE_T2 | Castillo-Marco, Nat Commun (2025) | plasma cfRNA | mixed cfRNA | T2 18–28 wk | no | 87 | yes |
| CastilloMarco2025_LOPE_T2 | Castillo-Marco, Nat Commun (2025) | plasma cfRNA | mixed cfRNA | T2 18–28 wk | no | 94 | partial |
| Hong2019_SLE_PE_P1_58probe | Hong, J Exp Med / PROMISSE (2019) | whole blood | mRNA | P1 <16 wk | no | 58 | yes |
| Hong2019_M5.15_neutrophil_module | Hong, J Exp Med (2019) | whole blood | mRNA | P1 <16 wk | no | 23 | yes |
| Mirzakhani2016_VDAART_348 | Mirzakhani, JCI / VDAART (2016) | whole blood | mRNA | 10–18 wk | no | 348 | yes |
| Yadama2020_VDAART_overlap238 | Yadama, Sci Rep / VDAART (2020) | whole blood | mRNA | 10–18 wk | no | 238 | yes |
| Mirzakhani2023_miRNA31 | Mirzakhani, Clin Transl Med (2023) | peripheral blood | miRNA | 10–18 wk | no | 31 | yes |
| Mirzakhani2023_miRNA16_PE_and_vitD | Mirzakhani, Clin Transl Med (2023) | peripheral blood | miRNA | 10–18 wk | no | 16 | yes |
| Chaiworapongsa2013_EOPE_atDx | Chaiworapongsa, J Perinat Med (2013) | whole blood | mRNA +H19 | at diagnosis (~31 wk) | **yes** | 43 | yes |
| Chaiworapongsa2013_LOPE_atDx | Chaiworapongsa, J Perinat Med (2013) | whole blood | mRNA +H19 | at diagnosis (~38 wk) | **yes** | 28 | yes |
| Tarca2021_untargeted16_wholeblood | Tarca, J Matern Fetal Neonatal Med (2021) | whole blood | mRNA +H19 | 9–34 wk longitudinal | no | 16 | yes |
| Tarca2021_LDA4_H19_FN1_TUBB6_FPR3 | Tarca (2021) | whole blood | mRNA +H19 | 22–28 wk | no | 4 | yes |
| Tarca2021_EVT20_wholeblood | Tarca (2021) | whole blood | mRNA | 11–17 wk | no | 20 | yes |
| Tarca2021_LDA3_MMP11_SLC6A2_AIF1L | Tarca (2021) | whole blood | mRNA | 11–17 wk | no | 3 | yes |
| Tarca2021_LDA3_MMP11_SLC6A2_IL18BP | Tarca (2021) | whole blood | mRNA | 11–17 wk | no | 3 | yes |
| DelVecchio2021_PreX_gHTN_LR5 | Del Vecchio, Epigenetics (2021) | plasma cfRNA | mRNA +HAGLR | 12–17 wk | no | 5 | yes |
| DelVecchio2021_APO_LR10 | Del Vecchio (2021) | plasma cfRNA | mRNA +HAGLR | 12–17 wk | no | 10 | yes |
| Timofeeva2023_core11 | Timofeeva, Int J Mol Sci (2023) | serum | miRNA | 11–14 wk | no | 11 | yes |
| Timofeeva2023_lPE_LR5a / _lPE_LR5b | Timofeeva (2023) | serum | miRNA | 11–14 wk | no | 5 / 5 | yes |
| Timofeeva2023_universal3_plus_clinical | Timofeeva (2023) | serum | miRNA | 11–14 wk | no | 3 | yes |
| Timofeeva2023_ePE_NGS_DE / _lPE_NGS_DE | Timofeeva (2023) | serum | miRNA | 11–14 wk | no | 2 / 10 | yes |
| Hedley2025_serum_miRNA7 | Hedley, BMC Pregnancy Childbirth (2025) | serum | miRNA | 10–14 wk | no | 7 | yes |
| Gal2024_EV_miRNA16 | Gál, Front Immunol (2024) | EV miRNA (exosomes) | miRNA | 10–13 wk | no | 16 | yes |
| Hromadnikova2022_CVDmiRNA_PE6 | Hromadnikova, Biomedicines (2022) | whole blood | miRNA | 10–13 wk | no | 6 | yes |
| Hromadnikova2019_exoC19MC_3 | Hromadnikova, Int J Mol Sci (2019) | EV miRNA (exosomes) | miRNA | 10–13 wk | no | 3 | yes |
| Hromadnikova2017_plasmaC19MC_3 | Hromadnikova, PLoS One (2017) | plasma | miRNA | 10–13 wk | no | 3 | yes |
| Gusar2022_EV_hypoxamiR2 | Gusar, Life (2022) | EV miRNA (exosomes) | miRNA | 11–14 wk | no | 2 | yes |
| Devor2020_EV_miRNA_T1 | Devor, J Matern Fetal Neonatal Med (2020) | EV miRNA (exosomes) | miRNA | first trimester | no | 4 | yes |
| Winger2018_buffycoat_miRNA8 | Winger, PLoS One (2018) | leukocytes (buffy coat) | miRNA | 12.7 wk | no | 8 | yes |
| Elovitz2025_placentalHDP_6gene | Elovitz, Nat Commun / Mirvie (2025) | plasma cfRNA | mRNA | 17w4d–22w0d | no | 6 | yes |
| Elovitz2025_placental_vs_immune_6gene | Elovitz (2025) | plasma cfRNA | mRNA | 17w4d–22w0d | no | 6 | yes |
| Elovitz2025_PE_vs_all_top2 | Elovitz (2025) | plasma cfRNA | mRNA | 17w4d–22w0d | no | 2 | partial |
| Gunaydin2025_LOPE_12gene | Gunaydin, Sci Rep / ITU (2025) | whole blood (PAXgene) | mRNA | 20–24 wk | no | 12 | yes |
| Lekva2021_PBMC_pilotarray_top20 | Lekva, Sci Rep / STORK (2021) | PBMC | mRNA/ncRNA | 22–24 wk | no | 20 | partial |
| Lekva2021_TERC_PBMC | Lekva (2021) | PBMC | lncRNA | 22–24 wk | no | 1 | yes |
| Myhrer2024_RPPH1_leukocytes | Myhrer/Frøystad, Biol Reprod (2024) | leukocytes | lncRNA | 22–24 wk | no | 1 | yes |
| Lekva2023_MT_TA_leukocytes | Lekva, BMC Med (2023) | leukocytes / plasma EV | mito tRNA/mRNA | 22–24 wk | no | 3 | yes |
| Froystad2026_PBMC_up10_22_24wk | Frøystad, Biol Reprod (2026) | PBMC | mRNA/lncRNA | 22–24 wk | no | 10 | yes |
| Froystad2026_NOTCH2NL_qPCR | Frøystad (2026) | PBMC | mRNA | 22–24 wk | no | 1 | yes |
| Weiner2022_EOP_NAMPT_APOA1 | Weiner, Diagnostics (2022) | plasma cfRNA | mRNA | 16.0–20.7 wk | no | 2 | yes |
| Weiner2026_EOP_NAMPT | Weiner, Diagnostics / FutureBIRTH (2026) | plasma cfRNA | mRNA | 12–20 wk | no | 1 | yes |
| Carter2025_inflammatory_cfRNA8 | Carter, Reprod Sci (2025) | plasma cfRNA | mRNA | first trimester | no | 8 | yes |
| Munchel2020_cfRNA_49 | Munchel, Sci Transl Med (2020) | plasma cfRNA | mRNA | at diagnosis | **yes** | 0 (placeholder) | partial |
| Nakano2025_LOPE_recurrent3 | Nakano, Healthcare (2025) | plasma cfRNA | mRNA | at diagnosis (re-analysis) | **yes** | 3 | partial |
| Zhou2023_cfRNA_13 | Zhou, Am J Obstet Gynecol (2023) | plasma cfRNA | mRNA+miRNA+lncRNA | before symptoms | no | 0 (placeholder) | partial |

## 2. Signatures marked `partial` — what is missing

| signature_id | what is missing / why |
|---|---|
| `Zhou2023_cfRNA_13` | Article is paywalled (AJOG; no PMC record; sciencedirect.com and ajog.org return HTTP 403). The 13 classifier cfRNAs and the 77 DE genes are in the article/appendix. **No feature name verified** — one placeholder row only. |
| `Munchel2020_cfRNA_49` | Sci Transl Med not open access; science.org blocked. The 30 DE transcripts and the 49-transcript classifier are in the paywalled figures/tables. **No gene name verified** — placeholder row only. The underlying data are GEO-free but the cohort is re-used by Nakano 2025. |
| `Nakano2025_LOPE_recurrent3` | Only the 3 genes recurring in >50% of CV repeats (HLA-G, IL17RB, KLRC4) are named in the text; the 52-gene and 87-gene elastic-net LO-PE signatures are in the MDPI supplementary files (not retrieved). Also `at_diagnosis = yes` (re-analysis of Munchel samples). |
| `Elovitz2025_PE_vs_all_top2` | Only the top two genes by effect size (PAPPA2, CD163) are named in the text; the full PE-vs-all DE list is not printed in the article. |
| `Lekva2021_PBMC_pilotarray_top20` | Supplementary Table S2 prints only the top 10 up- and top 10 down-regulated molecules of the pilot array (n = 5 PE / 4 controls); the full array result is not published. |
| `CastilloMarco2025_LOPE_T2` | Supplementary Table 5 prints **94** entries for this model while the main text says 92 cfRNAs. The list is recorded exactly as printed; the 2-entry discrepancy is unresolved. |

Additional caveats carried in the `notes` column rather than as `partial`:

- **`Mirzakhani2016_VDAART_348`** — direction is taken from the `Regulation` column of Supplemental
  Table 1. That column is consistent with the paper's text ("348 … 158 upregulated" in PE) and with
  Mirzakhani 2023 (which calls IGF1R, `Regulation = Down`, downregulated), but the neighbouring column
  header reads `FC.Control/Preeclampsia` and its values run the other way (Up genes have FC > 1).
  **Verify the sign in GEO GSE85307 before scoring.**
- **`Moufarrej2022_cfRNA_LR18`** — `direction_in_pe` is derived from the published odds ratio
  (OR < 1 → lower in PE), not from a separately stated DE direction.
- **`Rasmussen2022_cfRNA_7gene`**, **Castillo-Marco models**, **`Hong2019_SLE_PE_P1_58probe`**,
  **`Timofeeva2023_*`**, **`Carter2025_*`** — the papers give no per-feature direction.
- **`Hong2019_*`** — lupus (PROMISSE) pregnancies only; PE comparisons are SLE-PE vs SLE-NC.
- **`Yadama2020_VDAART_overlap238`** — intersection of the PE signature with asthma / vitamin D /
  BMI signatures, i.e. **not a pure PE signature**; the full 507-gene replicated PE list is not published.
- **`DelVecchio2021_*`** — "preeclampsia" model actually pools PE with gestational hypertension
  (PreX/gHTN); the APO model additionally pools GDM.
- **Window >24 wk**: `Rasmussen2022` (to 27 wk), `CastilloMarco2025_*_T2` (to 28 wk),
  `Tarca2021_LDA4` (22–28 wk), `Gunaydin2025` (20–24 wk), STORK signatures (22–24 wk).

## 3. Papers checked and excluded

| paper | reason for exclusion |
|---|---|
| Tarca et al., Cell Rep Med 2021, `10.1016/j.xcrm.2021.100323` (DREAM challenge) | Early-PE cases (n = 13) are only part of the gestational-age training set; the published lists are the 249-gene gestational-age transcriptome (Table S3) and the 402-RNA PPROM list (Table S7) — **no PE-vs-control gene list**. The same Detroit cohort's PE signatures are published in Tarca 2021 JMFNM (included) and Chaiworapongsa 2013 (included). |
| Tsang et al., PNAS 2017, `10.1073/pnas.1710470114` | Placental single-cell study; EVT signature is increased in plasma cfRNA **at early-PE diagnosis**, and the article is not open access (Europe PMC: "not open access"). The 20-gene EVT list is captured via Tarca 2021 JMFNM, which prints it in full. |
| Maric et al., Patterns 2022, `10.1016/j.patter.2022.100655` | Best early model is **urine metabolites**, not a blood transcriptomic signature. |
| Farina 2010 `10.1016/j.ajog.2010.07.043`; Sekizawa 2010 `10.1111/j.1471-0528.2010.02491.x`; Purwosunu 2009 `10.1016/j.ajog.2008.11.035`; Zanello 2014 `10.1159/000360854` | Targeted candidate mRNA panels (ENG/FLT1/PlGF/PP13/TGFB1/PLAC1 etc.) measured by qPCR, i.e. pre-specified assays rather than discovery signatures; abstracts only (no full text retrieved) so exact panel membership and directions could not be verified from the sources. |
| Salomon et al., JCEM 2017 `10.1210/jc.2017-00672` | Candidate exosomal miRNAs named only as hsa-miR-486-1-5p / -486-2-5p in the abstract; no full list retrieved (paywalled). |
| Ura et al., Taiwan J Obstet Gynecol 2014 `10.1016/j.tjog.2014.03.001` | Paywalled; only 4 validated miRNAs named in the abstract (miR-1233, miR-520, miR-210, miR-144) and "miR-520" is not a valid unique miRNA name — could not verify the printed list. |
| Timofeeva et al., Placenta 2018 `10.1016/j.placenta.2017.11.011` | Paywalled; plasma miRNAs named in the abstract only, superseded by the 2023 IJMS study (included). |
| Clément et al., Pregnancy Hypertens 2023 `10.1016/j.preghy.2023.10.011` | Paywalled; abstract names only miR-194-5p and miR-1278 within a model that also uses clinical factors; the 73-miRNA discovery list and the 5 replicated miRNAs were not retrievable. |
| Luque et al., Sci Rep 2014 `10.1038/srep04882` (PMC4013938) | **Negative study** — no first-trimester serum miRNA was confirmed as differentially abundant. |
| Winger et al., J Reprod Immunol 2015 `10.1016/j.jri.2015.03.005` | Paywalled; 7-miRNA panel members not named in the abstract. Superseded by Winger 2018 (included). |
| Wang et al., J Perinatol 2016 `10.1038/jp.2015.192` | Leucocyte miRNAs but samples taken **at diagnosis of severe PE**; the 4 validated miRNAs are not named in the abstract and the article is paywalled. |
| Hu et al., Clin Epigenetics 2018 `10.1186/s13148-018-0482-3` | Discovery in **placenta**; blood measurement is a single circRNA (hsa_circ_0036877) validated in whole blood, not a multi-feature early-pregnancy signature. |
| Liao et al., Pak J Med Sci 2024 `10.12669/pjms.40.11.10317` | n = 3 vs 3 circRNA microarray; no named feature list in the abstract, journal full text not retrieved. |
| Puvanendran & Brew, Placenta 2026 `10.1016/j.placenta.2026.09.008` | Re-analysis of GSE154377 + GSE192902 (already covered by Del Vecchio 2021 and Moufarrej 2022); paywalled; only 4 genes named in the abstract (MYZAP, CAV2, OXTR, GBP1) out of 2894 DE genes. |
| Elovitz-adjacent / Mirvie 2025 `10.1038/s41586-025-61931-7` (Gong et al., LEP + PAPPA2) | PE **with fetal growth restriction** phenotype; signature discovery is at **28 wk** (no DE genes at 12 or 20 wk) — outside the <24 wk window. Retained as context only. |
| cfDNA-based tests (Nat Med 2025 `10.1038/s41591-024-03509-2`; Nat Med 2023 methylome; Clin Epigenetics 2024; Nat Commun 2026 fragmentomics) | cell-free **DNA** (methylation/fragmentomics), not RNA. |
| Ghosh et al., Reprod Sci 2026 `10.1007/s43032-026-02133-4` | **Urinary** cfRNA, not maternal blood. |
| Hromadnikova first-trimester sPTB/FGR papers (e.g. `10.3390/ijms23073729`, `10.3390/biomedicines10102591`) | Outcome is preterm birth or FGR, not PE. Useful later as specificity comparators. |
| Placental-tissue transcriptome studies, HLA/NK-cell immunology papers, systematic reviews and meta-analyses returned by the searches | Not maternal-blood primary signatures. |

## 4. Notable cross-study overlap

- **Neutrophil/granulocyte degranulation genes dominate the whole-blood signatures.** `MMP8`,
  `MS4A3`, `DEFA4`, `BPI`, `LTF`, `CEACAM6`, `TCN1`, `CAMP`, `CTSG`, `OLFM4`, `LCN2`, `CRISP3` recur
  across Chaiworapongsa 2013 (at diagnosis), Mirzakhani 2016 / Yadama 2020 (10–18 wk) and Hong 2019's
  M5.15 neutrophil module — 12–15 shared genes per pair. This is the strongest convergence in the
  whole catalogue and is also the signal most likely to be shared with other inflammatory obstetric
  syndromes, so it is the key thing to test for syndrome specificity.
- **Direction conflicts on the same genes:** the neutrophil genes are *down* in whole blood at PE
  diagnosis (Chaiworapongsa 2013) but *down* (MMP8, LTF, LCN2, OLFM4) at 10–18 wk in VDAART, while the
  Hong 2019 module score is *up* in pre-PE lupus pregnancies at <16 wk. Do not assume a common sign.
- **Essentially no overlap between the cfRNA classifiers.** Moufarrej 2022 (18 genes), Rasmussen 2022
  (7 genes) and Castillo-Marco 2025 (36/19/87/94) share **zero** features with each other. The only
  cfRNA genes recurring across independent groups are **PAPPA2** (Rasmussen 2022, Elovitz 2025 ×3,
  and in the Tarca EVT-20 list) and **NAMPT** (Del Vecchio 2021 ×2, Weiner 2022/2026).
- **miRNA overlap is thin:** `miR-20a-5p` recurs in Hromadnikova 2022 (whole blood, up),
  Timofeeva 2023 (serum core list), Hedley 2025 (serum, **down**) and Gál 2024 (EV, up) — with
  conflicting directions. `miR-122-5p` and `miR-144-5p` are shared between Mirzakhani 2023 (whole
  blood, up/down) and Gál 2024 (EV, up). `miR-517-5p` is up in whole plasma (Hromadnikova 2017) but
  down in plasma exosomes (Hromadnikova 2019) in the same group's cohorts.
- **Within-group redundancy:** Mirzakhani 2016 and Yadama 2020 share 154 genes (same VDAART cohort);
  Tarca 2021 JMFNM's 16-gene untargeted signature is the Chaiworapongsa 2013 early-PE list re-tested
  prospectively. Treat these as one evidence stream, not independent replications.

## 5. Column definitions (`published_signatures.tsv`)

`signature_id`, `study`, `year`, `doi`, `compartment`, `molecule`, `sampling_window`,
`at_diagnosis` (yes/no), `pe_subtype`, `n_pe`, `n_control`, `feature` (symbol or miRNA name **exactly
as printed**), `feature_id_if_given` (Ensembl/Entrez/probe ID), `direction_in_pe` (up/down/NA),
`weight` (OR, log2FC, FC, AUC, Ct cut-off, MoM or importance — the kind is stated in the value or the
note), `source_location` (exact table/figure/file), `verified` (yes/partial), `notes`.

Two rows carry `feature = NOT_EXTRACTED` (Munchel 2020, Zhou 2023): placeholders recording that the
signature exists and where it lives, with no gene names claimed.
