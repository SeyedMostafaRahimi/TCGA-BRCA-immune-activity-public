# TCGA-BRCA Immune Activity and CCL19–IL10 Expression States

Status: manuscript in preparation

This repository contains the R analysis code, reproducibility checks, and selected derived outputs associated with a retrospective bulk-transcriptomic study of breast cancer. TCGA-BRCA was used as the discovery and tumour-microenvironment characterization cohort, and SCAN-B/GSE96058 was used for an independent evaluation of selected CCL19, IL10, and joint-expression survival findings.

The repository is a computational research and reproducibility record. It is not a substitute for the manuscript or a clinical prognostic tool.

## Study question

The study asked whether a transcriptome-derived immune-activity axis can characterize biologically distinct breast-cancer tumour microenvironments and whether moving from this broad axis to specific immune-network genes provides more informative survival associations.

The analysis therefore proceeded from broad transcriptomic characterization to focused candidate-gene evaluation:

1. construct an immune-activity axis from 17 prespecified ssGSEA signatures;
2. characterize the resulting ordered strata with complementary TME scores;
3. compare Immune-High and Immune-Low tumours by differential expression and pathway analysis;
4. prioritize highly connected genes in an immune-focused protein–protein interaction network;
5. evaluate survival associations for the prioritized hub genes;
6. examine CCL19 and IL10 as individual predictors and as four joint expression states; and
7. evaluate the selected prognostic findings in SCAN-B/GSE96058.

## Cohorts and data scope

| Cohort | Role in this project | Samples used | Expression scope in the supplied workflow |
|---|---|---:|---|
| TCGA-BRCA | Discovery, TME characterization, pathway and network analyses, and survival analyses | 1,106 primary tumours after sample selection; 1,083 patients in the valid overall-survival cohort; 150 deaths | Final matrix: 35,009 gene symbols × 1,106 primary tumours |
| SCAN-B/GSE96058 | Independent evaluation of CCL19, IL10, and their joint expression states | 3,069 tumours; 322 deaths; 2,911 complete cases and 303 deaths in the adjusted Cox model | The production script uses deterministic CCL19/IL10 expression values and phenotype data; the complete TCGA immune-axis and pathway-discovery workflow is not reconstructed in SCAN-B |

The initial TCGA import checks 1,226 matched expression/clinical samples. These comprise 1,106 Primary Tumour, 113 Solid Tissue Normal, and 7 Metastatic samples; only Primary Tumour samples are retained for the downstream analysis.

The supplied scripts do not perform FASTQ alignment, transcript quantification, or read-level quality control. They operate on archived expression matrices and clinical/phenotype records. The exact upstream transformation of the archived TCGA and SCAN-B expression inputs must be documented separately before claiming end-to-end RNA-seq reproducibility.

## Workflow

```text
Public TCGA-BRCA expression and clinical inputs
                    ↓
Sample matching and Primary Tumour selection
                    ↓
Ensembl-to-symbol reconstruction and zero-variance filtering
                    ↓
ESTIMATE, MCP-counter, and xCell TME profiling
                    ↓
17-signature ssGSEA → centred/scaled PCA → PC1 immune-activity axis
                    ↓
Immune-Low / Immune-Mid / Immune-High characterization
                    ↓
Immune-High versus Immune-Low limma differential expression
                    ↓
GO-BP and KEGG over-representation + Hallmark GSEA
                    ↓
Five-Hallmark leading-edge union → preserved STRING/Cytoscape network
                    ↓
CytoHubba MCC top-20 hub prioritization and hub/TME correlations
                    ↓
TCGA overall-survival and hub-gene Cox analyses
                    ↓
CCL19–IL10 joint-state survival and immune-contexture analyses
                    ↓
Independent CCL19–IL10 evaluation in SCAN-B/GSE96058
```

## Analysis methods

### TCGA preprocessing

- TCGA expression and clinical records are matched by sample identifier after replacing periods in expression sample identifiers with hyphens.
- Ensembl version suffixes are removed.
- The historical Ensembl-to-gene-symbol mapping is reconstructed from a frozen manifest using the historical `org.Hs.eg.db`/`AnnotationDbi::mapIds(..., multiVals = "first")` policy.
- Historically unmapped rows, duplicate symbols after the first occurrence, and 1,143 zero-variance genes are removed.
- The final matrix contains 35,009 genes across 1,106 primary tumours.
- The production preprocessing script does not add an expression transformation; it assumes that the archived matrix is already on the continuous scale used by the historical analysis.

### TME profiling and immune-activity axis

- ESTIMATE is applied to a frozen 9,859-gene common-gene input using `estimateScore(..., platform = "illumina")`; ImmuneScore, StromalScore, ESTIMATEScore, and estimated tumour purity are retained.
- Ten MCP-counter population scores and 67 xCell enrichment outputs are generated through the IOBR implementation layer (`IOBR::deconvo_mcpcounter()` and `IOBR::deconvo_xcell()`). xCell outputs are interpreted as enrichment scores, not direct cell counts.
- Seventeen prespecified signatures are scored with `IOBR::calculate_sig_score(..., method = "ssgsea")`. The signatures include immune, stromal, and epithelial–mesenchymal-transition programmes derived from the IOBR collection and the cited Ayers, Mariathasan, Jiang, and Bindea resources.
- The frozen signature names are `CD_8_T_effector`, `IFNG_signature_Ayers_et_al`, `T_cell_inflamed_GEP_Ayers_et_al`, `CAF_Peng_et_al`, `TAM_Peng_et_al`, `MDSC_Peng_et_al`, `Immune_Checkpoint`, `EMT1`, `EMT2`, `EMT3`, `Pan_F_TBRs`, `CD8_T_cells_Bindea_et_al`, `NK_cells_Bindea_et_al`, `Cytotoxic_cells_Bindea_et_al`, `Macrophages_Bindea_et_al`, `DC_Bindea_et_al`, and `B_cells_Bindea_et_al`.
- The 17 scores are centred and scaled, and PCA is performed with `stats::prcomp(center = TRUE, scale. = TRUE)`. PC1 is deterministically oriented so that higher values represent greater aggregate signature activity.
- PC1 is divided into type-7 tertiles, producing operational Immune-Low, Immune-Mid, and Immune-High strata. These are ordered strata of a continuous axis, not unsupervised molecular clusters.
- Age is compared with Kruskal–Wallis tests; pathological categories are compared with Pearson chi-square tests; TME features use Kruskal–Wallis and Dunn post-hoc tests with Benjamini–Hochberg (BH) correction. For xCell, BH correction is applied across all 67 omnibus tests and post-hoc tests are restricted to omnibus FDR-significant features.

### Differential expression and pathway analysis

- Differential expression compares Immune-High (n=369) with Immune-Low (n=369); Immune-Mid is excluded from this contrast.
- `limma::lmFit()`, `makeContrasts()`, `contrasts.fit()`, `eBayes()`, and `topTable(adjust.method = "BH")` are used on the 35,009-gene matrix.
- Differentially expressed genes are defined as BH-FDR < 0.05 and `|logFC| > 1`, with positive logFC denoting higher expression in Immune-High tumours.
- GO Biological Process and KEGG over-representation analyses use the combined significant DEG set and an explicit tested-gene background. These analyses are non-directional because positive and negative DEGs are analysed together.
- Hallmark GSEA uses the complete gene list ranked by Immune-High minus Immune-Low limma logFC. The frozen Hallmark collection is MSigDB 2026.1.Hs (50 gene sets; 7,322 pathway–gene pairs).
- The production enrichment module deliberately preserves the validated Hallmark GSEA result rather than rerunning a stochastic/backend-dependent analysis. It also preserves the corrected historical KEGG result because the historical online KEGG annotation has drifted; GO-BP is deterministically recomputed and numerically regression-tested.

### Protein-interaction network and hub prioritization

- Leading-edge genes from five selected immune-related Hallmarks are pooled and deduplicated to form a 455-gene immune-focused input.
- The preserved STRING network contains 389 represented nodes and 4,044 unique undirected interactions with combined scores from 0.700 to 0.999. One identifier reconciliation is recorded (`RIGI → DDX58`); 66 input genes are absent from the exported network.
- CytoHubba Maximal Clique Centrality (MCC) is used to retain the top 20 connected genes.
- The network was preserved from the validated historical artifact; the supplied production script does not query STRING again, rerun Cytoscape, or recalculate MCC. Hub status therefore means connectivity within this selected immune-focused network, not causal regulation or genome-wide master-regulator activity.
- Pearson correlations between the 20 hubs and 14 immune/TME features are evaluated across all 1,106 TCGA tumours, with BH correction across the complete 280-test family.

### Survival and CCL19–IL10 analyses

- Overall survival uses days to death for deceased patients and days to last follow-up for living patients. Valid positive survival time and recorded vital status are required.
- The valid TCGA survival cohort contains 1,083 patients, 150 deaths, and 933 censored observations. Immune-activity strata contain 361 patients each in this survival subset.
- Kaplan–Meier and log-rank analyses compare the three immune-activity strata. Each of the 20 hub genes is then dichotomized at its exact TCGA survival-cohort median, with values equal to the median assigned to High, and analysed in a univariate Cox model. BH correction is applied across all 20 hub tests.
- CCL19 and IL10 are evaluated together in a stage-stratified Cox model with age. The four joint states are High CCL19/Low IL10, High CCL19/High IL10, Low CCL19/Low IL10, and Low CCL19/High IL10.
- Joint states are compared using overall and pairwise log-rank tests, with BH correction across the six pairwise comparisons, followed by a stage-stratified Cox model. No `CCL19_group * IL10_group` interaction term is fitted; the four groups are joint expression states, not proof of statistical interaction.
- Six predefined TME features are compared across the four joint states in all 1,106 TCGA tumours. The state thresholds are derived from the 1,083-patient survival cohort and are reused for this descriptive/contexture analysis.

### SCAN-B/GSE96058 evaluation

- CCL19 and IL10 are classified using independent, cohort-specific exact medians calculated across all 3,069 SCAN-B tumours. TCGA numerical cutoffs are not transferred.
- Individual-gene and four-state Kaplan–Meier/log-rank analyses are performed, followed by BH-adjusted pairwise state comparisons.
- The adjusted Cox model includes joint state, age, lymph-node status, Nottingham histological grade, and tumour size; complete cases are used and proportional-hazards assumptions are checked with `survival::cox.zph()`.
- SCAN-B is treated as an independent evaluation of selected prognostic findings, not as validation of the complete TCGA immune-axis and pathway-discovery workflow.

## Key findings

1. The PCA-derived axis captured substantial TME biology. PC1 explained 48.1% of the variance in the 17 signature scores, and its tertiles showed coordinated differences in ESTIMATE, MCP-counter, xCell, and immune-signature measurements.
2. Immune-High versus Immune-Low tumours differed broadly at the transcriptomic level: 2,830 genes met the prespecified DEG criteria, and 19 Hallmark pathways were BH-FDR significant, all enriched toward Immune-High.
3. The broad immune-activity strata were not independently prognostic in TCGA. Overall survival did not differ significantly across the three strata (log-rank P=0.2407).
4. In the 20-hub survival screen, CCL19 and CXCL2 were favourable and IL10 was adverse after BH correction; the manuscript subsequently focuses on CCL19 and IL10. In the adjusted TCGA model, high CCL19 was associated with lower mortality (HR=0.592, P=0.00415) and high IL10 with higher mortality (HR=1.653, P=0.00591).
5. The TCGA Low-CCL19/High-IL10 state showed the least favourable survival in the four-state analysis (overall log-rank P=2.80×10⁻⁶; adjusted HR=2.445 versus High-CCL19/Low-IL10, P=0.00129). This is an associative, exploratory joint-state result; it is not a causal mechanism or a formal interaction test.
6. In SCAN-B, the favourable CCL19 association was reproduced (log-rank P=2.77×10⁻⁶). IL10 had the same adverse direction but was not statistically significant (P=0.0614). The unadjusted four-state pattern was reproduced, but the adjusted Low-CCL19/High-IL10 contrast was attenuated and non-significant (HR=1.319, P=0.165).

These findings support a cautious interpretation: broad immune activity is biologically informative but is not a validated prognostic classifier in this analysis. CCL19 is the most reproducible individual signal. IL10 and the CCL19–IL10 joint state remain context-dependent and require further independent, mechanistic, spatial, and prospective evaluation.

## Repository contents

The supplied code expects the following project-root layout. The `R/` directory name is intentional because the scripts call `source("R/00_setup.R")`.

```text
TCGA-BRCA-immune-activity/
├── README.md
├── R/
│   ├── 00_setup.R
│   ├── 01_import_tcga.R
│   ├── 02_preprocess_tcga.R
│   ├── 03_estimate.R
│   ├── 04_mcp_xcell.R
│   ├── 05_ssgsea_immune_axis.R
│   ├── 06_phenotype_characterization.R
│   ├── 07_differential_expression.R
│   ├── 08_functional_enrichment.R
│   ├── 09_ppi_hubs.R
│   ├── 10_tcga_survival.R
│   ├── 11_tcga_ccl19_il10.R
│   ├── 12_tcga_supporting_characterization.R
│   └── 13_scanb_validation.R
├── data_raw/                    # Source-derived inputs; normally not committed
├── data_processed/              # Generated RDS checkpoints; normally gitignored
├── external/                    # Frozen manifests and validated reference objects
│   ├── annotation/
│   ├── estimate/
│   ├── iobr/
│   ├── differential_expression/
│   ├── functional_enrichment/
│   ├── ppi/
│   ├── survival/
│   ├── tcga_support/
│   ├── scanb/
│   └── cytoscape/
├── results/
│   ├── figures/
│   ├── supplementary_figures/
│   ├── tables/
│   │   ├── phenotype_characterization/
│   │   ├── differential_expression/
│   │   ├── functional_enrichment/
│   │   ├── ppi_hubs/
│   │   ├── tcga_survival/
│   │   ├── tcga_ccl19_il10/
│   │   ├── tcga_supporting_characterization/
│   │   └── scanb_validation/
│   └── supplementary_tables/
├── logs/                        # QC summaries, hashes, and session information
├── docs/                        # Data provenance and run instructions
├── archive/                     # Optional historical materials; not required for production analysis
└── environment/                 # renv.lock, installation notes, and environment metadata
```

The generated `results/tables/` outputs should be preferred over uploading the complete manuscript tables as a Word document. Main and supplementary figures may be included as static outputs, but each figure should have a documented source table/object and generation script before it is described as reproducible.

## Reproducibility

### Required software and packages

The supplied scripts explicitly require or use:

- R (exact version not included in the supplied files);
- `IOBR` and `GSVA` for MCP-counter, xCell, and ssGSEA implementation;
- `estimate` for ESTIMATE scoring;
- `limma` for differential expression;
- `FSA` for Dunn post-hoc tests;
- `clusterProfiler` and `org.Hs.eg.db` for enrichment analysis;
- `survival` for Kaplan–Meier, log-rank, Cox, and proportional-hazards analyses;
- `survminer` and `dplyr` for SCAN-B analysis; and
- `Hmisc` for Pearson correlation calculations.

The manuscript also names `msigdbr` among the key packages, but the supplied production scripts preserve the Hallmark result from a frozen MSigDB 2026.1.Hs object rather than obtaining a fresh collection during the run. The historical role and version of `msigdbr` should therefore be documented explicitly.

### Run order

Run from the project root after placing the required inputs and frozen reference objects in their expected locations:

```r
source("R/00_setup.R")
source("R/01_import_tcga.R")
source("R/02_preprocess_tcga.R")
source("R/03_estimate.R")
source("R/04_mcp_xcell.R")
source("R/05_ssgsea_immune_axis.R")
source("R/06_phenotype_characterization.R")
source("R/07_differential_expression.R")
source("R/08_functional_enrichment.R")
source("R/09_ppi_hubs.R")
source("R/10_tcga_survival.R")
source("R/11_tcga_ccl19_il10.R")
source("R/12_tcga_supporting_characterization.R")
source("R/13_scanb_validation.R")
```

The scripts do not install packages automatically. They create and check `data_processed/`, `results/`, and `logs/` outputs and write `sessionInfo()` files at runtime. A clean public release should add an `renv.lock` file or an equivalent installation script and commit representative session information.

### What is regression-tested

The production code contains assertions and frozen-reference comparisons for sample matching, gene mapping, expression dimensions, ESTIMATE/MCP-counter/xCell outputs, ssGSEA values, PCA orientation and tertiles, DEG results, enrichment objects, PPI node/edge structure, survival cohorts, Cox models, and SCAN-B results. These checks are useful safeguards, but they depend on the reference files listed in the scripts.

### Current reproducibility boundary

The public repository documents the production workflow, validation logic, and selected reproducible outputs. Some restricted reference objects are intentionally not distributed. Where required, provenance information, checksums, and source locations are provided.The repository is designed as a transparent computational record rather than a fully self-contained reproduction package. Reproduction of specific analyses requires access to the documented source datasets and reference objects.

## Data sharing and privacy

The underlying cohorts are publicly available and de-identified, but this repository should not redistribute patient-level matrices by default. The public release should provide:

- cohort names and identifiers (`TCGA-BRCA` and `GSE96058`/SCAN-B);
- source and access instructions;
- the preprocessing and sample-selection code;
- frozen annotation/manifests and small, redistributable reference objects where their terms permit;
- summary tables, QC logs, and selected figures; and
- checksums or provenance records for inputs that are not stored in Git.

Complete expression matrices, clinical tables, raw FASTQ files, or any restricted/non-redistributable object should remain outside GitHub unless redistribution rights and privacy review have been documented. The supplied scripts currently require local copies of these files and do not download them automatically.

## Main outputs to include

When the repository is populated, the most useful public outputs are:

- main Figures 1–8 and selected supplementary figures;
- main Tables 1–4;
- CSV summaries for phenotype characterization, DEGs, enrichment, PPI hubs, TCGA survival, CCL19–IL10 states, supporting analyses, and SCAN-B evaluation;
- method-provenance and QC tables generated by modules 00–13; and
- session-information and checksum logs.

Figures should be referenced from the README only after their filenames, provenance, and generation scripts are stable. The complete manuscript PDF is not required for the initial computational repository and is intentionally not treated as a required input here.

## Limitations and interpretation

- This repository contains analysis scripts, selected processed outputs, validation scripts, and reproducibility documentation. Raw datasets and restricted reference objects are available from original repositories or upon request where applicable.
- This is a retrospective observational analysis; associations do not establish causality.
- Bulk expression cannot identify the cellular source or spatial organization of CCL19 or IL10.
- ESTIMATE, MCP-counter, xCell, and ssGSEA are transcriptome-derived computational measurements, not orthogonal cell-count assays.
- PC1 tertiles and gene-level median splits simplify continuous information and are not clinically established thresholds.
- The PPI network is restricted to selected Hallmark leading-edge genes and is not an unbiased whole-transcriptome regulatory network.
- The joint CCL19–IL10 groups are not a formal interaction test.
- SCAN-B used independent medians, so the external analysis evaluates direction and relative expression pattern rather than transportability of a fixed clinical cutoff.
- Discrimination, calibration, incremental predictive value, treatment interaction, and clinical utility were not evaluated.
- The adjusted joint-state contrast did not reach statistical significance in SCAN-B, and IL10 alone was also non-significant there.
- The supplied workflow analyses archived expression-level inputs; it is not an end-to-end raw-read RNA-seq pipeline.

## Manuscript status

The associated manuscript is in preparation. No journal publication, peer review, preprint DOI, or publication DOI should be inferred from this repository. This section can be updated later with a submission status, bioRxiv DOI, or journal DOI when those identifiers exist.

## Authors and acknowledgement

### Manuscript authors

- Seyed Mostafa Rahimi — Conceptualization, Methodology, Data Curation, Formal Analysis, Visualization, Investigation, Writing – Original Draft Preparation, Writing – Review & Editing
- Fatemeh Arab Mirrahmani — Methodological support, Writing – Review & Editing
- Mohsen Najafi — Supervision, Validation, Project oversight, Writing – Review & Editing

Affiliation: Department of Human Genetics, Negin Genetics & Pathobiology Laboratory, Negin Medical Complex, Sari, Iran.

The manuscript authors and GitHub collaborators are not assumed to be the same people. GitHub access should be managed separately and should not be inferred from manuscript authorship.

The authors acknowledge the TCGA and SCAN-B research groups for generating and making the datasets available.

## Citation

Please cite the associated manuscript when a citable version becomes available. Until then, this repository should be described as the computational analysis repository for a manuscript in preparation; no DOI is assigned here.

## License

The analysis code in this repository is released under the MIT License. Data and third-party reference artifacts remain subject to the terms and licenses of their original sources.
