############################################################
## 12_tcga_supporting_characterization.R
##
## PURPOSE
##
## Production implementation of the remaining validated
## TCGA-BRCA supporting analyses:
##
## A. Hub-gene × immune/TME correlations
## B. Immune-phenotype clinicopathological associations
## C. CCL19/IL10 joint-state immune-contexture analysis
##
##
## ANALYTICAL SCOPE
##
## A. HUB / IMMUNE CORRELATIONS
##
##   - 20 frozen MCC hub genes
##   - 14 immune/TME features
##   - 1,106 TCGA primary tumors
##   - Pearson correlation via Hmisc::rcorr
##   - 280 total tests
##   - BH correction globally across all 280 tests
##
##
## B. CLINICOPATHOLOGICAL ASSOCIATIONS
##
##   Age:
##     Kruskal-Wallis across Immune_Low / Mid / High
##
##   Stage, T, N, M:
##     Pearson chi-square tests
##
##
## C. JOINT-STATE IMMUNE CONTEXTURE
##
##   - all 1,106 primary tumors
##   - six immune-contexture features
##
## CRITICAL:
##
##   CCL19 and IL10 thresholds are NOT recalculated using
##   all 1,106 tumors.
##
##   The exact thresholds derived from the 1,083-patient
##   TCGA survival cohort in Module 11 are reused.
##
##   CCL19 = 8.5924570372680797
##   IL10  = 4.7548875021634682
##
##   Four states:
##
##     High_Low
##     High_High
##     Low_Low
##     Low_High
##
##   Omnibus:
##     Kruskal-Wallis for each of six features
##     BH across the six omnibus tests
##
##   Post-hoc:
##     FSA::dunnTest(..., method = "bh")
##     BH across the six pairwise contrasts WITHIN EACH
##     feature.
##
##
## INTERPRETATION GUARDRAILS
##
## - Hub/TME correlations are associations.
## - They do not establish causality.
##
## - Joint CCL19/IL10 groups are expression states.
## - No statistical interaction term is fitted here.
##
## - The 1,106-tumor contexture analysis is descriptive /
##   associative and distinct from the 1,083-patient
##   survival analysis.
############################################################


############################################################
## 12.1 — PROJECT SETUP
############################################################

if (!exists("PATHS")) {
  source("R/00_setup.R")
}


require_package(
  "Hmisc"
)


require_package(
  "FSA"
)


############################################################
## 12.2 — PRODUCTION INPUT FILES
############################################################

meta_file <- file.path(
  PATHS$processed,
  "05_tcga_metadata_immune_phenotype.rds"
)


expression_file <- file.path(
  PATHS$processed,
  "02_tcga_expression_final.rds"
)


hub_file <- file.path(
  PATHS$processed,
  "09_top20_MCC_hub_genes.rds"
)


survival_state_file <- file.path(
  PATHS$processed,
  "11_tcga_CCL19_IL10_joint_states.rds"
)


assert_file(
  meta_file
)


assert_file(
  expression_file
)


assert_file(
  hub_file
)


assert_file(
  survival_state_file
)


############################################################
## 12.3 — FROZEN REFERENCE FILES
##
## These were copied byte-for-byte from the historical
## FINAL_MANUSCRIPT_OUTPUTS directory during Pre-Module 12.
############################################################

reference_dir <- file.path(
  "external",
  "tcga_support"
)


S9_reference_file <- file.path(
  reference_dir,
  "TableS9_HubImmune_Correlations_FrozenObjects_FINAL.rds"
)


clinical_reference_file <- file.path(
  reference_dir,
  "TableS_TCGA_Clinical_Data_for_Table1.csv"
)


clinical_missing_reference_file <- file.path(
  reference_dir,
  "TableS_TCGA_Clinical_Missingness.csv"
)


S11D_reference_file <- file.path(
  reference_dir,
  "TableS11D_TCGA_JointState_Immune_Omnibus.csv"
)


S11E_reference_file <- file.path(
  reference_dir,
  "TableS11E_TCGA_JointState_Dunn_BH_AllComparisons.csv"
)


S11F_reference_file <- file.path(
  reference_dir,
  "TableS11F_TCGA_HighLow_vs_LowHigh_ImmuneComparison.csv"
)


S11G_reference_file <- file.path(
  reference_dir,
  "TableS11G_TCGA_JointState_ImmuneFeature_Medians.csv"
)


reference_files <- c(
  S9_reference_file,
  clinical_reference_file,
  clinical_missing_reference_file,
  S11D_reference_file,
  S11E_reference_file,
  S11F_reference_file,
  S11G_reference_file
)


for (f in reference_files) {
  assert_file(f)
}


############################################################
## 12.4 — VERIFY FROZEN REFERENCE MD5 VALUES
############################################################

expected_reference_md5 <- c(
  
  "TableS9_HubImmune_Correlations_FrozenObjects_FINAL.rds" =
    "beae6564a56bae2d1906eb4be11ce2e0",
  
  "TableS_TCGA_Clinical_Data_for_Table1.csv" =
    "78b00d7d862a8879186e0426a91847d8",
  
  "TableS_TCGA_Clinical_Missingness.csv" =
    "099342cf777ed52f606cc3080c2f4364",
  
  "TableS11D_TCGA_JointState_Immune_Omnibus.csv" =
    "56d4775635661ea33dc9989e3357b351",
  
  "TableS11E_TCGA_JointState_Dunn_BH_AllComparisons.csv" =
    "c851eb94e2d2f4d7f0ff40b939f55e11",
  
  "TableS11F_TCGA_HighLow_vs_LowHigh_ImmuneComparison.csv" =
    "33e8f8f057eb893d5b70405c6fb41058",
  
  "TableS11G_TCGA_JointState_ImmuneFeature_Medians.csv" =
    "9440c8da895972ef1a37e46de6f08966"
)


observed_reference_md5 <- unname(
  tools::md5sum(
    reference_files
  )
)


names(
  observed_reference_md5
) <- basename(
  reference_files
)


reference_md5_pass <- identical(
  observed_reference_md5[
    names(
      expected_reference_md5
    )
  ],
  expected_reference_md5
)


assert_true(
  reference_md5_pass,
  "One or more frozen Module-12 reference files failed MD5 verification."
)


############################################################
## 12.5 — LOAD PRODUCTION INPUTS
############################################################

meta05 <- readRDS(
  meta_file
)


expr_final <- readRDS(
  expression_file
)


hub09 <- readRDS(
  hub_file
)


tcga_state_survival <- readRDS(
  survival_state_file
)


############################################################
## 12.6 — LOAD FROZEN REFERENCES
############################################################

S9_ref <- readRDS(
  S9_reference_file
)


clinical_ref <- read.csv(
  clinical_reference_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


clinical_missing_ref <- read.csv(
  clinical_missing_reference_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


S11D_ref <- read.csv(
  S11D_reference_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


S11E_ref <- read.csv(
  S11E_reference_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


S11F_ref <- read.csv(
  S11F_reference_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


S11G_ref <- read.csv(
  S11G_reference_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


############################################################
## 12.7 — UPSTREAM INPUT QC
############################################################

assert_true(
  nrow(meta05) == 1106L,
  "Expected 1,106 TCGA primary tumors."
)


assert_true(
  nrow(expr_final) == 35009L,
  "Expected 35,009 genes in TCGA expression matrix."
)


assert_true(
  ncol(expr_final) == 1106L,
  "Expected 1,106 samples in TCGA expression matrix."
)


assert_true(
  nrow(hub09) == 20L,
  "Expected 20 Module-09 hub genes."
)


assert_true(
  nrow(tcga_state_survival) == 1083L,
  "Expected 1,083 patients in Module-11 survival-state object."
)


assert_unique(
  meta05$sample,
  "Module-05 samples"
)


assert_unique(
  colnames(expr_final),
  "TCGA expression samples"
)


assert_unique(
  rownames(expr_final),
  "TCGA expression genes"
)


assert_true(
  identical(
    colnames(expr_final),
    meta05$sample
  ),
  paste(
    "Module-02 expression and Module-05 metadata",
    "sample order differs."
  )
)


############################################################
############################################################
## PART A — HUB-GENE × IMMUNE/TME CORRELATIONS
############################################################
############################################################


############################################################
## 12.8 — EXACT HUB LIST
############################################################

S9_hub_genes <- as.character(
  hub09$Gene
)


assert_true(
  identical(
    S9_hub_genes,
    S9_ref$hub_genes
  ),
  "Module-09 hub list differs from frozen Table S9 hub list."
)


assert_true(
  length(S9_hub_genes) == 20L,
  "Exactly 20 hub genes are required."
)


assert_true(
  all(
    S9_hub_genes %in%
      rownames(expr_final)
  ),
  "One or more Table-S9 hub genes are absent from expression matrix."
)


############################################################
## 12.9 — EXACT 14 IMMUNE/TME FEATURES
############################################################

S9_immune_features <- c(
  "ImmuneScore",
  "StromalScore",
  "ESTIMATEScore",
  "TumorPurity",
  "CD_8_T_effector",
  "IFNG_signature_Ayers_et_al",
  "T_cell_inflamed_GEP_Ayers_et_al",
  "Immune_Checkpoint",
  "CAF_Peng_et_al",
  "TAM_Peng_et_al",
  "MDSC_Peng_et_al",
  "EMT1",
  "EMT2",
  "EMT3"
)


assert_true(
  identical(
    S9_immune_features,
    S9_ref$immune_features
  ),
  "Table-S9 immune/TME feature list differs from frozen reference."
)


assert_true(
  length(S9_immune_features) == 14L,
  "Exactly 14 immune/TME features are required."
)


assert_true(
  all(
    S9_immune_features %in%
      colnames(meta05)
  ),
  "One or more Table-S9 features are absent from Module-05 metadata."
)


############################################################
## 12.10 — BUILD ALIGNED CORRELATION INPUT
############################################################

S9_samples <- colnames(
  expr_final
)


S9_meta_index <- match(
  S9_samples,
  meta05$sample
)


assert_true(
  !anyNA(
    S9_meta_index
  ),
  "Unable to align Module-05 metadata to expression samples."
)


S9_hub_df <- as.data.frame(
  t(
    expr_final[
      S9_hub_genes,
      S9_samples,
      drop = FALSE
    ]
  )
)


S9_feature_df <- as.data.frame(
  meta05[
    S9_meta_index,
    S9_immune_features,
    drop = FALSE
  ]
)


rownames(
  S9_hub_df
) <- S9_samples


rownames(
  S9_feature_df
) <- S9_samples


assert_true(
  identical(
    rownames(S9_hub_df),
    rownames(S9_feature_df)
  ),
  "Hub-expression and immune-feature sample order differs."
)


############################################################
## 12.11 — PEARSON CORRELATIONS
##
## Historical method:
##
## Hmisc::rcorr(..., type = "pearson")
############################################################

S9_combined <- cbind(
  S9_hub_df,
  S9_feature_df
)


S9_corr <- Hmisc::rcorr(
  as.matrix(
    S9_combined
  ),
  type = "pearson"
)


S9_R <- S9_corr$r[
  S9_hub_genes,
  S9_immune_features,
  drop = FALSE
]


S9_P <- S9_corr$P[
  S9_hub_genes,
  S9_immune_features,
  drop = FALSE
]


S9_N <- S9_corr$n[
  S9_hub_genes,
  S9_immune_features,
  drop = FALSE
]


############################################################
## 12.12 — BH ACROSS ALL 280 TESTS
############################################################

S9_FDR <- matrix(
  
  stats::p.adjust(
    as.vector(
      S9_P
    ),
    method = "BH"
  ),
  
  nrow =
    nrow(S9_P),
  
  ncol =
    ncol(S9_P),
  
  dimnames =
    dimnames(S9_P)
)


assert_true(
  length(S9_R) == 280L,
  "Table S9 must contain 280 Pearson correlations."
)


############################################################
## 12.13 — REGRESS S9 AGAINST FROZEN MATRICES
############################################################

S9_R_max_diff <- max(
  abs(
    S9_R -
      S9_ref$R_matrix
  ),
  na.rm = TRUE
)


S9_P_max_diff <- max(
  abs(
    S9_P -
      S9_ref$raw_P_matrix
  ),
  na.rm = TRUE
)


S9_FDR_max_diff <- max(
  abs(
    S9_FDR -
      S9_ref$FDR_matrix
  ),
  na.rm = TRUE
)


S9_N_pass <- identical(
  S9_N,
  S9_ref$N_matrix
)


S9_regression_pass <- all(
  S9_R_max_diff < 1e-12,
  S9_P_max_diff < 1e-12,
  S9_FDR_max_diff < 1e-12,
  S9_N_pass
)


assert_true(
  S9_regression_pass,
  "Table S9 numerical regression failed."
)


############################################################
## 12.14 — CREATE COMPLETE S9 LONG TABLE
############################################################

S9_long <- expand.grid(
  
  Gene =
    rownames(
      S9_R
    ),
  
  Feature =
    colnames(
      S9_R
    ),
  
  KEEP.OUT.ATTRS =
    FALSE,
  
  stringsAsFactors =
    FALSE
)


S9_long$Pearson_r <- as.vector(
  S9_R
)


S9_long$Raw_P <- as.vector(
  S9_P
)


S9_long$FDR_BH <- as.vector(
  S9_FDR
)


S9_long$N <- as.vector(
  S9_N
)


S9_long$Direction <- ifelse(
  S9_long$Pearson_r > 0,
  "Positive",
  ifelse(
    S9_long$Pearson_r < 0,
    "Negative",
    "Zero"
  )
)


S9_long$FDR_significant_0_05 <-
  S9_long$FDR_BH < 0.05


############################################################
## 12.15 — S9 QC CHECKPOINTS
############################################################

S9_significant <- S9_FDR < 0.05


S9_positive_significant <- sum(
  S9_significant &
    S9_R > 0
)


S9_negative_significant <- sum(
  S9_significant &
    S9_R < 0
)


assert_true(
  sum(
    S9_significant
  ) == 277L,
  "Expected 277 FDR-significant hub/TME correlations."
)


assert_true(
  S9_positive_significant == 257L,
  "Expected 257 positive FDR-significant correlations."
)


assert_true(
  S9_negative_significant == 20L,
  "Expected 20 negative FDR-significant correlations."
)


EXPECTED_S9_MIN_R <-
  -0.753321578238781


EXPECTED_S9_MAX_R <-
  0.912999711568593


assert_true(
  abs(
    min(S9_R) -
      EXPECTED_S9_MIN_R
  ) < 1e-12,
  "Minimum frozen Pearson correlation differs."
)


assert_true(
  abs(
    max(S9_R) -
      EXPECTED_S9_MAX_R
  ) < 1e-12,
  "Maximum frozen Pearson correlation differs."
)


############################################################
############################################################
## PART B — CLINICOPATHOLOGICAL ASSOCIATIONS
############################################################
############################################################


############################################################
## 12.16 — BUILD CLINICAL TABLE
############################################################

clinical_table <- data.frame(
  
  sample =
    meta05$sample,
  
  ImmunePhenotype =
    factor(
      meta05$ImmunePhenotype,
      levels = c(
        "Immune_Low",
        "Immune_Mid",
        "Immune_High"
      )
    ),
  
  Age =
    meta05$
    age_at_earliest_diagnosis_in_years.diagnoses.xena_derived,
  
  stringsAsFactors =
    FALSE
)


############################################################
## 12.17 — CLEAN PATHOLOGICAL STAGE
############################################################

stage_raw <- as.character(
  meta05$ajcc_pathologic_stage.diagnoses
)


stage_clean <- rep(
  NA_character_,
  length(stage_raw)
)


stage_clean[
  !is.na(stage_raw) &
    grepl(
      "Stage IV",
      stage_raw,
      ignore.case = TRUE
    )
] <- "Stage IV"


stage_clean[
  is.na(stage_clean) &
    !is.na(stage_raw) &
    grepl(
      "Stage III",
      stage_raw,
      ignore.case = TRUE
    )
] <- "Stage III"


stage_clean[
  is.na(stage_clean) &
    !is.na(stage_raw) &
    grepl(
      "Stage II",
      stage_raw,
      ignore.case = TRUE
    )
] <- "Stage II"


stage_clean[
  is.na(stage_clean) &
    !is.na(stage_raw) &
    grepl(
      "Stage I",
      stage_raw,
      ignore.case = TRUE
    )
] <- "Stage I"


############################################################
## 12.18 — CLEAN T STAGE
############################################################

T_raw <- as.character(
  meta05$ajcc_pathologic_t.diagnoses
)


T_clean <- rep(
  NA_character_,
  length(T_raw)
)


for (lv in 1:4) {
  
  hit <- !is.na(T_raw) &
    grepl(
      paste0("^T", lv),
      T_raw,
      ignore.case = TRUE
    )
  
  T_clean[hit] <-
    paste0(
      "T",
      lv
    )
}


############################################################
## 12.19 — CLEAN N STAGE
############################################################

N_raw <- as.character(
  meta05$ajcc_pathologic_n.diagnoses
)


N_clean <- rep(
  NA_character_,
  length(N_raw)
)


for (lv in 0:3) {
  
  hit <- !is.na(N_raw) &
    grepl(
      paste0("^N", lv),
      N_raw,
      ignore.case = TRUE
    )
  
  N_clean[hit] <-
    paste0(
      "N",
      lv
    )
}


############################################################
## 12.20 — CLEAN M STAGE
############################################################

M_raw <- as.character(
  meta05$ajcc_pathologic_m.diagnoses
)


M_clean <- rep(
  NA_character_,
  length(M_raw)
)


M_clean[
  !is.na(M_raw) &
    grepl(
      "^M0",
      M_raw,
      ignore.case = TRUE
    )
] <- "M0"


M_clean[
  !is.na(M_raw) &
    grepl(
      "^M1",
      M_raw,
      ignore.case = TRUE
    )
] <- "M1"


############################################################
## 12.21 — ATTACH CLINICAL FACTORS
############################################################

clinical_table$Stage <- factor(
  stage_clean,
  levels = c(
    "Stage I",
    "Stage II",
    "Stage III",
    "Stage IV"
  )
)


clinical_table$Tstage <- factor(
  T_clean,
  levels = c(
    "T1",
    "T2",
    "T3",
    "T4"
  )
)


clinical_table$Nstage <- factor(
  N_clean,
  levels = c(
    "N0",
    "N1",
    "N2",
    "N3"
  )
)


clinical_table$Mstage <- factor(
  M_clean,
  levels = c(
    "M0",
    "M1"
  )
)


############################################################
## 12.22 — REGRESS CLINICAL TABLE AGAINST FROZEN REFERENCE
############################################################

clinical_reference_index <- match(
  clinical_table$sample,
  clinical_ref$sample
)


assert_true(
  !anyNA(
    clinical_reference_index
  ),
  "Unable to align clinical frozen reference."
)


clinical_ref_aligned <- clinical_ref[
  clinical_reference_index,
  ,
  drop = FALSE
]


clinical_sample_pass <- identical(
  clinical_table$sample,
  clinical_ref_aligned$sample
)


clinical_phenotype_pass <- identical(
  as.character(
    clinical_table$ImmunePhenotype
  ),
  as.character(
    clinical_ref_aligned$ImmunePhenotype
  )
)


clinical_age_pass <- all(
  is.na(
    clinical_table$Age
  ) ==
    is.na(
      clinical_ref_aligned$Age
    )
) &&
  max(
    abs(
      clinical_table$Age -
        clinical_ref_aligned$Age
    ),
    na.rm = TRUE
  ) < 1e-12


clinical_stage_pass <- identical(
  as.character(
    clinical_table$Stage
  ),
  as.character(
    clinical_ref_aligned$Stage
  )
)


clinical_T_pass <- identical(
  as.character(
    clinical_table$Tstage
  ),
  as.character(
    clinical_ref_aligned$Tstage
  )
)


clinical_N_pass <- identical(
  as.character(
    clinical_table$Nstage
  ),
  as.character(
    clinical_ref_aligned$Nstage
  )
)


clinical_M_pass <- identical(
  as.character(
    clinical_table$Mstage
  ),
  as.character(
    clinical_ref_aligned$Mstage
  )
)


clinical_data_regression_pass <- all(
  clinical_sample_pass,
  clinical_phenotype_pass,
  clinical_age_pass,
  clinical_stage_pass,
  clinical_T_pass,
  clinical_N_pass,
  clinical_M_pass
)


assert_true(
  clinical_data_regression_pass,
  "Clinical-data reconstruction differs from frozen reference."
)


############################################################
## 12.23 — CLINICAL MISSINGNESS
############################################################

clinical_missingness <- data.frame(
  
  Variable = c(
    "Age",
    "Stage",
    "Tstage",
    "Nstage",
    "Mstage"
  ),
  
  Missing_n = c(
    sum(
      is.na(
        clinical_table$Age
      )
    ),
    sum(
      is.na(
        clinical_table$Stage
      )
    ),
    sum(
      is.na(
        clinical_table$Tstage
      )
    ),
    sum(
      is.na(
        clinical_table$Nstage
      )
    ),
    sum(
      is.na(
        clinical_table$Mstage
      )
    )
  ),
  
  stringsAsFactors =
    FALSE
)


clinical_missingness$Missing_pct <-
  100 *
  clinical_missingness$Missing_n /
  nrow(
    clinical_table
  )


assert_true(
  identical(
    clinical_missingness$Missing_n,
    as.integer(
      clinical_missing_ref$Missing_n
    )
  ),
  "Clinical missingness counts differ from frozen reference."
)


assert_true(
  max(
    abs(
      clinical_missingness$Missing_pct -
        clinical_missing_ref$Missing_pct
    )
  ) < 1e-10,
  "Clinical missingness percentages differ from frozen reference."
)


############################################################
## 12.24 — CLINICOPATHOLOGICAL TESTS
############################################################

age_test <- stats::kruskal.test(
  Age ~ ImmunePhenotype,
  data =
    clinical_table
)


stage_test <- stats::chisq.test(
  table(
    clinical_table$ImmunePhenotype,
    clinical_table$Stage
  )
)


T_test <- stats::chisq.test(
  table(
    clinical_table$ImmunePhenotype,
    clinical_table$Tstage
  )
)


N_test <- stats::chisq.test(
  table(
    clinical_table$ImmunePhenotype,
    clinical_table$Nstage
  )
)


M_test <- stats::chisq.test(
  table(
    clinical_table$ImmunePhenotype,
    clinical_table$Mstage
  )
)


############################################################
## 12.25 — EXACT CLINICAL P-VALUE CHECKPOINTS
############################################################

EXPECTED_AGE_P <-
  0.0017490442794000139


EXPECTED_STAGE_P <-
  0.1462402476481697


EXPECTED_T_P <-
  0.12459611798652359


EXPECTED_N_P <-
  0.0038617759651983987


EXPECTED_M_P <-
  0.29112373558874149


assert_true(
  abs(
    age_test$p.value -
      EXPECTED_AGE_P
  ) < 1e-12,
  "Age Kruskal-Wallis P-value differs from validated result."
)


assert_true(
  abs(
    stage_test$p.value -
      EXPECTED_STAGE_P
  ) < 1e-12,
  "Stage chi-square P-value differs from validated result."
)


assert_true(
  abs(
    T_test$p.value -
      EXPECTED_T_P
  ) < 1e-12,
  "T-stage chi-square P-value differs from validated result."
)


assert_true(
  abs(
    N_test$p.value -
      EXPECTED_N_P
  ) < 1e-12,
  "N-stage chi-square P-value differs from validated result."
)


assert_true(
  abs(
    M_test$p.value -
      EXPECTED_M_P
  ) < 1e-12,
  "M-stage chi-square P-value differs from validated result."
)


############################################################
## 12.26 — CLINICAL TEST RESULTS TABLE
############################################################

clinical_associations <- data.frame(
  
  Variable = c(
    "Age",
    "Stage",
    "Tstage",
    "Nstage",
    "Mstage"
  ),
  
  Method = c(
    "Kruskal-Wallis",
    "Pearson chi-square",
    "Pearson chi-square",
    "Pearson chi-square",
    "Pearson chi-square"
  ),
  
  Statistic = c(
    as.numeric(
      age_test$statistic
    ),
    as.numeric(
      stage_test$statistic
    ),
    as.numeric(
      T_test$statistic
    ),
    as.numeric(
      N_test$statistic
    ),
    as.numeric(
      M_test$statistic
    )
  ),
  
  df = c(
    as.numeric(
      age_test$parameter
    ),
    as.numeric(
      stage_test$parameter
    ),
    as.numeric(
      T_test$parameter
    ),
    as.numeric(
      N_test$parameter
    ),
    as.numeric(
      M_test$parameter
    )
  ),
  
  P_value = c(
    age_test$p.value,
    stage_test$p.value,
    T_test$p.value,
    N_test$p.value,
    M_test$p.value
  ),
  
  Significant_0_05 = c(
    age_test$p.value < 0.05,
    stage_test$p.value < 0.05,
    T_test$p.value < 0.05,
    N_test$p.value < 0.05,
    M_test$p.value < 0.05
  ),
  
  stringsAsFactors =
    FALSE
)


############################################################
############################################################
## PART C — CCL19 / IL10 JOINT-STATE IMMUNE CONTEXTURE
############################################################
############################################################


############################################################
## 12.27 — RECOVER MODULE-11 FROZEN CUTOFFS
##
## These MUST remain based on the 1,083-patient survival
## cohort.
############################################################

CCL19_cutoff <- median(
  tcga_state_survival$CCL19_expr
)


IL10_cutoff <- median(
  tcga_state_survival$IL10_expr
)


EXPECTED_CCL19_CUTOFF <-
  8.5924570372680797


EXPECTED_IL10_CUTOFF <-
  4.7548875021634682


assert_true(
  abs(
    CCL19_cutoff -
      EXPECTED_CCL19_CUTOFF
  ) < 1e-12,
  "CCL19 frozen cutoff differs from Module-11 cutoff."
)


assert_true(
  abs(
    IL10_cutoff -
      EXPECTED_IL10_CUTOFF
  ) < 1e-12,
  "IL10 frozen cutoff differs from Module-11 cutoff."
)


############################################################
## 12.28 — APPLY FROZEN CUTOFFS TO ALL 1,106 TUMORS
############################################################

all_sample_index <- match(
  meta05$sample,
  colnames(
    expr_final
  )
)


assert_true(
  !anyNA(
    all_sample_index
  ),
  "Unable to align all 1,106 tumors to expression matrix."
)


tcga_state_all <- meta05


tcga_state_all$CCL19_expr <- as.numeric(
  expr_final[
    "CCL19",
    all_sample_index
  ]
)


tcga_state_all$IL10_expr <- as.numeric(
  expr_final[
    "IL10",
    all_sample_index
  ]
)


tcga_state_all$CCL19_group <- factor(
  
  ifelse(
    tcga_state_all$CCL19_expr >=
      CCL19_cutoff,
    "High",
    "Low"
  ),
  
  levels = c(
    "Low",
    "High"
  )
)


tcga_state_all$IL10_group <- factor(
  
  ifelse(
    tcga_state_all$IL10_expr >=
      IL10_cutoff,
    "High",
    "Low"
  ),
  
  levels = c(
    "Low",
    "High"
  )
)


tcga_state_all$ImmuneState <- factor(
  
  paste(
    tcga_state_all$CCL19_group,
    tcga_state_all$IL10_group,
    sep = "_"
  ),
  
  levels = c(
    "High_Low",
    "High_High",
    "Low_Low",
    "Low_High"
  )
)


############################################################
## 12.29 — ALL-1106 STATE COUNTS
############################################################

state_counts <- table(
  tcga_state_all$ImmuneState
)


assert_true(
  identical(
    as.integer(
      state_counts
    ),
    c(
      196L,
      354L,
      337L,
      219L
    )
  ),
  "All-1106 CCL19/IL10 state counts differ from validated result."
)


############################################################
## 12.30 — EXACT SIX IMMUNE-CONTEXTURE FEATURES
############################################################

immune_state_features <- c(
  "ImmuneScore",
  "CD8_T_cells_MCPcounter",
  "Cytotoxic_lymphocytes_MCPcounter",
  "Immune_Checkpoint",
  "CAF_Peng_et_al",
  "TAM_Peng_et_al"
)


feature_labels <- c(
  
  ImmuneScore =
    "ESTIMATE ImmuneScore",
  
  CD8_T_cells_MCPcounter =
    "CD8 T cells",
  
  Cytotoxic_lymphocytes_MCPcounter =
    "Cytotoxic lymphocytes",
  
  Immune_Checkpoint =
    "Immune checkpoint signature",
  
  CAF_Peng_et_al =
    "CAF signature",
  
  TAM_Peng_et_al =
    "TAM signature"
)


assert_true(
  all(
    immune_state_features %in%
      colnames(
        tcga_state_all
      )
  ),
  "One or more Figure-7C immune features are unavailable."
)


############################################################
## 12.31 — KRUSKAL-WALLIS OMNIBUS TESTS
############################################################

state_kw_list <- vector(
  "list",
  length(
    immune_state_features
  )
)


for (
  i in seq_along(
    immune_state_features
  )
) {
  
  feature <- immune_state_features[i]
  
  
  kw_test <- stats::kruskal.test(
    tcga_state_all[[feature]] ~
      tcga_state_all$ImmuneState
  )
  
  
  state_kw_list[[i]] <- data.frame(
    
    Internal_Feature =
      feature,
    
    Feature =
      unname(
        feature_labels[
          feature
        ]
      ),
    
    Kruskal_Wallis_Statistic =
      as.numeric(
        kw_test$statistic
      ),
    
    df =
      as.numeric(
        kw_test$parameter
      ),
    
    Kruskal_Wallis_P =
      kw_test$p.value,
    
    stringsAsFactors =
      FALSE
  )
}


state_kw <- do.call(
  rbind,
  state_kw_list
)


rownames(
  state_kw
) <- NULL


############################################################
## BH family = six omnibus feature tests
############################################################

state_kw$BH_FDR <- stats::p.adjust(
  state_kw$Kruskal_Wallis_P,
  method = "BH"
)


############################################################
## 12.32 — REGRESS OMNIBUS RESULTS
############################################################

S11D_index <- match(
  state_kw$Feature,
  S11D_ref$Feature
)


assert_true(
  !anyNA(
    S11D_index
  ),
  "Unable to align Figure-7C omnibus frozen reference."
)


S11D_aligned <- S11D_ref[
  S11D_index,
  ,
  drop = FALSE
]


kw_logP_max_diff <- max(
  abs(
    log10(
      state_kw$Kruskal_Wallis_P
    ) -
      log10(
        S11D_aligned$Kruskal_Wallis_P
      )
  )
)


kw_logFDR_max_diff <- max(
  abs(
    log10(
      state_kw$BH_FDR
    ) -
      log10(
        S11D_aligned$BH_FDR
      )
  )
)


kw_regression_pass <- all(
  kw_logP_max_diff < 1e-10,
  kw_logFDR_max_diff < 1e-10
)


assert_true(
  kw_regression_pass,
  "Figure-7C Kruskal-Wallis/BH regression failed."
)


############################################################
## 12.33 — JOINT-STATE MEDIANS
############################################################

state_levels <- c(
  "High_Low",
  "High_High",
  "Low_Low",
  "Low_High"
)


state_display_labels <- c(
  
  High_Low =
    "CCL19 High / IL10 Low",
  
  High_High =
    "CCL19 High / IL10 High",
  
  Low_Low =
    "CCL19 Low / IL10 Low",
  
  Low_High =
    "CCL19 Low / IL10 High"
)


state_medians <- data.frame(
  
  ImmuneState =
    unname(
      state_display_labels[
        state_levels
      ]
    ),
  
  N =
    as.integer(
      state_counts[
        state_levels
      ]
    ),
  
  stringsAsFactors =
    FALSE
)


median_output_names <- c(
  "ESTIMATE ImmuneScore",
  "CD8 T cells",
  "Cytotoxic lymphocytes",
  "Immune checkpoint signature",
  "CAF signature",
  "TAM signature"
)


for (
  j in seq_along(
    immune_state_features
  )
) {
  
  feature <- immune_state_features[j]
  
  
  values <- vapply(
    
    state_levels,
    
    function(state) {
      
      stats::median(
        tcga_state_all[
          tcga_state_all$ImmuneState ==
            state,
          feature
        ],
        na.rm = TRUE
      )
    },
    
    numeric(1)
  )
  
  
  state_medians[[median_output_names[j]]] <- values
}


############################################################
## 12.34 — REGRESS STATE MEDIANS
############################################################

S11G_index <- match(
  state_medians$ImmuneState,
  S11G_ref$ImmuneState
)


assert_true(
  !anyNA(
    S11G_index
  ),
  "Unable to align Figure-7C state-median frozen reference."
)


S11G_aligned <- S11G_ref[
  S11G_index,
  ,
  drop = FALSE
]


state_median_N_pass <- identical(
  as.integer(
    state_medians$N
  ),
  as.integer(
    S11G_aligned$N
  )
)


state_median_max_diff <- max(
  abs(
    as.matrix(
      state_medians[
        ,
        median_output_names,
        drop = FALSE
      ]
    ) -
      as.matrix(
        S11G_aligned[
          ,
          median_output_names,
          drop = FALSE
        ]
      )
  )
)


state_median_regression_pass <- all(
  state_median_N_pass,
  state_median_max_diff < 1e-10
)


assert_true(
  state_median_regression_pass,
  "Figure-7C state-median regression failed."
)


############################################################
## 12.35 — DUNN POST-HOC TESTS
##
## Historical method:
##
## FSA::dunnTest(..., method = "bh")
##
## BH is applied independently within each feature across
## that feature's six pairwise comparisons.
############################################################

dunn_list <- vector(
  "list",
  length(
    immune_state_features
  )
)


for (
  i in seq_along(
    immune_state_features
  )
) {
  
  feature <- immune_state_features[i]
  
  
  invisible(
    utils::capture.output(
      dunn_result <- FSA::dunnTest(
        tcga_state_all[[feature]] ~
          tcga_state_all$ImmuneState,
        method = "bh"
      )
    )
  )
  
  
  tmp <- dunn_result$res
  
  
  tmp$Feature <- unname(
    feature_labels[
      feature
    ]
  )
  
  
  dunn_list[[i]] <- tmp
}


dunn_results <- do.call(
  rbind,
  dunn_list
)


rownames(
  dunn_results
) <- NULL


############################################################
## Convert internal state labels to manuscript labels
############################################################

state_replace <- c(
  
  High_Low =
    "CCL19 High / IL10 Low",
  
  High_High =
    "CCL19 High / IL10 High",
  
  Low_Low =
    "CCL19 Low / IL10 Low",
  
  Low_High =
    "CCL19 Low / IL10 High"
)


for (
  nm in names(
    state_replace
  )
) {
  
  dunn_results$Comparison <- gsub(
    nm,
    state_replace[[nm]],
    dunn_results$Comparison,
    fixed = TRUE
  )
}


dunn_results <- data.frame(
  
  Comparison =
    as.character(
      dunn_results$Comparison
    ),
  
  Z_statistic =
    as.numeric(
      dunn_results$Z
    ),
  
  Unadjusted_P =
    as.numeric(
      dunn_results$P.unadj
    ),
  
  BH_adjusted_P =
    as.numeric(
      dunn_results$P.adj
    ),
  
  Feature =
    as.character(
      dunn_results$Feature
    ),
  
  stringsAsFactors =
    FALSE
)


assert_true(
  nrow(
    dunn_results
  ) == 36L,
  "Figure-7C Dunn output must contain 36 rows."
)


assert_true(
  all(
    table(
      dunn_results$Feature
    ) == 6L
  ),
  "Each Figure-7C feature must contain six Dunn contrasts."
)


############################################################
## 12.36 — REGRESS COMPLETE DUNN TABLE
############################################################

current_dunn_key <- paste(
  dunn_results$Feature,
  dunn_results$Comparison,
  sep = "|||"
)


reference_dunn_key <- paste(
  S11E_ref$Feature,
  S11E_ref$Comparison,
  sep = "|||"
)


assert_true(
  !anyDuplicated(
    current_dunn_key
  ),
  "Duplicate current Dunn keys detected."
)


assert_true(
  !anyDuplicated(
    reference_dunn_key
  ),
  "Duplicate frozen Dunn keys detected."
)


Dunn_index <- match(
  current_dunn_key,
  reference_dunn_key
)


assert_true(
  !anyNA(
    Dunn_index
  ),
  "One or more Dunn comparisons could not be aligned."
)


S11E_aligned <- S11E_ref[
  Dunn_index,
  ,
  drop = FALSE
]


Dunn_Z_max_diff <- max(
  abs(
    dunn_results$Z_statistic -
      S11E_aligned$Z_statistic
  )
)


Dunn_raw_logP_max_diff <- max(
  abs(
    log10(
      dunn_results$Unadjusted_P
    ) -
      log10(
        S11E_aligned$Unadjusted_P
      )
  )
)


Dunn_BH_logP_max_diff <- max(
  abs(
    log10(
      dunn_results$BH_adjusted_P
    ) -
      log10(
        S11E_aligned$BH_adjusted_P
      )
  )
)


Dunn_regression_pass <- all(
  Dunn_Z_max_diff < 1e-10,
  Dunn_raw_logP_max_diff < 1e-10,
  Dunn_BH_logP_max_diff < 1e-10
)


assert_true(
  Dunn_regression_pass,
  "Complete Figure-7C Dunn regression failed."
)


############################################################
## 12.37 — DIRECT HIGH_LOW vs LOW_HIGH COMPARISON
############################################################

direct_comparison_name <-
  "CCL19 High / IL10 Low - CCL19 Low / IL10 High"


direct_comparison <- dunn_results[
  dunn_results$Comparison ==
    direct_comparison_name,
  ,
  drop = FALSE
]


assert_true(
  nrow(
    direct_comparison
  ) == 6L,
  "Expected six direct High_Low vs Low_High feature comparisons."
)


direct_index <- match(
  direct_comparison$Feature,
  S11F_ref$Feature
)


assert_true(
  !anyNA(
    direct_index
  ),
  "Unable to align direct High_Low vs Low_High frozen reference."
)


S11F_aligned <- S11F_ref[
  direct_index,
  ,
  drop = FALSE
]


direct_Z_max_diff <- max(
  abs(
    direct_comparison$Z_statistic -
      S11F_aligned$Z_statistic
  )
)


direct_raw_P_max_diff <- max(
  abs(
    direct_comparison$Unadjusted_P -
      S11F_aligned$Unadjusted_P
  )
)


direct_BH_P_max_diff <- max(
  abs(
    direct_comparison$BH_adjusted_P -
      S11F_aligned$BH_adjusted_P
  )
)


direct_regression_pass <- all(
  direct_Z_max_diff < 1e-10,
  direct_raw_P_max_diff < 1e-10,
  direct_BH_P_max_diff < 1e-10
)


assert_true(
  direct_regression_pass,
  "Direct High_Low vs Low_High regression failed."
)


############################################################
## 12.38 — OUTPUT DIRECTORIES
############################################################

output_dir <- file.path(
  "results",
  "tables",
  "tcga_supporting_characterization"
)


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
## 12.39 — SAVE PROCESSED OBJECTS
############################################################

S9_production_object <- list(
  
  hub_genes =
    S9_hub_genes,
  
  immune_features =
    S9_immune_features,
  
  correlation_method =
    "Pearson",
  
  correlation_engine =
    "Hmisc::rcorr",
  
  FDR_method =
    "Benjamini-Hochberg",
  
  FDR_family =
    "All 280 hub-gene x immune/TME tests",
  
  R_matrix =
    S9_R,
  
  raw_P_matrix =
    S9_P,
  
  FDR_matrix =
    S9_FDR,
  
  N_matrix =
    S9_N,
  
  complete_long_table =
    S9_long
)


save_rds_checked(
  S9_production_object,
  file.path(
    PATHS$processed,
    "12_tcga_hub_immune_correlations.rds"
  )
)


save_rds_checked(
  clinical_table,
  file.path(
    PATHS$processed,
    "12_tcga_clinicopathological_data.rds"
  )
)


save_rds_checked(
  clinical_associations,
  file.path(
    PATHS$processed,
    "12_tcga_clinicopathological_associations.rds"
  )
)


save_rds_checked(
  tcga_state_all,
  file.path(
    PATHS$processed,
    "12_tcga_joint_state_all1106.rds"
  )
)


save_rds_checked(
  state_kw,
  file.path(
    PATHS$processed,
    "12_tcga_joint_state_immune_omnibus.rds"
  )
)


save_rds_checked(
  state_medians,
  file.path(
    PATHS$processed,
    "12_tcga_joint_state_immune_medians.rds"
  )
)


save_rds_checked(
  dunn_results,
  file.path(
    PATHS$processed,
    "12_tcga_joint_state_immune_Dunn_BH.rds"
  )
)


save_rds_checked(
  direct_comparison,
  file.path(
    PATHS$processed,
    "12_tcga_joint_state_HighLow_vs_LowHigh.rds"
  )
)


############################################################
## 12.40 — EXPORT TABLE S9
############################################################

write_csv_checked(
  S9_long,
  file.path(
    output_dir,
    "12_TCGA_HubImmune_Pearson_Correlations_All280.csv"
  )
)


############################################################
## 12.41 — EXPORT CLINICAL TABLES
############################################################

write_csv_checked(
  clinical_table,
  file.path(
    output_dir,
    "12_TCGA_clinicopathological_data.csv"
  )
)


write_csv_checked(
  clinical_missingness,
  file.path(
    output_dir,
    "12_TCGA_clinicopathological_missingness.csv"
  )
)


write_csv_checked(
  clinical_associations,
  file.path(
    output_dir,
    "12_TCGA_clinicopathological_associations.csv"
  )
)


############################################################
## 12.42 — EXPORT JOINT-STATE CONTEXTURE TABLES
############################################################

write_csv_checked(
  state_kw,
  file.path(
    output_dir,
    "12_TCGA_joint_state_immune_omnibus.csv"
  )
)


write_csv_checked(
  state_medians,
  file.path(
    output_dir,
    "12_TCGA_joint_state_immune_medians.csv"
  )
)


write_csv_checked(
  dunn_results,
  file.path(
    output_dir,
    "12_TCGA_joint_state_Dunn_BH_all36.csv"
  )
)


write_csv_checked(
  direct_comparison,
  file.path(
    output_dir,
    "12_TCGA_HighLow_vs_LowHigh_immune_comparison.csv"
  )
)


############################################################
## 12.43 — METHOD PROVENANCE
############################################################

method_provenance <- data.frame(
  
  Analysis = c(
    "Hub/TME correlations",
    "Hub/TME correlation multiplicity",
    "Age vs immune phenotype",
    "Stage vs immune phenotype",
    "T stage vs immune phenotype",
    "N stage vs immune phenotype",
    "M stage vs immune phenotype",
    "Joint-state threshold source",
    "Joint-state analysis cohort",
    "Joint-state omnibus",
    "Joint-state omnibus multiplicity",
    "Joint-state post-hoc",
    "Joint-state post-hoc multiplicity",
    "Statistical interaction",
    "Causal interpretation"
  ),
  
  Method = c(
    "Pearson correlation using Hmisc::rcorr",
    "BH across all 280 hub x feature tests",
    "Kruskal-Wallis",
    "Pearson chi-square",
    "Pearson chi-square",
    "Pearson chi-square",
    "Pearson chi-square",
    "Exact CCL19/IL10 cutoffs from 1083-patient TCGA survival cohort",
    "All 1106 TCGA primary tumors",
    "Kruskal-Wallis across four joint expression states",
    "BH across six immune-contexture features",
    "FSA::dunnTest(method='bh')",
    "BH across six state contrasts separately within each feature",
    "No interaction term fitted",
    "Associative only; no causal inference"
  ),
  
  stringsAsFactors =
    FALSE
)


write_csv_checked(
  method_provenance,
  file.path(
    output_dir,
    "12_TCGA_supporting_characterization_method_provenance.csv"
  )
)


############################################################
## 12.44 — REGRESSION LOG
############################################################

regression_log <- data.frame(
  
  Metric = c(
    "S9_R_max_difference",
    "S9_rawP_max_difference",
    "S9_FDR_max_difference",
    "S9_N_exact",
    "KW_log10P_max_difference",
    "KW_log10FDR_max_difference",
    "State_median_max_difference",
    "Dunn_Z_max_difference",
    "Dunn_raw_log10P_max_difference",
    "Dunn_BH_log10P_max_difference",
    "Direct_Z_max_difference",
    "Direct_rawP_max_difference",
    "Direct_BH_P_max_difference"
  ),
  
  Value = c(
    S9_R_max_diff,
    S9_P_max_diff,
    S9_FDR_max_diff,
    S9_N_pass,
    kw_logP_max_diff,
    kw_logFDR_max_diff,
    state_median_max_diff,
    Dunn_Z_max_diff,
    Dunn_raw_logP_max_diff,
    Dunn_BH_logP_max_diff,
    direct_Z_max_diff,
    direct_raw_P_max_diff,
    direct_BH_P_max_diff
  ),
  
  stringsAsFactors =
    FALSE
)


write_csv_checked(
  regression_log,
  file.path(
    PATHS$logs,
    "12_tcga_supporting_characterization_regression.csv"
  )
)


############################################################
## 12.45 — REFERENCE MD5 LOG
############################################################

reference_md5_table <- data.frame(
  
  File =
    names(
      observed_reference_md5
    ),
  
  Expected_MD5 =
    unname(
      expected_reference_md5[
        names(
          observed_reference_md5
        )
      ]
    ),
  
  Observed_MD5 =
    unname(
      observed_reference_md5
    ),
  
  Verified =
    unname(
      observed_reference_md5
    ) ==
    unname(
      expected_reference_md5[
        names(
          observed_reference_md5
        )
      ]
    ),
  
  stringsAsFactors =
    FALSE
)


write_csv_checked(
  reference_md5_table,
  file.path(
    PATHS$logs,
    "12_tcga_supporting_characterization_reference_MD5.csv"
  )
)


############################################################
## 12.46 — MODULE QC TABLE
############################################################

module12_qc <- data.frame(
  
  Metric = c(
    "TCGA tumors",
    "Hub genes",
    "Hub/TME features",
    "Hub/TME tests",
    "Hub/TME BH family",
    "Hub/TME FDR-significant",
    "Hub/TME significant positive",
    "Hub/TME significant negative",
    "Hub/TME minimum Pearson r",
    "Hub/TME maximum Pearson r",
    "Age P",
    "Stage P",
    "T-stage P",
    "N-stage P",
    "M-stage P",
    "Joint-state threshold source N",
    "CCL19 threshold",
    "IL10 threshold",
    "All-1106 High_Low",
    "All-1106 High_High",
    "All-1106 Low_Low",
    "All-1106 Low_High",
    "Joint-state contexture features",
    "Joint-state KW BH family",
    "Joint-state Dunn comparisons",
    "Dunn BH family",
    "Historical S9 regression",
    "Historical clinical regression",
    "Historical KW regression",
    "Historical median regression",
    "Historical Dunn regression",
    "Historical direct-comparison regression",
    "Frozen reference MD5 verification",
    "All-1106 medians used as thresholds",
    "SCAN-B validation performed"
  ),
  
  Value = c(
    1106,
    20,
    14,
    280,
    "All 280 tests",
    277,
    257,
    20,
    format(
      min(S9_R),
      digits = 17
    ),
    format(
      max(S9_R),
      digits = 17
    ),
    format(
      age_test$p.value,
      digits = 17
    ),
    format(
      stage_test$p.value,
      digits = 17
    ),
    format(
      T_test$p.value,
      digits = 17
    ),
    format(
      N_test$p.value,
      digits = 17
    ),
    format(
      M_test$p.value,
      digits = 17
    ),
    1083,
    format(
      CCL19_cutoff,
      digits = 17
    ),
    format(
      IL10_cutoff,
      digits = 17
    ),
    state_counts["High_Low"],
    state_counts["High_High"],
    state_counts["Low_Low"],
    state_counts["Low_High"],
    6,
    "Six omnibus feature tests",
    36,
    "Six contrasts separately within each feature",
    S9_regression_pass,
    clinical_data_regression_pass,
    kw_regression_pass,
    state_median_regression_pass,
    Dunn_regression_pass,
    direct_regression_pass,
    reference_md5_pass,
    FALSE,
    FALSE
  ),
  
  stringsAsFactors =
    FALSE
)


write_csv_checked(
  module12_qc,
  file.path(
    PATHS$logs,
    "12_tcga_supporting_characterization_QC.csv"
  )
)


############################################################
## 12.47 — SESSION INFORMATION
############################################################

capture.output(
  sessionInfo(),
  file = file.path(
    PATHS$logs,
    "12_tcga_supporting_characterization_sessionInfo.txt"
  )
)


############################################################
## 12.48 — FINAL STATUS
############################################################

cat(
  "\n=============================================\n"
)


cat(
  "12_tcga_supporting_characterization.R: PASS\n"
)


cat(
  "\nTable S9 hub / immune correlations\n"
)


cat(
  "TCGA tumors: 1106\n"
)


cat(
  "Hub genes: 20\n"
)


cat(
  "Immune/TME features: 14\n"
)


cat(
  "Pearson tests: 280\n"
)


cat(
  "BH family: all 280 tests\n"
)


cat(
  "FDR-significant: ",
  sum(
    S9_significant
  ),
  "\n",
  sep = ""
)


cat(
  "Positive significant: ",
  S9_positive_significant,
  "\n",
  sep = ""
)


cat(
  "Negative significant: ",
  S9_negative_significant,
  "\n",
  sep = ""
)


cat(
  "Pearson r range: ",
  format(
    min(S9_R),
    digits = 17
  ),
  " to ",
  format(
    max(S9_R),
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "Historical Table S9 regression: VERIFIED\n"
)


cat(
  "\nClinicopathological associations\n"
)


cat(
  "Age: Kruskal-Wallis, P = ",
  format(
    age_test$p.value,
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "Stage: Pearson chi-square, P = ",
  format(
    stage_test$p.value,
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "T stage: Pearson chi-square, P = ",
  format(
    T_test$p.value,
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "N stage: Pearson chi-square, P = ",
  format(
    N_test$p.value,
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "M stage: Pearson chi-square, P = ",
  format(
    M_test$p.value,
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "Historical clinical-data regression: VERIFIED\n"
)


cat(
  "\nJoint-state immune contexture\n"
)


cat(
  "Analysis cohort: all 1106 TCGA tumors\n"
)


cat(
  "Threshold source: Module-11 1083-patient survival cohort\n"
)


cat(
  "CCL19 cutoff: ",
  format(
    CCL19_cutoff,
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "IL10 cutoff: ",
  format(
    IL10_cutoff,
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "High_Low: ",
  state_counts["High_Low"],
  "\n",
  sep = ""
)


cat(
  "High_High: ",
  state_counts["High_High"],
  "\n",
  sep = ""
)


cat(
  "Low_Low: ",
  state_counts["Low_Low"],
  "\n",
  sep = ""
)


cat(
  "Low_High: ",
  state_counts["Low_High"],
  "\n",
  sep = ""
)


cat(
  "Immune-contexture features: 6\n"
)


cat(
  "Omnibus method: Kruskal-Wallis\n"
)


cat(
  "Omnibus BH family: six features\n"
)


cat(
  "Post-hoc method: Dunn\n"
)


cat(
  "Dunn BH family: six contrasts within each feature\n"
)


cat(
  "Historical omnibus regression: VERIFIED\n"
)


cat(
  "Historical state-median regression: VERIFIED\n"
)


cat(
  "Historical complete Dunn regression: VERIFIED\n"
)


cat(
  "Historical High_Low vs Low_High regression: VERIFIED\n"
)


cat(
  "\nProvenance\n"
)


cat(
  "Frozen reference MD5 verification: VERIFIED\n"
)


cat(
  "All-1106 medians used as joint-state thresholds: FALSE\n"
)


cat(
  "Statistical interaction fitted: FALSE\n"
)


cat(
  "Historical analysis changed: FALSE\n"
)


cat(
  "SCAN-B validation performed: FALSE\n"
)


cat(
  "=============================================\n"
)