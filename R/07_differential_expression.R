############################################################
## 07_differential_expression.R
##
## PURPOSE
##
##   Differential-expression analysis between the frozen
##   TCGA-BRCA Immune_High and Immune_Low phenotypes.
##
## HISTORICAL METHOD
##
##   Comparison:
##     Immune_High - Immune_Low
##
##   Sample sizes:
##     Immune_Low  = 369
##     Immune_High = 369
##
##   Genes:
##     35,009
##
##   limma workflow:
##
##     lmFit()
##     makeContrasts(
##       High_vs_Low = Immune_High - Immune_Low
##     )
##     contrasts.fit()
##     eBayes()
##     topTable(
##       coef = "High_vs_Low",
##       number = Inf,
##       adjust.method = "BH"
##     )
##
## SIGNIFICANT DEG DEFINITION
##
##   adj.P.Val < 0.05
##   AND
##   abs(logFC) > 1
##
## INTERPRETATION
##
##   Positive logFC:
##     higher expression in Immune_High
##
##   Negative logFC:
##     higher expression in Immune_Low
##
## IMPORTANT
##
##   Immune_Mid is excluded from DEG analysis.
##
##   Expression preprocessing and phenotype construction
##   are NOT recomputed here.
############################################################


############################################################
## 7.1 — REQUIRE SETUP
############################################################

if (!exists("PATHS")) {
  source("R/00_setup.R")
}


############################################################
## 7.2 — REQUIRED PACKAGE
############################################################

require_package("limma")


############################################################
## 7.3 — INPUT FILES
############################################################

expression_file <- file.path(
  PATHS$processed,
  "02_tcga_expression_final.rds"
)


metadata_file <- file.path(
  PATHS$processed,
  "05_tcga_metadata_immune_phenotype.rds"
)


historical_reference_file <- file.path(
  "external",
  "differential_expression",
  "TCGA_BRCA_historical_DEG_results_reference.csv"
)


deg_full_reference_file <- file.path(
  "external",
  "differential_expression",
  "TCGA_BRCA_complete_DEG_High_vs_Low_full_results.csv"
)


assert_file(expression_file)
assert_file(metadata_file)
assert_file(historical_reference_file)
assert_file(deg_full_reference_file)


############################################################
## 7.4 — OUTPUT DIRECTORY
############################################################

deg_output_dir <- file.path(
  "results",
  "tables",
  "differential_expression"
)


dir.create(
  deg_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
## 7.5 — LOAD INPUT DATA
############################################################

expr_final <- readRDS(
  expression_file
)


meta05 <- readRDS(
  metadata_file
)


deg_hist <- read.csv(
  historical_reference_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


deg_full_reference <- read.csv(
  deg_full_reference_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


############################################################
## 7.6 — INPUT QC
############################################################

assert_true(
  nrow(expr_final) == 35009L,
  paste0(
    "Expected 35,009 expression genes; observed ",
    nrow(expr_final),
    "."
  )
)


assert_true(
  ncol(expr_final) == 1106L,
  paste0(
    "Expected 1,106 TCGA tumors; observed ",
    ncol(expr_final),
    "."
  )
)


assert_true(
  nrow(meta05) == 1106L,
  "Module-05 metadata must contain 1,106 tumors."
)


assert_true(
  "sample" %in% colnames(meta05),
  "Module-05 metadata is missing the sample column."
)


assert_true(
  "ImmunePhenotype" %in% colnames(meta05),
  "Module-05 metadata is missing ImmunePhenotype."
)


assert_identical_order(
  colnames(expr_final),
  meta05$sample,
  "expression samples",
  "Module-05 metadata samples"
)


assert_unique(
  rownames(expr_final),
  "expression gene symbols"
)


assert_unique(
  colnames(expr_final),
  "expression sample IDs"
)


assert_true(
  sum(is.na(expr_final)) == 0L,
  "Missing values detected in Module-02 expression matrix."
)


############################################################
## 7.7 — HISTORICAL REFERENCE QC
############################################################

assert_true(
  nrow(deg_hist) == 35009L,
  "Historical DEG reference must contain 35,009 genes."
)


assert_true(
  "Gene" %in% colnames(deg_hist),
  "Historical DEG reference lacks Gene column."
)


assert_true(
  all(
    c(
      "logFC",
      "AveExpr",
      "t",
      "P.Value",
      "adj.P.Val",
      "B"
    ) %in% colnames(deg_hist)
  ),
  "Historical DEG reference lacks required limma columns."
)


assert_unique(
  deg_hist$Gene,
  "historical DEG genes"
)


assert_true(
  nrow(deg_full_reference) == 35009L,
  "Complete DEG reference must contain 35,009 genes."
)


assert_true(
  "Gene" %in% colnames(deg_full_reference),
  "Complete DEG reference lacks Gene column."
)


assert_true(
  all(
    c(
      "logFC",
      "AveExpr",
      "t",
      "P_value",
      "FDR_BH",
      "B",
      "Direction",
      "DEG_status"
    ) %in% colnames(deg_full_reference)
  ),
  "Complete DEG reference lacks required columns."
)


assert_unique(
  deg_full_reference$Gene,
  "Complete DEG reference genes"
)


############################################################
## 7.8 — SELECT EXACT HISTORICAL DEG COHORT
############################################################

keep <- as.character(
  meta05$ImmunePhenotype
) %in% c(
  "Immune_Low",
  "Immune_High"
)


meta_deg <- meta05[
  keep,
  ,
  drop = FALSE
]


expr_deg <- expr_final[
  ,
  meta_deg$sample,
  drop = FALSE
]


assert_true(
  nrow(meta_deg) == 738L,
  paste0(
    "High/Low DEG cohort must contain 738 tumors; observed ",
    nrow(meta_deg),
    "."
  )
)


assert_true(
  ncol(expr_deg) == 738L,
  "DEG expression matrix must contain 738 tumors."
)


assert_true(
  nrow(expr_deg) == 35009L,
  "DEG expression matrix must contain 35,009 genes."
)


assert_identical_order(
  colnames(expr_deg),
  meta_deg$sample,
  "DEG expression samples",
  "DEG metadata samples"
)


############################################################
## 7.9 — GROUP FACTOR
############################################################

group <- factor(
  as.character(
    meta_deg$ImmunePhenotype
  ),
  levels = c(
    "Immune_Low",
    "Immune_High"
  )
)


group_counts <- table(
  group
)


assert_true(
  identical(
    levels(group),
    c(
      "Immune_Low",
      "Immune_High"
    )
  ),
  "DEG group-factor level order is incorrect."
)


assert_true(
  identical(
    as.integer(group_counts),
    c(
      369L,
      369L
    )
  ),
  "Expected 369 Immune-Low and 369 Immune-High tumors."
)


############################################################
## 7.10 — HISTORICAL DESIGN MATRIX
############################################################

design <- model.matrix(
  ~0 + group
)


colnames(design) <- c(
  "Immune_Low",
  "Immune_High"
)


rownames(design) <- meta_deg$sample


assert_true(
  identical(
    dim(design),
    c(
      738L,
      2L
    )
  ),
  "DEG design matrix must be 738 x 2."
)


assert_identical_order(
  rownames(design),
  colnames(expr_deg),
  "design samples",
  "DEG expression samples"
)


assert_true(
  sum(
    design[, "Immune_Low"]
  ) == 369L,
  "Immune-Low design count is incorrect."
)


assert_true(
  sum(
    design[, "Immune_High"]
  ) == 369L,
  "Immune-High design count is incorrect."
)


############################################################
## 7.11 — LIMMA FIT
############################################################

fit <- limma::lmFit(
  expr_deg,
  design
)


############################################################
## 7.12 — HISTORICAL CONTRAST
############################################################

contrast_matrix <- limma::makeContrasts(
  High_vs_Low =
    Immune_High -
    Immune_Low,
  levels = design
)


assert_true(
  contrast_matrix[
    "Immune_Low",
    "High_vs_Low"
  ] == -1,
  "Immune-Low contrast coefficient must be -1."
)


assert_true(
  contrast_matrix[
    "Immune_High",
    "High_vs_Low"
  ] == 1,
  "Immune-High contrast coefficient must be +1."
)


############################################################
## 7.13 — CONTRAST FIT
############################################################

fit2 <- limma::contrasts.fit(
  fit,
  contrast_matrix
)


############################################################
## 7.14 — EMPIRICAL BAYES
############################################################

## The historical analysis may produce the warning:
##
## "Zero sample variances detected, have been offset away
##  from zero"
##
## This is expected and was observed during validated
## historical regression.

fit2 <- limma::eBayes(
  fit2
)


############################################################
## 7.15 — HISTORICAL topTable CALL
############################################################

deg_core <- limma::topTable(
  fit2,
  coef = "High_vs_Low",
  number = Inf,
  adjust.method = "BH"
)


############################################################
## 7.16 — CORE DEG QC
############################################################

assert_true(
  nrow(deg_core) == 35009L,
  "limma result must contain 35,009 genes."
)


assert_true(
  ncol(deg_core) == 6L,
  "Core limma table must contain six statistical columns."
)


assert_true(
  all(
    c(
      "logFC",
      "AveExpr",
      "t",
      "P.Value",
      "adj.P.Val",
      "B"
    ) %in% colnames(deg_core)
  ),
  "Expected limma columns are absent."
)


assert_unique(
  rownames(deg_core),
  "limma result genes"
)


assert_true(
  sum(
    is.na(deg_core)
  ) == 0L,
  "Missing values detected in limma result."
)


assert_true(
  setequal(
    rownames(deg_core),
    rownames(expr_final)
  ),
  "limma gene set differs from Module-02 expression matrix."
)


############################################################
## 7.17 — SIGNIFICANCE THRESHOLDS
############################################################

FDR_THRESHOLD <- 0.05

LOGFC_THRESHOLD <- 1


significant_mask <-
  deg_core$adj.P.Val <
  FDR_THRESHOLD &
  abs(
    deg_core$logFC
  ) >
  LOGFC_THRESHOLD


############################################################
## 7.18 — CLEAN DIRECTION LABEL
############################################################

direction <- ifelse(
  deg_core$logFC > 0,
  "Higher_in_Immune_High",
  ifelse(
    deg_core$logFC < 0,
    "Higher_in_Immune_Low",
    "No_change"
  )
)


############################################################
## 7.19 — DEG STATUS
############################################################

deg_status <- ifelse(
  significant_mask,
  "Significant",
  "Not_significant"
)


############################################################
## 7.20 — HISTORICAL SIGNIFICANCE LABEL
############################################################

historical_significance <- ifelse(
  significant_mask &
    deg_core$logFC > 0,
  "Upregulated",
  ifelse(
    significant_mask &
      deg_core$logFC < 0,
    "Downregulated",
    "Not Significant"
  )
)


############################################################
## 7.21 — COMPLETE PRODUCTION DEG TABLE
############################################################

deg_complete <- data.frame(
  Gene =
    rownames(deg_core),
  
  logFC =
    deg_core$logFC,
  
  AveExpr =
    deg_core$AveExpr,
  
  t =
    deg_core$t,
  
  P_value =
    deg_core$P.Value,
  
  FDR_BH =
    deg_core$adj.P.Val,
  
  B =
    deg_core$B,
  
  Direction =
    direction,
  
  DEG_status =
    deg_status,
  
  stringsAsFactors = FALSE,
  check.names = FALSE
)


assert_true(
  nrow(deg_complete) == 35009L,
  "Complete DEG table must contain 35,009 genes."
)


assert_unique(
  deg_complete$Gene,
  "production DEG genes"
)


############################################################
## 7.22 — HISTORICAL-COMPATIBLE DEG OBJECT
############################################################

deg_results <- deg_core


deg_results$Significance <-
  historical_significance


## Historical Label was used for figure annotation.
## It is not part of the statistical DEG calculation.
##
## Figure-specific labels will be handled in the figure
## production/synchronization stage.

deg_results$Label <- ""


############################################################
## 7.23 — SIGNIFICANT DEG TABLE
############################################################

deg_significant <- deg_complete[
  deg_complete$DEG_status ==
    "Significant",
  ,
  drop = FALSE
]


assert_true(
  nrow(deg_significant) == 2830L,
  paste0(
    "Expected 2,830 significant DEGs; observed ",
    nrow(deg_significant),
    "."
  )
)


n_higher_high <- sum(
  deg_significant$Direction ==
    "Higher_in_Immune_High"
)


n_higher_low <- sum(
  deg_significant$Direction ==
    "Higher_in_Immune_Low"
)


assert_true(
  n_higher_high == 2463L,
  paste0(
    "Expected 2,463 DEGs higher in Immune-High; observed ",
    n_higher_high,
    "."
  )
)


assert_true(
  n_higher_low == 367L,
  paste0(
    "Expected 367 DEGs higher in Immune-Low; observed ",
    n_higher_low,
    "."
  )
)


############################################################
## 7.24 — DIRECT HIGH-MINUS-LOW ORIENTATION CHECK
############################################################

phen_deg <- as.character(
  meta_deg$ImmunePhenotype
)


high_index <-
  phen_deg ==
  "Immune_High"


low_index <-
  phen_deg ==
  "Immune_Low"


assert_true(
  sum(high_index) == 369L,
  "Direct-orientation check expected 369 Immune-High tumors."
)


assert_true(
  sum(low_index) == 369L,
  "Direct-orientation check expected 369 Immune-Low tumors."
)


mean_high <- rowMeans(
  expr_deg[
    ,
    high_index,
    drop = FALSE
  ]
)


mean_low <- rowMeans(
  expr_deg[
    ,
    low_index,
    drop = FALSE
  ]
)


direct_high_minus_low <-
  mean_high -
  mean_low


direct_aligned <-
  direct_high_minus_low[
    rownames(deg_core)
  ]


assert_true(
  !anyNA(direct_aligned),
  "Direct High-Low expression differences could not be aligned."
)


orientation_correlation <- cor(
  deg_core$logFC,
  direct_aligned
)


orientation_max_difference <- max(
  abs(
    deg_core$logFC -
      direct_aligned
  )
)


assert_true(
  abs(
    orientation_correlation -
      1
  ) <
    1e-12,
  "limma logFC does not agree with direct High-Low expression difference."
)


assert_true(
  orientation_max_difference <
    1e-10,
  "limma logFC differs unexpectedly from direct High-Low difference."
)


############################################################
## 7.25 — ALIGN HISTORICAL DEG REFERENCE
############################################################

hist_index <- match(
  deg_complete$Gene,
  deg_hist$Gene
)


assert_true(
  !anyNA(hist_index),
  "Historical DEG reference is missing production genes."
)


deg_hist_aligned <- deg_hist[
  hist_index,
  ,
  drop = FALSE
]


assert_identical_order(
  deg_complete$Gene,
  deg_hist_aligned$Gene,
  "production DEG genes",
  "historical DEG genes"
)


############################################################
## 7.26 — HISTORICAL NUMERICAL REGRESSION
############################################################

limma_metrics <- c(
  "logFC",
  "AveExpr",
  "t",
  "P.Value",
  "adj.P.Val",
  "B"
)


historical_regression <- data.frame(
  Metric =
    limma_metrics,
  
  Max_abs_difference =
    NA_real_,
  
  Mean_abs_difference =
    NA_real_,
  
  Pearson_r =
    NA_real_,
  
  stringsAsFactors = FALSE
)


for (
  i in seq_along(
    limma_metrics
  )
) {
  
  metric <-
    limma_metrics[i]
  
  
  current_values <-
    as.numeric(
      deg_core[[metric]]
    )
  
  
  historical_values <-
    as.numeric(
      deg_hist_aligned[[metric]]
    )
  
  
  historical_regression$Max_abs_difference[i] <-
    max(
      abs(
        current_values -
          historical_values
      )
    )
  
  
  historical_regression$Mean_abs_difference[i] <-
    mean(
      abs(
        current_values -
          historical_values
      )
    )
  
  
  historical_regression$Pearson_r[i] <-
    suppressWarnings(
      cor(
        current_values,
        historical_values
      )
    )
}


HISTORICAL_TOLERANCE <- 1e-10


historical_regression_pass <-
  max(
    historical_regression$Max_abs_difference
  ) <
  HISTORICAL_TOLERANCE


assert_true(
  historical_regression_pass,
  paste0(
    "Historical 35,009-gene limma regression failed. Max difference = ",
    format(
      max(
        historical_regression$Max_abs_difference
      ),
      scientific = TRUE,
      digits = 17
    ),
    "."
  )
)


############################################################
## 7.27 — HISTORICAL SIGNIFICANT GENE-SET REGRESSION
############################################################

historical_sig_mask <-
  deg_hist$adj.P.Val < 0.05 &
  abs(
    deg_hist$logFC
  ) > 1


historical_sig_genes <-
  deg_hist$Gene[
    historical_sig_mask
  ]


significant_set_regression_pass <-
  setequal(
    deg_significant$Gene,
    historical_sig_genes
  )


assert_true(
  significant_set_regression_pass,
  "Significant DEG gene set differs from historical result."
)


############################################################
## 7.28 — ALIGN Complete DEG reference
############################################################

deg_full_reference_index <- match(
  deg_complete$Gene,
  deg_full_reference$Gene
)


assert_true(
  !anyNA(deg_full_reference_index),
  "Complete DEG reference is missing production genes."
)


deg_full_reference_aligned <- deg_full_reference[
  deg_full_reference_index,
  ,
  drop = FALSE
]


assert_identical_order(
  deg_complete$Gene,
  deg_full_reference_aligned$Gene,
  "production DEG genes",
  "Complete DEG reference genes"
)


############################################################
## 7.29 — COMPLETE DEG REFERENCE NUMERICAL REGRESSION
############################################################

deg_reference_mapping <- c(
  logFC =
    "logFC",
  
  AveExpr =
    "AveExpr",
  
  t =
    "t",
  
  P_value =
    "P_value",
  
  FDR_BH =
    "FDR_BH",
  
  B =
    "B"
)


deg_reference_regression <- data.frame(
  Metric =
    names(
      deg_reference_mapping
    ),
  
  Max_abs_difference =
    NA_real_,
  
  Mean_abs_difference =
    NA_real_,
  
  stringsAsFactors = FALSE
)


for (
  i in seq_along(
    deg_reference_mapping
  )
) {
  
  production_column <-
    names(
      deg_reference_mapping
    )[i]
  
  
  frozen_column <-
    unname(
      deg_reference_mapping[i]
    )
  
  
  ##########################################################
  ## CORRECTED SYNTAX:
  ## use [[column_name]] for programmatic column extraction.
  ##########################################################
  
  current_values <-
    as.numeric(
      deg_complete[[production_column]]
    )
  
  
  frozen_values <-
    as.numeric(
      deg_full_reference_aligned[[frozen_column]]
    )
  
  
  deg_reference_regression$Max_abs_difference[i] <-
    max(
      abs(
        current_values -
          frozen_values
      )
    )
  
  
  deg_reference_regression$Mean_abs_difference[i] <-
    mean(
      abs(
        current_values -
          frozen_values
      )
    )
}


DEG_REFERENCE_TOLERANCE <- 1e-10


deg_reference_regression_pass <-
  max(
    deg_reference_regression$Max_abs_difference
  ) <
  DEG_REFERENCE_TOLERANCE


assert_true(
  deg_reference_regression_pass,
  paste0(
    "Complete DEG reference numerical regression failed. Max difference = ",
    format(
      max(
        deg_reference_regression$Max_abs_difference
      ),
      scientific = TRUE,
      digits = 17
    ),
    "."
  )
)


############################################################
## 7.30 — COMPLETE DEG REFERENCE DIRECTION / STATUS REGRESSION
############################################################

direction_regression_pass <-
  identical(
    as.character(
      deg_complete$Direction
    ),
    as.character(
      deg_full_reference_aligned$Direction
    )
  )


status_regression_pass <-
  identical(
    as.character(
      deg_complete$DEG_status
    ),
    as.character(
      deg_full_reference_aligned$DEG_status
    )
  )


assert_true(
  direction_regression_pass,
  "Frozen DEG direction labels differ."
)


assert_true(
  status_regression_pass,
  "Frozen DEG status labels differ."
)


############################################################
## 7.31 — SUMMARY TABLE
############################################################

deg_summary <- data.frame(
  Metric = c(
    "Genes tested",
    "Immune Low samples",
    "Immune High samples",
    "Significant DEGs",
    "Higher in Immune High",
    "Higher in Immune Low",
    "FDR threshold",
    "Absolute logFC threshold",
    "Contrast"
  ),
  
  Value = c(
    nrow(deg_complete),
    as.integer(
      group_counts["Immune_Low"]
    ),
    as.integer(
      group_counts["Immune_High"]
    ),
    nrow(deg_significant),
    n_higher_high,
    n_higher_low,
    FDR_THRESHOLD,
    LOGFC_THRESHOLD,
    "Immune_High - Immune_Low"
  ),
  
  stringsAsFactors = FALSE
)


############################################################
## 7.32 — TOP DEG TABLE
############################################################

top_n <- 25L


top_high <- deg_complete[
  deg_complete$DEG_status ==
    "Significant" &
    deg_complete$Direction ==
    "Higher_in_Immune_High",
  ,
  drop = FALSE
]


top_high <- top_high[
  order(
    top_high$logFC,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]


top_high <- head(
  top_high,
  top_n
)


top_high$Rank_group <-
  seq_len(
    nrow(top_high)
  )


top_high$Group <-
  "Higher_in_Immune_High"


top_low <- deg_complete[
  deg_complete$DEG_status ==
    "Significant" &
    deg_complete$Direction ==
    "Higher_in_Immune_Low",
  ,
  drop = FALSE
]


top_low <- top_low[
  order(
    top_low$logFC,
    decreasing = FALSE
  ),
  ,
  drop = FALSE
]


top_low <- head(
  top_low,
  top_n
)


top_low$Rank_group <-
  seq_len(
    nrow(top_low)
  )


top_low$Group <-
  "Higher_in_Immune_Low"


top_deg_table <- rbind(
  top_high,
  top_low
)


rownames(top_deg_table) <- NULL


############################################################
## 7.33 — SAVE LIMMA MODEL OBJECTS
############################################################

deg_model <- list(
  group =
    group,
  
  design =
    design,
  
  contrast_matrix =
    contrast_matrix,
  
  fit =
    fit,
  
  fit2 =
    fit2,
  
  contrast =
    "Immune_High - Immune_Low",
  
  FDR_threshold =
    FDR_THRESHOLD,
  
  absolute_logFC_threshold =
    LOGFC_THRESHOLD
)


save_rds_checked(
  deg_model,
  file.path(
    PATHS$processed,
    "07_tcga_limma_model.rds"
  )
)


############################################################
## 7.34 — SAVE HISTORICAL-COMPATIBLE DEG OBJECT
############################################################

save_rds_checked(
  deg_results,
  file.path(
    PATHS$processed,
    "07_tcga_deg_results.rds"
  )
)


############################################################
## 7.35 — SAVE CLEAN COMPLETE TABLE
############################################################

save_rds_checked(
  deg_complete,
  file.path(
    PATHS$processed,
    "07_tcga_deg_complete.rds"
  )
)


write_csv_checked(
  deg_complete,
  file.path(
    deg_output_dir,
    "07_TCGA_complete_DEG_High_vs_Low.csv"
  )
)


############################################################
## 7.36 — SAVE SIGNIFICANT DEG TABLE
############################################################

save_rds_checked(
  deg_significant,
  file.path(
    PATHS$processed,
    "07_tcga_deg_significant.rds"
  )
)


write_csv_checked(
  deg_significant,
  file.path(
    deg_output_dir,
    "07_TCGA_significant_DEG_High_vs_Low.csv"
  )
)


############################################################
## 7.37 — SAVE SUMMARY + TOP GENES
############################################################

write_csv_checked(
  deg_summary,
  file.path(
    deg_output_dir,
    "07_TCGA_DEG_summary.csv"
  )
)


write_csv_checked(
  top_deg_table,
  file.path(
    deg_output_dir,
    "07_TCGA_top_DEGs.csv"
  )
)


############################################################
## 7.38 — SAVE HISTORICAL REGRESSION
############################################################

write_csv_checked(
  historical_regression,
  file.path(
    PATHS$logs,
    "07_limma_historical_regression.csv"
  )
)


write_csv_checked(
  deg_reference_regression,
  file.path(
    PATHS$logs,
    "07_complete_DEG_reference_regression.csv"
  )
)


############################################################
## 7.39 — SAVE ORIENTATION QC
############################################################

orientation_qc <- data.frame(
  Metric = c(
    "Genes tested",
    "Immune High samples",
    "Immune Low samples",
    "logFC versus direct High-Low correlation",
    "Maximum direct High-Low discrepancy"
  ),
  
  Value = c(
    nrow(deg_core),
    sum(high_index),
    sum(low_index),
    orientation_correlation,
    orientation_max_difference
  ),
  
  stringsAsFactors = FALSE
)


write_csv_checked(
  orientation_qc,
  file.path(
    PATHS$logs,
    "07_DEG_orientation_QC.csv"
  )
)


############################################################
## 7.40 — FILE HASHES
############################################################

historical_reference_md5 <- unname(
  tools::md5sum(
    historical_reference_file
  )
)


deg_reference_md5 <- unname(
  tools::md5sum(
    deg_full_reference_file
  )
)


############################################################
## 7.41 — MODULE QC
############################################################

module07_qc <- data.frame(
  Metric = c(
    "Genes tested",
    "DEG cohort tumors",
    "Immune Low tumors",
    "Immune High tumors",
    "Significant DEGs",
    "Higher in Immune High",
    "Higher in Immune Low",
    "FDR threshold",
    "Absolute logFC threshold",
    "Contrast",
    "Historical limma max absolute difference",
    "Historical limma regression",
    "Complete DEG reference max absolute difference",
    "Complete DEG reference numerical regression",
    "Significant gene-set regression",
    "Direction regression",
    "DEG-status regression",
    "logFC direct-difference correlation",
    "logFC maximum direct-difference discrepancy",
    "limma version",
    "Historical DEG reference MD5",
    "Complete DEG reference MD5"
  ),
  
  Value = c(
    nrow(deg_complete),
    
    nrow(meta_deg),
    
    as.integer(
      group_counts["Immune_Low"]
    ),
    
    as.integer(
      group_counts["Immune_High"]
    ),
    
    nrow(deg_significant),
    
    n_higher_high,
    
    n_higher_low,
    
    FDR_THRESHOLD,
    
    LOGFC_THRESHOLD,
    
    "Immune_High - Immune_Low",
    
    format(
      max(
        historical_regression$Max_abs_difference
      ),
      scientific = TRUE,
      digits = 17
    ),
    
    historical_regression_pass,
    
    format(
      max(
        deg_reference_regression$Max_abs_difference
      ),
      scientific = TRUE,
      digits = 17
    ),
    
    deg_reference_regression_pass,
    
    significant_set_regression_pass,
    
    direction_regression_pass,
    
    status_regression_pass,
    
    format(
      orientation_correlation,
      digits = 17
    ),
    
    format(
      orientation_max_difference,
      scientific = TRUE,
      digits = 17
    ),
    
    as.character(
      packageVersion("limma")
    ),
    
    historical_reference_md5,
    
    deg_reference_md5
  ),
  
  stringsAsFactors = FALSE
)


write_csv_checked(
  module07_qc,
  file.path(
    PATHS$logs,
    "07_differential_expression_QC.csv"
  )
)


############################################################
## 7.42 — SESSION INFORMATION
############################################################

capture.output(
  sessionInfo(),
  file = file.path(
    PATHS$logs,
    "07_differential_expression_sessionInfo.txt"
  )
)


############################################################
## 7.43 — FINAL STATUS
############################################################

cat(
  "\n=============================================\n"
)


cat(
  "07_differential_expression.R: PASS\n"
)


cat(
  "Genes tested: ",
  nrow(deg_complete),
  "\n",
  sep = ""
)


cat(
  "DEG cohort: ",
  nrow(meta_deg),
  " tumors\n",
  sep = ""
)


cat(
  "Immune Low / High: ",
  group_counts["Immune_Low"],
  " / ",
  group_counts["Immune_High"],
  "\n",
  sep = ""
)


cat(
  "Contrast: Immune_High - Immune_Low\n"
)


cat(
  "Significance criterion: BH-FDR < 0.05 AND |logFC| > 1\n"
)


cat(
  "Significant DEGs: ",
  nrow(deg_significant),
  "\n",
  sep = ""
)


cat(
  "Higher in Immune High: ",
  n_higher_high,
  "\n",
  sep = ""
)


cat(
  "Higher in Immune Low: ",
  n_higher_low,
  "\n",
  sep = ""
)


cat(
  "Historical limma max difference: ",
  format(
    max(
      historical_regression$Max_abs_difference
    ),
    scientific = TRUE,
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "Complete DEG reference max difference: ",
  format(
    max(
      deg_reference_regression$Max_abs_difference
    ),
    scientific = TRUE,
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "Historical limma numerical regression: VERIFIED\n"
)


cat(
  "Complete DEG reference numerical regression: VERIFIED\n"
)


cat(
  "Significant DEG gene set: VERIFIED\n"
)


cat(
  "Direction/status labels: VERIFIED\n"
)


cat(
  "logFC orientation: VERIFIED AS HIGH MINUS LOW\n"
)


cat(
  "Expression preprocessing recomputed: FALSE\n"
)


cat(
  "Immune phenotype recomputed: FALSE\n"
)


cat(
  "=============================================\n"
)