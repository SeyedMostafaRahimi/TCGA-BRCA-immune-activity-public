############################################################
## 05_ssgsea_immune_axis.R
##
## PURPOSE
##   Reproduce the historical 17-signature ssGSEA immune
##   axis and TCGA-BRCA Immune_Low / Immune_Mid /
##   Immune_High phenotypes.
##
## HISTORICAL METHOD
##
##   IOBR::calculate_sig_score(
##     eset = expr_final,
##     signature = selected_signatures2,
##     method = "ssgsea"
##   )
##
##   followed by:
##
##   prcomp(
##     ssGSEA_matrix,
##     center = TRUE,
##     scale. = TRUE
##   )
##
##   PC1 tertiles:
##     quantile(
##       PC1,
##       probs = c(1/3, 2/3),
##       type = 7
##     )
##
##   Phenotypes:
##     Immune_Low
##     Immune_Mid
##     Immune_High
##
## REPRODUCIBILITY PROTECTION
##
##   - Exact historical 17 gene sets are frozen locally.
##   - IOBR remains the implementation layer.
##   - Samples are aligned explicitly by TCGA sample ID.
##   - PCA PC1 orientation is deterministic:
##       higher PC1 = higher aggregate signature activity.
##   - Historical ssGSEA, PCA, PC1 and phenotype assignments
##     are regression-tested.
############################################################


############################################################
## 5.1 — REQUIRE SETUP
############################################################

if (!exists("PATHS")) {
  source("R/00_setup.R")
}


############################################################
## 5.2 — REQUIRED PACKAGES
############################################################

require_package("IOBR")
require_package("GSVA")


############################################################
## 5.3 — EXPLICIT INPUT FILES
############################################################

expression_file <- file.path(
  PATHS$processed,
  "02_tcga_expression_final.rds"
)

metadata_file <- file.path(
  PATHS$processed,
  "04_tcga_metadata_tme.rds"
)

signature_file <- file.path(
  "external",
  "iobr",
  "TCGA_BRCA_historical_17_signatures.rds"
)

ssgsea_reference_file <- file.path(
  "external",
  "iobr",
  "TCGA_BRCA_historical_ssGSEA17_reference.csv"
)

phenotype_reference_file <- file.path(
  "external",
  "iobr",
  "TCGA_BRCA_historical_PC1_phenotype_reference.csv"
)

pca_reference_file <- file.path(
  "external",
  "iobr",
  "TCGA_BRCA_historical_ssGSEA_PCA_reference.rds"
)


assert_file(expression_file)
assert_file(metadata_file)
assert_file(signature_file)
assert_file(ssgsea_reference_file)
assert_file(phenotype_reference_file)
assert_file(pca_reference_file)


############################################################
## 5.4 — LOAD INPUTS
############################################################

expr_final <- readRDS(
  expression_file
)

meta_input <- readRDS(
  metadata_file
)

selected_signatures <- readRDS(
  signature_file
)

ssgsea_hist <- read.csv(
  ssgsea_reference_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

phenotype_hist <- read.csv(
  phenotype_reference_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

pca_hist <- readRDS(
  pca_reference_file
)


############################################################
## 5.5 — EXPRESSION / METADATA QC
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
    "Expected 1,106 tumors; observed ",
    ncol(expr_final),
    "."
  )
)

assert_true(
  nrow(meta_input) == 1106L,
  paste0(
    "Expected 1,106 metadata rows; observed ",
    nrow(meta_input),
    "."
  )
)

assert_identical_order(
  colnames(expr_final),
  meta_input$sample,
  "expression samples",
  "Module-04 metadata samples"
)

assert_unique(
  colnames(expr_final),
  "expression sample IDs"
)

assert_unique(
  rownames(expr_final),
  "expression gene symbols"
)


############################################################
## 5.6 — SIGNATURE QC
############################################################

assert_true(
  is.list(selected_signatures),
  "Historical signature object is not a list."
)

assert_true(
  length(selected_signatures) == 17L,
  paste0(
    "Expected 17 historical signatures; observed ",
    length(selected_signatures),
    "."
  )
)

assert_unique(
  names(selected_signatures),
  "historical signature names"
)


ssgsea_features <- setdiff(
  colnames(ssgsea_hist),
  "sample"
)


assert_true(
  length(ssgsea_features) == 17L,
  "Historical ssGSEA reference does not contain 17 signatures."
)

assert_true(
  identical(
    names(selected_signatures),
    ssgsea_features
  ),
  "Historical signature order differs from frozen ssGSEA reference."
)


############################################################
## 5.7 — HISTORICAL REFERENCE QC
############################################################

assert_true(
  nrow(ssgsea_hist) == 1106L,
  "Historical ssGSEA reference does not contain 1,106 tumors."
)

assert_true(
  nrow(phenotype_hist) == 1106L,
  "Historical phenotype reference does not contain 1,106 tumors."
)

assert_unique(
  ssgsea_hist$sample,
  "historical ssGSEA sample IDs"
)

assert_unique(
  phenotype_hist$sample,
  "historical phenotype sample IDs"
)


assert_true(
  setequal(
    ssgsea_hist$sample,
    meta_input$sample
  ),
  "Historical and production ssGSEA sample sets differ."
)

assert_true(
  setequal(
    phenotype_hist$sample,
    meta_input$sample
  ),
  "Historical and production phenotype sample sets differ."
)


############################################################
## 5.8 — ALIGN HISTORICAL REFERENCES TO PRODUCTION ORDER
############################################################

ssgsea_hist <- ssgsea_hist[
  match(
    meta_input$sample,
    ssgsea_hist$sample
  ),
  ,
  drop = FALSE
]

phenotype_hist <- phenotype_hist[
  match(
    meta_input$sample,
    phenotype_hist$sample
  ),
  ,
  drop = FALSE
]


assert_identical_order(
  ssgsea_hist$sample,
  meta_input$sample,
  "historical ssGSEA samples",
  "production samples"
)

assert_identical_order(
  phenotype_hist$sample,
  meta_input$sample,
  "historical phenotype samples",
  "production samples"
)


############################################################
## 5.9 — VERIFY HISTORICAL IOBR FUNCTION
############################################################

assert_true(
  exists(
    "calculate_sig_score",
    envir = asNamespace("IOBR"),
    inherits = FALSE
  ),
  "IOBR::calculate_sig_score() is unavailable."
)


############################################################
## 5.10 — RUN HISTORICAL ssGSEA METHOD
############################################################

cat(
  "\nRunning historical IOBR 17-signature ssGSEA implementation...\n"
)


ssgsea_raw <- IOBR::calculate_sig_score(
  eset = expr_final,
  signature = selected_signatures,
  method = "ssgsea"
)


############################################################
## 5.11 — RAW ssGSEA QC
############################################################

ssgsea_raw <- as.data.frame(
  ssgsea_raw,
  check.names = FALSE
)


assert_true(
  nrow(ssgsea_raw) == 1106L,
  paste0(
    "Expected 1,106 ssGSEA rows; observed ",
    nrow(ssgsea_raw),
    "."
  )
)

assert_true(
  "ID" %in% colnames(ssgsea_raw),
  "IOBR ssGSEA output is missing ID column."
)

assert_true(
  all(
    ssgsea_features %in%
      colnames(ssgsea_raw)
  ),
  "IOBR ssGSEA output is missing historical signatures."
)

assert_unique(
  ssgsea_raw$ID,
  "IOBR ssGSEA sample IDs"
)

assert_true(
  setequal(
    ssgsea_raw$ID,
    meta_input$sample
  ),
  "IOBR ssGSEA sample set differs from production TCGA samples."
)


############################################################
## 5.12 — ALIGN IOBR OUTPUT BY SAMPLE ID
############################################################

ssgsea_raw <- ssgsea_raw[
  match(
    meta_input$sample,
    ssgsea_raw$ID
  ),
  ,
  drop = FALSE
]


assert_identical_order(
  ssgsea_raw$ID,
  meta_input$sample,
  "aligned ssGSEA samples",
  "production samples"
)


############################################################
## 5.13 — BUILD CLEAN ssGSEA TABLE
############################################################

ssgsea_scores <- data.frame(
  
  sample =
    ssgsea_raw$ID,
  
  ssgsea_raw[
    ,
    ssgsea_features,
    drop = FALSE
  ],
  
  check.names = FALSE,
  stringsAsFactors = FALSE
)


assert_true(
  nrow(ssgsea_scores) == 1106L,
  "Clean ssGSEA table does not contain 1,106 tumors."
)

assert_true(
  ncol(ssgsea_scores) == 18L,
  "Clean ssGSEA table does not contain sample + 17 signatures."
)

assert_no_missing(
  ssgsea_scores[
    ,
    ssgsea_features,
    drop = FALSE
  ],
  "17-signature ssGSEA scores"
)


############################################################
## 5.14 — ssGSEA HISTORICAL REGRESSION
############################################################

ssgsea_regression <- data.frame(
  
  Signature =
    ssgsea_features,
  
  Max_abs_difference =
    NA_real_,
  
  Mean_abs_difference =
    NA_real_,
  
  Pearson_r =
    NA_real_,
  
  stringsAsFactors = FALSE
)


for (i in seq_along(ssgsea_features)) {
  
  feature <- ssgsea_features[i]
  
  production_values <-
    as.numeric(
      ssgsea_scores[[feature]]
    )
  
  historical_values <-
    as.numeric(
      ssgsea_hist[[feature]]
    )
  
  
  ssgsea_regression$Max_abs_difference[i] <-
    max(
      abs(
        production_values -
          historical_values
      )
    )
  
  
  ssgsea_regression$Mean_abs_difference[i] <-
    mean(
      abs(
        production_values -
          historical_values
      )
    )
  
  
  ssgsea_regression$Pearson_r[i] <-
    suppressWarnings(
      cor(
        production_values,
        historical_values,
        method = "pearson"
      )
    )
}


ssgsea_global_max <- max(
  ssgsea_regression$Max_abs_difference
)


SSGSEA_TOLERANCE <- 1e-12


ssgsea_regression_pass <-
  ssgsea_global_max <
  SSGSEA_TOLERANCE


assert_true(
  ssgsea_regression_pass,
  paste0(
    "Historical ssGSEA regression failed. ",
    "Maximum absolute difference = ",
    format(
      ssgsea_global_max,
      scientific = TRUE,
      digits = 17
    ),
    "."
  )
)


############################################################
## 5.15 — CREATE PCA INPUT MATRIX
############################################################

ss_mat <- as.matrix(
  ssgsea_scores[
    ,
    ssgsea_features,
    drop = FALSE
  ]
)

storage.mode(
  ss_mat
) <- "double"

rownames(
  ss_mat
) <- ssgsea_scores$sample


assert_true(
  identical(
    dim(ss_mat),
    c(1106L, 17L)
  ),
  "PCA input matrix must be 1,106 x 17."
)

assert_true(
  !anyNA(ss_mat),
  "PCA input contains missing values."
)


############################################################
## 5.16 — HISTORICAL PCA METHOD
############################################################

ss_pca <- prcomp(
  ss_mat,
  center = TRUE,
  scale. = TRUE
)


############################################################
## 5.17 — DETERMINISTIC PC1 ORIENTATION
############################################################

## Principal-component signs are mathematically arbitrary.
##
## The historical/frozen PC1 had all 17 positive loadings.
## Production therefore orients PC1 so the sum of its
## 17 loadings is positive.
##
## This prevents an arbitrary PCA sign reversal from
## switching the biological Low/High interpretation.

pc1_flipped <- FALSE


if (
  sum(
    ss_pca$rotation[
      ,
      "PC1"
    ]
  ) < 0
) {
  
  ss_pca$rotation[
    ,
    "PC1"
  ] <-
    -ss_pca$rotation[
      ,
      "PC1"
    ]
  
  ss_pca$x[
    ,
    "PC1"
  ] <-
    -ss_pca$x[
      ,
      "PC1"
    ]
  
  pc1_flipped <- TRUE
}


############################################################
## 5.18 — VERIFY PC1 DIRECTION
############################################################

pc1_loadings <- ss_pca$rotation[
  ssgsea_features,
  "PC1"
]


assert_true(
  all(
    pc1_loadings > 0
  ),
  "Final PC1 does not have the historical positive-loading orientation."
)


############################################################
## 5.19 — PCA VARIANCE EXPLAINED
############################################################

pca_variance <- ss_pca$sdev^2 /
  sum(
    ss_pca$sdev^2
  )


pc1_variance_percent <-
  100 *
  pca_variance[1]


############################################################
## 5.20 — PCA REGRESSION AGAINST HISTORICAL REFERENCE
############################################################

loading_max_difference <- max(
  abs(
    pc1_loadings -
      pca_hist$rotation[
        ssgsea_features,
        "PC1"
      ]
  )
)


center_max_difference <- max(
  abs(
    ss_pca$center[
      ssgsea_features
    ] -
      pca_hist$center[
        ssgsea_features
      ]
  )
)


scale_max_difference <- max(
  abs(
    ss_pca$scale[
      ssgsea_features
    ] -
      pca_hist$scale[
        ssgsea_features
      ]
  )
)


variance_difference <- abs(
  pc1_variance_percent -
    pca_hist$PC1_variance_percent
)


PCA_COMPONENT_TOLERANCE <- 1e-12


assert_true(
  loading_max_difference <
    PCA_COMPONENT_TOLERANCE,
  "Historical PC1 loading regression failed."
)

assert_true(
  center_max_difference <
    PCA_COMPONENT_TOLERANCE,
  "Historical PCA centering regression failed."
)

assert_true(
  scale_max_difference <
    PCA_COMPONENT_TOLERANCE,
  "Historical PCA scaling regression failed."
)

assert_true(
  variance_difference <
    PCA_COMPONENT_TOLERANCE,
  "Historical PC1 variance regression failed."
)


############################################################
## 5.21 — SAMPLE-LEVEL PCA TABLE
############################################################

pca_scores <- data.frame(
  
  sample =
    rownames(ss_pca$x),
  
  PC1 =
    as.numeric(
      ss_pca$x[
        ,
        "PC1"
      ]
    ),
  
  PC2 =
    as.numeric(
      ss_pca$x[
        ,
        "PC2"
      ]
    ),
  
  PC3 =
    as.numeric(
      ss_pca$x[
        ,
        "PC3"
      ]
    ),
  
  stringsAsFactors = FALSE
)


assert_identical_order(
  pca_scores$sample,
  meta_input$sample,
  "PCA samples",
  "production samples"
)


############################################################
## 5.22 — SAMPLE-LEVEL PC1 HISTORICAL REGRESSION
############################################################

pc1_max_difference <- max(
  abs(
    pca_scores$PC1 -
      phenotype_hist$PC1
  )
)


pc1_mean_difference <- mean(
  abs(
    pca_scores$PC1 -
      phenotype_hist$PC1
  )
)


pc1_correlation <- cor(
  pca_scores$PC1,
  phenotype_hist$PC1,
  method = "pearson"
)


PC1_TOLERANCE <- 1e-10


assert_true(
  pc1_max_difference <
    PC1_TOLERANCE,
  paste0(
    "Historical sample-level PC1 regression failed. ",
    "Maximum difference = ",
    format(
      pc1_max_difference,
      scientific = TRUE,
      digits = 17
    ),
    "."
  )
)


############################################################
## 5.23 — TYPE-7 PC1 TERTILES
############################################################

pc1_tertiles <- quantile(
  
  pca_scores$PC1,
  
  probs = c(
    1 / 3,
    2 / 3
  ),
  
  type = 7,
  
  names = FALSE
)


lower_tertile <- pc1_tertiles[1]
upper_tertile <- pc1_tertiles[2]


############################################################
## 5.24 — HISTORICAL TERTILE REGRESSION
############################################################

historical_tertiles <-
  pca_hist$PC1_tertiles_type7


TERTILE_TOLERANCE <- 1e-10


assert_true(
  max(
    abs(
      pc1_tertiles -
        historical_tertiles
    )
  ) <
    TERTILE_TOLERANCE,
  "Historical PC1 tertile regression failed."
)


############################################################
## 5.25 — CONSTRUCT IMMUNE PHENOTYPE
############################################################

immune_phenotype <- cut(
  
  pca_scores$PC1,
  
  breaks = c(
    -Inf,
    lower_tertile,
    upper_tertile,
    Inf
  ),
  
  labels = c(
    "Immune_Low",
    "Immune_Mid",
    "Immune_High"
  ),
  
  right = TRUE,
  
  include.lowest = TRUE,
  
  ordered_result = TRUE
)


immune_phenotype <- factor(
  
  immune_phenotype,
  
  levels = c(
    "Immune_Low",
    "Immune_Mid",
    "Immune_High"
  ),
  
  ordered = TRUE
)


############################################################
## 5.26 — PHENOTYPE COUNTS
############################################################

phenotype_counts <- table(
  immune_phenotype
)


assert_true(
  identical(
    as.integer(
      phenotype_counts
    ),
    c(
      369L,
      368L,
      369L
    )
  ),
  paste0(
    "Unexpected phenotype counts: ",
    paste(
      as.integer(
        phenotype_counts
      ),
      collapse = "/"
    ),
    "."
  )
)


############################################################
## 5.27 — EXACT HISTORICAL PHENOTYPE REGRESSION
############################################################

phenotype_matches <-
  as.character(
    immune_phenotype
  ) ==
  as.character(
    phenotype_hist$ImmunePhenotype
  )


phenotype_mismatches <-
  sum(
    !phenotype_matches
  )


assert_true(
  phenotype_mismatches == 0L,
  paste0(
    "Historical immune-phenotype regression failed: ",
    phenotype_mismatches,
    " samples differ."
  )
)


############################################################
## 5.28 — BUILD IMMUNE-AXIS TABLE
############################################################

immune_axis <- data.frame(
  
  sample =
    pca_scores$sample,
  
  PC1 =
    pca_scores$PC1,
  
  PC2 =
    pca_scores$PC2,
  
  PC3 =
    pca_scores$PC3,
  
  ImmunePhenotype =
    immune_phenotype,
  
  stringsAsFactors = FALSE
)


############################################################
## 5.29 — ADD ssGSEA + PCA + PHENOTYPE TO METADATA
############################################################

meta_immune <- meta_input


for (feature in ssgsea_features) {
  
  assert_true(
    !feature %in%
      colnames(meta_immune),
    paste0(
      "ssGSEA column already exists in metadata: ",
      feature
    )
  )
  
  meta_immune[[feature]] <-
    ssgsea_scores[[feature]]
}


meta_immune$PC1 <-
  pca_scores$PC1

meta_immune$PC2 <-
  pca_scores$PC2

meta_immune$PC3 <-
  pca_scores$PC3

meta_immune$ImmunePhenotype <-
  immune_phenotype


assert_identical_order(
  meta_immune$sample,
  colnames(expr_final),
  "immune-phenotype metadata",
  "expression samples"
)


############################################################
## 5.30 — FINAL MISSINGNESS
############################################################

assert_true(
  sum(
    is.na(
      meta_immune[
        ,
        ssgsea_features,
        drop = FALSE
      ]
    )
  ) == 0L,
  "Missing ssGSEA values detected in final metadata."
)

assert_true(
  sum(
    is.na(
      meta_immune$PC1
    )
  ) == 0L,
  "Missing PC1 values detected."
)

assert_true(
  sum(
    is.na(
      meta_immune$ImmunePhenotype
    )
  ) == 0L,
  "Missing immune phenotypes detected."
)


############################################################
## 5.31 — SAVE PRODUCTION OUTPUTS
############################################################

save_rds_checked(
  ssgsea_scores,
  file.path(
    PATHS$processed,
    "05_tcga_ssgsea17_scores.rds"
  )
)


save_rds_checked(
  ss_pca,
  file.path(
    PATHS$processed,
    "05_tcga_ssgsea_pca.rds"
  )
)


save_rds_checked(
  immune_axis,
  file.path(
    PATHS$processed,
    "05_tcga_immune_axis.rds"
  )
)


save_rds_checked(
  meta_immune,
  file.path(
    PATHS$processed,
    "05_tcga_metadata_immune_phenotype.rds"
  )
)


write_csv_checked(
  ssgsea_scores,
  file.path(
    PATHS$processed,
    "05_tcga_ssgsea17_scores.csv"
  )
)


write_csv_checked(
  immune_axis,
  file.path(
    PATHS$processed,
    "05_tcga_immune_axis.csv"
  )
)


############################################################
## 5.32 — SAVE PCA LOADINGS
############################################################

pc1_loading_table <- data.frame(
  
  Signature =
    names(pc1_loadings),
  
  PC1_loading =
    as.numeric(pc1_loadings),
  
  stringsAsFactors = FALSE
)


write_csv_checked(
  pc1_loading_table,
  file.path(
    PATHS$processed,
    "05_tcga_PC1_loadings.csv"
  )
)


############################################################
## 5.33 — SAVE ssGSEA HISTORICAL REGRESSION
############################################################

write_csv_checked(
  ssgsea_regression,
  file.path(
    PATHS$logs,
    "05_ssGSEA17_historical_regression.csv"
  )
)


############################################################
## 5.34 — SAVE PCA / PHENOTYPE REGRESSION
############################################################

pca_regression_summary <- data.frame(
  
  Metric = c(
    "PC1 variance percent",
    "Historical PC1 variance percent",
    "PC1 variance absolute difference",
    "PC1 loading max absolute difference",
    "PCA center max absolute difference",
    "PCA scale max absolute difference",
    "PC1 sample max absolute difference",
    "PC1 sample mean absolute difference",
    "PC1 sample Pearson correlation",
    "Lower tertile type 7",
    "Upper tertile type 7",
    "Immune Low count",
    "Immune Mid count",
    "Immune High count",
    "Phenotype exact matches",
    "Phenotype mismatches",
    "PC1 sign flipped during deterministic orientation"
  ),
  
  Value = c(
    pc1_variance_percent,
    pca_hist$PC1_variance_percent,
    variance_difference,
    loading_max_difference,
    center_max_difference,
    scale_max_difference,
    pc1_max_difference,
    pc1_mean_difference,
    pc1_correlation,
    lower_tertile,
    upper_tertile,
    phenotype_counts["Immune_Low"],
    phenotype_counts["Immune_Mid"],
    phenotype_counts["Immune_High"],
    sum(phenotype_matches),
    phenotype_mismatches,
    pc1_flipped
  ),
  
  stringsAsFactors = FALSE
)


write_csv_checked(
  pca_regression_summary,
  file.path(
    PATHS$logs,
    "05_PCA_phenotype_historical_regression.csv"
  )
)


############################################################
## 5.35 — FILE / SOFTWARE PROVENANCE
############################################################

signature_md5 <- unname(
  tools::md5sum(
    signature_file
  )
)

ssgsea_reference_md5 <- unname(
  tools::md5sum(
    ssgsea_reference_file
  )
)

phenotype_reference_md5 <- unname(
  tools::md5sum(
    phenotype_reference_file
  )
)

pca_reference_md5 <- unname(
  tools::md5sum(
    pca_reference_file
  )
)


############################################################
## 5.36 — MODULE QC SUMMARY
############################################################

module05_qc <- data.frame(
  
  Metric = c(
    "TCGA tumors",
    "Input expression genes",
    "ssGSEA signatures",
    "ssGSEA missing values",
    "ssGSEA historical max absolute difference",
    "ssGSEA regression tolerance",
    "ssGSEA regression passed",
    "PCA centered",
    "PCA scaled",
    "Positive PC1 loadings",
    "PC1 variance explained percent",
    "PC1 historical max absolute difference",
    "PC1 regression tolerance",
    "Lower PC1 tertile type 7",
    "Upper PC1 tertile type 7",
    "Immune Low tumors",
    "Immune Mid tumors",
    "Immune High tumors",
    "Historical phenotype exact matches",
    "Historical phenotype mismatches",
    "IOBR version",
    "GSVA version",
    "Historical signature RDS MD5",
    "Historical ssGSEA reference MD5",
    "Historical phenotype reference MD5",
    "Historical PCA reference MD5",
    "ssGSEA implementation",
    "PCA implementation",
    "PC1 orientation rule"
  ),
  
  Value = c(
    nrow(meta_immune),
    nrow(expr_final),
    length(ssgsea_features),
    sum(
      is.na(
        ssgsea_scores[
          ,
          ssgsea_features,
          drop = FALSE
        ]
      )
    ),
    format(
      ssgsea_global_max,
      scientific = TRUE,
      digits = 17
    ),
    SSGSEA_TOLERANCE,
    ssgsea_regression_pass,
    TRUE,
    TRUE,
    sum(pc1_loadings > 0),
    format(
      pc1_variance_percent,
      digits = 17
    ),
    format(
      pc1_max_difference,
      scientific = TRUE,
      digits = 17
    ),
    PC1_TOLERANCE,
    format(
      lower_tertile,
      digits = 17
    ),
    format(
      upper_tertile,
      digits = 17
    ),
    phenotype_counts["Immune_Low"],
    phenotype_counts["Immune_Mid"],
    phenotype_counts["Immune_High"],
    sum(phenotype_matches),
    phenotype_mismatches,
    as.character(
      packageVersion("IOBR")
    ),
    as.character(
      packageVersion("GSVA")
    ),
    signature_md5,
    ssgsea_reference_md5,
    phenotype_reference_md5,
    pca_reference_md5,
    "IOBR::calculate_sig_score(method = 'ssgsea')",
    "stats::prcomp(center = TRUE, scale. = TRUE)",
    "PC1 loading sum forced positive"
  ),
  
  stringsAsFactors = FALSE
)


write_csv_checked(
  module05_qc,
  file.path(
    PATHS$logs,
    "05_ssgsea_immune_axis_QC.csv"
  )
)


############################################################
## 5.37 — SESSION INFORMATION
############################################################

capture.output(
  sessionInfo(),
  file = file.path(
    PATHS$logs,
    "05_ssgsea_immune_axis_sessionInfo.txt"
  )
)


############################################################
## 5.38 — STATUS
############################################################

cat(
  "\n=============================================\n"
)

cat(
  "05_ssgsea_immune_axis.R: PASS\n"
)

cat(
  "TCGA tumors: ",
  nrow(meta_immune),
  "\n",
  sep = ""
)

cat(
  "Historical ssGSEA signatures: ",
  length(ssgsea_features),
  "\n",
  sep = ""
)

cat(
  "ssGSEA values evaluated: ",
  1106L * 17L,
  "\n",
  sep = ""
)

cat(
  "ssGSEA historical max difference: ",
  format(
    ssgsea_global_max,
    scientific = TRUE,
    digits = 17
  ),
  "\n",
  sep = ""
)

cat(
  "PC1 variance explained: ",
  format(
    pc1_variance_percent,
    digits = 17
  ),
  "%\n",
  sep = ""
)

cat(
  "Positive PC1 loadings: ",
  sum(pc1_loadings > 0),
  " / 17\n",
  sep = ""
)

cat(
  "PC1 historical max difference: ",
  format(
    pc1_max_difference,
    scientific = TRUE,
    digits = 17
  ),
  "\n",
  sep = ""
)

cat(
  "PC1 tertiles: ",
  format(
    lower_tertile,
    digits = 17
  ),
  " / ",
  format(
    upper_tertile,
    digits = 17
  ),
  "\n",
  sep = ""
)

cat(
  "Phenotype counts: ",
  phenotype_counts["Immune_Low"],
  " / ",
  phenotype_counts["Immune_Mid"],
  " / ",
  phenotype_counts["Immune_High"],
  "\n",
  sep = ""
)

cat(
  "Historical phenotype matches: ",
  sum(phenotype_matches),
  " / 1106\n",
  sep = ""
)

cat(
  "IOBR ssGSEA implementation: PRESERVED\n"
)

cat(
  "PCA center/scale method: PRESERVED\n"
)

cat(
  "PC1 orientation: DETERMINISTIC, HIGHER = GREATER SIGNATURE ACTIVITY\n"
)

cat(
  "=============================================\n"
)