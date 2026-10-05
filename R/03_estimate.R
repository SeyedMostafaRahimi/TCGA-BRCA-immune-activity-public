############################################################
## 03_estimate.R
##
## PURPOSE
##   Reproduce the historical TCGA-BRCA ESTIMATE analysis
##   from the frozen Module-02 expression matrix.
##
## HISTORICAL METHOD
##   - Historical common-gene input: 9,859 genes
##   - estimateScore(..., platform = "illumina")
##   - ImmuneScore
##   - StromalScore
##   - ESTIMATEScore
##   - TumorPurity =
##       cos(0.6049872018 +
##           0.0001467884 * ESTIMATEScore)
##
## REPRODUCIBILITY REPAIRS
##   1. Historical 9,859-gene ESTIMATE input is frozen as
##      an explicit gene manifest.
##   2. Current estimate 1.0.13 does not automatically expose
##      SI_geneset to estimateScore(), so the original package
##      data are explicitly supplied to the unchanged function.
##   3. Historical accidental Description.1 pseudo-column is
##      not treated as a biological sample.
##   4. Final scores are regression-tested against the exact
##      historical 1,106-sample score reference.
##
## No change to the ESTIMATE analytical method is introduced.
############################################################


############################################################
## 3.1 — REQUIRE SETUP
############################################################

if (!exists("PATHS")) {
  source("R/00_setup.R")
}


############################################################
## 3.2 — REQUIRED PACKAGE
############################################################

require_package("estimate")


############################################################
## 3.3 — EXPLICIT INPUT FILES
############################################################

expression_file <- file.path(
  PATHS$processed,
  "02_tcga_expression_final.rds"
)

metadata_file <- file.path(
  PATHS$processed,
  "02_tcga_metadata_final.rds"
)

common_manifest_file <- file.path(
  "external",
  "estimate",
  "TCGA_BRCA_historical_ESTIMATE_common_genes.csv"
)

historical_reference_file <- file.path(
  "external",
  "estimate",
  "TCGA_BRCA_historical_ESTIMATE_scores_reference.csv"
)


assert_file(expression_file)
assert_file(metadata_file)
assert_file(common_manifest_file)
assert_file(historical_reference_file)


############################################################
## 3.4 — LOAD MODULE-02 INPUTS
############################################################

expr_final <- readRDS(
  expression_file
)

meta_final <- readRDS(
  metadata_file
)


############################################################
## 3.5 — MODULE-02 INPUT QC
############################################################

assert_true(
  nrow(expr_final) == 35009L,
  paste0(
    "Expected 35,009 genes; observed ",
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
  nrow(meta_final) == 1106L,
  paste0(
    "Expected 1,106 metadata rows; observed ",
    nrow(meta_final),
    "."
  )
)

assert_identical_order(
  colnames(expr_final),
  meta_final$sample,
  "expression samples",
  "clinical samples"
)

assert_unique(
  rownames(expr_final),
  "Module-02 gene symbols"
)

assert_unique(
  colnames(expr_final),
  "Module-02 sample IDs"
)


############################################################
## 3.6 — LOAD FROZEN HISTORICAL ESTIMATE GENE MANIFEST
############################################################

common_manifest <- read.csv(
  common_manifest_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

assert_true(
  all(
    c(
      "Historical_ESTIMATE_order",
      "Gene"
    ) %in% colnames(common_manifest)
  ),
  "ESTIMATE common-gene manifest is missing required columns."
)

assert_true(
  nrow(common_manifest) == 9859L,
  paste0(
    "Expected 9,859 historical ESTIMATE genes; observed ",
    nrow(common_manifest),
    "."
  )
)

assert_unique(
  common_manifest$Gene,
  "historical ESTIMATE common genes"
)

assert_true(
  identical(
    common_manifest$Historical_ESTIMATE_order,
    seq_len(9859L)
  ),
  "Historical ESTIMATE gene order is inconsistent."
)

assert_true(
  all(
    common_manifest$Gene %in%
      rownames(expr_final)
  ),
  "Not all historical ESTIMATE genes are present in expr_final."
)


############################################################
## 3.7 — RECONSTRUCT EXACT HISTORICAL ESTIMATE INPUT
############################################################

estimate_genes <- common_manifest$Gene

estimate_expr <- expr_final[
  estimate_genes,
  ,
  drop = FALSE
]


assert_true(
  nrow(estimate_expr) == 9859L,
  "ESTIMATE input does not contain 9,859 genes."
)

assert_true(
  ncol(estimate_expr) == 1106L,
  "ESTIMATE input does not contain 1,106 tumors."
)

assert_true(
  identical(
    rownames(estimate_expr),
    estimate_genes
  ),
  "ESTIMATE gene ordering is inconsistent."
)

assert_true(
  identical(
    colnames(estimate_expr),
    meta_final$sample
  ),
  "ESTIMATE sample ordering is inconsistent."
)

assert_true(
  !anyNA(estimate_expr),
  "ESTIMATE expression input contains missing values."
)


############################################################
## 3.8 — LOAD ORIGINAL ESTIMATE SI_geneset PACKAGE DATA
############################################################

estimate_data_env <- new.env(
  parent = emptyenv()
)

si_load <- try(
  utils::data(
    "SI_geneset",
    package = "estimate",
    envir = estimate_data_env
  ),
  silent = TRUE
)


assert_true(
  !inherits(si_load, "try-error"),
  "Could not load ESTIMATE SI_geneset package data."
)

assert_true(
  exists(
    "SI_geneset",
    envir = estimate_data_env,
    inherits = FALSE
  ),
  "ESTIMATE SI_geneset was not recovered."
)


SI_geneset_repro <- get(
  "SI_geneset",
  envir = estimate_data_env
)


assert_true(
  is.data.frame(SI_geneset_repro),
  "SI_geneset is not a data.frame."
)

assert_true(
  identical(
    dim(SI_geneset_repro),
    c(2L, 142L)
  ),
  paste0(
    "Unexpected SI_geneset dimensions: ",
    paste(dim(SI_geneset_repro), collapse = " x "),
    "."
  )
)


############################################################
## 3.9 — VERIFY SIGNATURE OVERLAP
############################################################

stromal_signature <- as.character(
  SI_geneset_repro[
    1,
    -1,
    drop = TRUE
  ]
)

immune_signature <- as.character(
  SI_geneset_repro[
    2,
    -1,
    drop = TRUE
  ]
)


stromal_overlap <- length(
  intersect(
    stromal_signature,
    estimate_genes
  )
)

immune_overlap <- length(
  intersect(
    immune_signature,
    estimate_genes
  )
)


assert_true(
  stromal_overlap == 136L,
  paste0(
    "Expected StromalSignature overlap = 136; observed ",
    stromal_overlap,
    "."
  )
)

assert_true(
  immune_overlap == 139L,
  paste0(
    "Expected ImmuneSignature overlap = 139; observed ",
    immune_overlap,
    "."
  )
)


############################################################
## 3.10 — CREATE COMPATIBILITY VERSION OF estimateScore()
############################################################

## IMPORTANT:
## The estimateScore() function body is NOT modified.
##
## Current estimate 1.0.13 does not automatically expose
## SI_geneset to the function. We provide the original
## package-data object through its evaluation environment.

estimateScore_repro <- estimate::estimateScore

estimate_score_env <- new.env(
  parent = environment(
    estimate::estimateScore
  )
)

estimate_score_env$SI_geneset <-
  SI_geneset_repro

environment(
  estimateScore_repro
) <- estimate_score_env


assert_true(
  identical(
    body(estimateScore_repro),
    body(estimate::estimateScore)
  ),
  "estimateScore function body was unexpectedly altered."
)


############################################################
## 3.11 — CREATE RUNTIME DIRECTORY
############################################################

runtime_dir <- file.path(
  PATHS$logs,
  "estimate_runtime"
)

dir.create(
  runtime_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


estimate_input_gct <- file.path(
  runtime_dir,
  "03_TCGA_ESTIMATE_commonGenes_input.gct"
)

estimate_output_gct <- file.path(
  runtime_dir,
  "03_TCGA_ESTIMATE_scores.gct"
)


if (file.exists(estimate_input_gct)) {
  file.remove(estimate_input_gct)
}

if (file.exists(estimate_output_gct)) {
  file.remove(estimate_output_gct)
}


############################################################
## 3.12 — WRITE CLEAN ESTIMATE GCT INPUT
############################################################

## Historical analysis used outputGCT().
## We retain the same GCT-writing function.
##
## The historical filtering workaround produced an accidental
## Description.1 pseudo-column. It was subsequently removed
## before downstream analysis.
##
## Production writes only the 9,859 real genes × 1,106 real
## tumor samples. This removes a file-format artifact, not an
## analytical variable.

estimate::outputGCT(
  estimate_expr,
  estimate_input_gct
)


assert_file(
  estimate_input_gct
)


############################################################
## 3.13 — VERIFY GCT DIMENSION HEADER
############################################################

gct_header <- readLines(
  estimate_input_gct,
  n = 3
)

assert_true(
  length(gct_header) >= 3L,
  "ESTIMATE input GCT header is incomplete."
)


gct_dimension_fields <- strsplit(
  gct_header[2],
  "\t",
  fixed = TRUE
)[[1]]


gct_gene_count <- suppressWarnings(
  as.integer(
    gct_dimension_fields[1]
  )
)

gct_sample_count <- suppressWarnings(
  as.integer(
    gct_dimension_fields[2]
  )
)


assert_true(
  gct_gene_count == 9859L,
  paste0(
    "GCT reports ",
    gct_gene_count,
    " genes instead of 9,859."
  )
)

assert_true(
  gct_sample_count == 1106L,
  paste0(
    "GCT reports ",
    gct_sample_count,
    " samples instead of 1,106."
  )
)


############################################################
## 3.14 — RUN HISTORICAL ESTIMATE METHOD
############################################################

estimateScore_repro(
  input.ds = estimate_input_gct,
  output.ds = estimate_output_gct,
  platform = "illumina"
)


assert_file(
  estimate_output_gct
)


############################################################
## 3.15 — READ ESTIMATE OUTPUT
############################################################

scores_raw <- read.table(
  estimate_output_gct,
  skip = 2,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)


assert_true(
  nrow(scores_raw) == 3L,
  paste0(
    "Expected 3 ESTIMATE score rows; observed ",
    nrow(scores_raw),
    "."
  )
)


############################################################
## 3.16 — IDENTIFY REAL SAMPLE COLUMNS
############################################################

technical_columns <- c(
  "NAME",
  "Name",
  "Description",
  "Description.1"
)


score_sample_columns <- setdiff(
  colnames(scores_raw),
  technical_columns
)


score_sample_ids <- gsub(
  "\\.",
  "-",
  score_sample_columns
)


assert_true(
  length(score_sample_ids) == 1106L,
  paste0(
    "Expected 1,106 ESTIMATE sample columns; observed ",
    length(score_sample_ids),
    "."
  )
)

assert_same_set(
  score_sample_ids,
  meta_final$sample,
  "ESTIMATE score samples",
  "Module-02 samples"
)

assert_identical_order(
  score_sample_ids,
  meta_final$sample,
  "ESTIMATE score samples",
  "Module-02 samples"
)


############################################################
## 3.17 — BUILD NUMERIC SCORE MATRIX
############################################################

score_matrix <- as.matrix(
  scores_raw[
    ,
    score_sample_columns,
    drop = FALSE
  ]
)

storage.mode(
  score_matrix
) <- "double"


rownames(score_matrix) <-
  scores_raw$NAME

colnames(score_matrix) <-
  score_sample_ids


required_scores <- c(
  "ImmuneScore",
  "StromalScore",
  "ESTIMATEScore"
)


assert_true(
  all(
    required_scores %in%
      rownames(score_matrix)
  ),
  "Required ESTIMATE score rows are missing."
)


score_matrix <- score_matrix[
  required_scores,
  meta_final$sample,
  drop = FALSE
]


############################################################
## 3.18 — CREATE FINAL ESTIMATE TABLE
############################################################

estimate_scores <- data.frame(
  
  sample =
    meta_final$sample,
  
  ImmuneScore =
    as.numeric(
      score_matrix[
        "ImmuneScore",
      ]
    ),
  
  StromalScore =
    as.numeric(
      score_matrix[
        "StromalScore",
      ]
    ),
  
  ESTIMATEScore =
    as.numeric(
      score_matrix[
        "ESTIMATEScore",
      ]
    ),
  
  stringsAsFactors = FALSE
)


############################################################
## 3.19 — HISTORICAL TUMOR-PURITY FORMULA
############################################################

estimate_scores$TumorPurity <-
  cos(
    0.6049872018 +
      0.0001467884 *
      estimate_scores$ESTIMATEScore
  )


############################################################
## 3.20 — FINAL SCORE INTEGRITY
############################################################

assert_true(
  nrow(estimate_scores) == 1106L,
  "Final ESTIMATE table does not contain 1,106 samples."
)

assert_unique(
  estimate_scores$sample,
  "ESTIMATE sample IDs"
)

assert_identical_order(
  estimate_scores$sample,
  meta_final$sample,
  "ESTIMATE samples",
  "clinical samples"
)

assert_no_missing(
  estimate_scores[
    ,
    c(
      "ImmuneScore",
      "StromalScore",
      "ESTIMATEScore",
      "TumorPurity"
    )
  ],
  "ESTIMATE features"
)


############################################################
## 3.21 — LOAD HISTORICAL SCORE REFERENCE
############################################################

historical_reference <- read.csv(
  historical_reference_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


assert_true(
  nrow(historical_reference) == 1106L,
  paste0(
    "Expected 1,106 historical reference samples; observed ",
    nrow(historical_reference),
    "."
  )
)

assert_identical_order(
  estimate_scores$sample,
  historical_reference$sample,
  "production ESTIMATE samples",
  "historical ESTIMATE samples"
)


############################################################
## 3.22 — HISTORICAL SCORE REGRESSION TEST
############################################################

features <- c(
  "ImmuneScore",
  "StromalScore",
  "ESTIMATEScore",
  "TumorPurity"
)


regression_summary <- data.frame(
  
  Feature = features,
  
  Max_abs_difference =
    NA_real_,
  
  Mean_abs_difference =
    NA_real_,
  
  Pearson_r =
    NA_real_,
  
  stringsAsFactors = FALSE
)


for (i in seq_along(features)) {
  
  feature <- features[i]
  
  production_values <-
    estimate_scores[[feature]]
  
  historical_values <-
    historical_reference[[feature]]
  
  
  regression_summary$Max_abs_difference[i] <-
    max(
      abs(
        production_values -
          historical_values
      )
    )
  
  
  regression_summary$Mean_abs_difference[i] <-
    mean(
      abs(
        production_values -
          historical_values
      )
    )
  
  
  regression_summary$Pearson_r[i] <-
    cor(
      production_values,
      historical_values,
      method = "pearson"
    )
}


global_max_difference <- max(
  regression_summary$Max_abs_difference
)


ESTIMATE_SCORE_TOLERANCE <- 1e-8


regression_pass <-
  global_max_difference <
  ESTIMATE_SCORE_TOLERANCE


assert_true(
  regression_pass,
  paste0(
    "Historical ESTIMATE regression failed. ",
    "Maximum absolute difference = ",
    format(
      global_max_difference,
      scientific = TRUE,
      digits = 17
    ),
    "."
  )
)


############################################################
## 3.23 — ADD ESTIMATE FEATURES TO METADATA
############################################################

## Preserve Module-02 ordering explicitly.
## No merge() is needed because sample order has already been
## verified as identical.

meta_estimate <- meta_final

meta_estimate$ImmuneScore <-
  estimate_scores$ImmuneScore

meta_estimate$StromalScore <-
  estimate_scores$StromalScore

meta_estimate$ESTIMATEScore <-
  estimate_scores$ESTIMATEScore

meta_estimate$TumorPurity <-
  estimate_scores$TumorPurity


assert_identical_order(
  meta_estimate$sample,
  colnames(expr_final),
  "ESTIMATE-augmented metadata",
  "expression samples"
)


############################################################
## 3.24 — SAVE PRODUCTION OUTPUTS
############################################################

save_rds_checked(
  estimate_scores,
  file.path(
    PATHS$processed,
    "03_tcga_estimate_scores.rds"
  )
)

save_rds_checked(
  meta_estimate,
  file.path(
    PATHS$processed,
    "03_tcga_metadata_estimate.rds"
  )
)


write_csv_checked(
  estimate_scores,
  file.path(
    PATHS$processed,
    "03_tcga_estimate_scores.csv"
  )
)


############################################################
## 3.25 — SAVE REGRESSION RESULTS
############################################################

write_csv_checked(
  regression_summary,
  file.path(
    PATHS$logs,
    "03_ESTIMATE_historical_regression.csv"
  )
)


############################################################
## 3.26 — PACKAGE / SIGNATURE PROVENANCE
############################################################

si_gmt_file <- system.file(
  "extdata",
  "SI_geneset.gmt",
  package = "estimate"
)


si_gmt_md5 <- if (
  nzchar(si_gmt_file) &&
  file.exists(si_gmt_file)
) {
  
  unname(
    tools::md5sum(
      si_gmt_file
    )
  )
  
} else {
  
  NA_character_
}


common_manifest_md5 <- unname(
  tools::md5sum(
    common_manifest_file
  )
)

historical_reference_md5 <- unname(
  tools::md5sum(
    historical_reference_file
  )
)


############################################################
## 3.27 — MODULE QC SUMMARY
############################################################

module03_qc <- data.frame(
  
  Metric = c(
    "TCGA tumors",
    "Input expression genes",
    "Historical ESTIMATE common genes",
    "StromalSignature overlap",
    "ImmuneSignature overlap",
    "Missing ImmuneScore",
    "Missing StromalScore",
    "Missing ESTIMATEScore",
    "Missing TumorPurity",
    "Historical regression max absolute difference",
    "Historical regression tolerance",
    "Historical regression passed",
    "estimate package version",
    "SI_geneset dimensions",
    "SI_geneset GMT MD5",
    "Historical common-gene manifest MD5",
    "Historical score-reference MD5",
    "estimateScore function body modified",
    "Historical Description.1 pseudo-column retained as sample"
  ),
  
  Value = c(
    nrow(estimate_scores),
    nrow(expr_final),
    length(estimate_genes),
    stromal_overlap,
    immune_overlap,
    sum(is.na(estimate_scores$ImmuneScore)),
    sum(is.na(estimate_scores$StromalScore)),
    sum(is.na(estimate_scores$ESTIMATEScore)),
    sum(is.na(estimate_scores$TumorPurity)),
    format(
      global_max_difference,
      scientific = TRUE,
      digits = 17
    ),
    ESTIMATE_SCORE_TOLERANCE,
    regression_pass,
    as.character(
      packageVersion("estimate")
    ),
    paste(
      dim(SI_geneset_repro),
      collapse = " x "
    ),
    si_gmt_md5,
    common_manifest_md5,
    historical_reference_md5,
    "NO",
    "NO"
  ),
  
  stringsAsFactors = FALSE
)


write_csv_checked(
  module03_qc,
  file.path(
    PATHS$logs,
    "03_estimate_QC.csv"
  )
)


############################################################
## 3.28 — SESSION INFORMATION
############################################################

capture.output(
  sessionInfo(),
  file = file.path(
    PATHS$logs,
    "03_estimate_sessionInfo.txt"
  )
)


############################################################
## 3.29 — REMOVE LARGE TEMPORARY INPUT AFTER SUCCESS
############################################################

## The GCT input is deterministically recreated every run
## and does not need to remain as a permanent 100+ MB file.

if (file.exists(estimate_input_gct)) {
  file.remove(estimate_input_gct)
}


############################################################
## 3.30 — STATUS
############################################################

cat(
  "\n=============================================\n"
)

cat(
  "03_estimate.R: PASS\n"
)

cat(
  "TCGA tumors: ",
  nrow(estimate_scores),
  "\n",
  sep = ""
)

cat(
  "Historical ESTIMATE genes: ",
  length(estimate_genes),
  "\n",
  sep = ""
)

cat(
  "StromalSignature overlap: ",
  stromal_overlap,
  "\n",
  sep = ""
)

cat(
  "ImmuneSignature overlap: ",
  immune_overlap,
  "\n",
  sep = ""
)

cat(
  "Missing ESTIMATE values: 0\n"
)

cat(
  "Historical regression max difference: ",
  format(
    global_max_difference,
    scientific = TRUE,
    digits = 17
  ),
  "\n",
  sep = ""
)

cat(
  "Historical score regression: VERIFIED\n"
)

cat(
  "estimateScore(platform = 'illumina'): PRESERVED\n"
)

cat(
  "TumorPurity formula: PRESERVED\n"
)

cat(
  "Description.1 pseudo-column: EXCLUDED FROM BIOLOGICAL SAMPLES\n"
)

cat(
  "=============================================\n"
)