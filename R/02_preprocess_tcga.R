############################################################
## 02_preprocess_tcga.R
##
## PURPOSE
##   Reproduce the historical TCGA-BRCA gene preprocessing
##   exactly from the Module-01 Primary Tumor matrix.
##
## HISTORICAL ANALYSIS POLICY
##   1. Ensembl version suffix removal
##   2. Ensembl -> SYMBOL mapping
##   3. Remove historically unmapped genes
##   4. Retain first occurrence of duplicated SYMBOL
##   5. Remove zero-variance genes
##
## REPRODUCIBILITY NOTE
##   The original mapping was performed with org.Hs.eg.db
##   using AnnotationDbi::mapIds(..., multiVals = "first").
##
##   Gene annotation databases change over time. Re-running
##   the same mapping policy with a later OrgDb release was
##   demonstrated to alter the gene universe.
##
##   Therefore the exact historical mapping output has been
##   frozen as an explicit Ensembl-to-SYMBOL manifest.
##
##   The variance filter is still recomputed directly from
##   the expression matrix.
##
##   No expression transformation is added here.
############################################################


############################################################
## 2.1 — REQUIRE SETUP
############################################################

if (!exists("PATHS")) {
  source("R/00_setup.R")
}


############################################################
## 2.2 — INPUT FILES
############################################################

expression_input <- file.path(
  PATHS$processed,
  "01_tcga_primary_tumor_expression_raw.rds"
)

clinical_input <- file.path(
  PATHS$processed,
  "01_tcga_primary_tumor_clinical_raw.rds"
)

manifest_file <- file.path(
  "external",
  "annotation",
  "TCGA_BRCA_historical_Ensembl_SYMBOL_manifest.csv"
)


assert_file(expression_input)
assert_file(clinical_input)
assert_file(manifest_file)


############################################################
## 2.3 — LOAD EXPLICIT INPUTS
############################################################

expr_tumor_raw <- readRDS(
  expression_input
)

meta_tumor_raw <- readRDS(
  clinical_input
)

manifest <- read.csv(
  manifest_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


############################################################
## 2.4 — MODULE-01 INPUT QC
############################################################

assert_true(
  nrow(expr_tumor_raw) == 60660L,
  paste0(
    "Expected 60,660 raw Ensembl rows; observed ",
    nrow(expr_tumor_raw),
    "."
  )
)

assert_true(
  ncol(expr_tumor_raw) == 1106L,
  paste0(
    "Expected 1,106 Primary Tumor samples; observed ",
    ncol(expr_tumor_raw),
    "."
  )
)

assert_true(
  nrow(meta_tumor_raw) == 1106L,
  paste0(
    "Expected 1,106 clinical rows; observed ",
    nrow(meta_tumor_raw),
    "."
  )
)

assert_identical_order(
  colnames(expr_tumor_raw),
  meta_tumor_raw$sample,
  "expression samples",
  "clinical samples"
)

assert_unique(
  rownames(expr_tumor_raw),
  "raw Ensembl IDs"
)

assert_unique(
  colnames(expr_tumor_raw),
  "Primary Tumor sample IDs"
)


############################################################
## 2.5 — MANIFEST STRUCTURE QC
############################################################

required_manifest_columns <- c(
  "Raw_row",
  "Ensembl_original",
  "Ensembl_clean",
  "Historical_SYMBOL",
  "Mapped_historical",
  "Kept_first_symbol",
  "Prevariance_order",
  "Kept_final",
  "Final_order"
)

assert_true(
  all(
    required_manifest_columns %in%
      colnames(manifest)
  ),
  "Historical annotation manifest is missing required columns."
)

assert_true(
  nrow(manifest) == 60660L,
  paste0(
    "Expected 60,660 manifest rows; observed ",
    nrow(manifest),
    "."
  )
)

assert_true(
  identical(
    manifest$Raw_row,
    seq_len(60660L)
  ),
  "Manifest Raw_row does not equal 1:60660."
)


############################################################
## 2.6 — VERIFY MANIFEST AGAINST RAW EXPRESSION
############################################################

raw_ensembl_original <- rownames(
  expr_tumor_raw
)

raw_ensembl_clean <- sub(
  "\\..*",
  "",
  raw_ensembl_original
)


assert_true(
  identical(
    manifest$Ensembl_original,
    raw_ensembl_original
  ),
  paste0(
    "Raw Ensembl IDs do not exactly match the ",
    "historical annotation manifest."
  )
)

assert_true(
  identical(
    manifest$Ensembl_clean,
    raw_ensembl_clean
  ),
  paste0(
    "Version-stripped Ensembl IDs do not exactly ",
    "match the historical manifest."
  )
)


############################################################
## 2.7 — HISTORICAL MAPPING COUNTS
############################################################

assert_true(
  sum(manifest$Mapped_historical) == 36257L,
  paste0(
    "Expected 36,257 historically mapped rows; observed ",
    sum(manifest$Mapped_historical),
    "."
  )
)

assert_true(
  sum(!manifest$Mapped_historical) == 24403L,
  paste0(
    "Expected 24,403 historically unmapped rows; observed ",
    sum(!manifest$Mapped_historical),
    "."
  )
)


############################################################
## 2.8 — VERIFY FIRST-DUPLICATE POLICY
############################################################

mapped_idx <- which(
  manifest$Mapped_historical
)

mapped_symbols <- manifest$Historical_SYMBOL[
  mapped_idx
]

assert_true(
  !anyNA(mapped_symbols),
  "Mapped historical rows contain missing SYMBOL values."
)


derived_keep_first <- !duplicated(
  mapped_symbols
)

manifest_keep_first_mapped <-
  manifest$Kept_first_symbol[
    mapped_idx
  ]


assert_true(
  identical(
    derived_keep_first,
    manifest_keep_first_mapped
  ),
  paste0(
    "Historical first-occurrence duplicate-symbol ",
    "policy cannot be reproduced from the manifest."
  )
)


assert_true(
  sum(!derived_keep_first) == 105L,
  paste0(
    "Expected 105 duplicate-symbol rows; observed ",
    sum(!derived_keep_first),
    "."
  )
)

assert_true(
  sum(derived_keep_first) == 36152L,
  paste0(
    "Expected 36,152 unique symbols before variance ",
    "filtering; observed ",
    sum(derived_keep_first),
    "."
  )
)


############################################################
## 2.9 — RECONSTRUCT HISTORICAL UNIQUE-SYMBOL MATRIX
############################################################

prevariance_raw_rows <- mapped_idx[
  derived_keep_first
]

prevariance_symbols <- mapped_symbols[
  derived_keep_first
]


assert_true(
  length(prevariance_raw_rows) == 36152L,
  "Unexpected pre-variance row count."
)

assert_unique(
  prevariance_symbols,
  "historical gene symbols before variance filtering"
)


## Verify recorded pre-variance ordering.

assert_true(
  identical(
    manifest$Prevariance_order[
      prevariance_raw_rows
    ],
    seq_len(36152L)
  ),
  "Manifest pre-variance ordering is inconsistent."
)


expr_symbol <- expr_tumor_raw[
  prevariance_raw_rows,
  ,
  drop = FALSE
]

rownames(expr_symbol) <- prevariance_symbols


############################################################
## 2.10 — EXPRESSION INTEGRITY BEFORE VARIANCE FILTER
############################################################

assert_true(
  !anyNA(expr_symbol),
  "Expression matrix contains missing values."
)

finite_by_sample <- vapply(
  expr_symbol,
  function(x) {
    all(is.finite(x))
  },
  logical(1)
)

assert_true(
  all(finite_by_sample),
  "Expression matrix contains non-finite values."
)


############################################################
## 2.11 — HISTORICAL ZERO-VARIANCE FILTER
############################################################

## Preserve original analytical implementation:
##
## gene_var <- apply(expr_final, 1, var)
## expr_final <- expr_final[gene_var > 0, ]

gene_var <- apply(
  expr_symbol,
  1,
  var
)

assert_true(
  !anyNA(gene_var),
  "Variance calculation produced missing values."
)


keep_variance <- gene_var > 0

n_zero_variance <- sum(
  !keep_variance
)


assert_true(
  n_zero_variance == 1143L,
  paste0(
    "Expected 1,143 zero-variance genes; observed ",
    n_zero_variance,
    "."
  )
)

assert_true(
  sum(keep_variance) == 35009L,
  paste0(
    "Expected 35,009 genes after variance filtering; observed ",
    sum(keep_variance),
    "."
  )
)


############################################################
## 2.12 — VERIFY VARIANCE FILTER AGAINST FROZEN MANIFEST
############################################################

derived_final_raw_rows <- prevariance_raw_rows[
  keep_variance
]

manifest_final_raw_rows <- which(
  manifest$Kept_final
)


assert_true(
  identical(
    derived_final_raw_rows,
    manifest_final_raw_rows
  ),
  paste0(
    "Recomputed zero-variance filtering does not reproduce ",
    "the frozen historical final-row selection."
  )
)


assert_true(
  identical(
    manifest$Final_order[
      manifest_final_raw_rows
    ],
    seq_len(35009L)
  ),
  "Manifest final-gene ordering is inconsistent."
)


############################################################
## 2.13 — CREATE FINAL EXPRESSION MATRIX
############################################################

expr_final <- expr_symbol[
  keep_variance,
  ,
  drop = FALSE
]

meta_final <- meta_tumor_raw


############################################################
## 2.14 — FINAL FROZEN CHECKPOINTS
############################################################

assert_true(
  nrow(expr_final) == 35009L,
  paste0(
    "FROZEN CHECK FAILED: expected 35,009 genes; observed ",
    nrow(expr_final),
    "."
  )
)

assert_true(
  ncol(expr_final) == 1106L,
  paste0(
    "FROZEN CHECK FAILED: expected 1,106 tumors; observed ",
    ncol(expr_final),
    "."
  )
)


assert_unique(
  rownames(expr_final),
  "final historical gene symbols"
)

assert_unique(
  colnames(expr_final),
  "final tumor sample IDs"
)


assert_identical_order(
  colnames(expr_final),
  meta_final$sample,
  "final expression samples",
  "final clinical samples"
)


assert_true(
  identical(
    rownames(expr_final),
    manifest$Historical_SYMBOL[
      manifest$Kept_final
    ]
  ),
  paste0(
    "Final expression gene order does not exactly match ",
    "the frozen historical annotation manifest."
  )
)


assert_true(
  !anyNA(expr_final),
  "Final expression matrix contains missing values."
)


############################################################
## 2.15 — SAVE PRODUCTION OUTPUTS
############################################################

save_rds_checked(
  expr_final,
  file.path(
    PATHS$processed,
    "02_tcga_expression_final.rds"
  )
)

save_rds_checked(
  meta_final,
  file.path(
    PATHS$processed,
    "02_tcga_metadata_final.rds"
  )
)


############################################################
## 2.16 — SAMPLE-LEVEL EXPRESSION QC
############################################################

sample_qc <- data.frame(
  
  sample =
    colnames(expr_final),
  
  median_expression =
    vapply(
      expr_final,
      median,
      numeric(1),
      na.rm = TRUE
    ),
  
  expression_IQR =
    vapply(
      expr_final,
      IQR,
      numeric(1),
      na.rm = TRUE
    ),
  
  stringsAsFactors = FALSE
)


write_csv_checked(
  sample_qc,
  file.path(
    PATHS$processed,
    "02_tcga_sample_expression_QC.csv"
  )
)


############################################################
## 2.17 — MODULE QC SUMMARY
############################################################

manifest_md5 <- unname(
  tools::md5sum(
    manifest_file
  )
)


module02_qc <- data.frame(
  
  Metric = c(
    "Raw Ensembl rows",
    "Historical mapped rows",
    "Historical unmapped rows",
    "Duplicate-symbol rows removed",
    "Unique symbols before variance filter",
    "Zero-variance genes removed",
    "Final genes",
    "Final tumors",
    "Final gene symbols unique",
    "Expression-clinical sample order identical",
    "Historical manifest MD5",
    "Live annotation database used to define final genes"
  ),
  
  Value = c(
    nrow(expr_tumor_raw),
    sum(manifest$Mapped_historical),
    sum(!manifest$Mapped_historical),
    sum(!derived_keep_first),
    nrow(expr_symbol),
    n_zero_variance,
    nrow(expr_final),
    ncol(expr_final),
    !anyDuplicated(rownames(expr_final)),
    identical(
      colnames(expr_final),
      meta_final$sample
    ),
    manifest_md5,
    "NO"
  ),
  
  stringsAsFactors = FALSE
)


write_csv_checked(
  module02_qc,
  file.path(
    PATHS$logs,
    "02_preprocess_tcga_QC.csv"
  )
)


############################################################
## 2.18 — SESSION INFORMATION
############################################################

capture.output(
  sessionInfo(),
  file = file.path(
    PATHS$logs,
    "02_preprocess_tcga_sessionInfo.txt"
  )
)


############################################################
## 2.19 — STATUS
############################################################

cat(
  "\n=============================================\n"
)

cat(
  "02_preprocess_tcga.R: PASS\n"
)

cat(
  "Raw Ensembl rows: ",
  nrow(expr_tumor_raw),
  "\n",
  sep = ""
)

cat(
  "Historically mapped rows: ",
  sum(manifest$Mapped_historical),
  "\n",
  sep = ""
)

cat(
  "Historically unmapped rows removed: ",
  sum(!manifest$Mapped_historical),
  "\n",
  sep = ""
)

cat(
  "Duplicate-symbol rows removed: ",
  sum(!derived_keep_first),
  "\n",
  sep = ""
)

cat(
  "Genes before variance filter: ",
  nrow(expr_symbol),
  "\n",
  sep = ""
)

cat(
  "Zero-variance genes removed: ",
  n_zero_variance,
  "\n",
  sep = ""
)

cat(
  "Final expression matrix: ",
  nrow(expr_final),
  " genes x ",
  ncol(expr_final),
  " tumors\n",
  sep = ""
)

cat(
  "Historical annotation manifest: VERIFIED\n"
)

cat(
  "Historical duplicate-symbol policy: PRESERVED\n"
)

cat(
  "Historical zero-variance filter: RECOMPUTED AND VERIFIED\n"
)

cat(
  "Additional expression transformation: NONE\n"
)

cat(
  "=============================================\n"
)