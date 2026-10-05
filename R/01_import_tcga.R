############################################################
## 01_import_tcga.R
##
## PURPOSE
##   Import the exact TCGA-BRCA expression and clinical
##   files used in the exploratory analysis, align samples,
##   and retain Primary Tumor samples.
##
## IMPORTANT
##   This module preserves the historical import and
##   sample-selection logic. It does NOT perform gene
##   annotation or downstream preprocessing.
############################################################


############################################################
## 1.1 — REQUIRE SETUP
############################################################

if (!exists("PATHS")) {
  source("R/00_setup.R")
}


############################################################
## 1.2 — INPUT FILES
############################################################

expression_file <- file.path(
  PATHS$raw,
  "TCGA-BRCA.star_counts.tsv"
)

clinical_file <- file.path(
  PATHS$raw,
  "TCGA-BRCA.clinical.tsv"
)

assert_file(expression_file)
assert_file(clinical_file)


############################################################
## 1.3 — VERIFY INPUT FINGERPRINTS
############################################################

EXPECTED_EXPR_MD5 <-
  "7c6f02ee1c80277c2a3f8c07d9db6210"

EXPECTED_CLINICAL_MD5 <-
  "b66ea0f7089920e951444a4e71861e09"

observed_expr_md5 <- unname(
  tools::md5sum(expression_file)
)

observed_clinical_md5 <- unname(
  tools::md5sum(clinical_file)
)

assert_true(
  identical(
    observed_expr_md5,
    EXPECTED_EXPR_MD5
  ),
  "Expression input MD5 does not match the frozen source file."
)

assert_true(
  identical(
    observed_clinical_md5,
    EXPECTED_CLINICAL_MD5
  ),
  "Clinical input MD5 does not match the frozen source file."
)


############################################################
## 1.4 — HISTORICAL IMPORT METHOD
############################################################

## Preserve the original read.table() approach.
## The exploratory notebook used:
##
## expr <- read.table(
##   file.choose(),
##   header = TRUE,
##   row.names = 1,
##   sep = "\t"
## )
##
## meta <- read.table(
##   file.choose(),
##   header = TRUE,
##   sep = "\t"
## )
##
## Only file.choose() is replaced by deterministic paths.

expr <- read.table(
  expression_file,
  header = TRUE,
  row.names = 1,
  sep = "\t"
)

meta <- read.table(
  clinical_file,
  header = TRUE,
  sep = "\t"
)


############################################################
## 1.5 — RAW IMPORT QC
############################################################

assert_true(
  ncol(expr) == 1226L,
  paste0(
    "Expected 1226 expression samples; observed ",
    ncol(expr),
    "."
  )
)

assert_true(
  nrow(meta) == 1255L,
  paste0(
    "Expected 1255 clinical rows; observed ",
    nrow(meta),
    "."
  )
)

assert_true(
  "sample" %in% colnames(meta),
  "Clinical table is missing the 'sample' column."
)

assert_true(
  "sample_type.samples" %in% colnames(meta),
  "Clinical table is missing 'sample_type.samples'."
)

assert_unique(
  meta$sample,
  "clinical sample IDs"
)


############################################################
## 1.6 — PRESERVE HISTORICAL SAMPLE-ID NORMALIZATION
############################################################

expr_ids <- colnames(expr)

expr_ids_clean <- gsub(
  "\\.",
  "-",
  expr_ids
)

meta_ids <- meta$sample

common <- intersect(
  expr_ids_clean,
  meta_ids
)

assert_true(
  length(common) == 1226L,
  paste0(
    "Expected 1226 common expression/clinical samples; observed ",
    length(common),
    "."
  )
)


############################################################
## 1.7 — ALIGN EXPRESSION AND CLINICAL DATA
############################################################

colnames(expr) <- expr_ids_clean

expr <- expr[
  ,
  common,
  drop = FALSE
]

meta <- meta[
  meta$sample %in% common,
  ,
  drop = FALSE
]

meta <- meta[
  match(
    colnames(expr),
    meta$sample
  ),
  ,
  drop = FALSE
]


assert_identical_order(
  colnames(expr),
  meta$sample,
  "expression samples",
  "clinical samples"
)


############################################################
## 1.8 — SAMPLE-TYPE QC
############################################################

sample_type_counts <- table(
  meta$sample_type.samples,
  useNA = "ifany"
)

assert_true(
  unname(
    sample_type_counts["Primary Tumor"]
  ) == 1106L,
  "Primary Tumor count is not 1106."
)

assert_true(
  unname(
    sample_type_counts["Solid Tissue Normal"]
  ) == 113L,
  "Solid Tissue Normal count is not 113."
)

assert_true(
  unname(
    sample_type_counts["Metastatic"]
  ) == 7L,
  "Metastatic count is not 7."
)


############################################################
## 1.9 — HISTORICAL PRIMARY-TUMOR SELECTION
############################################################

tumor_idx <-
  meta$sample_type.samples == "Primary Tumor"

expr_tumor_raw <- expr[
  ,
  tumor_idx,
  drop = FALSE
]

meta_tumor_raw <- meta[
  tumor_idx,
  ,
  drop = FALSE
]


############################################################
## 1.10 — FINAL MODULE QC
############################################################

assert_true(
  ncol(expr_tumor_raw) == 1106L,
  paste0(
    "Expected 1106 Primary Tumor expression samples; observed ",
    ncol(expr_tumor_raw),
    "."
  )
)

assert_true(
  nrow(meta_tumor_raw) == 1106L,
  paste0(
    "Expected 1106 Primary Tumor clinical rows; observed ",
    nrow(meta_tumor_raw),
    "."
  )
)

assert_identical_order(
  colnames(expr_tumor_raw),
  meta_tumor_raw$sample,
  "tumor expression samples",
  "tumor clinical samples"
)

assert_unique(
  colnames(expr_tumor_raw),
  "Primary Tumor expression IDs"
)

assert_unique(
  meta_tumor_raw$sample,
  "Primary Tumor clinical IDs"
)


############################################################
## 1.11 — SAVE MODULE OUTPUTS
############################################################

save_rds_checked(
  expr_tumor_raw,
  file.path(
    PATHS$processed,
    "01_tcga_primary_tumor_expression_raw.rds"
  )
)

save_rds_checked(
  meta_tumor_raw,
  file.path(
    PATHS$processed,
    "01_tcga_primary_tumor_clinical_raw.rds"
  )
)


############################################################
## 1.12 — SAVE QC SUMMARY
############################################################

module01_qc <- data.frame(
  Metric = c(
    "Expression input MD5",
    "Clinical input MD5",
    "Raw expression samples",
    "Raw clinical rows",
    "Matched expression-clinical samples",
    "Primary Tumor",
    "Solid Tissue Normal",
    "Metastatic",
    "Final Primary Tumor expression samples",
    "Final Primary Tumor clinical rows",
    "Expression-clinical order identical"
  ),
  
  Value = c(
    observed_expr_md5,
    observed_clinical_md5,
    ncol(expr),
    nrow(meta),
    length(common),
    unname(sample_type_counts["Primary Tumor"]),
    unname(sample_type_counts["Solid Tissue Normal"]),
    unname(sample_type_counts["Metastatic"]),
    ncol(expr_tumor_raw),
    nrow(meta_tumor_raw),
    identical(
      colnames(expr_tumor_raw),
      meta_tumor_raw$sample
    )
  ),
  
  stringsAsFactors = FALSE
)

write_csv_checked(
  module01_qc,
  file.path(
    PATHS$logs,
    "01_import_tcga_QC.csv"
  )
)


############################################################
## 1.13 — STATUS
############################################################

cat(
  "\n=============================================\n"
)

cat(
  "01_import_tcga.R: PASS\n"
)

cat(
  "Raw expression samples: ",
  ncol(expr),
  "\n",
  sep = ""
)

cat(
  "Matched samples: ",
  length(common),
  "\n",
  sep = ""
)

cat(
  "Primary Tumor samples: ",
  ncol(expr_tumor_raw),
  "\n",
  sep = ""
)

cat(
  "Expression/clinical ordering: VERIFIED\n"
)

cat(
  "Historical sample-selection method: PRESERVED\n"
)

cat(
  "=============================================\n"
)