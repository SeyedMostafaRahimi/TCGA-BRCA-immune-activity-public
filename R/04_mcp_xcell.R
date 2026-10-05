############################################################
## 04_mcp_xcell.R
##
## PURPOSE
##   Reproduce the historical TCGA-BRCA MCP-counter and
##   xCell analyses using IOBR.
##
## HISTORICAL IMPLEMENTATION
##
## MCP-counter:
##   IOBR::deconvo_mcpcounter(
##     eset = expr_final,
##     project = "TCGA-BRCA"
##   )
##
## xCell:
##   IOBR::deconvo_xcell(
##     eset = expr_final,
##     arrays = FALSE
##   )
##
## IMPORTANT
##   IOBR is preserved as the implementation layer.
##   The analysis is NOT replaced by standalone MCPcounter
##   or standalone xCell calls.
##
## Historical regression references:
##   - 10 MCP-counter features
##   - 67 xCell features
##   - 1,106 TCGA Primary Tumor samples
############################################################


############################################################
## 4.1 — REQUIRE SETUP
############################################################

if (!exists("PATHS")) {
  source("R/00_setup.R")
}


############################################################
## 4.2 — REQUIRED PACKAGES
############################################################

require_package("IOBR")

## xCell is executed through IOBR but currently uses GSVA
## internally. Record its availability/version explicitly.

require_package("GSVA")


############################################################
## 4.3 — EXPLICIT INPUT FILES
############################################################

expression_file <- file.path(
  PATHS$processed,
  "02_tcga_expression_final.rds"
)

metadata_file <- file.path(
  PATHS$processed,
  "03_tcga_metadata_estimate.rds"
)

mcp_reference_file <- file.path(
  "external",
  "iobr",
  "TCGA_BRCA_historical_MCPcounter_reference.csv"
)

xcell_reference_file <- file.path(
  "external",
  "iobr",
  "TCGA_BRCA_historical_xCell_reference.csv"
)


assert_file(expression_file)
assert_file(metadata_file)
assert_file(mcp_reference_file)
assert_file(xcell_reference_file)


############################################################
## 4.4 — LOAD INPUTS
############################################################

expr_final <- readRDS(
  expression_file
)

meta_input <- readRDS(
  metadata_file
)

mcp_hist <- read.csv(
  mcp_reference_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

xcell_hist <- read.csv(
  xcell_reference_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


############################################################
## 4.5 — INPUT QC
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
  "metadata samples"
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
## 4.6 — HISTORICAL REFERENCE QC
############################################################

mcp_features <- setdiff(
  colnames(mcp_hist),
  "sample"
)

xcell_features <- setdiff(
  colnames(xcell_hist),
  "sample"
)


assert_true(
  length(mcp_features) == 10L,
  paste0(
    "Expected 10 MCP-counter reference features; observed ",
    length(mcp_features),
    "."
  )
)

assert_true(
  length(xcell_features) == 67L,
  paste0(
    "Expected 67 xCell reference features; observed ",
    length(xcell_features),
    "."
  )
)

assert_true(
  nrow(mcp_hist) == 1106L,
  "Historical MCP-counter reference does not contain 1,106 samples."
)

assert_true(
  nrow(xcell_hist) == 1106L,
  "Historical xCell reference does not contain 1,106 samples."
)

assert_unique(
  mcp_hist$sample,
  "historical MCP-counter sample IDs"
)

assert_unique(
  xcell_hist$sample,
  "historical xCell sample IDs"
)

assert_identical_order(
  mcp_hist$sample,
  meta_input$sample,
  "historical MCP-counter samples",
  "production samples"
)

assert_identical_order(
  xcell_hist$sample,
  meta_input$sample,
  "historical xCell samples",
  "production samples"
)


############################################################
## 4.7 — VERIFY HISTORICAL FUNCTIONS ARE AVAILABLE
############################################################

assert_true(
  exists(
    "deconvo_mcpcounter",
    envir = asNamespace("IOBR"),
    inherits = FALSE
  ),
  "IOBR::deconvo_mcpcounter() is unavailable."
)

assert_true(
  exists(
    "deconvo_xcell",
    envir = asNamespace("IOBR"),
    inherits = FALSE
  ),
  "IOBR::deconvo_xcell() is unavailable."
)


############################################################
## 4.8 — RUN HISTORICAL MCP-COUNTER METHOD
############################################################

cat(
  "\nRunning historical IOBR MCP-counter implementation...\n"
)

mcp_raw <- IOBR::deconvo_mcpcounter(
  eset = expr_final,
  project = "TCGA-BRCA"
)


############################################################
## 4.9 — MCP STRUCTURAL QC
############################################################

assert_true(
  is.data.frame(mcp_raw),
  "IOBR MCP-counter output is not a data.frame."
)

assert_true(
  nrow(mcp_raw) == 1106L,
  paste0(
    "Expected 1,106 MCP-counter rows; observed ",
    nrow(mcp_raw),
    "."
  )
)

assert_true(
  "ID" %in% colnames(mcp_raw),
  "MCP-counter output is missing ID column."
)

assert_true(
  all(
    mcp_features %in%
      colnames(mcp_raw)
  ),
  "MCP-counter output is missing historical features."
)

assert_identical_order(
  mcp_raw$ID,
  meta_input$sample,
  "MCP-counter samples",
  "production samples"
)


############################################################
## 4.10 — BUILD CLEAN MCP-COUNTER TABLE
############################################################

mcp_scores <- data.frame(
  sample = mcp_raw$ID,
  mcp_raw[
    ,
    mcp_features,
    drop = FALSE
  ],
  check.names = FALSE,
  stringsAsFactors = FALSE
)


assert_true(
  nrow(mcp_scores) == 1106L,
  "Clean MCP-counter table does not contain 1,106 samples."
)

assert_true(
  ncol(mcp_scores) == 11L,
  "Clean MCP-counter table does not contain sample + 10 features."
)

assert_no_missing(
  mcp_scores[
    ,
    mcp_features,
    drop = FALSE
  ],
  "MCP-counter scores"
)


############################################################
## 4.11 — MCP HISTORICAL REGRESSION
############################################################

mcp_regression <- data.frame(
  
  Feature =
    mcp_features,
  
  Max_abs_difference =
    NA_real_,
  
  Mean_abs_difference =
    NA_real_,
  
  Pearson_r =
    NA_real_,
  
  stringsAsFactors = FALSE
)


for (i in seq_along(mcp_features)) {
  
  feature <- mcp_features[i]
  
  production_values <-
    as.numeric(
      mcp_scores[[feature]]
    )
  
  historical_values <-
    as.numeric(
      mcp_hist[[feature]]
    )
  
  
  mcp_regression$Max_abs_difference[i] <-
    max(
      abs(
        production_values -
          historical_values
      )
    )
  
  
  mcp_regression$Mean_abs_difference[i] <-
    mean(
      abs(
        production_values -
          historical_values
      )
    )
  
  
  mcp_regression$Pearson_r[i] <-
    cor(
      production_values,
      historical_values,
      method = "pearson"
    )
}


mcp_global_max <- max(
  mcp_regression$Max_abs_difference
)

MCP_TOLERANCE <- 1e-10

mcp_regression_pass <-
  mcp_global_max <
  MCP_TOLERANCE


assert_true(
  mcp_regression_pass,
  paste0(
    "Historical MCP-counter regression failed. ",
    "Maximum absolute difference = ",
    format(
      mcp_global_max,
      scientific = TRUE,
      digits = 17
    ),
    "."
  )
)


############################################################
## 4.12 — RUN HISTORICAL xCELL METHOD
############################################################

cat(
  "\nRunning historical IOBR xCell implementation...\n"
)

xcell_raw <- IOBR::deconvo_xcell(
  eset = expr_final,
  arrays = FALSE
)


############################################################
## 4.13 — xCELL STRUCTURAL QC
############################################################

assert_true(
  is.data.frame(xcell_raw),
  "IOBR xCell output is not a data.frame."
)

assert_true(
  nrow(xcell_raw) == 1106L,
  paste0(
    "Expected 1,106 xCell rows; observed ",
    nrow(xcell_raw),
    "."
  )
)

assert_true(
  "ID" %in% colnames(xcell_raw),
  "xCell output is missing ID column."
)

assert_true(
  all(
    xcell_features %in%
      colnames(xcell_raw)
  ),
  "xCell output is missing historical features."
)

assert_identical_order(
  xcell_raw$ID,
  meta_input$sample,
  "xCell samples",
  "production samples"
)


############################################################
## 4.14 — BUILD CLEAN xCELL TABLE
############################################################

xcell_scores <- data.frame(
  sample = xcell_raw$ID,
  xcell_raw[
    ,
    xcell_features,
    drop = FALSE
  ],
  check.names = FALSE,
  stringsAsFactors = FALSE
)


assert_true(
  nrow(xcell_scores) == 1106L,
  "Clean xCell table does not contain 1,106 samples."
)

assert_true(
  ncol(xcell_scores) == 68L,
  "Clean xCell table does not contain sample + 67 features."
)

assert_no_missing(
  xcell_scores[
    ,
    xcell_features,
    drop = FALSE
  ],
  "xCell scores"
)


############################################################
## 4.15 — xCELL HISTORICAL REGRESSION
############################################################

xcell_regression <- data.frame(
  
  Feature =
    xcell_features,
  
  Max_abs_difference =
    NA_real_,
  
  Mean_abs_difference =
    NA_real_,
  
  Pearson_r =
    NA_real_,
  
  stringsAsFactors = FALSE
)


for (i in seq_along(xcell_features)) {
  
  feature <- xcell_features[i]
  
  production_values <-
    as.numeric(
      xcell_scores[[feature]]
    )
  
  historical_values <-
    as.numeric(
      xcell_hist[[feature]]
    )
  
  
  xcell_regression$Max_abs_difference[i] <-
    max(
      abs(
        production_values -
          historical_values
      )
    )
  
  
  xcell_regression$Mean_abs_difference[i] <-
    mean(
      abs(
        production_values -
          historical_values
      )
    )
  
  
  xcell_regression$Pearson_r[i] <-
    suppressWarnings(
      cor(
        production_values,
        historical_values,
        method = "pearson"
      )
    )
}


xcell_global_max <- max(
  xcell_regression$Max_abs_difference
)


## Historical/current observed discrepancy was approximately
## 4.47e-19. 1e-12 therefore remains an extremely strict
## numerical regression threshold while accommodating normal
## floating-point serialization/platform effects.

XCELL_TOLERANCE <- 1e-12


xcell_regression_pass <-
  xcell_global_max <
  XCELL_TOLERANCE


assert_true(
  xcell_regression_pass,
  paste0(
    "Historical xCell regression failed. ",
    "Maximum absolute difference = ",
    format(
      xcell_global_max,
      scientific = TRUE,
      digits = 17
    ),
    "."
  )
)


############################################################
## 4.16 — FEATURE ORDER CHECKPOINTS
############################################################

assert_true(
  identical(
    colnames(mcp_scores)[-1],
    mcp_features
  ),
  "MCP-counter feature order changed."
)

assert_true(
  identical(
    colnames(xcell_scores)[-1],
    xcell_features
  ),
  "xCell feature order changed."
)


############################################################
## 4.17 — ADD SCORES TO METADATA
############################################################

## Existing metadata already contains ESTIMATE features from
## Module 03.
##
## Sample ordering is identical, so no merge() is necessary.

assert_true(
  !any(
    mcp_features %in%
      colnames(meta_input)
  ),
  "MCP-counter columns already exist in input metadata."
)

assert_true(
  !any(
    xcell_features %in%
      colnames(meta_input)
  ),
  "xCell columns already exist in input metadata."
)


meta_tme <- meta_input


for (feature in mcp_features) {
  
  meta_tme[[feature]] <-
    mcp_scores[[feature]]
}


for (feature in xcell_features) {
  
  meta_tme[[feature]] <-
    xcell_scores[[feature]]
}


assert_identical_order(
  meta_tme$sample,
  colnames(expr_final),
  "TME metadata samples",
  "expression samples"
)


############################################################
## 4.18 — FINAL COMPLETENESS QC
############################################################

assert_true(
  sum(
    is.na(
      meta_tme[
        ,
        mcp_features,
        drop = FALSE
      ]
    )
  ) == 0L,
  "Missing MCP-counter values detected."
)


assert_true(
  sum(
    is.na(
      meta_tme[
        ,
        xcell_features,
        drop = FALSE
      ]
    )
  ) == 0L,
  "Missing xCell values detected."
)


############################################################
## 4.19 — SAVE PRODUCTION OUTPUTS
############################################################

save_rds_checked(
  mcp_scores,
  file.path(
    PATHS$processed,
    "04_tcga_mcpcounter_scores.rds"
  )
)

save_rds_checked(
  xcell_scores,
  file.path(
    PATHS$processed,
    "04_tcga_xcell_scores.rds"
  )
)

save_rds_checked(
  meta_tme,
  file.path(
    PATHS$processed,
    "04_tcga_metadata_tme.rds"
  )
)


write_csv_checked(
  mcp_scores,
  file.path(
    PATHS$processed,
    "04_tcga_mcpcounter_scores.csv"
  )
)

write_csv_checked(
  xcell_scores,
  file.path(
    PATHS$processed,
    "04_tcga_xcell_scores.csv"
  )
)


############################################################
## 4.20 — SAVE HISTORICAL REGRESSION AUDITS
############################################################

write_csv_checked(
  mcp_regression,
  file.path(
    PATHS$logs,
    "04_MCPcounter_historical_regression.csv"
  )
)

write_csv_checked(
  xcell_regression,
  file.path(
    PATHS$logs,
    "04_xCell_historical_regression.csv"
  )
)


############################################################
## 4.21 — REFERENCE FILE FINGERPRINTS
############################################################

mcp_reference_md5 <- unname(
  tools::md5sum(
    mcp_reference_file
  )
)

xcell_reference_md5 <- unname(
  tools::md5sum(
    xcell_reference_file
  )
)


############################################################
## 4.22 — MODULE QC SUMMARY
############################################################

module04_qc <- data.frame(
  
  Metric = c(
    "TCGA tumors",
    "Input expression genes",
    "MCP-counter features",
    "xCell features",
    "MCP-counter missing values",
    "xCell missing values",
    "MCP-counter historical max absolute difference",
    "MCP-counter regression tolerance",
    "MCP-counter regression passed",
    "xCell historical max absolute difference",
    "xCell regression tolerance",
    "xCell regression passed",
    "IOBR version",
    "GSVA version",
    "MCP historical reference MD5",
    "xCell historical reference MD5",
    "MCP implementation",
    "xCell implementation"
  ),
  
  Value = c(
    nrow(meta_tme),
    nrow(expr_final),
    length(mcp_features),
    length(xcell_features),
    sum(
      is.na(
        mcp_scores[
          ,
          mcp_features,
          drop = FALSE
        ]
      )
    ),
    sum(
      is.na(
        xcell_scores[
          ,
          xcell_features,
          drop = FALSE
        ]
      )
    ),
    format(
      mcp_global_max,
      scientific = TRUE,
      digits = 17
    ),
    MCP_TOLERANCE,
    mcp_regression_pass,
    format(
      xcell_global_max,
      scientific = TRUE,
      digits = 17
    ),
    XCELL_TOLERANCE,
    xcell_regression_pass,
    as.character(
      packageVersion("IOBR")
    ),
    as.character(
      packageVersion("GSVA")
    ),
    mcp_reference_md5,
    xcell_reference_md5,
    "IOBR::deconvo_mcpcounter",
    "IOBR::deconvo_xcell"
  ),
  
  stringsAsFactors = FALSE
)


write_csv_checked(
  module04_qc,
  file.path(
    PATHS$logs,
    "04_mcp_xcell_QC.csv"
  )
)


############################################################
## 4.23 — SESSION INFORMATION
############################################################

capture.output(
  sessionInfo(),
  file = file.path(
    PATHS$logs,
    "04_mcp_xcell_sessionInfo.txt"
  )
)


############################################################
## 4.24 — STATUS
############################################################

cat(
  "\n=============================================\n"
)

cat(
  "04_mcp_xcell.R: PASS\n"
)

cat(
  "TCGA tumors: ",
  nrow(meta_tme),
  "\n",
  sep = ""
)

cat(
  "MCP-counter features: ",
  length(mcp_features),
  "\n",
  sep = ""
)

cat(
  "xCell features: ",
  length(xcell_features),
  "\n",
  sep = ""
)

cat(
  "MCP-counter missing values: 0\n"
)

cat(
  "xCell missing values: 0\n"
)

cat(
  "MCP historical max difference: ",
  format(
    mcp_global_max,
    scientific = TRUE,
    digits = 17
  ),
  "\n",
  sep = ""
)

cat(
  "xCell historical max difference: ",
  format(
    xcell_global_max,
    scientific = TRUE,
    digits = 17
  ),
  "\n",
  sep = ""
)

cat(
  "MCP historical regression: VERIFIED\n"
)

cat(
  "xCell historical regression: VERIFIED\n"
)

cat(
  "MCP implementation: IOBR PRESERVED\n"
)

cat(
  "xCell implementation: IOBR PRESERVED\n"
)

cat(
  "=============================================\n"
)