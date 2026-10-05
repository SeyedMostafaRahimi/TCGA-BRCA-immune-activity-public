############################################################
## 00_setup.R
##
## TCGA-BRCA immune transcriptomics project
## Production pipeline
##
## PURPOSE
##   - establish reproducible project paths
##   - define immutable project constants
##   - define QC/assertion helpers
##   - establish phenotype terminology and colors
##   - establish deterministic random seed
##
## THIS SCRIPT DOES NOT PERFORM BIOLOGICAL ANALYSIS.
############################################################


############################################################
## 0.1 — PROJECT ROOT CHECK
############################################################

PROJECT_SENTINEL <- "PROJECT_ROOT"

if (!file.exists(PROJECT_SENTINEL)) {
  
  stop(
    paste0(
      "\nPROJECT ROOT NOT DETECTED.\n\n",
      "Open TCGA_BRCA_Immune_Project.Rproj in RStudio ",
      "before running the production pipeline.\n\n",
      "Current working directory:\n",
      getwd(),
      "\n"
    ),
    call. = FALSE
  )
}


PROJECT_ROOT <- normalizePath(
  ".",
  winslash = "/",
  mustWork = TRUE
)


############################################################
## 0.2 — REPRODUCIBILITY OPTIONS
############################################################

options(
  stringsAsFactors = FALSE,
  scipen = 999,
  warn = 1
)

GLOBAL_SEED <- 12345L

set.seed(
  GLOBAL_SEED
)


############################################################
## 0.3 — CANONICAL PHENOTYPE DEFINITION
############################################################

PHENOTYPE_LEVELS <- c(
  "Immune_Low",
  "Immune_Mid",
  "Immune_High"
)


PHENOTYPE_COLORS <- c(
  Immune_Low  = "#1F78B4",
  Immune_Mid  = "#A6A6A6",
  Immune_High = "#E31A1C"
)


PHENOTYPE_DISPLAY <- c(
  Immune_Low  = "Immune-Low",
  Immune_Mid  = "Immune-Mid",
  Immune_High = "Immune-High"
)


############################################################
## 0.4 — EXACT FINAL 17 ssGSEA SIGNATURES
############################################################

IMMUNE_SIGNATURES_17 <- c(
  
  "CD_8_T_effector",
  
  "IFNG_signature_Ayers_et_al",
  
  "T_cell_inflamed_GEP_Ayers_et_al",
  
  "CAF_Peng_et_al",
  
  "TAM_Peng_et_al",
  
  "MDSC_Peng_et_al",
  
  "Immune_Checkpoint",
  
  "EMT1",
  
  "EMT2",
  
  "EMT3",
  
  "Pan_F_TBRs",
  
  "CD8_T_cells_Bindea_et_al",
  
  "NK_cells_Bindea_et_al",
  
  "Cytotoxic_cells_Bindea_et_al",
  
  "Macrophages_Bindea_et_al",
  
  "DC_Bindea_et_al",
  
  "B_cells_Bindea_et_al"
)


stopifnot(
  length(IMMUNE_SIGNATURES_17) == 17L,
  !anyDuplicated(IMMUNE_SIGNATURES_17)
)


############################################################
## 0.5 — PROJECT PATHS
############################################################

PATHS <- list(
  
  raw =
    "data_raw",
  
  processed =
    "data_processed",
  
  results =
    "results",
  
  tables =
    file.path(
      "results",
      "tables"
    ),
  
  figures =
    file.path(
      "results",
      "figures"
    ),
  
  supplementary_tables =
    file.path(
      "results",
      "supplementary_tables"
    ),
  
  supplementary_figures =
    file.path(
      "results",
      "supplementary_figures"
    ),
  
  external =
    "external",
  
  cytoscape =
    file.path(
      "external",
      "cytoscape"
    ),
  
  docs =
    "docs",
  
  logs =
    "logs"
)


for (p in unname(PATHS)) {
  
  dir.create(
    p,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


############################################################
## 0.6 — CORE ASSERTION HELPERS
############################################################

assert_true <- function(
    condition,
    message
) {
  
  if (
    length(condition) != 1L ||
    is.na(condition) ||
    !isTRUE(condition)
  ) {
    
    stop(
      paste0(
        "\nQC FAILURE: ",
        message,
        "\n"
      ),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}


assert_file <- function(
    path,
    allow_empty = FALSE
) {
  
  assert_true(
    file.exists(path),
    paste0(
      "Required file does not exist: ",
      path
    )
  )
  
  if (!allow_empty) {
    
    assert_true(
      file.info(path)$size > 0,
      paste0(
        "File exists but is empty: ",
        path
      )
    )
  }
  
  invisible(TRUE)
}


assert_no_missing <- function(
    x,
    object_name = deparse(
      substitute(x)
    )
) {
  
  n_missing <- sum(
    is.na(x)
  )
  
  assert_true(
    n_missing == 0L,
    paste0(
      object_name,
      " contains ",
      n_missing,
      " missing value(s)."
    )
  )
  
  invisible(TRUE)
}


assert_unique <- function(
    x,
    object_name = deparse(
      substitute(x)
    )
) {
  
  n_dup <- sum(
    duplicated(x)
  )
  
  assert_true(
    n_dup == 0L,
    paste0(
      object_name,
      " contains ",
      n_dup,
      " duplicated value(s)."
    )
  )
  
  invisible(TRUE)
}


assert_same_set <- function(
    x,
    y,
    x_name = "x",
    y_name = "y"
) {
  
  missing_from_y <- setdiff(
    x,
    y
  )
  
  missing_from_x <- setdiff(
    y,
    x
  )
  
  assert_true(
    length(missing_from_y) == 0L &&
      length(missing_from_x) == 0L,
    paste0(
      x_name,
      " and ",
      y_name,
      " do not contain identical ID sets."
    )
  )
  
  invisible(TRUE)
}


assert_identical_order <- function(
    x,
    y,
    x_name = "x",
    y_name = "y"
) {
  
  assert_true(
    identical(
      as.character(x),
      as.character(y)
    ),
    paste0(
      x_name,
      " and ",
      y_name,
      " are not in identical order."
    )
  )
  
  invisible(TRUE)
}


############################################################
## 0.7 — SAFE SAMPLE MATCHING HELPER
############################################################

match_samples <- function(
    target_ids,
    source_ids,
    target_name = "target",
    source_name = "source"
) {
  
  assert_unique(
    target_ids,
    paste0(
      target_name,
      " IDs"
    )
  )
  
  assert_unique(
    source_ids,
    paste0(
      source_name,
      " IDs"
    )
  )
  
  idx <- match(
    target_ids,
    source_ids
  )
  
  assert_true(
    !anyNA(idx),
    paste0(
      "Some ",
      target_name,
      " IDs were not found in ",
      source_name,
      "."
    )
  )
  
  assert_identical_order(
    target_ids,
    source_ids[idx],
    target_name,
    paste0(
      source_name,
      "[matched]"
    )
  )
  
  idx
}


############################################################
## 0.8 — PACKAGE CHECK HELPER
############################################################

require_package <- function(
    package
) {
  
  if (
    !requireNamespace(
      package,
      quietly = TRUE
    )
  ) {
    
    stop(
      paste0(
        "\nRequired package is not installed: ",
        package,
        "\n\n",
        "Production scripts never install packages automatically.\n",
        "Install the package explicitly before continuing.\n"
      ),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}


############################################################
## 0.9 — SAFE OUTPUT HELPERS
############################################################

save_rds_checked <- function(
    object,
    path,
    compress = FALSE
) {
  
  dir.create(
    dirname(path),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  saveRDS(
    object,
    file = path,
    compress = compress
  )
  
  assert_file(
    path
  )
  
  invisible(path)
}


write_csv_checked <- function(
    object,
    path,
    row.names = FALSE
) {
  
  dir.create(
    dirname(path),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  write.csv(
    object,
    file = path,
    row.names = row.names
  )
  
  assert_file(
    path
  )
  
  invisible(path)
}


############################################################
## 0.10 — SESSION INFORMATION HELPER
############################################################

write_session_info <- function(
    filename = "sessionInfo.txt"
) {
  
  path <- file.path(
    PATHS$logs,
    filename
  )
  
  capture.output(
    sessionInfo(),
    file = path
  )
  
  assert_file(
    path
  )
  
  invisible(path)
}


############################################################
## 0.11 — VERIFY PROJECT STRUCTURE
############################################################

required_dirs <- unname(
  PATHS
)

for (p in required_dirs) {
  
  assert_true(
    dir.exists(p),
    paste0(
      "Required project directory missing: ",
      p
    )
  )
}


############################################################
## 0.12 — RECORD INITIAL SESSION
############################################################

write_session_info(
  "00_setup_sessionInfo.txt"
)


############################################################
## 0.13 — SETUP QC SUMMARY
############################################################

setup_qc <- data.frame(
  
  Check = c(
    "Project root detected",
    "Phenotype levels",
    "Final ssGSEA signatures",
    "Global seed",
    "Required directories",
    "Session-info log"
  ),
  
  Result = c(
    PROJECT_ROOT,
    paste(
      PHENOTYPE_LEVELS,
      collapse = " | "
    ),
    length(
      IMMUNE_SIGNATURES_17
    ),
    GLOBAL_SEED,
    length(
      required_dirs
    ),
    file.path(
      PATHS$logs,
      "00_setup_sessionInfo.txt"
    )
  ),
  
  stringsAsFactors = FALSE
)


write_csv_checked(
  setup_qc,
  file.path(
    PATHS$logs,
    "00_setup_QC.csv"
  )
)


############################################################
## 0.14 — FINAL STATUS
############################################################

cat(
  "\n=============================================\n"
)

cat(
  "00_setup.R: PASS\n"
)

cat(
  "Project root: ",
  PROJECT_ROOT,
  "\n",
  sep = ""
)

cat(
  "Final phenotype levels: ",
  paste(
    PHENOTYPE_LEVELS,
    collapse = " -> "
  ),
  "\n",
  sep = ""
)

cat(
  "Final ssGSEA signatures: ",
  length(
    IMMUNE_SIGNATURES_17
  ),
  "\n",
  sep = ""
)

cat(
  "Global seed: ",
  GLOBAL_SEED,
  "\n",
  sep = ""
)

cat(
  "=============================================\n"
)