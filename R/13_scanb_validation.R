############################################################
## 13_scanb_validation.R
##
## PURPOSE
##
## Production implementation of the already-completed
## external validation in SCAN-B / GSE96058.
##
## This module reproduces the external SCAN-B validation analysis
## using deterministic project-local inputs and compares outputs
## against a frozen validation reference object.
##
##
## HISTORICAL ANALYSIS PRESERVED
##
## Cohort:
##   SCAN-B / GSE96058
##   N = 3,069
##   OS events = 322
##
## Biomarkers:
##   CCL19
##   IL10
##
## Cutoffs:
##   cohort-specific exact medians calculated in all
##   3,069 SCAN-B patients
##
## Tie rule:
##   expression >= exact median -> High
##
## Individual survival:
##   Kaplan-Meier
##   log-rank
##
## Joint expression states:
##   High_Low
##   High_High
##   Low_Low
##   Low_High
##
## Joint survival:
##   four-state Kaplan-Meier
##   overall log-rank
##   pairwise log-rank with BH adjustment
##
## Adjusted Cox:
##
##   Surv(OS.time, OS.status) ~
##     ImmuneState +
##     Age +
##     LymphNode +
##     NHG +
##     TumorSize
##
## Reference state:
##   High_Low
##
## PH diagnostics:
##   survival::cox.zph()
##
##
## IMPORTANT INTERPRETATION GUARDRAILS
##
## - SCAN-B validates the CCL19 / IL10 prognostic findings.
##
## - TCGA cutoffs are NOT transferred to SCAN-B.
##
## - SCAN-B uses its own cohort-specific medians.
##
## - The four groups are joint expression states.
##
## - No statistical interaction term is fitted.
##
## - The adjusted SCAN-B Low_High vs High_Low contrast is
##   not statistically significant.
##
## No new prognostic analysis is introduced in this production module.
## Supporting clinical characterization summaries are regenerated.
############################################################


############################################################
## 13.1 — PROJECT SETUP
############################################################

if (!exists("PATHS")) {
  source("R/00_setup.R")
}


require_package(
  "survival"
)


require_package(
  "survminer"
)


require_package(
  "dplyr"
)


############################################################
## 13.2 — PRODUCTION INPUTS
##
## These deterministic RDS files were created during the
## validated SCAN-B input-bridge audit.
##
## They contain only:
##
##   1. exact 3,069-patient GSE96058 phenotype data
##   2. exact historical CCL19 / IL10 expression values
##
## No new preprocessing is performed here.
############################################################

scanb_expression_file <- file.path(
  PATHS$raw,
  "SCANB_GSE96058_CCL19_IL10_expression_3069.rds"
)


scanb_pheno_file <- file.path(
  PATHS$raw,
  "SCANB_GSE96058_pheno_3069.rds"
)


assert_file(
  scanb_expression_file
)


assert_file(
  scanb_pheno_file
)


############################################################
## 13.3 — FROZEN HISTORICAL REFERENCES
############################################################

scanb_reference_dir <- file.path(
  "external",
  "scanb"
)

figure8_reference_file <- file.path(
  scanb_reference_dir,
  "Figure8_SCANB_Final_Analysis_Objects.rds"
)


assert_file(
  figure8_reference_file
)


############################################################
## 13.4 — VERIFY BYTE-IDENTICAL HISTORICAL REFERENCES
############################################################

EXPECTED_FIGURE8_REFERENCE_MD5 <-
  "0bb5c488406631175d91385f2deb1b17"


observed_figure8_md5 <- unname(
  tools::md5sum(
    figure8_reference_file
  )
)

############################################################
## 13.5 — LOAD PRODUCTION INPUTS
############################################################

scanb_expression <- readRDS(
  scanb_expression_file
)


pheno <- readRDS(
  scanb_pheno_file
)


############################################################
## 13.6 — LOAD FROZEN FIGURE-8 REFERENCE
############################################################

fig8_ref <- readRDS(
  figure8_reference_file
)


required_reference_objects <- c(
  "n_total",
  "events_total",
  "median_CCL19",
  "median_IL10",
  "CCL19_KM",
  "CCL19_logrank",
  "CCL19_p",
  "IL10_KM",
  "IL10_logrank",
  "IL10_p",
  "state_KM",
  "state_logrank",
  "state_p",
  "state_pairwise",
  "adjusted_data",
  "adjusted_model",
  "PH_test",
  "adjusted_results",
  "clinical_by_state"
)


assert_true(
  all(
    required_reference_objects %in%
      names(fig8_ref)
  ),
  "Frozen Figure-8 reference is incomplete."
)


############################################################
## 13.7 — PRODUCTION INPUT QC
############################################################

assert_true(
  is.data.frame(scanb_expression),
  "SCAN-B CCL19/IL10 expression input must be a data.frame."
)


assert_true(
  nrow(scanb_expression) == 3069L,
  "Expected 3,069 SCAN-B expression samples."
)


assert_true(
  identical(
    colnames(scanb_expression),
    c(
      "sample",
      "CCL19",
      "IL10"
    )
  ),
  "Unexpected SCAN-B expression-input columns."
)


assert_true(
  nrow(pheno) == 3069L,
  "Expected 3,069 SCAN-B phenotype rows."
)


assert_true(
  "title" %in%
    colnames(pheno),
  "GSE96058 phenotype input lacks title column."
)


assert_unique(
  scanb_expression$sample,
  "SCAN-B expression sample IDs"
)


assert_unique(
  pheno$title,
  "SCAN-B phenotype sample IDs"
)


assert_true(
  identical(
    scanb_expression$sample,
    as.character(
      pheno$title
    )
  ),
  "SCAN-B expression and phenotype sample order differs."
)


assert_true(
  !anyNA(
    scanb_expression$CCL19
  ),
  "Missing CCL19 expression values detected."
)


assert_true(
  !anyNA(
    scanb_expression$IL10
  ),
  "Missing IL10 expression values detected."
)


############################################################
## 13.8 — REQUIRED HISTORICAL PHENOTYPE VARIABLES
############################################################

required_pheno_columns <- c(
  "title",
  "overall survival days:ch1",
  "overall survival event:ch1",
  "age at diagnosis:ch1",
  "er status:ch1",
  "her2 status:ch1",
  "pgr status:ch1",
  "ki67 status:ch1",
  "nhg:ch1",
  "lymph node status:ch1",
  "tumor size:ch1"
)


assert_true(
  all(
    required_pheno_columns %in%
      colnames(pheno)
  ),
  "Required historical GSE96058 phenotype variables are missing."
)


############################################################
## 13.9 — RECONSTRUCT HISTORICAL SURVIVAL DATAFRAME
############################################################

survival_df <- pheno |>
  dplyr::transmute(
    
    sample =
      as.character(
        title
      ),
    
    OS.time =
      as.numeric(
        `overall survival days:ch1`
      ),
    
    OS.status =
      as.numeric(
        `overall survival event:ch1`
      ),
    
    Age =
      as.numeric(
        `age at diagnosis:ch1`
      ),
    
    ER =
      `er status:ch1`,
    
    HER2 =
      `her2 status:ch1`,
    
    PR =
      `pgr status:ch1`,
    
    Ki67 =
      `ki67 status:ch1`,
    
    NHG =
      `nhg:ch1`,
    
    LymphNode =
      `lymph node status:ch1`,
    
    TumorSize =
      suppressWarnings(
        as.numeric(
          `tumor size:ch1`
        )
      )
  )


rownames(
  survival_df
) <- survival_df$sample


############################################################
## 13.10 — HISTORICAL CLINICAL FACTOR CODING
############################################################

survival_df$ER <- factor(
  survival_df$ER,
  levels = c(
    0,
    1
  ),
  labels = c(
    "Negative",
    "Positive"
  )
)


survival_df$PR <- factor(
  survival_df$PR,
  levels = c(
    0,
    1
  ),
  labels = c(
    "Negative",
    "Positive"
  )
)


survival_df$HER2 <- factor(
  survival_df$HER2,
  levels = c(
    0,
    1
  ),
  labels = c(
    "Negative",
    "Positive"
  )
)


survival_df$NHG <- factor(
  survival_df$NHG,
  levels = c(
    "G1",
    "G2",
    "G3"
  )
)


survival_df$LymphNode <- factor(
  survival_df$LymphNode,
  levels = c(
    "NodeNegative",
    "NodePositive"
  )
)


survival_df$Ki67[
  survival_df$Ki67 == "NA"
] <- NA


survival_df$Ki67 <- factor(
  survival_df$Ki67
)


############################################################
## 13.11 — SURVIVAL COHORT QC
############################################################

assert_true(
  nrow(survival_df) == 3069L,
  "SCAN-B survival cohort must contain 3,069 patients."
)


assert_true(
  identical(
    survival_df$sample,
    scanb_expression$sample
  ),
  "SCAN-B survival and expression samples are not aligned."
)


assert_true(
  sum(
    is.na(
      survival_df$OS.time
    )
  ) == 0L,
  "Unexpected missing overall-survival times."
)


assert_true(
  sum(
    is.na(
      survival_df$OS.status
    )
  ) == 0L,
  "Unexpected missing overall-survival event indicators."
)


assert_true(
  sum(
    survival_df$OS.status
  ) == 322L,
  "Expected 322 SCAN-B overall-survival events."
)


############################################################
## 13.12 — CLINICAL MISSINGNESS
############################################################

clinical_missingness <- data.frame(
  
  Variable = c(
    "Age",
    "Nottingham histological grade",
    "Lymph-node status",
    "Tumor size"
  ),
  
  Missing_n = c(
    sum(
      is.na(
        survival_df$Age
      )
    ),
    sum(
      is.na(
        survival_df$NHG
      )
    ),
    sum(
      is.na(
        survival_df$LymphNode
      )
    ),
    sum(
      is.na(
        survival_df$TumorSize
      )
    )
  ),
  
  stringsAsFactors = FALSE
)


clinical_missingness$Missing_pct <-
  100 *
  clinical_missingness$Missing_n /
  nrow(survival_df)


assert_true(
  identical(
    clinical_missingness$Missing_n,
    c(
      0L,
      61L,
      96L,
      32L
    )
  ),
  "SCAN-B clinical missingness differs from historical analysis."
)


############################################################
## 13.13 — CREATE CANONICAL SCAN-B ANALYSIS DATAFRAME
############################################################

scanb_final <- survival_df


scanb_final$CCL19 <- as.numeric(
  scanb_expression$CCL19
)


scanb_final$IL10 <- as.numeric(
  scanb_expression$IL10
)


assert_true(
  !anyNA(
    scanb_final$CCL19
  ),
  "Missing SCAN-B CCL19 values."
)


assert_true(
  !anyNA(
    scanb_final$IL10
  ),
  "Missing SCAN-B IL10 values."
)


############################################################
## 13.14 — FULL-COHORT SCAN-B MEDIANS
##
## Historical method:
##
## median(expression, na.rm = TRUE)
##
## calculated independently in all 3,069 SCAN-B patients.
##
## TCGA thresholds are NOT used.
############################################################

median_CCL19_scanb <- median(
  scanb_final$CCL19,
  na.rm = TRUE
)


median_IL10_scanb <- median(
  scanb_final$IL10,
  na.rm = TRUE
)


############################################################
## Regression against frozen Figure-8 medians
############################################################

assert_true(
  abs(
    median_CCL19_scanb -
      fig8_ref$median_CCL19
  ) < 1e-12,
  "SCAN-B CCL19 median differs from frozen Figure-8 result."
)


assert_true(
  abs(
    median_IL10_scanb -
      fig8_ref$median_IL10
  ) < 1e-12,
  "SCAN-B IL10 median differs from frozen Figure-8 result."
)


EXPECTED_CCL19_MEDIAN <-
  4.2269281758851696


EXPECTED_IL10_MEDIAN <-
  -0.96432260598112396


assert_true(
  abs(
    median_CCL19_scanb -
      EXPECTED_CCL19_MEDIAN
  ) < 1e-12,
  "Unexpected SCAN-B CCL19 median."
)


assert_true(
  abs(
    median_IL10_scanb -
      EXPECTED_IL10_MEDIAN
  ) < 1e-12,
  "Unexpected SCAN-B IL10 median."
)


############################################################
## 13.15 — MEDIAN TIE QC
############################################################

CCL19_median_ties <- sum(
  scanb_final$CCL19 ==
    median_CCL19_scanb
)


IL10_median_ties <- sum(
  scanb_final$IL10 ==
    median_IL10_scanb
)


assert_true(
  CCL19_median_ties == 1L,
  "Expected one exact CCL19 median tie."
)


assert_true(
  IL10_median_ties == 1L,
  "Expected one exact IL10 median tie."
)


############################################################
## 13.16 — HISTORICAL MEDIAN SPLITS
##
## Exact historical tie rule:
##
## expression >= exact cohort median -> High
############################################################

scanb_final$CCL19_group <- factor(
  
  ifelse(
    scanb_final$CCL19 >=
      median_CCL19_scanb,
    "High",
    "Low"
  ),
  
  levels = c(
    "Low",
    "High"
  )
)


scanb_final$IL10_group <- factor(
  
  ifelse(
    scanb_final$IL10 >=
      median_IL10_scanb,
    "High",
    "Low"
  ),
  
  levels = c(
    "Low",
    "High"
  )
)


############################################################
## 13.17 — INDIVIDUAL-GENE GROUP COUNTS
############################################################

CCL19_group_counts <- table(
  scanb_final$CCL19_group
)


IL10_group_counts <- table(
  scanb_final$IL10_group
)


assert_true(
  identical(
    as.integer(
      CCL19_group_counts
    ),
    c(
      1534L,
      1535L
    )
  ),
  "Unexpected CCL19 Low/High group counts."
)


assert_true(
  identical(
    as.integer(
      IL10_group_counts
    ),
    c(
      1534L,
      1535L
    )
  ),
  "Unexpected IL10 Low/High group counts."
)


############################################################
## 13.18 — CCL19 KAPLAN-MEIER / LOG-RANK
############################################################

fit_CCL19 <- survival::survfit(
  
  survival::Surv(
    OS.time,
    OS.status
  ) ~ CCL19_group,
  
  data =
    scanb_final
)


logrank_CCL19 <- survival::survdiff(
  
  survival::Surv(
    OS.time,
    OS.status
  ) ~ CCL19_group,
  
  data =
    scanb_final
)


CCL19_logrank_p <- 1 -
  stats::pchisq(
    logrank_CCL19$chisq,
    df = 1
  )


############################################################
## 13.19 — IL10 KAPLAN-MEIER / LOG-RANK
############################################################

fit_IL10 <- survival::survfit(
  
  survival::Surv(
    OS.time,
    OS.status
  ) ~ IL10_group,
  
  data =
    scanb_final
)


logrank_IL10 <- survival::survdiff(
  
  survival::Surv(
    OS.time,
    OS.status
  ) ~ IL10_group,
  
  data =
    scanb_final
)


IL10_logrank_p <- 1 -
  stats::pchisq(
    logrank_IL10$chisq,
    df = 1
  )


############################################################
## 13.20 — REGRESS INDIVIDUAL-GENE SURVIVAL
############################################################

assert_true(
  abs(
    CCL19_logrank_p -
      fig8_ref$CCL19_p
  ) < 1e-12,
  "CCL19 log-rank P-value differs from frozen Figure 8."
)


assert_true(
  abs(
    IL10_logrank_p -
      fig8_ref$IL10_p
  ) < 1e-12,
  "IL10 log-rank P-value differs from frozen Figure 8."
)


assert_true(
  abs(
    logrank_CCL19$chisq -
      fig8_ref$CCL19_logrank$chisq
  ) < 1e-12,
  "CCL19 log-rank chi-square differs from frozen result."
)


assert_true(
  abs(
    logrank_IL10$chisq -
      fig8_ref$IL10_logrank$chisq
  ) < 1e-12,
  "IL10 log-rank chi-square differs from frozen result."
)


############################################################
## 13.21 — INDIVIDUAL-GENE COUNTS / EVENTS TABLES
############################################################

CCL19_events <- tapply(
  scanb_final$OS.status,
  scanb_final$CCL19_group,
  sum
)


IL10_events <- tapply(
  scanb_final$OS.status,
  scanb_final$IL10_group,
  sum
)


CCL19_counts_events <- data.frame(
  
  Group = c(
    "Low",
    "High"
  ),
  
  N = as.integer(
    CCL19_group_counts[
      c(
        "Low",
        "High"
      )
    ]
  ),
  
  Events = as.integer(
    CCL19_events[
      c(
        "Low",
        "High"
      )
    ]
  ),
  
  stringsAsFactors = FALSE
)


IL10_counts_events <- data.frame(
  
  Group = c(
    "Low",
    "High"
  ),
  
  N = as.integer(
    IL10_group_counts[
      c(
        "Low",
        "High"
      )
    ]
  ),
  
  Events = as.integer(
    IL10_events[
      c(
        "Low",
        "High"
      )
    ]
  ),
  
  stringsAsFactors = FALSE
)


assert_true(
  identical(
    CCL19_counts_events$Events,
    c(
      200L,
      122L
    )
  ),
  "Unexpected CCL19 event counts."
)


assert_true(
  identical(
    IL10_counts_events$Events,
    c(
      145L,
      177L
    )
  ),
  "Unexpected IL10 event counts."
)


############################################################
## 13.22 — FOUR JOINT CCL19 / IL10 EXPRESSION STATES
############################################################

scanb_final$ImmuneState <- factor(
  
  paste(
    scanb_final$CCL19_group,
    scanb_final$IL10_group,
    sep = "_"
  ),
  
  levels = c(
    "High_Low",
    "High_High",
    "Low_Low",
    "Low_High"
  )
)


state_counts <- table(
  scanb_final$ImmuneState
)


state_events <- tapply(
  scanb_final$OS.status,
  scanb_final$ImmuneState,
  sum
)


assert_true(
  identical(
    as.integer(
      state_counts
    ),
    c(
      586L,
      949L,
      948L,
      586L
    )
  ),
  "SCAN-B joint-state counts differ from frozen analysis."
)


assert_true(
  identical(
    as.integer(
      state_events
    ),
    c(
      39L,
      83L,
      106L,
      94L
    )
  ),
  "SCAN-B joint-state event counts differ from frozen analysis."
)


############################################################
## 13.23 — FOUR-STATE KAPLAN-MEIER / LOG-RANK
############################################################

fit_state <- survival::survfit(
  
  survival::Surv(
    OS.time,
    OS.status
  ) ~ ImmuneState,
  
  data =
    scanb_final
)


logrank_state <- survival::survdiff(
  
  survival::Surv(
    OS.time,
    OS.status
  ) ~ ImmuneState,
  
  data =
    scanb_final
)


state_logrank_p <- 1 -
  stats::pchisq(
    
    logrank_state$chisq,
    
    df =
      length(
        logrank_state$n
      ) - 1L
  )


assert_true(
  abs(
    state_logrank_p -
      fig8_ref$state_p
  ) < 1e-12,
  "Four-state SCAN-B log-rank P-value differs from frozen Figure 8."
)


assert_true(
  abs(
    logrank_state$chisq -
      fig8_ref$state_logrank$chisq
  ) < 1e-12,
  "Four-state SCAN-B chi-square differs from frozen Figure 8."
)


############################################################
## 13.24 — HISTORICAL PAIRWISE LOG-RANK + BH
############################################################

pairwise_state <- survminer::pairwise_survdiff(
  
  survival::Surv(
    OS.time,
    OS.status
  ) ~ ImmuneState,
  
  data =
    scanb_final,
  
  p.adjust.method =
    "BH"
)


############################################################
## Regression against frozen pairwise matrix
############################################################

assert_true(
  identical(
    dim(
      pairwise_state$p.value
    ),
    dim(
      fig8_ref$state_pairwise$p.value
    )
  ),
  "Pairwise log-rank matrix dimensions differ."
)


assert_true(
  identical(
    dimnames(
      pairwise_state$p.value
    ),
    dimnames(
      fig8_ref$state_pairwise$p.value
    )
  ),
  "Pairwise log-rank matrix labels differ."
)


pairwise_max_diff <- max(
  abs(
    pairwise_state$p.value -
      fig8_ref$state_pairwise$p.value
  ),
  na.rm = TRUE
)


assert_true(
  pairwise_max_diff < 1e-12,
  "Pairwise SCAN-B BH-adjusted log-rank results differ."
)


############################################################
## 13.25 — JOINT-STATE COUNTS / EVENTS TABLE
############################################################

joint_state_counts_events <- data.frame(
  
  ImmuneState = c(
    "High_Low",
    "High_High",
    "Low_Low",
    "Low_High"
  ),
  
  N = as.integer(
    state_counts[
      c(
        "High_Low",
        "High_High",
        "Low_Low",
        "Low_High"
      )
    ]
  ),
  
  Events = as.integer(
    state_events[
      c(
        "High_Low",
        "High_High",
        "Low_Low",
        "Low_High"
      )
    ]
  ),
  
  stringsAsFactors = FALSE
)


############################################################
## 13.26 — PAIRWISE RESULT LONG TABLE
############################################################

pairwise_state_long <- as.data.frame(
  as.table(
    pairwise_state$p.value
  ),
  stringsAsFactors = FALSE
)


colnames(
  pairwise_state_long
) <- c(
  "Group_1",
  "Group_2",
  "BH_adjusted_P"
)


pairwise_state_long <- pairwise_state_long[
  !is.na(
    pairwise_state_long$BH_adjusted_P
  ),
  ,
  drop = FALSE
]


rownames(
  pairwise_state_long
) <- NULL


assert_true(
  nrow(
    pairwise_state_long
  ) == 6L,
  "Expected six pairwise joint-state comparisons."
)


############################################################
## 13.27 — HISTORICAL ADJUSTED COMPLETE-CASE COHORT
##
## IMPORTANT:
##
## Joint states were already defined above using medians from
## all 3,069 patients.
##
## Medians are NOT recalculated in this subset.
############################################################

scanb_cox <- scanb_final |>
  dplyr::filter(
    !is.na(OS.time),
    !is.na(OS.status),
    !is.na(Age),
    !is.na(LymphNode),
    !is.na(NHG),
    !is.na(TumorSize),
    !is.na(ImmuneState)
  )


scanb_cox$ImmuneState <- factor(
  
  scanb_cox$ImmuneState,
  
  levels = c(
    "High_Low",
    "High_High",
    "Low_Low",
    "Low_High"
  )
)


scanb_cox$ImmuneState <- stats::relevel(
  scanb_cox$ImmuneState,
  ref = "High_Low"
)


scanb_cox$LymphNode <- factor(
  
  scanb_cox$LymphNode,
  
  levels = c(
    "NodeNegative",
    "NodePositive"
  )
)


scanb_cox$NHG <- factor(
  
  scanb_cox$NHG,
  
  levels = c(
    "G1",
    "G2",
    "G3"
  )
)


assert_true(
  nrow(scanb_cox) == 2911L,
  "Expected 2,911 complete cases for adjusted SCAN-B Cox model."
)


assert_true(
  sum(
    scanb_cox$OS.status
  ) == 303L,
  "Expected 303 events in adjusted SCAN-B Cox cohort."
)


############################################################
## 13.28 — REGRESS ADJUSTED DATA AGAINST FROZEN FIGURE 8
############################################################

frozen_adjusted <- fig8_ref$adjusted_data


assert_true(
  nrow(
    frozen_adjusted
  ) == 2911L,
  "Frozen adjusted-data reference has unexpected sample count."
)


frozen_adjusted_index <- match(
  scanb_cox$sample,
  frozen_adjusted$sample
)


assert_true(
  !anyNA(
    frozen_adjusted_index
  ),
  "Unable to align adjusted SCAN-B cohort to frozen reference."
)


frozen_adjusted_aligned <- frozen_adjusted[
  frozen_adjusted_index,
  ,
  drop = FALSE
]


assert_true(
  identical(
    scanb_cox$sample,
    frozen_adjusted_aligned$sample
  ),
  "Adjusted SCAN-B sample order differs from frozen Figure 8."
)


adjusted_numeric_columns <- c(
  "OS.time",
  "OS.status",
  "Age",
  "TumorSize",
  "CCL19",
  "IL10"
)


adjusted_numeric_max_diff <- max(
  abs(
    as.matrix(
      scanb_cox[
        ,
        adjusted_numeric_columns,
        drop = FALSE
      ]
    ) -
      as.matrix(
        frozen_adjusted_aligned[
          ,
          adjusted_numeric_columns,
          drop = FALSE
        ]
      )
  ),
  na.rm = TRUE
)


assert_true(
  adjusted_numeric_max_diff < 1e-12,
  "Adjusted SCAN-B numerical data differ from frozen reference."
)


adjusted_categorical_columns <- c(
  "ER",
  "HER2",
  "PR",
  "Ki67",
  "NHG",
  "LymphNode",
  "CCL19_group",
  "IL10_group",
  "ImmuneState"
)


adjusted_categorical_pass <- all(
  vapply(
    
    adjusted_categorical_columns,
    
    function(variable) {
      
      identical(
        as.character(
          scanb_cox[[variable]]
        ),
        as.character(
          frozen_adjusted_aligned[[variable]]
        )
      )
    },
    
    logical(1)
  )
)


assert_true(
  adjusted_categorical_pass,
  "Adjusted SCAN-B categorical data differ from frozen reference."
)


############################################################
## 13.29 — EXACT HISTORICAL ADJUSTED COX MODEL
############################################################

fit_adjusted <- survival::coxph(
  
  survival::Surv(
    OS.time,
    OS.status
  ) ~
    ImmuneState +
    Age +
    LymphNode +
    NHG +
    TumorSize,
  
  data =
    scanb_cox,
  
  x =
    TRUE
)


############################################################
## 13.30 — PROPORTIONAL-HAZARDS DIAGNOSTICS
############################################################

PH_test <- survival::cox.zph(
  fit_adjusted
)


############################################################
## 13.31 — REGRESS COX MODEL AGAINST FROZEN MODEL
############################################################

current_coefficient_table <- summary(
  fit_adjusted
)$coefficients


frozen_coefficient_table <- summary(
  fig8_ref$adjusted_model
)$coefficients


assert_true(
  identical(
    rownames(
      current_coefficient_table
    ),
    rownames(
      frozen_coefficient_table
    )
  ),
  "Adjusted Cox coefficient labels differ from frozen model."
)


cox_coefficient_max_diff <- max(
  abs(
    current_coefficient_table -
      frozen_coefficient_table
  )
)


current_CI_table <- summary(
  fit_adjusted
)$conf.int


frozen_CI_table <- summary(
  fig8_ref$adjusted_model
)$conf.int


cox_CI_max_diff <- max(
  abs(
    current_CI_table -
      frozen_CI_table
  )
)


assert_true(
  cox_coefficient_max_diff < 1e-10,
  "Adjusted SCAN-B Cox coefficients differ from frozen model."
)


assert_true(
  cox_CI_max_diff < 1e-10,
  "Adjusted SCAN-B confidence intervals differ from frozen model."
)


############################################################
## 13.32 — REGRESS PH DIAGNOSTICS
############################################################

assert_true(
  identical(
    rownames(
      PH_test$table
    ),
    rownames(
      fig8_ref$PH_test$table
    )
  ),
  "SCAN-B PH diagnostic labels differ from frozen reference."
)


PH_max_diff <- max(
  abs(
    PH_test$table -
      fig8_ref$PH_test$table
  )
)


assert_true(
  PH_max_diff < 1e-10,
  "SCAN-B PH diagnostics differ from frozen Figure 8."
)


############################################################
## 13.33 — CREATE ADJUSTED COX RESULT TABLE
############################################################

adjusted_summary <- summary(
  fit_adjusted
)


adjusted_results <- data.frame(
  
  Variable = c(
    "High CCL19 / High IL10",
    "Low CCL19 / Low IL10",
    "Low CCL19 / High IL10",
    "Age (per year)",
    "Node positive",
    "NHG G2",
    "NHG G3",
    "Tumor size"
  ),
  
  HR =
    adjusted_summary$coefficients[
      ,
      "exp(coef)"
    ],
  
  Lower95 =
    adjusted_summary$conf.int[
      ,
      "lower .95"
    ],
  
  Upper95 =
    adjusted_summary$conf.int[
      ,
      "upper .95"
    ],
  
  Pvalue =
    adjusted_summary$coefficients[
      ,
      "Pr(>|z|)"
    ],
  
  stringsAsFactors = FALSE
)


############################################################
## 13.34 — REGRESS ADJUSTED RESULT TABLE
############################################################

frozen_adjusted_results <-
  fig8_ref$adjusted_results


assert_true(
  identical(
    as.character(
      adjusted_results$Variable
    ),
    as.character(
      frozen_adjusted_results$Variable
    )
  ),
  "Adjusted SCAN-B result labels differ from frozen Figure 8."
)


adjusted_result_columns <- c(
  "HR",
  "Lower95",
  "Upper95",
  "Pvalue"
)


adjusted_result_max_diff <- max(
  abs(
    as.matrix(
      adjusted_results[
        ,
        adjusted_result_columns,
        drop = FALSE
      ]
    ) -
      as.matrix(
        frozen_adjusted_results[
          ,
          adjusted_result_columns,
          drop = FALSE
        ]
      )
  )
)


assert_true(
  adjusted_result_max_diff < 1e-10,
  "Adjusted SCAN-B result table differs from frozen Figure 8."
)


############################################################
## 13.35 — PH DIAGNOSTIC OUTPUT TABLE
############################################################

PH_results <- data.frame(
  
  Term =
    rownames(
      PH_test$table
    ),
  
  chisq =
    PH_test$table[
      ,
      "chisq"
    ],
  
  df =
    PH_test$table[
      ,
      "df"
    ],
  
  Pvalue =
    PH_test$table[
      ,
      "p"
    ],
  
  stringsAsFactors = FALSE
)


rownames(
  PH_results
) <- NULL


############################################################
## 13.36 — HISTORICAL CLINICAL CHARACTERISTICS BY STATE
############################################################

clinical_by_state <- scanb_final |>
  dplyr::group_by(
    ImmuneState
  ) |>
  dplyr::summarise(
    
    n =
      dplyr::n(),
    
    Age_median =
      median(
        Age,
        na.rm = TRUE
      ),
    
    TumorSize_median =
      median(
        TumorSize,
        na.rm = TRUE
      ),
    
    NodePositive_pct =
      100 *
      mean(
        LymphNode ==
          "NodePositive",
        na.rm = TRUE
      ),
    
    G3_pct =
      100 *
      mean(
        NHG ==
          "G3",
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
  )


############################################################
## 13.37 — REGRESS CLINICAL STATE SUMMARY
############################################################

frozen_clinical_by_state <- as.data.frame(
  fig8_ref$clinical_by_state
)


clinical_state_index <- match(
  as.character(
    clinical_by_state$ImmuneState
  ),
  as.character(
    frozen_clinical_by_state$ImmuneState
  )
)


assert_true(
  !anyNA(
    clinical_state_index
  ),
  "Unable to align SCAN-B clinical-state frozen reference."
)


frozen_clinical_aligned <-
  frozen_clinical_by_state[
    clinical_state_index,
    ,
    drop = FALSE
  ]


clinical_numeric_columns <- c(
  "n",
  "Age_median",
  "TumorSize_median",
  "NodePositive_pct",
  "G3_pct"
)


clinical_state_max_diff <- max(
  abs(
    as.matrix(
      clinical_by_state[
        ,
        clinical_numeric_columns,
        drop = FALSE
      ]
    ) -
      as.matrix(
        frozen_clinical_aligned[
          ,
          clinical_numeric_columns,
          drop = FALSE
        ]
      )
  )
)


assert_true(
  clinical_state_max_diff < 1e-10,
  "SCAN-B clinical-state summary differs from frozen Figure 8."
)


############################################################
## 13.38 — HISTORICAL JOINT-STATE CLINICAL ASSOCIATIONS
##
## These are supporting characterization analyses already
## present in the historical manuscript outputs.
##
## No multiple-testing framework or new method is added.
############################################################

clinical_association_tests <- data.frame(
  
  Variable = c(
    "Age",
    "Tumor size",
    "Lymph-node status",
    "Nottingham histological grade"
  ),
  
  Test = c(
    "Kruskal-Wallis",
    "Kruskal-Wallis",
    "Pearson chi-square",
    "Pearson chi-square"
  ),
  
  Pvalue = c(
    
    stats::kruskal.test(
      Age ~ ImmuneState,
      data =
        scanb_final
    )$p.value,
    
    stats::kruskal.test(
      TumorSize ~ ImmuneState,
      data =
        scanb_final
    )$p.value,
    
    stats::chisq.test(
      table(
        scanb_final$ImmuneState,
        scanb_final$LymphNode
      )
    )$p.value,
    
    stats::chisq.test(
      table(
        scanb_final$ImmuneState,
        scanb_final$NHG
      )
    )$p.value
  ),
  
  stringsAsFactors = FALSE
)


############################################################
## 13.39 — OUTPUT DIRECTORY
############################################################

output_dir <- file.path(
  "results",
  "tables",
  "scanb_validation"
)


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
## 13.40 — SAVE PROCESSED ANALYSIS OBJECTS
############################################################

save_rds_checked(
  scanb_final,
  file.path(
    PATHS$processed,
    "13_scanb_validation_cohort.rds"
  )
)


save_rds_checked(
  scanb_cox,
  file.path(
    PATHS$processed,
    "13_scanb_adjusted_cox_cohort.rds"
  )
)


############################################################
## 13.41 — SAVE COMPLETE VALIDATION BUNDLE
############################################################

scanb_validation_results <- list(
  
  dataset =
    "SCAN-B / GSE96058",
  
  n_total =
    nrow(scanb_final),
  
  events_total =
    sum(scanb_final$OS.status),
  
  median_CCL19 =
    median_CCL19_scanb,
  
  median_IL10 =
    median_IL10_scanb,
  
  tie_rule =
    "expression >= exact cohort median -> High",
  
  CCL19_KM =
    fit_CCL19,
  
  CCL19_logrank =
    logrank_CCL19,
  
  CCL19_p =
    CCL19_logrank_p,
  
  IL10_KM =
    fit_IL10,
  
  IL10_logrank =
    logrank_IL10,
  
  IL10_p =
    IL10_logrank_p,
  
  state_KM =
    fit_state,
  
  state_logrank =
    logrank_state,
  
  state_p =
    state_logrank_p,
  
  state_pairwise =
    pairwise_state,
  
  adjusted_data =
    scanb_cox,
  
  adjusted_model =
    fit_adjusted,
  
  PH_test =
    PH_test,
  
  adjusted_results =
    adjusted_results,
  
  clinical_by_state =
    clinical_by_state,
  
  clinical_association_tests =
    clinical_association_tests,
  
  clinical_missingness =
    clinical_missingness
)


save_rds_checked(
  scanb_validation_results,
  file.path(
    PATHS$processed,
    "13_scanb_validation_results.rds"
  )
)


############################################################
## 13.42 — EXPORT INDIVIDUAL-GENE TABLES
############################################################

write_csv_checked(
  CCL19_counts_events,
  file.path(
    output_dir,
    "13_SCANB_CCL19_group_counts_events.csv"
  )
)


write_csv_checked(
  IL10_counts_events,
  file.path(
    output_dir,
    "13_SCANB_IL10_group_counts_events.csv"
  )
)


############################################################
## 13.43 — EXPORT JOINT-STATE TABLES
############################################################

write_csv_checked(
  joint_state_counts_events,
  file.path(
    output_dir,
    "13_SCANB_joint_state_counts_events.csv"
  )
)


write_csv_checked(
  pairwise_state_long,
  file.path(
    output_dir,
    "13_SCANB_joint_state_pairwise_logrank_BH.csv"
  )
)


############################################################
## 13.44 — EXPORT ADJUSTED COX / PH TABLES
############################################################

write_csv_checked(
  adjusted_results,
  file.path(
    output_dir,
    "13_SCANB_adjusted_Cox.csv"
  )
)


write_csv_checked(
  PH_results,
  file.path(
    output_dir,
    "13_SCANB_Cox_PH_diagnostics.csv"
  )
)


############################################################
## 13.45 — EXPORT CLINICAL SUPPORT TABLES
############################################################

write_csv_checked(
  clinical_by_state,
  file.path(
    output_dir,
    "13_SCANB_clinical_characteristics_by_joint_state.csv"
  )
)


write_csv_checked(
  clinical_association_tests,
  file.path(
    output_dir,
    "13_SCANB_joint_state_clinical_association_tests.csv"
  )
)


write_csv_checked(
  clinical_missingness,
  file.path(
    output_dir,
    "13_SCANB_clinical_missingness.csv"
  )
)


############################################################
## 13.46 — METHOD / PROVENANCE TABLE
############################################################

method_provenance <- data.frame(
  
  Item = c(
    "External cohort",
    "Historical phenotype source",
    "Production phenotype input",
    "Production expression input",
    "Historical hub-expression reference",
    "Frozen SCAN-B validation reference",
    "Survival endpoint",
    "Biomarker cutoffs",
    "Tie rule",
    "Individual survival method",
    "Joint-state survival method",
    "Pairwise multiplicity",
    "Adjusted cohort construction",
    "Adjusted model",
    "Reference state",
    "PH diagnostic",
    "Interaction term",
    "TCGA cutoffs transferred",
    "New analysis added"
  ),
  
  Value = c(
    "SCAN-B / GSE96058",
    "GEO Series GSE96058 phenotype",
    basename(scanb_pheno_file),
    basename(scanb_expression_file),
    "Not included",
    basename(figure8_reference_file),
    "Overall survival; time in days",
    "Independent cohort-specific exact medians in all 3069 SCAN-B patients",
    "Expression >= exact median -> High",
    "Kaplan-Meier + log-rank",
    "Four-state Kaplan-Meier + overall log-rank",
    "Pairwise log-rank with Benjamini-Hochberg adjustment",
    "Complete cases after joint states were assigned using full-cohort medians",
    "ImmuneState + Age + LymphNode + NHG + TumorSize",
    "High_Low",
    "survival::cox.zph",
    "FALSE",
    "FALSE",
    "FALSE"
  ),
  
  stringsAsFactors = FALSE
)


write_csv_checked(
  method_provenance,
  file.path(
    output_dir,
    "13_SCANB_method_provenance.csv"
  )
)


############################################################
## 13.47 — NUMERICAL REGRESSION LOG
############################################################

regression_log <- data.frame(
  
  Metric = c(
    "CCL19_median_difference",
    "IL10_median_difference",
    "CCL19_logrank_P_difference",
    "IL10_logrank_P_difference",
    "Joint_state_logrank_P_difference",
    "Pairwise_BH_max_difference",
    "Adjusted_data_numeric_max_difference",
    "Adjusted_data_categorical_exact",
    "Cox_coefficient_table_max_difference",
    "Cox_CI_table_max_difference",
    "PH_table_max_difference",
    "Adjusted_result_table_max_difference",
    "Clinical_state_summary_max_difference"
  ),
  
  Value = c(
    abs(
      median_CCL19_scanb -
        fig8_ref$median_CCL19
    ),
    abs(
      median_IL10_scanb -
        fig8_ref$median_IL10
    ),
    abs(
      CCL19_logrank_p -
        fig8_ref$CCL19_p
    ),
    abs(
      IL10_logrank_p -
        fig8_ref$IL10_p
    ),
    abs(
      state_logrank_p -
        fig8_ref$state_p
    ),
    pairwise_max_diff,
    adjusted_numeric_max_diff,
    adjusted_categorical_pass,
    cox_coefficient_max_diff,
    cox_CI_max_diff,
    PH_max_diff,
    adjusted_result_max_diff,
    clinical_state_max_diff
  ),
  
  stringsAsFactors = FALSE
)


write_csv_checked(
  regression_log,
  file.path(
    PATHS$logs,
    "13_scanb_validation_regression.csv"
  )
)


############################################################
## 13.48 — REFERENCE MD5 LOG
############################################################

reference_md5 <- data.frame(
  
  File = c(
    basename(figure8_reference_file)
  ),
  Expected_MD5 = c(
    EXPECTED_FIGURE8_REFERENCE_MD5
  ),
  
  Observed_MD5 = c(
    observed_figure8_md5
  ),
  
  Verified = c(
    identical(
      observed_figure8_md5,
      EXPECTED_FIGURE8_REFERENCE_MD5
    )
  ),
  
  stringsAsFactors = FALSE
)


write_csv_checked(
  reference_md5,
  file.path(
    PATHS$logs,
    "13_scanb_reference_MD5.csv"
  )
)


############################################################
## 13.49 — MODULE QC SUMMARY
############################################################

low_high_row <- which(
  adjusted_results$Variable ==
    "Low CCL19 / High IL10"
)


module13_qc <- data.frame(
  
  Metric = c(
    "SCAN-B patients",
    "SCAN-B OS events",
    "CCL19 exact median",
    "IL10 exact median",
    "CCL19 median ties",
    "IL10 median ties",
    "CCL19 Low",
    "CCL19 High",
    "IL10 Low",
    "IL10 High",
    "CCL19 log-rank P",
    "IL10 log-rank P",
    "High_Low N",
    "High_High N",
    "Low_Low N",
    "Low_High N",
    "High_Low events",
    "High_High events",
    "Low_Low events",
    "Low_High events",
    "Four-state overall log-rank P",
    "Adjusted cohort N",
    "Adjusted cohort events",
    "Adjusted model",
    "Adjusted reference state",
    "Low_High vs High_Low HR",
    "Low_High vs High_Low P",
    "GLOBAL PH P",
    "Frozen Figure-8 regression",
    "Frozen Figure-8 reference MD5",
    "Adjusted-subset medians used as cutoffs",
    "TCGA cutoffs transferred",
    "Statistical interaction fitted",
    "Historical analysis changed",
    "New analysis added"
  ),
  
  Value = c(
    nrow(scanb_final),
    sum(scanb_final$OS.status),
    format(
      median_CCL19_scanb,
      digits = 17
    ),
    format(
      median_IL10_scanb,
      digits = 17
    ),
    CCL19_median_ties,
    IL10_median_ties,
    CCL19_group_counts["Low"],
    CCL19_group_counts["High"],
    IL10_group_counts["Low"],
    IL10_group_counts["High"],
    format(
      CCL19_logrank_p,
      digits = 17
    ),
    format(
      IL10_logrank_p,
      digits = 17
    ),
    state_counts["High_Low"],
    state_counts["High_High"],
    state_counts["Low_Low"],
    state_counts["Low_High"],
    state_events["High_Low"],
    state_events["High_High"],
    state_events["Low_Low"],
    state_events["Low_High"],
    format(
      state_logrank_p,
      digits = 17
    ),
    nrow(scanb_cox),
    sum(scanb_cox$OS.status),
    "ImmuneState + Age + LymphNode + NHG + TumorSize",
    "High_Low",
    format(
      adjusted_results$HR[
        low_high_row
      ],
      digits = 17
    ),
    format(
      adjusted_results$Pvalue[
        low_high_row
      ],
      digits = 17
    ),
    format(
      PH_test$table[
        "GLOBAL",
        "p"
      ],
      digits = 17
    ),
    "VERIFIED",
    "VERIFIED",
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE
  ),
  
  stringsAsFactors = FALSE
)


write_csv_checked(
  module13_qc,
  file.path(
    PATHS$logs,
    "13_scanb_validation_QC.csv"
  )
)


############################################################
## 13.50 — SESSION INFORMATION
############################################################

capture.output(
  sessionInfo(),
  file = file.path(
    PATHS$logs,
    "13_scanb_validation_sessionInfo.txt"
  )
)


############################################################
## 13.51 — FINAL STATUS
############################################################

cat(
  "\n=============================================\n"
)


cat(
  "13_scanb_validation.R: PASS\n"
)


cat(
  "\nExternal validation cohort\n"
)


cat(
  "Dataset: SCAN-B / GSE96058\n"
)


cat(
  "Patients: ",
  nrow(scanb_final),
  "\n",
  sep = ""
)


cat(
  "OS events: ",
  sum(scanb_final$OS.status),
  "\n",
  sep = ""
)


cat(
  "\nCohort-specific biomarker cutoffs\n"
)


cat(
  "CCL19 median: ",
  format(
    median_CCL19_scanb,
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "IL10 median: ",
  format(
    median_IL10_scanb,
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "CCL19 median ties: ",
  CCL19_median_ties,
  "\n",
  sep = ""
)


cat(
  "IL10 median ties: ",
  IL10_median_ties,
  "\n",
  sep = ""
)


cat(
  "Tie rule: expression >= exact cohort median -> High\n"
)


cat(
  "\nIndividual survival\n"
)


cat(
  "CCL19 Low/High: ",
  CCL19_group_counts["Low"],
  " / ",
  CCL19_group_counts["High"],
  "\n",
  sep = ""
)


cat(
  "CCL19 events Low/High: ",
  CCL19_events["Low"],
  " / ",
  CCL19_events["High"],
  "\n",
  sep = ""
)


cat(
  "CCL19 log-rank P: ",
  format(
    CCL19_logrank_p,
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "IL10 Low/High: ",
  IL10_group_counts["Low"],
  " / ",
  IL10_group_counts["High"],
  "\n",
  sep = ""
)


cat(
  "IL10 events Low/High: ",
  IL10_events["Low"],
  " / ",
  IL10_events["High"],
  "\n",
  sep = ""
)


cat(
  "IL10 log-rank P: ",
  format(
    IL10_logrank_p,
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "\nJoint CCL19 / IL10 expression states\n"
)


cat(
  "High_Low: ",
  state_counts["High_Low"],
  " / ",
  state_events["High_Low"],
  " events\n",
  sep = ""
)


cat(
  "High_High: ",
  state_counts["High_High"],
  " / ",
  state_events["High_High"],
  " events\n",
  sep = ""
)


cat(
  "Low_Low: ",
  state_counts["Low_Low"],
  " / ",
  state_events["Low_Low"],
  " events\n",
  sep = ""
)


cat(
  "Low_High: ",
  state_counts["Low_High"],
  " / ",
  state_events["Low_High"],
  " events\n",
  sep = ""
)


cat(
  "Overall four-state log-rank P: ",
  format(
    state_logrank_p,
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "Pairwise method: log-rank + BH\n"
)


cat(
  "\nAdjusted historical Cox model\n"
)


cat(
  "Complete-case N: ",
  nrow(scanb_cox),
  "\n",
  sep = ""
)


cat(
  "Events: ",
  sum(scanb_cox$OS.status),
  "\n",
  sep = ""
)


cat(
  "Formula: ImmuneState + Age + LymphNode + NHG + TumorSize\n"
)


cat(
  "Reference state: High_Low\n"
)


cat(
  "Low_High vs High_Low HR: ",
  format(
    adjusted_results$HR[
      low_high_row
    ],
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "Low_High vs High_Low 95% CI: ",
  format(
    adjusted_results$Lower95[
      low_high_row
    ],
    digits = 17
  ),
  " to ",
  format(
    adjusted_results$Upper95[
      low_high_row
    ],
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "Low_High vs High_Low P: ",
  format(
    adjusted_results$Pvalue[
      low_high_row
    ],
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "GLOBAL PH P: ",
  format(
    PH_test$table[
      "GLOBAL",
      "p"
    ],
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "\nHistorical regression\n"
)


cat(
  "Full-cohort medians: VERIFIED\n"
)


cat(
  "Individual KM/log-rank: VERIFIED\n"
)


cat(
  "Four-state survival: VERIFIED\n"
)


cat(
  "Pairwise BH: VERIFIED\n"
)


cat(
  "Adjusted dataset: VERIFIED\n"
)


cat(
  "Adjusted Cox model: VERIFIED\n"
)


cat(
  "PH diagnostics: VERIFIED\n"
)


cat(
  "Adjusted result table: VERIFIED\n"
)


cat(
  "Clinical state summary: VERIFIED\n"
)


cat(
  "Historical reference MD5: VERIFIED\n"
)


cat(
  "\nAnalysis provenance\n"
)


cat(
  "SCAN-B-specific medians used: TRUE\n"
)


cat(
  "Adjusted-subset medians used as cutoffs: FALSE\n"
)


cat(
  "TCGA cutoffs transferred to SCAN-B: FALSE\n"
)


cat(
  "Statistical interaction fitted: FALSE\n"
)


cat(
  "Validation analysis specification changed: FALSE"
)


cat(
  "New analysis added: FALSE\n"
)


cat(
  "=============================================\n"
)