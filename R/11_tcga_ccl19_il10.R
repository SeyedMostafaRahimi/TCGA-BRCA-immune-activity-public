############################################################
## 11_tcga_ccl19_il10.R
##
## PURPOSE
##
## Production implementation of the validated focused
## CCL19 / IL10 TCGA-BRCA survival analysis.
##
## ANALYSES
##
##   1. Exact cohort-specific CCL19 median split
##   2. Exact cohort-specific IL10 median split
##   3. Four CCL19 / IL10 joint expression states
##   4. Overall four-state Kaplan-Meier / log-rank analysis
##   5. Pairwise four-state log-rank comparisons + BH
##   6. Final adjusted CCL19 + IL10 Cox model
##   7. Final adjusted four-state Cox model
##   8. Proportional-hazards diagnostics
##
##
## IMPORTANT ANALYTICAL GUARDRAILS
##
## - Medians are calculated on the complete 1,083-patient
##   TCGA survival cohort.
##
## - Median ties are assigned to High:
##
##       expression >= exact median -> High
##
## - The four groups represent JOINT EXPRESSION STATES.
##
## - They are NOT evidence of a statistical interaction
##   because no CCL19_group * IL10_group interaction term
##   is fitted.
##
## - Stage is handled by stratification in the final models.
##
## - Final main model:
##
##   Surv(OS.time, OS.status) ~
##       CCL19_group +
##       IL10_group +
##       Age +
##       strata(Stage)
##
## - Final joint-state model:
##
##   Surv(OS.time, OS.status) ~
##       ImmuneState +
##       Age +
##       strata(Stage)
##
## - Reference joint state:
##
##       High CCL19 / Low IL10
##
############################################################


############################################################
## 11.1 — PROJECT SETUP
############################################################

if (!exists("PATHS")) {
  source("R/00_setup.R")
}


require_package(
  "survival"
)


############################################################
## 11.2 — INPUT FILES
############################################################

survival_file <- file.path(
  PATHS$processed,
  "10_tcga_survival_cohort.rds"
)


expression_file <- file.path(
  PATHS$processed,
  "02_tcga_expression_final.rds"
)


cox_reference_file <- file.path(
  "external",
  "survival",
  "TCGA_BRCA_frozen_final_multivariable_Cox_reference.rds"
)


joint_reference_file <- file.path(
  "external",
  "survival",
  "TCGA_BRCA_frozen_CCL19_IL10_joint_state_reference.rds"
)


assert_file(
  survival_file
)


assert_file(
  expression_file
)


assert_file(
  cox_reference_file
)


assert_file(
  joint_reference_file
)


############################################################
## 11.3 — OUTPUT DIRECTORY
############################################################

module11_output_dir <- file.path(
  "results",
  "tables",
  "tcga_ccl19_il10"
)


dir.create(
  module11_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
## 11.4 — LOAD INPUTS
############################################################

tcga_surv <- readRDS(
  survival_file
)


expr_final <- readRDS(
  expression_file
)


cox_reference <- readRDS(
  cox_reference_file
)


joint_reference <- readRDS(
  joint_reference_file
)


############################################################
## 11.5 — INPUT QC
############################################################

assert_true(
  nrow(tcga_surv) == 1083L,
  "Module-10 survival cohort must contain 1,083 patients."
)


assert_true(
  sum(
    tcga_surv$OS.status == 1
  ) == 150L,
  "Module-10 survival cohort must contain 150 events."
)


assert_unique(
  tcga_surv$sample,
  "TCGA survival samples"
)


assert_true(
  nrow(expr_final) == 35009L,
  "Expression matrix must contain 35,009 genes."
)


assert_unique(
  rownames(expr_final),
  "TCGA expression genes"
)


assert_unique(
  colnames(expr_final),
  "TCGA expression samples"
)


assert_true(
  all(
    tcga_surv$sample %in%
      colnames(expr_final)
  ),
  "One or more TCGA survival samples are absent from expression matrix."
)


assert_true(
  all(
    c(
      "CCL19",
      "IL10"
    ) %in%
      rownames(expr_final)
  ),
  "CCL19 and/or IL10 are absent from expression matrix."
)


############################################################
## 11.6 — EXPRESSION EXTRACTION HELPER
############################################################

get_tcga_gene <- function(
    gene,
    samples
) {
  
  assert_true(
    gene %in%
      rownames(expr_final),
    paste0(
      "Gene absent from expression matrix: ",
      gene
    )
  )
  
  
  sample_index <- match(
    samples,
    colnames(expr_final)
  )
  
  
  assert_true(
    !anyNA(sample_index),
    "Unable to match all requested samples to expression matrix."
  )
  
  
  as.numeric(
    expr_final[
      gene,
      sample_index
    ]
  )
}


############################################################
## 11.7 — EXTRACT CCL19 / IL10 EXPRESSION
############################################################

CCL19_expression <- get_tcga_gene(
  "CCL19",
  tcga_surv$sample
)


IL10_expression <- get_tcga_gene(
  "IL10",
  tcga_surv$sample
)


assert_true(
  length(CCL19_expression) == 1083L,
  "CCL19 expression vector must contain 1,083 values."
)


assert_true(
  length(IL10_expression) == 1083L,
  "IL10 expression vector must contain 1,083 values."
)


assert_true(
  !anyNA(CCL19_expression),
  "Missing CCL19 expression detected."
)


assert_true(
  !anyNA(IL10_expression),
  "Missing IL10 expression detected."
)


############################################################
## 11.8 — EXACT COHORT MEDIANS
############################################################

CCL19_median <- median(
  CCL19_expression,
  na.rm = TRUE
)


IL10_median <- median(
  IL10_expression,
  na.rm = TRUE
)


CCL19_ties <- sum(
  CCL19_expression ==
    CCL19_median
)


IL10_ties <- sum(
  IL10_expression ==
    IL10_median
)


############################################################
## Validated exact checkpoints
############################################################

EXPECTED_CCL19_MEDIAN <-
  8.5924570372680797


EXPECTED_IL10_MEDIAN <-
  4.7548875021634682


assert_true(
  abs(
    CCL19_median -
      EXPECTED_CCL19_MEDIAN
  ) < 1e-12,
  "CCL19 cohort median differs from validated cutoff."
)


assert_true(
  abs(
    IL10_median -
      EXPECTED_IL10_MEDIAN
  ) < 1e-12,
  "IL10 cohort median differs from validated cutoff."
)


assert_true(
  CCL19_ties == 3L,
  "Expected 3 CCL19 observations equal to exact median."
)


assert_true(
  IL10_ties == 19L,
  "Expected 19 IL10 observations equal to exact median."
)


############################################################
## 11.9 — CREATE GENE GROUPS
##
## Frozen tie rule:
##
## expression >= exact median -> High
############################################################

tcga_state <- tcga_surv


tcga_state$CCL19_expr <-
  CCL19_expression


tcga_state$IL10_expr <-
  IL10_expression


tcga_state$CCL19_group <- factor(
  
  ifelse(
    tcga_state$CCL19_expr >=
      CCL19_median,
    "High",
    "Low"
  ),
  
  levels = c(
    "Low",
    "High"
  )
)


tcga_state$IL10_group <- factor(
  
  ifelse(
    tcga_state$IL10_expr >=
      IL10_median,
    "High",
    "Low"
  ),
  
  levels = c(
    "Low",
    "High"
  )
)


############################################################
## 11.10 — GENE-GROUP QC
############################################################

CCL19_counts <- table(
  tcga_state$CCL19_group
)


IL10_counts <- table(
  tcga_state$IL10_group
)


assert_true(
  identical(
    as.integer(
      CCL19_counts
    ),
    c(
      540L,
      543L
    )
  ),
  "CCL19 Low/High group counts differ from validated result."
)


assert_true(
  identical(
    as.integer(
      IL10_counts
    ),
    c(
      524L,
      559L
    )
  ),
  "IL10 Low/High group counts differ from validated result."
)


############################################################
## 11.11 — DEFINE FOUR JOINT EXPRESSION STATES
############################################################

tcga_state$ImmuneState <- paste(
  tcga_state$CCL19_group,
  tcga_state$IL10_group,
  sep = "_"
)


tcga_state$ImmuneState <- factor(
  
  tcga_state$ImmuneState,
  
  levels = c(
    "High_Low",
    "High_High",
    "Low_Low",
    "Low_High"
  )
)


############################################################
## 11.12 — FOUR-STATE COUNTS
############################################################

state_counts <- table(
  tcga_state$ImmuneState
)


state_events <- with(
  tcga_state,
  tapply(
    OS.status,
    ImmuneState,
    sum
  )
)


assert_true(
  identical(
    as.integer(state_counts),
    c(
      196L,
      347L,
      328L,
      212L
    )
  ),
  "Joint-state sample counts differ from validated result."
)


assert_true(
  identical(
    as.integer(state_events),
    c(
      21L,
      42L,
      41L,
      46L
    )
  ),
  "Joint-state event counts differ from validated result."
)


state_summary <- data.frame(
  
  ImmuneState = c(
    "High_Low",
    "High_High",
    "Low_Low",
    "Low_High"
  ),
  
  CCL19 = c(
    "High",
    "High",
    "Low",
    "Low"
  ),
  
  IL10 = c(
    "Low",
    "High",
    "Low",
    "High"
  ),
  
  N =
    as.integer(
      state_counts
    ),
  
  Events =
    as.integer(
      state_events
    ),
  
  stringsAsFactors =
    FALSE
)


############################################################
## 11.13 — FOUR-STATE KAPLAN-MEIER
############################################################

joint_survfit <- survival::survfit(
  
  survival::Surv(
    OS.time,
    OS.status
  ) ~ ImmuneState,
  
  data =
    tcga_state
)


############################################################
## 11.14 — OVERALL FOUR-STATE LOG-RANK
############################################################

joint_logrank <- survival::survdiff(
  
  survival::Surv(
    OS.time,
    OS.status
  ) ~ ImmuneState,
  
  data =
    tcga_state
)


joint_logrank_df <-
  length(
    joint_logrank$n
  ) - 1L


joint_logrank_p <- stats::pchisq(
  
  joint_logrank$chisq,
  
  df =
    joint_logrank_df,
  
  lower.tail =
    FALSE
)


############################################################
## 11.15 — LOG-RANK CHECKPOINTS
############################################################

EXPECTED_JOINT_CHISQ <-
  28.535057701090778


EXPECTED_JOINT_P <-
  0.0000028042089441644108


assert_true(
  abs(
    joint_logrank$chisq -
      EXPECTED_JOINT_CHISQ
  ) < 1e-10,
  "Joint-state log-rank chi-square differs from validated result."
)


assert_true(
  joint_logrank_df == 3L,
  "Joint-state overall log-rank must have df = 3."
)


assert_true(
  abs(
    joint_logrank_p -
      EXPECTED_JOINT_P
  ) < 1e-12,
  "Joint-state overall log-rank P-value differs from validated result."
)


joint_logrank_table <- data.frame(
  
  N =
    nrow(
      tcga_state
    ),
  
  Events =
    sum(
      tcga_state$OS.status
    ),
  
  Chi_square =
    as.numeric(
      joint_logrank$chisq
    ),
  
  df =
    joint_logrank_df,
  
  P_value =
    joint_logrank_p,
  
  stringsAsFactors =
    FALSE
)


############################################################
## 11.16 — PAIRWISE LOG-RANK COMPARISONS
##
## Six comparisons.
## BH correction across these six pairwise tests.
############################################################

state_levels <- levels(
  tcga_state$ImmuneState
)


state_pairs <- combn(
  state_levels,
  2,
  simplify = FALSE
)


state_label <- c(
  High_Low =
    "CCL19 High / IL10 Low",
  
  High_High =
    "CCL19 High / IL10 High",
  
  Low_Low =
    "CCL19 Low / IL10 Low",
  
  Low_High =
    "CCL19 Low / IL10 High"
)


pairwise_list <- vector(
  "list",
  length(
    state_pairs
  )
)


for (
  i in seq_along(
    state_pairs
  )
) {
  
  pair <- state_pairs[[i]]
  
  
  state_1 <- pair[1]
  state_2 <- pair[2]
  
  
  dat_pair <- tcga_state[
    tcga_state$ImmuneState %in%
      pair,
    ,
    drop = FALSE
  ]
  
  
  dat_pair$ImmuneState <-
    droplevels(
      dat_pair$ImmuneState
    )
  
  
  pair_test <- survival::survdiff(
    
    survival::Surv(
      OS.time,
      OS.status
    ) ~ ImmuneState,
    
    data =
      dat_pair
  )
  
  
  pair_p <- stats::pchisq(
    
    pair_test$chisq,
    
    df = 1,
    
    lower.tail =
      FALSE
  )
  
  
  pairwise_list[[i]] <- data.frame(
    
    State_1 =
      state_1,
    
    State_2 =
      state_2,
    
    Comparison =
      paste(
        state_label[
          state_2
        ],
        "vs",
        state_label[
          state_1
        ]
      ),
    
    N_state_1 =
      sum(
        tcga_state$ImmuneState ==
          state_1
      ),
    
    N_state_2 =
      sum(
        tcga_state$ImmuneState ==
          state_2
      ),
    
    Events_state_1 =
      sum(
        tcga_state$OS.status[
          tcga_state$ImmuneState ==
            state_1
        ]
      ),
    
    Events_state_2 =
      sum(
        tcga_state$OS.status[
          tcga_state$ImmuneState ==
            state_2
        ]
      ),
    
    Chi_square =
      as.numeric(
        pair_test$chisq
      ),
    
    P_unadjusted =
      pair_p,
    
    stringsAsFactors =
      FALSE
  )
}


pairwise_logrank <- do.call(
  rbind,
  pairwise_list
)


rownames(
  pairwise_logrank
) <- NULL


pairwise_logrank$P_BH <- stats::p.adjust(
  pairwise_logrank$P_unadjusted,
  method = "BH"
)


############################################################
## 11.17 — PAIRWISE CHECKPOINTS
############################################################

assert_true(
  nrow(
    pairwise_logrank
  ) == 6L,
  "Expected exactly six pairwise joint-state comparisons."
)


pairwise_BH_recheck <- stats::p.adjust(
  pairwise_logrank$P_unadjusted,
  method = "BH"
)


assert_true(
  max(
    abs(
      pairwise_BH_recheck -
        pairwise_logrank$P_BH
    )
  ) < 1e-15,
  "Pairwise joint-state P-values were not BH adjusted across all six tests."
)


############################################################
## Key validated pairwise comparisons
############################################################

pair_HL_LH <- pairwise_logrank[
  pairwise_logrank$State_1 ==
    "High_Low" &
    pairwise_logrank$State_2 ==
    "Low_High",
  ,
  drop = FALSE
]


assert_true(
  nrow(
    pair_HL_LH
  ) == 1L,
  "High_Low vs Low_High comparison missing."
)


assert_true(
  abs(
    pair_HL_LH$P_BH -
      0.00001540144
  ) < 1e-10,
  "High_Low vs Low_High pairwise BH P-value differs from validated result."
)


############################################################
## 11.18 — FINAL ADJUSTED COHORT
##
## Complete Age and Stage.
##
## Cutoffs remain those from the full 1,083-patient
## survival cohort and are NOT recalculated here.
############################################################

tcga_cox <- tcga_state[
  !is.na(
    tcga_state$Age
  ) &
    !is.na(
      tcga_state$Stage
    ),
  ,
  drop = FALSE
]


assert_true(
  nrow(
    tcga_cox
  ) == 1046L,
  "Final adjusted-model cohort must contain 1,046 patients."
)


assert_true(
  sum(
    tcga_cox$OS.status
  ) == 140L,
  "Final adjusted-model cohort must contain 140 events."
)


adjusted_stage_counts <- table(
  tcga_cox$Stage
)


assert_true(
  identical(
    as.integer(
      adjusted_stage_counts
    ),
    c(
      183L,
      602L,
      241L,
      20L
    )
  ),
  "Adjusted-model Stage counts differ from validated cohort."
)


############################################################
## 11.19 — FINAL CCL19 + IL10 COX MODEL
##
## Stage handled by stratification.
############################################################

fit6B_final <- survival::coxph(
  
  survival::Surv(
    OS.time,
    OS.status
  ) ~
    CCL19_group +
    IL10_group +
    Age +
    strata(Stage),
  
  data =
    tcga_cox,
  
  x =
    TRUE
)


############################################################
## 11.20 — MAIN MODEL PH DIAGNOSTICS
############################################################

ph6B_final <- survival::cox.zph(
  fit6B_final
)


############################################################
## 11.21 — FINAL FOUR-STATE ADJUSTED COX MODEL
############################################################

tcga_state_cox <- tcga_cox


tcga_state_cox$ImmuneState <- factor(
  
  tcga_state_cox$ImmuneState,
  
  levels = c(
    "High_Low",
    "High_High",
    "Low_Low",
    "Low_High"
  )
)


tcga_state_cox$ImmuneState <- stats::relevel(
  tcga_state_cox$ImmuneState,
  ref = "High_Low"
)


fit7B <- survival::coxph(
  
  survival::Surv(
    OS.time,
    OS.status
  ) ~
    ImmuneState +
    Age +
    strata(Stage),
  
  data =
    tcga_state_cox,
  
  x =
    TRUE
)


############################################################
## 11.22 — JOINT-STATE MODEL PH DIAGNOSTICS
############################################################

ph7B_final <- survival::cox.zph(
  fit7B
)


############################################################
## 11.23 — MODEL REGRESSION HELPER
############################################################

compare_cox_models <- function(
    production_fit,
    reference_fit
) {
  
  production_summary <- summary(
    production_fit
  )
  
  
  reference_summary <- summary(
    reference_fit
  )
  
  
  assert_true(
    identical(
      rownames(
        production_summary$coefficients
      ),
      rownames(
        reference_summary$coefficients
      )
    ),
    "Cox coefficient names differ from frozen reference."
  )
  
  
  coefficient_difference <- max(
    abs(
      production_summary$coefficients -
        reference_summary$coefficients
    ),
    na.rm = TRUE
  )
  
  
  confidence_difference <- max(
    abs(
      production_summary$conf.int -
        reference_summary$conf.int
    ),
    na.rm = TRUE
  )
  
  
  c(
    coefficient_max_difference =
      coefficient_difference,
    
    CI_max_difference =
      confidence_difference
  )
}


############################################################
## 11.24 — REGRESS FINAL MODELS AGAINST FROZEN REFERENCES
############################################################

fit6_regression <- compare_cox_models(
  fit6B_final,
  cox_reference$fit6B_final
)


fit7_regression <- compare_cox_models(
  fit7B,
  cox_reference$fit7B
)


assert_true(
  all(
    fit6_regression <
      1e-10
  ),
  "Final CCL19 + IL10 Cox model differs from frozen reference."
)


assert_true(
  all(
    fit7_regression <
      1e-10
  ),
  "Final joint-state Cox model differs from frozen reference."
)


############################################################
## 11.25 — PH REGRESSION HELPER
############################################################

compare_ph_models <- function(
    production_ph,
    reference_ph
) {
  
  production_table <- as.data.frame(
    production_ph$table
  )
  
  
  reference_table <- as.data.frame(
    reference_ph$table
  )
  
  
  assert_true(
    identical(
      rownames(
        production_table
      ),
      rownames(
        reference_table
      )
    ),
    "PH diagnostic terms differ from frozen reference."
  )
  
  
  max(
    abs(
      as.matrix(
        production_table
      ) -
        as.matrix(
          reference_table
        )
    ),
    na.rm = TRUE
  )
}


############################################################
## 11.26 — PH DIAGNOSTIC REGRESSION
############################################################

ph6_max_difference <- compare_ph_models(
  ph6B_final,
  cox_reference$ph6B_final
)


ph7_max_difference <- compare_ph_models(
  ph7B_final,
  cox_reference$ph7B_final
)


assert_true(
  ph6_max_difference <
    1e-10,
  "Main-model PH diagnostics differ from frozen reference."
)


assert_true(
  ph7_max_difference <
    1e-10,
  "Joint-state PH diagnostics differ from frozen reference."
)


############################################################
## 11.27 — MAIN MODEL RESULTS TABLE
############################################################

s6 <- summary(
  fit6B_final
)


main_model_table <- data.frame(
  
  Model =
    "CCL19 + IL10 model",
  
  Predictor = c(
    "CCL19",
    "IL10",
    "Age"
  ),
  
  Comparison = c(
    "High vs Low",
    "High vs Low",
    "Per 1-year increase"
  ),
  
  HR = c(
    s6$coefficients[
      "CCL19_groupHigh",
      "exp(coef)"
    ],
    s6$coefficients[
      "IL10_groupHigh",
      "exp(coef)"
    ],
    s6$coefficients[
      "Age",
      "exp(coef)"
    ]
  ),
  
  Lower95 = c(
    s6$conf.int[
      "CCL19_groupHigh",
      "lower .95"
    ],
    s6$conf.int[
      "IL10_groupHigh",
      "lower .95"
    ],
    s6$conf.int[
      "Age",
      "lower .95"
    ]
  ),
  
  Upper95 = c(
    s6$conf.int[
      "CCL19_groupHigh",
      "upper .95"
    ],
    s6$conf.int[
      "IL10_groupHigh",
      "upper .95"
    ],
    s6$conf.int[
      "Age",
      "upper .95"
    ]
  ),
  
  Pvalue = c(
    s6$coefficients[
      "CCL19_groupHigh",
      "Pr(>|z|)"
    ],
    s6$coefficients[
      "IL10_groupHigh",
      "Pr(>|z|)"
    ],
    s6$coefficients[
      "Age",
      "Pr(>|z|)"
    ]
  ),
  
  stringsAsFactors =
    FALSE
)


############################################################
## 11.28 — JOINT-STATE MODEL RESULTS TABLE
############################################################

s7 <- summary(
  fit7B
)


joint_model_table <- data.frame(
  
  Model =
    "CCL19-IL10 joint-state model",
  
  Predictor = c(
    "High CCL19 / High IL10",
    "Low CCL19 / Low IL10",
    "Low CCL19 / High IL10",
    "Age"
  ),
  
  Comparison = c(
    "vs High CCL19 / Low IL10",
    "vs High CCL19 / Low IL10",
    "vs High CCL19 / Low IL10",
    "Per 1-year increase"
  ),
  
  HR = c(
    s7$coefficients[
      "ImmuneStateHigh_High",
      "exp(coef)"
    ],
    s7$coefficients[
      "ImmuneStateLow_Low",
      "exp(coef)"
    ],
    s7$coefficients[
      "ImmuneStateLow_High",
      "exp(coef)"
    ],
    s7$coefficients[
      "Age",
      "exp(coef)"
    ]
  ),
  
  Lower95 = c(
    s7$conf.int[
      "ImmuneStateHigh_High",
      "lower .95"
    ],
    s7$conf.int[
      "ImmuneStateLow_Low",
      "lower .95"
    ],
    s7$conf.int[
      "ImmuneStateLow_High",
      "lower .95"
    ],
    s7$conf.int[
      "Age",
      "lower .95"
    ]
  ),
  
  Upper95 = c(
    s7$conf.int[
      "ImmuneStateHigh_High",
      "upper .95"
    ],
    s7$conf.int[
      "ImmuneStateLow_Low",
      "upper .95"
    ],
    s7$conf.int[
      "ImmuneStateLow_High",
      "upper .95"
    ],
    s7$conf.int[
      "Age",
      "upper .95"
    ]
  ),
  
  Pvalue = c(
    s7$coefficients[
      "ImmuneStateHigh_High",
      "Pr(>|z|)"
    ],
    s7$coefficients[
      "ImmuneStateLow_Low",
      "Pr(>|z|)"
    ],
    s7$coefficients[
      "ImmuneStateLow_High",
      "Pr(>|z|)"
    ],
    s7$coefficients[
      "Age",
      "Pr(>|z|)"
    ]
  ),
  
  stringsAsFactors =
    FALSE
)


############################################################
## 11.29 — COMBINED FINAL MODEL TABLE
############################################################

multivariable_results <- rbind(
  main_model_table,
  joint_model_table
)


############################################################
## 11.30 — NUMERICAL RESULT CHECKPOINTS
############################################################

assert_true(
  abs(
    main_model_table$HR[
      main_model_table$Predictor ==
        "CCL19"
    ] -
      0.59229769392161646
  ) < 1e-10,
  "CCL19 adjusted HR differs from validated result."
)


assert_true(
  abs(
    main_model_table$HR[
      main_model_table$Predictor ==
        "IL10"
    ] -
      1.65290411040879492
  ) < 1e-10,
  "IL10 adjusted HR differs from validated result."
)


assert_true(
  abs(
    joint_model_table$HR[
      joint_model_table$Predictor ==
        "Low CCL19 / High IL10"
    ] -
      2.4444647241026751
  ) < 1e-10,
  "Low CCL19 / High IL10 adjusted HR differs from validated result."
)


############################################################
## 11.31 — PH DIAGNOSTICS TABLE
############################################################

extract_ph_table <- function(
    ph_object,
    model_name
) {
  
  x <- as.data.frame(
    ph_object$table
  )
  
  
  x$Term <- rownames(
    x
  )
  
  
  rownames(
    x
  ) <- NULL
  
  
  data.frame(
    
    Model =
      model_name,
    
    Term =
      x$Term,
    
    chisq =
      x$chisq,
    
    df =
      x$df,
    
    p =
      x$p,
    
    stringsAsFactors =
      FALSE
  )
}


ph_main_table <- extract_ph_table(
  ph6B_final,
  "CCL19 + IL10 model"
)


ph_joint_table <- extract_ph_table(
  ph7B_final,
  "CCL19-IL10 joint-state model"
)


ph_diagnostics <- rbind(
  ph_main_table,
  ph_joint_table
)


############################################################
## 11.32 — GLOBAL PH CHECKPOINTS
############################################################

ph6_global_p <- ph6B_final$table[
  "GLOBAL",
  "p"
]


ph7_global_p <- ph7B_final$table[
  "GLOBAL",
  "p"
]


assert_true(
  abs(
    ph6_global_p -
      0.34935589443222409
  ) < 1e-10,
  "Main-model global PH P-value differs from validated result."
)


assert_true(
  abs(
    ph7_global_p -
      0.33904654159774278
  ) < 1e-10,
  "Joint-state global PH P-value differs from validated result."
)


############################################################
## 11.33 — CUTOFF / TIE AUDIT TABLE
############################################################

cutoff_audit <- data.frame(
  
  Gene = c(
    "CCL19",
    "IL10"
  ),
  
  Exact_median = c(
    CCL19_median,
    IL10_median
  ),
  
  Exact_median_full_precision = c(
    format(
      CCL19_median,
      digits = 17
    ),
    format(
      IL10_median,
      digits = 17
    )
  ),
  
  N_equal_to_exact_median = c(
    CCL19_ties,
    IL10_ties
  ),
  
  Rule =
    "High if expression >= exact cohort median",
  
  stringsAsFactors =
    FALSE
)


############################################################
## 11.34 — MODEL COHORT SUMMARY
############################################################

model_cohort_summary <- data.frame(
  
  Metric = c(
    "Full_survival_cohort",
    "Full_survival_events",
    "Adjusted_model_cohort",
    "Adjusted_model_events",
    "Adjusted_Stage_I",
    "Adjusted_Stage_II",
    "Adjusted_Stage_III",
    "Adjusted_Stage_IV",
    "CCL19_Low",
    "CCL19_High",
    "IL10_Low",
    "IL10_High",
    "Joint_High_Low",
    "Joint_High_High",
    "Joint_Low_Low",
    "Joint_Low_High"
  ),
  
  Value = c(
    nrow(tcga_state),
    sum(tcga_state$OS.status),
    nrow(tcga_cox),
    sum(tcga_cox$OS.status),
    adjusted_stage_counts[
      "Stage I"
    ],
    adjusted_stage_counts[
      "Stage II"
    ],
    adjusted_stage_counts[
      "Stage III"
    ],
    adjusted_stage_counts[
      "Stage IV"
    ],
    CCL19_counts[
      "Low"
    ],
    CCL19_counts[
      "High"
    ],
    IL10_counts[
      "Low"
    ],
    IL10_counts[
      "High"
    ],
    state_counts[
      "High_Low"
    ],
    state_counts[
      "High_High"
    ],
    state_counts[
      "Low_Low"
    ],
    state_counts[
      "Low_High"
    ]
  ),
  
  stringsAsFactors =
    FALSE
)


############################################################
## 11.35 — METHOD PROVENANCE
############################################################

method_provenance <- data.frame(
  
  Item = c(
    "CCL19_cutoff_source",
    "IL10_cutoff_source",
    "Median_tie_rule",
    "Joint_state_definition",
    "Overall_joint_state_survival_test",
    "Pairwise_joint_state_test",
    "Pairwise_multiple_testing",
    "Adjusted_model_cohort",
    "Main_Cox_formula",
    "Joint_state_Cox_formula",
    "Stage_handling",
    "Main_gene_reference_level",
    "Joint_state_reference",
    "PH_assumption_test",
    "Statistical_interaction_term",
    "Interpretation"
  ),
  
  Value = c(
    "Exact median in complete 1083-patient TCGA survival cohort",
    "Exact median in complete 1083-patient TCGA survival cohort",
    "High if expression >= exact median",
    "Four combinations of CCL19 High/Low and IL10 High/Low",
    "Kaplan-Meier with overall log-rank test",
    "Pairwise log-rank tests",
    "BH across six pairwise comparisons",
    "Complete Age and pathological Stage",
    "Surv(OS.time, OS.status) ~ CCL19_group + IL10_group + Age + strata(Stage)",
    "Surv(OS.time, OS.status) ~ ImmuneState + Age + strata(Stage)",
    "Stratified baseline hazards by pathological Stage",
    "Low",
    "High CCL19 / Low IL10",
    "survival::cox.zph",
    "Not fitted",
    "Joint expression-state model; associative and non-causal"
  ),
  
  stringsAsFactors =
    FALSE
)


############################################################
## 11.36 — SAVE PROCESSED OBJECTS
############################################################

save_rds_checked(
  tcga_state,
  file.path(
    PATHS$processed,
    "11_tcga_CCL19_IL10_joint_states.rds"
  )
)


save_rds_checked(
  tcga_cox,
  file.path(
    PATHS$processed,
    "11_tcga_CCL19_IL10_adjusted_cohort.rds"
  )
)


save_rds_checked(
  joint_survfit,
  file.path(
    PATHS$processed,
    "11_tcga_joint_state_survfit.rds"
  )
)


save_rds_checked(
  joint_logrank,
  file.path(
    PATHS$processed,
    "11_tcga_joint_state_logrank.rds"
  )
)


save_rds_checked(
  pairwise_logrank,
  file.path(
    PATHS$processed,
    "11_tcga_joint_state_pairwise_logrank.rds"
  )
)


save_rds_checked(
  fit6B_final,
  file.path(
    PATHS$processed,
    "11_tcga_CCL19_IL10_final_cox.rds"
  )
)


save_rds_checked(
  fit7B,
  file.path(
    PATHS$processed,
    "11_tcga_joint_state_final_cox.rds"
  )
)


save_rds_checked(
  ph6B_final,
  file.path(
    PATHS$processed,
    "11_tcga_CCL19_IL10_final_cox_PH.rds"
  )
)


save_rds_checked(
  ph7B_final,
  file.path(
    PATHS$processed,
    "11_tcga_joint_state_final_cox_PH.rds"
  )
)


############################################################
## 11.37 — EXPORT TABLES
############################################################

write_csv_checked(
  cutoff_audit,
  file.path(
    module11_output_dir,
    "11_TCGA_CCL19_IL10_cutoff_tie_audit.csv"
  )
)


write_csv_checked(
  state_summary,
  file.path(
    module11_output_dir,
    "11_TCGA_joint_state_survival_counts.csv"
  )
)


write_csv_checked(
  joint_logrank_table,
  file.path(
    module11_output_dir,
    "11_TCGA_joint_state_overall_logrank.csv"
  )
)


write_csv_checked(
  pairwise_logrank,
  file.path(
    module11_output_dir,
    "11_TCGA_joint_state_pairwise_logrank_BH.csv"
  )
)


write_csv_checked(
  multivariable_results,
  file.path(
    module11_output_dir,
    "11_TCGA_final_multivariable_Cox.csv"
  )
)


write_csv_checked(
  ph_diagnostics,
  file.path(
    module11_output_dir,
    "11_TCGA_final_Cox_PH_diagnostics.csv"
  )
)


write_csv_checked(
  model_cohort_summary,
  file.path(
    module11_output_dir,
    "11_TCGA_CCL19_IL10_cohort_summary.csv"
  )
)


write_csv_checked(
  method_provenance,
  file.path(
    module11_output_dir,
    "11_TCGA_CCL19_IL10_method_provenance.csv"
  )
)


############################################################
## 11.38 — REFERENCE MD5
############################################################

reference_hashes <- tools::md5sum(
  c(
    cox_reference_file,
    joint_reference_file
  )
)


reference_MD5 <- data.frame(
  
  File =
    names(
      reference_hashes
    ),
  
  MD5 =
    unname(
      reference_hashes
    ),
  
  stringsAsFactors =
    FALSE
)


write_csv_checked(
  reference_MD5,
  file.path(
    PATHS$logs,
    "11_TCGA_CCL19_IL10_reference_MD5.csv"
  )
)


############################################################
## 11.39 — REGRESSION LOG
############################################################

model_regression_log <- data.frame(
  
  Metric = c(
    "fit6B_coefficient_max_difference",
    "fit6B_CI_max_difference",
    "fit7B_coefficient_max_difference",
    "fit7B_CI_max_difference",
    "fit6B_PH_max_difference",
    "fit7B_PH_max_difference"
  ),
  
  Value = c(
    fit6_regression[
      "coefficient_max_difference"
    ],
    fit6_regression[
      "CI_max_difference"
    ],
    fit7_regression[
      "coefficient_max_difference"
    ],
    fit7_regression[
      "CI_max_difference"
    ],
    ph6_max_difference,
    ph7_max_difference
  ),
  
  stringsAsFactors =
    FALSE
)


write_csv_checked(
  model_regression_log,
  file.path(
    PATHS$logs,
    "11_TCGA_CCL19_IL10_model_regression.csv"
  )
)


############################################################
## 11.40 — MODULE QC
############################################################

module11_qc <- data.frame(
  
  Metric = c(
    "Full survival cohort",
    "Full events",
    "CCL19 exact median",
    "CCL19 ties",
    "IL10 exact median",
    "IL10 ties",
    "Median tie rule",
    "High_Low N",
    "High_Low events",
    "High_High N",
    "High_High events",
    "Low_Low N",
    "Low_Low events",
    "Low_High N",
    "Low_High events",
    "Overall joint-state chi-square",
    "Overall joint-state df",
    "Overall joint-state P",
    "Pairwise comparisons",
    "Pairwise BH family",
    "Adjusted model N",
    "Adjusted model events",
    "Stage treatment",
    "Main model",
    "Joint-state model",
    "Joint-state reference",
    "Main-model global PH P",
    "Joint-state global PH P",
    "Statistical interaction term fitted",
    "SCAN-B validation performed"
  ),
  
  Value = c(
    nrow(tcga_state),
    sum(tcga_state$OS.status),
    format(
      CCL19_median,
      digits = 17
    ),
    CCL19_ties,
    format(
      IL10_median,
      digits = 17
    ),
    IL10_ties,
    "expression >= exact median -> High",
    state_counts[
      "High_Low"
    ],
    state_events[
      "High_Low"
    ],
    state_counts[
      "High_High"
    ],
    state_events[
      "High_High"
    ],
    state_counts[
      "Low_Low"
    ],
    state_events[
      "Low_Low"
    ],
    state_counts[
      "Low_High"
    ],
    state_events[
      "Low_High"
    ],
    format(
      joint_logrank$chisq,
      digits = 17
    ),
    joint_logrank_df,
    format(
      joint_logrank_p,
      digits = 17
    ),
    nrow(
      pairwise_logrank
    ),
    "BH across six tests",
    nrow(
      tcga_cox
    ),
    sum(
      tcga_cox$OS.status
    ),
    "strata(Stage)",
    "CCL19_group + IL10_group + Age + strata(Stage)",
    "ImmuneState + Age + strata(Stage)",
    "High_Low",
    format(
      ph6_global_p,
      digits = 17
    ),
    format(
      ph7_global_p,
      digits = 17
    ),
    FALSE,
    FALSE
  ),
  
  stringsAsFactors =
    FALSE
)


write_csv_checked(
  module11_qc,
  file.path(
    PATHS$logs,
    "11_tcga_ccl19_il10_QC.csv"
  )
)


############################################################
## 11.41 — SESSION INFORMATION
############################################################

capture.output(
  sessionInfo(),
  file = file.path(
    PATHS$logs,
    "11_tcga_ccl19_il10_sessionInfo.txt"
  )
)


############################################################
## 11.42 — FINAL STATUS
############################################################

cat(
  "\n=============================================\n"
)


cat(
  "11_tcga_ccl19_il10.R: PASS\n"
)


cat(
  "\nExact TCGA cutoffs\n"
)


cat(
  "CCL19 median: ",
  format(
    CCL19_median,
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "CCL19 median ties: ",
  CCL19_ties,
  "\n",
  sep = ""
)


cat(
  "IL10 median: ",
  format(
    IL10_median,
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "IL10 median ties: ",
  IL10_ties,
  "\n",
  sep = ""
)


cat(
  "Tie rule: expression >= exact median -> High\n"
)


cat(
  "\nJoint expression states\n"
)


cat(
  "High_Low: ",
  state_counts[
    "High_Low"
  ],
  " / ",
  state_events[
    "High_Low"
  ],
  " events\n",
  sep = ""
)


cat(
  "High_High: ",
  state_counts[
    "High_High"
  ],
  " / ",
  state_events[
    "High_High"
  ],
  " events\n",
  sep = ""
)


cat(
  "Low_Low: ",
  state_counts[
    "Low_Low"
  ],
  " / ",
  state_events[
    "Low_Low"
  ],
  " events\n",
  sep = ""
)


cat(
  "Low_High: ",
  state_counts[
    "Low_High"
  ],
  " / ",
  state_events[
    "Low_High"
  ],
  " events\n",
  sep = ""
)


cat(
  "Overall log-rank chi-square: ",
  format(
    joint_logrank$chisq,
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "Overall log-rank df: ",
  joint_logrank_df,
  "\n",
  sep = ""
)


cat(
  "Overall log-rank P: ",
  format(
    joint_logrank_p,
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "Pairwise log-rank comparisons: 6\n"
)


cat(
  "Pairwise correction: BH across six tests\n"
)


cat(
  "\nFinal adjusted models\n"
)


cat(
  "Adjusted cohort: ",
  nrow(
    tcga_cox
  ),
  "\n",
  sep = ""
)


cat(
  "Adjusted events: ",
  sum(
    tcga_cox$OS.status
  ),
  "\n",
  sep = ""
)


cat(
  "Stage treatment: STRATIFIED\n"
)


cat(
  "Main model: CCL19_group + IL10_group + Age + strata(Stage)\n"
)


cat(
  "CCL19 HR: ",
  format(
    main_model_table$HR[
      main_model_table$Predictor ==
        "CCL19"
    ],
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "IL10 HR: ",
  format(
    main_model_table$HR[
      main_model_table$Predictor ==
        "IL10"
    ],
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "Main-model GLOBAL PH P: ",
  format(
    ph6_global_p,
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "\nJoint-state adjusted model\n"
)


cat(
  "Reference: High CCL19 / Low IL10\n"
)


cat(
  "Low CCL19 / High IL10 HR: ",
  format(
    joint_model_table$HR[
      joint_model_table$Predictor ==
        "Low CCL19 / High IL10"
    ],
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "Low CCL19 / High IL10 P: ",
  format(
    joint_model_table$Pvalue[
      joint_model_table$Predictor ==
        "Low CCL19 / High IL10"
    ],
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "Joint-state GLOBAL PH P: ",
  format(
    ph7_global_p,
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "\nInterpretation\n"
)


cat(
  "Joint-state model is an expression-state model, not a statistical interaction model.\n"
)


cat(
  "Associations are prognostic associations; causal regulation is not inferred.\n"
)


cat(
  "Historical final models regression: VERIFIED\n"
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