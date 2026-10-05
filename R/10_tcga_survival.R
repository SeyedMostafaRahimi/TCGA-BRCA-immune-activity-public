############################################################
## 10_tcga_survival.R
##
## PURPOSE
##
## Production/refactoring module for the validated TCGA-BRCA
## overall-survival analyses.
##
## THIS MODULE PRESERVES THE HISTORICAL ANALYSIS:
##
##   1. Construct overall survival
##   2. Define the valid survival cohort
##   3. Kaplan-Meier / log-rank analysis by immune phenotype
##   4. Univariate Cox screening of the 20 MCC hub genes
##   5. BH correction across all 20 hub-gene tests
##
##
## THIS MODULE DOES NOT:
##
##   - redefine the survival endpoint
##   - change the survival cohort
##   - change the median-split rule
##   - change the Cox model
##   - change the multiple-testing family
##   - perform adjusted CCL19 / IL10 models
##   - perform the CCL19 × IL10 four-state analysis
##
## Those focused CCL19 / IL10 analyses are handled later.
##
##
## HISTORICAL SURVIVAL DEFINITION
##
## Dead:
##   OS.time = days_to_death.demographic
##
## Alive:
##   OS.time = days_to_last_follow_up.diagnoses
##
## OS.status:
##   Dead  = 1
##   Alive = 0
##
## Include only:
##
##   !is.na(OS.time)
##   OS.time > 0
##   vital_status.demographic %in% c("Alive", "Dead")
##
##
## EXPECTED FINAL COHORT
##
##   N       = 1083
##   Events  = 150
##   Censored= 933
##
## Immune phenotype:
##
##   Immune_Low  = 361
##   Immune_Mid  = 361
##   Immune_High = 361
##
##
## HUB SURVIVAL SCREEN
##
## For each of the 20 MCC hubs:
##
##   - expression explicitly matched to survival sample IDs
##   - exact cohort-specific median calculated
##   - High if expression >= median
##   - Low otherwise
##   - Low used as reference level
##
## Model:
##
##   Cox PH:
##   Surv(OS.time, OS.status) ~ GeneGroup
##
## Multiple testing:
##
##   Benjamini-Hochberg across ALL 20 hub tests
##
############################################################


############################################################
## 10.1 — REQUIRE PROJECT SETUP
############################################################

if (!exists("PATHS")) {
  source("R/00_setup.R")
}


############################################################
## 10.2 — REQUIRED PACKAGE
############################################################

require_package(
  "survival"
)


############################################################
## 10.3 — INPUT FILES
############################################################

metadata_file <- file.path(
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


survival_reference_file <- file.path(
  "external",
  "survival",
  "TCGA_BRCA_historical_survival_reference.rds"
)


assert_file(
  metadata_file
)


assert_file(
  expression_file
)


assert_file(
  hub_file
)


assert_file(
  survival_reference_file
)


############################################################
## 10.4 — OUTPUT DIRECTORY
############################################################

survival_output_dir <- file.path(
  "results",
  "tables",
  "tcga_survival"
)


dir.create(
  survival_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
## 10.5 — LOAD INPUTS
############################################################

meta05 <- readRDS(
  metadata_file
)


expr_final <- readRDS(
  expression_file
)


hub09 <- readRDS(
  hub_file
)


tcga_surv_reference <- readRDS(
  survival_reference_file
)


############################################################
## 10.6 — UPSTREAM INPUT QC
############################################################

assert_true(
  nrow(meta05) == 1106L,
  "Module-05 metadata must contain 1,106 primary tumors."
)


assert_true(
  nrow(expr_final) == 35009L,
  "Expression matrix must contain 35,009 genes."
)


assert_true(
  ncol(expr_final) == 1106L,
  "Expression matrix must contain 1,106 tumor samples."
)


assert_true(
  nrow(hub09) == 20L,
  "Module-09 hub table must contain 20 genes."
)


assert_unique(
  meta05$sample,
  "Module-05 samples"
)


assert_unique(
  colnames(expr_final),
  "Expression samples"
)


assert_unique(
  rownames(expr_final),
  "Expression genes"
)


assert_true(
  identical(
    colnames(expr_final),
    meta05$sample
  ),
  paste(
    "Module-05 metadata and Module-02 expression",
    "sample order must be identical."
  )
)


############################################################
## 10.7 — HUB INPUT QC
############################################################

expected_hubs <- c(
  "CCL2",
  "CXCL10",
  "CXCL8",
  "CCL5",
  "CCL4",
  "CCL22",
  "CCR5",
  "CCR2",
  "CCR7",
  "CXCR3",
  "CXCL11",
  "CCL19",
  "IL6",
  "CCL20",
  "CXCL9",
  "CXCL1",
  "CX3CL1",
  "CXCL2",
  "IL10",
  "IFNG"
)


hub_genes <- as.character(
  hub09$Gene
)


assert_true(
  identical(
    hub_genes,
    expected_hubs
  ),
  "Module-09 top-20 hub order differs from frozen MCC ranking."
)


assert_true(
  all(
    hub_genes %in%
      rownames(expr_final)
  ),
  "One or more hub genes are absent from expression matrix."
)


############################################################
## 10.8 — CONSTRUCT OVERALL SURVIVAL
##
## Exact historical definition.
############################################################

tcga_surv <- meta05


tcga_surv$OS.time <- ifelse(
  
  tcga_surv$vital_status.demographic ==
    "Dead",
  
  tcga_surv$days_to_death.demographic,
  
  tcga_surv$days_to_last_follow_up.diagnoses
)


tcga_surv$OS.status <- ifelse(
  
  tcga_surv$vital_status.demographic ==
    "Dead",
  
  1L,
  
  0L
)


############################################################
## 10.9 — CANONICAL AGE
##
## Preserved for downstream focused Cox models.
############################################################

tcga_surv$Age <-
  tcga_surv$
  age_at_earliest_diagnosis_in_years.diagnoses.xena_derived


############################################################
## 10.10 — CANONICAL PATHOLOGICAL STAGE
##
## Stage substages are collapsed to Stage I-IV exactly as
## used in the validated downstream Cox workflow.
############################################################

stage_raw <- as.character(
  tcga_surv$ajcc_pathologic_stage.diagnoses
)


stage_clean <- rep(
  NA_character_,
  length(stage_raw)
)


stage_clean[
  grepl(
    "Stage IV",
    stage_raw
  )
] <- "Stage IV"


stage_clean[
  is.na(stage_clean) &
    grepl(
      "Stage III",
      stage_raw
    )
] <- "Stage III"


stage_clean[
  is.na(stage_clean) &
    grepl(
      "Stage II",
      stage_raw
    )
] <- "Stage II"


stage_clean[
  is.na(stage_clean) &
    grepl(
      "Stage I",
      stage_raw
    )
] <- "Stage I"


tcga_surv$Stage <- stage_clean


############################################################
## 10.11 — APPLY EXACT HISTORICAL SURVIVAL FILTER
############################################################

survival_keep <- (
  !is.na(
    tcga_surv$OS.time
  ) &
    tcga_surv$OS.time > 0 &
    tcga_surv$vital_status.demographic %in%
    c(
      "Alive",
      "Dead"
    )
)


tcga_surv <- tcga_surv[
  survival_keep,
  ,
  drop = FALSE
]


############################################################
## 10.12 — FACTOR ORDER
############################################################

tcga_surv$ImmunePhenotype <- factor(
  
  tcga_surv$ImmunePhenotype,
  
  levels = c(
    "Immune_Low",
    "Immune_Mid",
    "Immune_High"
  )
)


tcga_surv$Stage <- factor(
  
  tcga_surv$Stage,
  
  levels = c(
    "Stage I",
    "Stage II",
    "Stage III",
    "Stage IV"
  )
)


############################################################
## 10.13 — SURVIVAL COHORT QC
############################################################

assert_true(
  nrow(tcga_surv) == 1083L,
  paste0(
    "Expected 1,083 survival cases; observed ",
    nrow(tcga_surv),
    "."
  )
)


assert_true(
  sum(
    tcga_surv$OS.status == 1
  ) == 150L,
  "Expected 150 deaths/events."
)


assert_true(
  sum(
    tcga_surv$OS.status == 0
  ) == 933L,
  "Expected 933 censored observations."
)


assert_true(
  all(
    tcga_surv$OS.time > 0
  ),
  "Non-positive survival times remain in TCGA cohort."
)


assert_true(
  sum(
    is.na(
      tcga_surv$OS.time
    )
  ) == 0L,
  "Missing OS.time values remain after survival filtering."
)


assert_true(
  sum(
    is.na(
      tcga_surv$OS.status
    )
  ) == 0L,
  "Missing OS.status values remain after survival filtering."
)


phenotype_counts_internal <- table(
  tcga_surv$ImmunePhenotype
)


assert_true(
  identical(
    as.integer(
      phenotype_counts_internal
    ),
    c(
      361L,
      361L,
      361L
    )
  ),
  "Expected 361 survival cases in each immune phenotype."
)


############################################################
## 10.14 — HISTORICAL SURVIVAL-COHORT REGRESSION
##
## Historical reference row order differs from the current
## production metadata order.
##
## Reference is therefore aligned ONLY for regression.
##
## The production survival cohort itself is not reordered.
############################################################

assert_true(
  nrow(
    tcga_surv_reference
  ) == 1083L,
  "Historical survival reference must contain 1,083 samples."
)


assert_unique(
  tcga_surv_reference$sample,
  "Historical survival-reference samples"
)


survival_sample_set_pass <- setequal(
  tcga_surv$sample,
  tcga_surv_reference$sample
)


assert_true(
  survival_sample_set_pass,
  "Production and historical survival sample sets differ."
)


reference_index <- match(
  tcga_surv$sample,
  tcga_surv_reference$sample
)


assert_true(
  !anyNA(
    reference_index
  ),
  "Unable to align historical survival reference by sample ID."
)


tcga_surv_reference_aligned <-
  tcga_surv_reference[
    reference_index,
    ,
    drop = FALSE
  ]


sample_alignment_pass <- identical(
  
  as.character(
    tcga_surv$sample
  ),
  
  as.character(
    tcga_surv_reference_aligned$sample
  )
)


os_time_max_difference <- max(
  abs(
    as.numeric(
      tcga_surv$OS.time
    ) -
      as.numeric(
        tcga_surv_reference_aligned$OS.time
      )
  )
)


os_status_pass <- identical(
  
  as.integer(
    tcga_surv$OS.status
  ),
  
  as.integer(
    tcga_surv_reference_aligned$OS.status
  )
)


phenotype_reference_pass <- identical(
  
  as.character(
    tcga_surv$ImmunePhenotype
  ),
  
  as.character(
    tcga_surv_reference_aligned$ImmunePhenotype
  )
)


age_missingness_pass <- identical(
  
  is.na(
    tcga_surv$Age
  ),
  
  is.na(
    tcga_surv_reference_aligned$Age
  )
)


age_max_difference <- max(
  abs(
    as.numeric(
      tcga_surv$Age
    ) -
      as.numeric(
        tcga_surv_reference_aligned$Age
      )
  ),
  na.rm = TRUE
)


stage_reference_pass <- identical(
  
  as.character(
    tcga_surv$Stage
  ),
  
  as.character(
    tcga_surv_reference_aligned$Stage
  )
)


survival_reference_pass <- all(
  sample_alignment_pass,
  os_time_max_difference < 1e-12,
  os_status_pass,
  phenotype_reference_pass,
  age_missingness_pass,
  age_max_difference < 1e-12,
  stage_reference_pass
)


assert_true(
  survival_reference_pass,
  "Production survival cohort differs from historical reference."
)


############################################################
## 10.15 — IMMUNE-PHENOTYPE KAPLAN-MEIER
############################################################

phenotype_survfit <- survival::survfit(
  
  survival::Surv(
    OS.time,
    OS.status
  ) ~ ImmunePhenotype,
  
  data =
    tcga_surv
)


############################################################
## 10.16 — IMMUNE-PHENOTYPE LOG-RANK TEST
############################################################

phenotype_logrank <- survival::survdiff(
  
  survival::Surv(
    OS.time,
    OS.status
  ) ~ ImmunePhenotype,
  
  data =
    tcga_surv
)


phenotype_logrank_df <-
  length(
    phenotype_logrank$n
  ) - 1L


phenotype_logrank_p <- stats::pchisq(
  
  phenotype_logrank$chisq,
  
  df =
    phenotype_logrank_df,
  
  lower.tail =
    FALSE
)


############################################################
## 10.17 — LOG-RANK CHECKPOINTS
############################################################

EXPECTED_LOGRANK_CHISQ <-
  2.8485864102032799


EXPECTED_LOGRANK_P <-
  0.24067851345139299


assert_true(
  abs(
    phenotype_logrank$chisq -
      EXPECTED_LOGRANK_CHISQ
  ) < 1e-12,
  "Immune-phenotype log-rank chi-square differs from validated result."
)


assert_true(
  abs(
    phenotype_logrank_p -
      EXPECTED_LOGRANK_P
  ) < 1e-12,
  "Immune-phenotype log-rank P-value differs from validated result."
)


assert_true(
  phenotype_logrank_df == 2L,
  "Immune-phenotype log-rank test must have df = 2."
)


############################################################
## 10.18 — PHENOTYPE SURVIVAL TABLES
############################################################

phenotype_labels <- c(
  "Immune Low",
  "Immune Mid",
  "Immune High"
)


phenotype_survival_counts <- data.frame(
  
  ImmunePhenotype =
    phenotype_labels,
  
  N =
    as.integer(
      phenotype_logrank$n
    ),
  
  Observed_deaths =
    as.integer(
      phenotype_logrank$obs
    ),
  
  Expected_deaths =
    as.numeric(
      phenotype_logrank$exp
    ),
  
  stringsAsFactors =
    FALSE
)


phenotype_logrank_table <- data.frame(
  
  Analysis =
    "Immune activity strata",
  
  N =
    nrow(
      tcga_surv
    ),
  
  Events =
    sum(
      tcga_surv$OS.status
    ),
  
  Chi_square =
    as.numeric(
      phenotype_logrank$chisq
    ),
  
  df =
    phenotype_logrank_df,
  
  Pvalue =
    phenotype_logrank_p,
  
  stringsAsFactors =
    FALSE
)


############################################################
## 10.19 — EXPRESSION EXTRACTION HELPER
##
## Explicit sample matching preserves historical behavior.
############################################################

get_tcga_gene <- function(
    gene,
    samples
) {
  
  assert_true(
    gene %in%
      rownames(
        expr_final
      ),
    paste0(
      "Gene not found in TCGA expression matrix: ",
      gene
    )
  )
  
  
  sample_index <- match(
    samples,
    colnames(
      expr_final
    )
  )
  
  
  assert_true(
    !anyNA(
      sample_index
    ),
    "One or more survival samples are absent from expression matrix."
  )
  
  
  as.numeric(
    expr_final[
      gene,
      sample_index
    ]
  )
}


############################################################
## 10.20 — 20-HUB UNIVARIATE COX SCREEN
############################################################

hub_cox_list <- vector(
  "list",
  length(
    hub_genes
  )
)


hub_cutoff_list <- vector(
  "list",
  length(
    hub_genes
  )
)


for (
  i in seq_along(
    hub_genes
  )
) {
  
  gene <- hub_genes[i]
  
  
  gene_expression <- get_tcga_gene(
    gene,
    tcga_surv$sample
  )
  
  
  gene_median <- median(
    gene_expression,
    na.rm = TRUE
  )
  
  
  gene_group <- factor(
    
    ifelse(
      gene_expression >=
        gene_median,
      "High",
      "Low"
    ),
    
    levels = c(
      "Low",
      "High"
    )
  )
  
  
  gene_survival_data <- data.frame(
    
    OS.time =
      tcga_surv$OS.time,
    
    OS.status =
      tcga_surv$OS.status,
    
    GeneGroup =
      gene_group
  )
  
  
  gene_cox <- survival::coxph(
    
    survival::Surv(
      OS.time,
      OS.status
    ) ~ GeneGroup,
    
    data =
      gene_survival_data
  )
  
  
  gene_cox_summary <- summary(
    gene_cox
  )
  
  
  hub_cox_list[[i]] <- data.frame(
    
    Gene =
      gene,
    
    HR =
      gene_cox_summary$coefficients[
        "GeneGroupHigh",
        "exp(coef)"
      ],
    
    Lower95 =
      gene_cox_summary$conf.int[
        "GeneGroupHigh",
        "lower .95"
      ],
    
    Upper95 =
      gene_cox_summary$conf.int[
        "GeneGroupHigh",
        "upper .95"
      ],
    
    Pvalue =
      gene_cox_summary$coefficients[
        "GeneGroupHigh",
        "Pr(>|z|)"
      ],
    
    stringsAsFactors =
      FALSE
  )
  
  
  hub_cutoff_list[[i]] <- data.frame(
    
    Gene =
      gene,
    
    Median =
      gene_median,
    
    N_equal_to_median =
      sum(
        gene_expression ==
          gene_median,
        na.rm = TRUE
      ),
    
    N_Low =
      sum(
        gene_group ==
          "Low"
      ),
    
    N_High =
      sum(
        gene_group ==
          "High"
      ),
    
    Rule =
      "High if expression >= exact cohort median",
    
    stringsAsFactors =
      FALSE
  )
}


hub_cox_results <- do.call(
  rbind,
  hub_cox_list
)


hub_median_cutoffs <- do.call(
  rbind,
  hub_cutoff_list
)


rownames(
  hub_cox_results
) <- NULL


rownames(
  hub_median_cutoffs
) <- NULL


############################################################
## 10.21 — BH CORRECTION ACROSS ALL 20 HUB TESTS
############################################################

hub_cox_results$FDR <- stats::p.adjust(
  
  hub_cox_results$Pvalue,
  
  method =
    "BH"
)


############################################################
## Preserve historical final ordering:
## ascending raw P-value.
############################################################

hub_cox_results <- hub_cox_results[
  order(
    hub_cox_results$Pvalue
  ),
  ,
  drop = FALSE
]


rownames(
  hub_cox_results
) <- NULL


############################################################
## 10.22 — HUB-SCREEN CHECKPOINTS
############################################################

assert_true(
  nrow(
    hub_cox_results
  ) == 20L,
  "Hub survival screen must contain 20 tests."
)


assert_true(
  setequal(
    hub_cox_results$Gene,
    hub_genes
  ),
  "Hub survival result genes differ from Module-09 hub set."
)


############################################################
## BH correction-family validation
############################################################

hub_FDR_recheck <- stats::p.adjust(
  hub_cox_results$Pvalue,
  method = "BH"
)


assert_true(
  max(
    abs(
      hub_FDR_recheck -
        hub_cox_results$FDR
    )
  ) < 1e-12,
  "Hub FDR column is not BH correction across all 20 tests."
)


############################################################
## 10.23 — FDR-SIGNIFICANT HUBS
############################################################

significant_hubs <- hub_cox_results[
  hub_cox_results$FDR <
    0.05,
  ,
  drop = FALSE
]


expected_significant_hubs <- c(
  "CCL19",
  "CXCL2",
  "IL10"
)


assert_true(
  identical(
    as.character(
      significant_hubs$Gene
    ),
    expected_significant_hubs
  ),
  "FDR-significant hub set differs from validated result."
)


assert_true(
  all(
    significant_hubs$HR[
      significant_hubs$Gene %in%
        c(
          "CCL19",
          "CXCL2"
        )
    ] < 1
  ),
  "CCL19 and CXCL2 should show favorable univariate associations."
)


assert_true(
  significant_hubs$HR[
    significant_hubs$Gene ==
      "IL10"
  ] > 1,
  "IL10 should show an adverse univariate association."
)


############################################################
## 10.24 — VALIDATED NUMERICAL HUB CHECKPOINTS
##
## These confirm the deterministic production result has not
## drifted away from the previously validated analysis.
############################################################

expected_hub_results <- data.frame(
  
  Gene = c(
    "CCL19",
    "CXCL2",
    "IL10"
  ),
  
  HR = c(
    0.59513746192959305,
    0.62426657884083003,
    1.58428285503148003
  ),
  
  FDR = c(
    0.034630720718148898,
    0.041925571038161302,
    0.041925571038161302
  ),
  
  stringsAsFactors =
    FALSE
)


expected_index <- match(
  expected_hub_results$Gene,
  hub_cox_results$Gene
)


assert_true(
  !anyNA(
    expected_index
  ),
  "Unable to find validated significant hubs in production result."
)


expected_observed <- hub_cox_results[
  expected_index,
  ,
  drop = FALSE
]


assert_true(
  max(
    abs(
      expected_observed$HR -
        expected_hub_results$HR
    )
  ) < 1e-10,
  "Validated significant-hub HR regression failed."
)


assert_true(
  max(
    abs(
      expected_observed$FDR -
        expected_hub_results$FDR
    )
  ) < 1e-10,
  "Validated significant-hub FDR regression failed."
)


############################################################
## 10.25 — SURVIVAL COHORT SUMMARY
############################################################

survival_cohort_summary <- data.frame(
  
  Metric = c(
    "Module05_tumors",
    "Survival_cohort",
    "Excluded_from_survival",
    "Events",
    "Censored",
    "Minimum_OS_days",
    "Maximum_OS_days",
    "Immune_Low",
    "Immune_Mid",
    "Immune_High",
    "Missing_Age",
    "Missing_Stage"
  ),
  
  Value = c(
    nrow(meta05),
    nrow(tcga_surv),
    nrow(meta05) -
      nrow(tcga_surv),
    sum(
      tcga_surv$OS.status == 1
    ),
    sum(
      tcga_surv$OS.status == 0
    ),
    min(
      tcga_surv$OS.time
    ),
    max(
      tcga_surv$OS.time
    ),
    as.integer(
      phenotype_counts_internal[
        "Immune_Low"
      ]
    ),
    as.integer(
      phenotype_counts_internal[
        "Immune_Mid"
      ]
    ),
    as.integer(
      phenotype_counts_internal[
        "Immune_High"
      ]
    ),
    sum(
      is.na(
        tcga_surv$Age
      )
    ),
    sum(
      is.na(
        tcga_surv$Stage
      )
    )
  ),
  
  stringsAsFactors =
    FALSE
)


############################################################
## 10.26 — ANALYSIS PROVENANCE
############################################################

survival_provenance <- data.frame(
  
  Item = c(
    "OS_time_dead",
    "OS_time_alive",
    "OS_status_dead",
    "OS_status_alive",
    "Valid_survival_filter",
    "Phenotype_survival_method",
    "Phenotype_logrank_method",
    "Hub_expression_matching",
    "Hub_grouping",
    "Median_tie_rule",
    "Hub_survival_model",
    "Hub_reference_level",
    "Multiple_testing",
    "FDR_family",
    "CCL19_IL10_joint_state_in_this_module"
  ),
  
  Value = c(
    "days_to_death.demographic",
    "days_to_last_follow_up.diagnoses",
    "1",
    "0",
    "!is.na(OS.time) & OS.time > 0 & vital status Alive/Dead",
    "Kaplan-Meier",
    "survival::survdiff",
    "Explicit sample-ID matching",
    "Cohort-specific median split",
    "High if expression >= exact median",
    "coxph(Surv(OS.time, OS.status) ~ GeneGroup)",
    "Low",
    "Benjamini-Hochberg",
    "All 20 MCC hub-gene Cox tests",
    "FALSE"
  ),
  
  stringsAsFactors =
    FALSE
)


############################################################
## 10.27 — SAVE PROCESSED SURVIVAL OBJECTS
############################################################

save_rds_checked(
  tcga_surv,
  file.path(
    PATHS$processed,
    "10_tcga_survival_cohort.rds"
  )
)


save_rds_checked(
  phenotype_survfit,
  file.path(
    PATHS$processed,
    "10_tcga_immune_phenotype_survfit.rds"
  )
)


save_rds_checked(
  phenotype_logrank,
  file.path(
    PATHS$processed,
    "10_tcga_immune_phenotype_logrank.rds"
  )
)


save_rds_checked(
  hub_cox_results,
  file.path(
    PATHS$processed,
    "10_tcga_hub_univariate_cox.rds"
  )
)


save_rds_checked(
  hub_median_cutoffs,
  file.path(
    PATHS$processed,
    "10_tcga_hub_median_cutoffs.rds"
  )
)


save_rds_checked(
  significant_hubs,
  file.path(
    PATHS$processed,
    "10_tcga_FDR_significant_survival_hubs.rds"
  )
)


############################################################
## 10.28 — EXPORT SURVIVAL TABLES
############################################################

write_csv_checked(
  survival_cohort_summary,
  file.path(
    survival_output_dir,
    "10_TCGA_survival_cohort_summary.csv"
  )
)


write_csv_checked(
  phenotype_logrank_table,
  file.path(
    survival_output_dir,
    "10_TCGA_immune_phenotype_logrank.csv"
  )
)


write_csv_checked(
  phenotype_survival_counts,
  file.path(
    survival_output_dir,
    "10_TCGA_immune_phenotype_survival_counts.csv"
  )
)


write_csv_checked(
  hub_cox_results,
  file.path(
    survival_output_dir,
    "10_TCGA_20HubGenes_univariate_Cox.csv"
  )
)


write_csv_checked(
  hub_median_cutoffs,
  file.path(
    survival_output_dir,
    "10_TCGA_20HubGenes_median_cutoffs.csv"
  )
)


write_csv_checked(
  significant_hubs,
  file.path(
    survival_output_dir,
    "10_TCGA_FDR_significant_survival_hubs.csv"
  )
)


write_csv_checked(
  survival_provenance,
  file.path(
    survival_output_dir,
    "10_TCGA_survival_method_provenance.csv"
  )
)


############################################################
## 10.29 — HISTORICAL REFERENCE MD5
############################################################

survival_reference_md5 <- tools::md5sum(
  survival_reference_file
)


survival_reference_hash <- data.frame(
  
  File =
    names(
      survival_reference_md5
    ),
  
  MD5 =
    unname(
      survival_reference_md5
    ),
  
  stringsAsFactors =
    FALSE
)


write_csv_checked(
  survival_reference_hash,
  file.path(
    PATHS$logs,
    "10_TCGA_survival_reference_MD5.csv"
  )
)


############################################################
## 10.30 — MODULE QC
############################################################

module10_qc <- data.frame(
  
  Metric = c(
    "Module05 tumors",
    "Survival cohort",
    "Excluded from survival cohort",
    "Events",
    "Censored",
    "Immune Low",
    "Immune Mid",
    "Immune High",
    "Historical sample-set regression",
    "Historical sample-alignment regression",
    "OS-time max difference",
    "OS-status regression",
    "Immune-phenotype regression",
    "Age regression",
    "Stage regression",
    "Phenotype log-rank chi-square",
    "Phenotype log-rank df",
    "Phenotype log-rank P",
    "Hub genes tested",
    "Hub median rule",
    "Hub Cox model",
    "BH correction family",
    "FDR-significant hubs",
    "CCL19 direction",
    "CXCL2 direction",
    "IL10 direction",
    "survival package version",
    "CCL19/IL10 joint-state model performed"
  ),
  
  Value = c(
    nrow(meta05),
    nrow(tcga_surv),
    nrow(meta05) -
      nrow(tcga_surv),
    sum(
      tcga_surv$OS.status == 1
    ),
    sum(
      tcga_surv$OS.status == 0
    ),
    phenotype_counts_internal[
      "Immune_Low"
    ],
    phenotype_counts_internal[
      "Immune_Mid"
    ],
    phenotype_counts_internal[
      "Immune_High"
    ],
    survival_sample_set_pass,
    sample_alignment_pass,
    format(
      os_time_max_difference,
      scientific = TRUE,
      digits = 17
    ),
    os_status_pass,
    phenotype_reference_pass,
    age_max_difference <
      1e-12,
    stage_reference_pass,
    format(
      phenotype_logrank$chisq,
      digits = 17
    ),
    phenotype_logrank_df,
    format(
      phenotype_logrank_p,
      digits = 17
    ),
    nrow(
      hub_cox_results
    ),
    "High if expression >= exact cohort median",
    "Surv(OS.time, OS.status) ~ GeneGroup",
    "BH across all 20 hub tests",
    paste(
      significant_hubs$Gene,
      collapse = ";"
    ),
    "Favorable",
    "Favorable",
    "Adverse",
    as.character(
      packageVersion(
        "survival"
      )
    ),
    FALSE
  ),
  
  stringsAsFactors =
    FALSE
)


write_csv_checked(
  module10_qc,
  file.path(
    PATHS$logs,
    "10_tcga_survival_QC.csv"
  )
)


############################################################
## 10.31 — SESSION INFORMATION
############################################################

capture.output(
  sessionInfo(),
  file = file.path(
    PATHS$logs,
    "10_tcga_survival_sessionInfo.txt"
  )
)


############################################################
## 10.32 — FINAL STATUS
############################################################

cat(
  "\n=============================================\n"
)


cat(
  "10_tcga_survival.R: PASS\n"
)


cat(
  "\nTCGA survival cohort\n"
)


cat(
  "Primary tumors upstream: ",
  nrow(meta05),
  "\n",
  sep = ""
)


cat(
  "Survival cohort: ",
  nrow(tcga_surv),
  "\n",
  sep = ""
)


cat(
  "Excluded: ",
  nrow(meta05) -
    nrow(tcga_surv),
  "\n",
  sep = ""
)


cat(
  "Events: ",
  sum(
    tcga_surv$OS.status == 1
  ),
  "\n",
  sep = ""
)


cat(
  "Censored: ",
  sum(
    tcga_surv$OS.status == 0
  ),
  "\n",
  sep = ""
)


cat(
  "Historical survival cohort regression: VERIFIED\n"
)


cat(
  "\nImmune phenotype survival\n"
)


cat(
  "Immune Low: ",
  phenotype_counts_internal[
    "Immune_Low"
  ],
  "\n",
  sep = ""
)


cat(
  "Immune Mid: ",
  phenotype_counts_internal[
    "Immune_Mid"
  ],
  "\n",
  sep = ""
)


cat(
  "Immune High: ",
  phenotype_counts_internal[
    "Immune_High"
  ],
  "\n",
  sep = ""
)


cat(
  "Log-rank chi-square: ",
  format(
    phenotype_logrank$chisq,
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "Log-rank df: ",
  phenotype_logrank_df,
  "\n",
  sep = ""
)


cat(
  "Log-rank P: ",
  format(
    phenotype_logrank_p,
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "Immune phenotype OS association: NOT SIGNIFICANT\n"
)


cat(
  "\n20-hub survival screen\n"
)


cat(
  "Hub genes tested: ",
  nrow(
    hub_cox_results
  ),
  "\n",
  sep = ""
)


cat(
  "Grouping: cohort-specific median split\n"
)


cat(
  "Tie rule: expression >= median -> High\n"
)


cat(
  "Reference group: Low\n"
)


cat(
  "BH correction family: all 20 hub-gene tests\n"
)


cat(
  "FDR-significant hubs: ",
  paste(
    significant_hubs$Gene,
    collapse = " / "
  ),
  "\n",
  sep = ""
)


cat(
  "CCL19: favorable\n"
)


cat(
  "CXCL2: favorable\n"
)


cat(
  "IL10: adverse\n"
)


cat(
  "\nScope\n"
)


cat(
  "Adjusted CCL19/IL10 models performed: FALSE\n"
)


cat(
  "CCL19/IL10 joint-state model performed: FALSE\n"
)


cat(
  "Historical survival method changed: FALSE\n"
)


cat(
  "=============================================\n"
)