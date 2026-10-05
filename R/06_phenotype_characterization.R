############################################################
## 06_phenotype_characterization.R
##
## PURPOSE
##
##   Characterize the already-frozen TCGA-BRCA immune
##   phenotypes using:
##
##     - ESTIMATE
##     - MCP-counter
##     - xCell
##
## IMPORTANT
##
##   This module DOES NOT recompute:
##
##     ESTIMATE
##     MCP-counter deconvolution
##     xCell deconvolution
##     ssGSEA
##     PCA
##     ImmunePhenotype
##
##   It reads the frozen output of Module 05.
##
## HISTORICAL INFERENCE
##
##   Omnibus:
##     stats::kruskal.test()
##
##   Multiple-testing correction:
##     stats::p.adjust(method = "BH")
##
##   Post-hoc:
##     FSA::dunnTest(method = "bh")
##
##   xCell Dunn testing:
##     only features with omnibus BH-FDR < 0.05
##
## REPRODUCIBILITY
##
##   All descriptive, Kruskal-Wallis, BH-FDR and Dunn-BH
##   results are regression-tested against frozen Table S3
##   historical references.
############################################################


############################################################
## 6.1 — REQUIRE SETUP
############################################################

if (!exists("PATHS")) {
  source("R/00_setup.R")
}


############################################################
## 6.2 — REQUIRED PACKAGE
############################################################

require_package("FSA")


############################################################
## 6.3 — EXPLICIT INPUT FILES
############################################################

metadata_file <- file.path(
  PATHS$processed,
  "05_tcga_metadata_immune_phenotype.rds"
)


reference_file <- file.path(
  "external",
  "tme_statistics",
  "TCGA_BRCA_historical_TME_statistics_reference.rds"
)


assert_file(metadata_file)
assert_file(reference_file)


############################################################
## 6.4 — OUTPUT DIRECTORY
############################################################

table_output_dir <- file.path(
  "results",
  "tables",
  "phenotype_characterization"
)


dir.create(
  table_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
## 6.5 — LOAD INPUTS
############################################################

meta05 <- readRDS(
  metadata_file
)


ref06 <- readRDS(
  reference_file
)


############################################################
## 6.6 — BASIC INPUT QC
############################################################

assert_true(
  is.data.frame(meta05),
  "Module-05 metadata is not a data.frame."
)


assert_true(
  nrow(meta05) == 1106L,
  paste0(
    "Expected 1,106 TCGA tumors; observed ",
    nrow(meta05),
    "."
  )
)


assert_true(
  "sample" %in% colnames(meta05),
  "sample column is missing."
)


assert_true(
  "ImmunePhenotype" %in% colnames(meta05),
  "ImmunePhenotype column is missing."
)


assert_unique(
  meta05$sample,
  "TCGA sample IDs"
)


############################################################
## 6.7 — FIX PHENOTYPE LEVEL ORDER
############################################################

phenotype_levels <- c(
  "Immune_Low",
  "Immune_Mid",
  "Immune_High"
)


## Ordinary factor is deliberate.
##
## FSA::dunnTest() does not require an ordered factor and
## the historical comparison order is fully preserved by
## explicitly fixing the levels.

meta05$ImmunePhenotype <- factor(
  as.character(
    meta05$ImmunePhenotype
  ),
  levels = phenotype_levels
)


phenotype_counts <- table(
  meta05$ImmunePhenotype
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
  "Immune phenotype counts differ from the frozen definition."
)


############################################################
## 6.8 — RECOVER EXACT HISTORICAL FEATURE ORDERS
############################################################

estimate_features <- unique(
  as.character(
    ref06$ESTIMATE_descriptive$Feature
  )
)


mcp_features <- unique(
  as.character(
    ref06$MCP_descriptive$Feature
  )
)


xcell_features <- unique(
  as.character(
    ref06$xCell_descriptive$Feature
  )
)


assert_true(
  length(estimate_features) == 4L,
  "Expected four ESTIMATE features."
)


assert_true(
  length(mcp_features) == 10L,
  "Expected ten MCP-counter populations."
)


assert_true(
  length(xcell_features) == 67L,
  "Expected 67 xCell outputs."
)


assert_true(
  all(
    estimate_features %in%
      colnames(meta05)
  ),
  "One or more ESTIMATE features are absent."
)


assert_true(
  all(
    mcp_features %in%
      colnames(meta05)
  ),
  "One or more MCP-counter features are absent."
)


assert_true(
  all(
    xcell_features %in%
      colnames(meta05)
  ),
  "One or more xCell features are absent."
)


############################################################
## 6.9 — FEATURE LABEL MAPS
############################################################

estimate_map <- unique(
  ref06$ESTIMATE_KW[
    ,
    c(
      "Feature",
      "Label"
    ),
    drop = FALSE
  ]
)


mcp_map <- unique(
  ref06$MCP_KW[
    ,
    c(
      "Feature",
      "Label"
    ),
    drop = FALSE
  ]
)


xcell_map <- unique(
  ref06$xCell_KW[
    ,
    c(
      "Feature",
      "Label",
      "Feature_type"
    ),
    drop = FALSE
  ]
)


############################################################
## 6.10 — GLOBAL MISSINGNESS QC
############################################################

all_tme_features <- c(
  estimate_features,
  mcp_features,
  xcell_features
)


method_missingness <- data.frame(
  
  Method = c(
    rep(
      "ESTIMATE",
      length(estimate_features)
    ),
    
    rep(
      "MCP-counter",
      length(mcp_features)
    ),
    
    rep(
      "xCell",
      length(xcell_features)
    )
  ),
  
  Feature =
    all_tme_features,
  
  Missing_n =
    vapply(
      meta05[
        ,
        all_tme_features,
        drop = FALSE
      ],
      function(x) {
        sum(is.na(x))
      },
      integer(1)
    ),
  
  stringsAsFactors = FALSE
)


method_missingness$Missing_pct <-
  100 *
  method_missingness$Missing_n /
  nrow(meta05)


assert_true(
  sum(
    method_missingness$Missing_n
  ) == 0L,
  "Missing TME characterization values detected."
)


############################################################
## 6.11 — DESCRIPTIVE SUMMARY HELPER
############################################################

summarize_features <- function(
    data,
    features,
    phenotype_levels
) {
  
  out <- lapply(
    
    features,
    
    function(feature) {
      
      do.call(
        
        rbind,
        
        lapply(
          
          phenotype_levels,
          
          function(group) {
            
            x <- data[
              data$ImmunePhenotype == group,
              feature
            ]
            
            
            data.frame(
              
              Feature =
                feature,
              
              ImmunePhenotype =
                group,
              
              N_nonmissing =
                sum(
                  !is.na(x)
                ),
              
              Median =
                median(
                  x,
                  na.rm = TRUE
                ),
              
              Q1 =
                as.numeric(
                  quantile(
                    x,
                    0.25,
                    na.rm = TRUE
                  )
                ),
              
              Q3 =
                as.numeric(
                  quantile(
                    x,
                    0.75,
                    na.rm = TRUE
                  )
                ),
              
              stringsAsFactors = FALSE
            )
          }
        )
      )
    }
  )
  
  
  out <- do.call(
    rbind,
    out
  )
  
  
  rownames(out) <- NULL
  
  out
}


############################################################
## 6.12 — xCELL DESCRIPTIVE SUMMARY HELPER
############################################################

summarize_xcell <- function(
    data,
    features,
    phenotype_levels
) {
  
  out <- lapply(
    
    features,
    
    function(feature) {
      
      do.call(
        
        rbind,
        
        lapply(
          
          phenotype_levels,
          
          function(group) {
            
            x <- data[
              data$ImmunePhenotype == group,
              feature
            ]
            
            
            data.frame(
              
              Feature =
                feature,
              
              ImmunePhenotype =
                group,
              
              N_nonmissing =
                sum(
                  !is.na(x)
                ),
              
              Zero_n =
                sum(
                  x == 0,
                  na.rm = TRUE
                ),
              
              Zero_pct =
                100 *
                mean(
                  x == 0,
                  na.rm = TRUE
                ),
              
              Median =
                median(
                  x,
                  na.rm = TRUE
                ),
              
              Q1 =
                as.numeric(
                  quantile(
                    x,
                    0.25,
                    na.rm = TRUE
                  )
                ),
              
              Q3 =
                as.numeric(
                  quantile(
                    x,
                    0.75,
                    na.rm = TRUE
                  )
                ),
              
              stringsAsFactors = FALSE
            )
          }
        )
      )
    }
  )
  
  
  out <- do.call(
    rbind,
    out
  )
  
  
  rownames(out) <- NULL
  
  out
}


############################################################
## 6.13 — KRUSKAL-WALLIS HELPER
############################################################

run_kw <- function(
    data,
    features
) {
  
  out <- lapply(
    
    features,
    
    function(feature) {
      
      keep <-
        !is.na(
          data[[feature]]
        ) &
        !is.na(
          data$ImmunePhenotype
        )
      
      
      kt <- stats::kruskal.test(
        
        x =
          data[[feature]][keep],
        
        g =
          factor(
            data$ImmunePhenotype[keep],
            levels = phenotype_levels
          )
      )
      
      
      data.frame(
        
        Feature =
          feature,
        
        KW_chisq =
          unname(
            kt$statistic
          ),
        
        df =
          unname(
            kt$parameter
          ),
        
        Pvalue =
          kt$p.value,
        
        stringsAsFactors = FALSE
      )
    }
  )
  
  
  out <- do.call(
    rbind,
    out
  )
  
  
  rownames(out) <- NULL
  
  
  out$BH_FDR <- stats::p.adjust(
    out$Pvalue,
    method = "BH"
  )
  
  
  out
}


############################################################
## 6.14 — DUNN-BH HELPER
############################################################

run_dunn_BH <- function(
    x,
    group
) {
  
  keep <-
    !is.na(x) &
    !is.na(group)
  
  
  x <- x[keep]
  
  
  group <- factor(
    as.character(
      group[keep]
    ),
    levels = phenotype_levels
  )
  
  
  group <- droplevels(
    group
  )
  
  
  dunn_result <- NULL
  
  
  ## Suppress the long console tables printed by the
  ## dunn.test backend. The returned statistics are unchanged.
  
  invisible(
    capture.output(
      
      dunn_result <-
        FSA::dunnTest(
          x,
          group,
          method = "bh"
        )$res
    )
  )
  
  
  data.frame(
    
    Comparison =
      as.character(
        dunn_result$Comparison
      ),
    
    Z =
      as.numeric(
        dunn_result$Z
      ),
    
    P.unadj =
      as.numeric(
        dunn_result$P.unadj
      ),
    
    P.adj =
      as.numeric(
        dunn_result$P.adj
      ),
    
    stringsAsFactors = FALSE
  )
}


############################################################
## 6.15 — ESTIMATE DESCRIPTIVE TABLE
############################################################

estimate_summary <- summarize_features(
  data = meta05,
  features = estimate_features,
  phenotype_levels = phenotype_levels
)


estimate_summary$Label <-
  estimate_map$Label[
    match(
      estimate_summary$Feature,
      estimate_map$Feature
    )
  ]


estimate_summary <- estimate_summary[
  ,
  c(
    "Feature",
    "ImmunePhenotype",
    "N_nonmissing",
    "Median",
    "Q1",
    "Q3",
    "Label"
  ),
  drop = FALSE
]


assert_true(
  nrow(estimate_summary) == 12L,
  "ESTIMATE descriptive table should contain 12 rows."
)


############################################################
## 6.16 — ESTIMATE KRUSKAL-WALLIS
############################################################

estimate_kw <- run_kw(
  data = meta05,
  features = estimate_features
)


estimate_kw$Label <-
  estimate_map$Label[
    match(
      estimate_kw$Feature,
      estimate_map$Feature
    )
  ]


estimate_kw <- estimate_kw[
  ,
  c(
    "Feature",
    "KW_chisq",
    "df",
    "Pvalue",
    "BH_FDR",
    "Label"
  ),
  drop = FALSE
]


assert_true(
  nrow(estimate_kw) == 4L,
  "ESTIMATE omnibus table should contain four rows."
)


############################################################
## 6.17 — ESTIMATE DUNN-BH
############################################################

estimate_dunn <- do.call(
  
  rbind,
  
  lapply(
    
    estimate_features,
    
    function(feature) {
      
      tmp <- run_dunn_BH(
        x = meta05[[feature]],
        group = meta05$ImmunePhenotype
      )
      
      
      tmp$Feature <- feature
      
      tmp
    }
  )
)


rownames(
  estimate_dunn
) <- NULL


estimate_dunn$Label <-
  estimate_map$Label[
    match(
      estimate_dunn$Feature,
      estimate_map$Feature
    )
  ]


estimate_dunn <- estimate_dunn[
  ,
  c(
    "Feature",
    "Label",
    "Comparison",
    "Z",
    "P.unadj",
    "P.adj"
  ),
  drop = FALSE
]


assert_true(
  nrow(estimate_dunn) == 12L,
  "ESTIMATE Dunn table should contain 12 rows."
)


############################################################
## 6.18 — MCP-COUNTER DESCRIPTIVE TABLE
############################################################

mcp_summary <- summarize_features(
  data = meta05,
  features = mcp_features,
  phenotype_levels = phenotype_levels
)


mcp_summary$Label <-
  mcp_map$Label[
    match(
      mcp_summary$Feature,
      mcp_map$Feature
    )
  ]


mcp_summary <- mcp_summary[
  ,
  c(
    "Feature",
    "ImmunePhenotype",
    "N_nonmissing",
    "Median",
    "Q1",
    "Q3",
    "Label"
  ),
  drop = FALSE
]


assert_true(
  nrow(mcp_summary) == 30L,
  "MCP-counter descriptive table should contain 30 rows."
)


############################################################
## 6.19 — MCP-COUNTER KRUSKAL-WALLIS
############################################################

mcp_kw <- run_kw(
  data = meta05,
  features = mcp_features
)


mcp_kw$Label <-
  mcp_map$Label[
    match(
      mcp_kw$Feature,
      mcp_map$Feature
    )
  ]


mcp_kw <- mcp_kw[
  ,
  c(
    "Feature",
    "KW_chisq",
    "df",
    "Pvalue",
    "BH_FDR",
    "Label"
  ),
  drop = FALSE
]


assert_true(
  nrow(mcp_kw) == 10L,
  "MCP-counter omnibus table should contain ten rows."
)


############################################################
## 6.20 — MCP-COUNTER DUNN-BH
############################################################

mcp_dunn <- do.call(
  
  rbind,
  
  lapply(
    
    mcp_features,
    
    function(feature) {
      
      tmp <- run_dunn_BH(
        x = meta05[[feature]],
        group = meta05$ImmunePhenotype
      )
      
      
      tmp$Feature <- feature
      
      tmp
    }
  )
)


rownames(
  mcp_dunn
) <- NULL


mcp_dunn$Label <-
  mcp_map$Label[
    match(
      mcp_dunn$Feature,
      mcp_map$Feature
    )
  ]


mcp_dunn <- mcp_dunn[
  ,
  c(
    "Feature",
    "Label",
    "Comparison",
    "Z",
    "P.unadj",
    "P.adj"
  ),
  drop = FALSE
]


assert_true(
  nrow(mcp_dunn) == 30L,
  "MCP-counter Dunn table should contain 30 rows."
)


############################################################
## 6.21 — xCELL DESCRIPTIVE TABLE
############################################################

xcell_summary <- summarize_xcell(
  data = meta05,
  features = xcell_features,
  phenotype_levels = phenotype_levels
)


xcell_summary$Label <-
  xcell_map$Label[
    match(
      xcell_summary$Feature,
      xcell_map$Feature
    )
  ]


xcell_summary$Feature_type <-
  xcell_map$Feature_type[
    match(
      xcell_summary$Feature,
      xcell_map$Feature
    )
  ]


xcell_summary <- xcell_summary[
  ,
  c(
    "Feature",
    "ImmunePhenotype",
    "N_nonmissing",
    "Zero_n",
    "Zero_pct",
    "Median",
    "Q1",
    "Q3",
    "Label",
    "Feature_type"
  ),
  drop = FALSE
]


assert_true(
  nrow(xcell_summary) == 201L,
  "xCell descriptive table should contain 201 rows."
)


############################################################
## 6.22 — xCELL GLOBAL DISTRIBUTION QC
############################################################

xcell_distribution_QC <- do.call(
  
  rbind,
  
  lapply(
    
    xcell_features,
    
    function(feature) {
      
      x <- meta05[[feature]]
      
      
      data.frame(
        
        Feature =
          feature,
        
        N =
          length(x),
        
        Missing_n =
          sum(
            is.na(x)
          ),
        
        Zero_n =
          sum(
            x == 0,
            na.rm = TRUE
          ),
        
        Zero_pct =
          100 *
          mean(
            x == 0,
            na.rm = TRUE
          ),
        
        Min =
          min(
            x,
            na.rm = TRUE
          ),
        
        Q1 =
          as.numeric(
            quantile(
              x,
              0.25,
              na.rm = TRUE
            )
          ),
        
        Median =
          median(
            x,
            na.rm = TRUE
          ),
        
        Q3 =
          as.numeric(
            quantile(
              x,
              0.75,
              na.rm = TRUE
            )
          ),
        
        Max =
          max(
            x,
            na.rm = TRUE
          ),
        
        stringsAsFactors = FALSE
      )
    }
  )
)


rownames(
  xcell_distribution_QC
) <- NULL


############################################################
## 6.23 — xCELL KRUSKAL-WALLIS
############################################################

xcell_kw <- run_kw(
  data = meta05,
  features = xcell_features
)


## Historical table was ordered by BH-FDR and then raw P.

xcell_kw <- xcell_kw[
  order(
    xcell_kw$BH_FDR,
    xcell_kw$Pvalue
  ),
  ,
  drop = FALSE
]


rownames(
  xcell_kw
) <- NULL


xcell_kw$Label <-
  xcell_map$Label[
    match(
      xcell_kw$Feature,
      xcell_map$Feature
    )
  ]


xcell_kw$Feature_type <-
  xcell_map$Feature_type[
    match(
      xcell_kw$Feature,
      xcell_map$Feature
    )
  ]


xcell_kw$Zero_n <-
  xcell_distribution_QC$Zero_n[
    match(
      xcell_kw$Feature,
      xcell_distribution_QC$Feature
    )
  ]


xcell_kw$Zero_pct <-
  xcell_distribution_QC$Zero_pct[
    match(
      xcell_kw$Feature,
      xcell_distribution_QC$Feature
    )
  ]


xcell_kw$Significant_FDR05 <-
  xcell_kw$BH_FDR < 0.05


xcell_kw$Significant_FDR01 <-
  xcell_kw$BH_FDR < 0.01


xcell_kw <- xcell_kw[
  ,
  c(
    "Feature",
    "KW_chisq",
    "df",
    "Pvalue",
    "BH_FDR",
    "Label",
    "Feature_type",
    "Zero_n",
    "Zero_pct",
    "Significant_FDR05",
    "Significant_FDR01"
  ),
  drop = FALSE
]


assert_true(
  nrow(xcell_kw) == 67L,
  "xCell omnibus table should contain 67 rows."
)


############################################################
## 6.24 — xCELL SIGNIFICANT FEATURE SET
############################################################

xcell_sig_features <-
  xcell_kw$Feature[
    xcell_kw$Significant_FDR05
  ]


assert_true(
  length(
    xcell_sig_features
  ) == 25L,
  paste0(
    "Expected 25 xCell FDR-significant features; observed ",
    length(xcell_sig_features),
    "."
  )
)


############################################################
## 6.25 — xCELL DUNN-BH
##
## Only omnibus BH-FDR < 0.05 features are tested.
############################################################

xcell_dunn <- do.call(
  
  rbind,
  
  lapply(
    
    xcell_sig_features,
    
    function(feature) {
      
      tmp <- run_dunn_BH(
        x = meta05[[feature]],
        group = meta05$ImmunePhenotype
      )
      
      
      tmp$Feature <- feature
      
      tmp
    }
  )
)


rownames(
  xcell_dunn
) <- NULL


xcell_dunn$Label <-
  xcell_map$Label[
    match(
      xcell_dunn$Feature,
      xcell_map$Feature
    )
  ]


xcell_dunn$Feature_type <-
  xcell_map$Feature_type[
    match(
      xcell_dunn$Feature,
      xcell_map$Feature
    )
  ]


xcell_dunn <- xcell_dunn[
  ,
  c(
    "Feature",
    "Label",
    "Feature_type",
    "Comparison",
    "Z",
    "P.unadj",
    "P.adj"
  ),
  drop = FALSE
]


assert_true(
  nrow(xcell_dunn) == 75L,
  "xCell Dunn table should contain 75 rows."
)


############################################################
## 6.26 — xCELL LOW-vs-HIGH TABLE
############################################################

xcell_low_high <- xcell_dunn[
  
  grepl(
    "Immune_Low",
    xcell_dunn$Comparison
  ) &
    
    grepl(
      "Immune_High",
      xcell_dunn$Comparison
    ),
  
  ,
  drop = FALSE
]


assert_true(
  nrow(xcell_low_high) == 25L,
  "xCell Low-vs-High table should contain 25 rows."
)


############################################################
## 6.27 — DESCRIPTIVE xCELL DIRECTION TABLE
############################################################

xcell_direction <- do.call(
  
  rbind,
  
  lapply(
    
    xcell_features,
    
    function(feature) {
      
      ss <- xcell_summary[
        xcell_summary$Feature == feature,
        ,
        drop = FALSE
      ]
      
      
      low_med <-
        ss$Median[
          ss$ImmunePhenotype ==
            "Immune_Low"
        ]
      
      
      mid_med <-
        ss$Median[
          ss$ImmunePhenotype ==
            "Immune_Mid"
        ]
      
      
      high_med <-
        ss$Median[
          ss$ImmunePhenotype ==
            "Immune_High"
        ]
      
      
      data.frame(
        
        Feature =
          feature,
        
        Low_median =
          low_med,
        
        Mid_median =
          mid_med,
        
        High_median =
          high_med,
        
        Direction_Low_to_High =
          ifelse(
            high_med > low_med,
            "Higher in Immune High",
            ifelse(
              high_med < low_med,
              "Lower in Immune High",
              "Equal median"
            )
          ),
        
        stringsAsFactors = FALSE
      )
    }
  )
)


xcell_direction$Label <-
  xcell_map$Label[
    match(
      xcell_direction$Feature,
      xcell_map$Feature
    )
  ]


xcell_direction$Feature_type <-
  xcell_map$Feature_type[
    match(
      xcell_direction$Feature,
      xcell_map$Feature
    )
  ]


xcell_direction$BH_FDR <-
  xcell_kw$BH_FDR[
    match(
      xcell_direction$Feature,
      xcell_kw$Feature
    )
  ]


xcell_direction$Significant_FDR05 <-
  xcell_direction$BH_FDR < 0.05


############################################################
## 6.28 — HISTORICAL TABLE REGRESSION HELPER
############################################################

compare_to_reference <- function(
    current,
    historical,
    key_columns,
    numeric_columns
) {
  
  current_key <- do.call(
    paste,
    c(
      current[
        ,
        key_columns,
        drop = FALSE
      ],
      sep = "|||"
    )
  )
  
  
  historical_key <- do.call(
    paste,
    c(
      historical[
        ,
        key_columns,
        drop = FALSE
      ],
      sep = "|||"
    )
  )
  
  
  assert_true(
    !anyDuplicated(current_key),
    "Duplicate production regression keys detected."
  )
  
  
  assert_true(
    !anyDuplicated(historical_key),
    "Duplicate historical regression keys detected."
  )
  
  
  assert_true(
    setequal(
      current_key,
      historical_key
    ),
    "Production and historical table keys differ."
  )
  
  
  historical_aligned <- historical[
    match(
      current_key,
      historical_key
    ),
    ,
    drop = FALSE
  ]
  
  
  numeric_differences <- vapply(
    
    numeric_columns,
    
    function(column) {
      
      max(
        abs(
          as.numeric(
            current[[column]]
          ) -
            as.numeric(
              historical_aligned[[column]]
            )
        ),
        na.rm = TRUE
      )
    },
    
    numeric(1)
  )
  
  
  nonnumeric_columns <- setdiff(
    intersect(
      colnames(current),
      colnames(historical_aligned)
    ),
    numeric_columns
  )
  
  
  exact_matches <- vapply(
    
    nonnumeric_columns,
    
    function(column) {
      
      identical(
        as.character(
          current[[column]]
        ),
        as.character(
          historical_aligned[[column]]
        )
      )
    },
    
    logical(1)
  )
  
  
  list(
    
    max_numeric_difference =
      max(
        numeric_differences
      ),
    
    numeric_differences =
      numeric_differences,
    
    exact_nonnumeric_match =
      all(
        exact_matches
      ),
    
    nonnumeric_matches =
      exact_matches
  )
}


############################################################
## 6.29 — DESCRIPTIVE HISTORICAL REGRESSION
############################################################

estimate_desc_reg <- compare_to_reference(
  
  current =
    estimate_summary,
  
  historical =
    ref06$ESTIMATE_descriptive,
  
  key_columns = c(
    "Feature",
    "ImmunePhenotype"
  ),
  
  numeric_columns = c(
    "N_nonmissing",
    "Median",
    "Q1",
    "Q3"
  )
)


mcp_desc_reg <- compare_to_reference(
  
  current =
    mcp_summary,
  
  historical =
    ref06$MCP_descriptive,
  
  key_columns = c(
    "Feature",
    "ImmunePhenotype"
  ),
  
  numeric_columns = c(
    "N_nonmissing",
    "Median",
    "Q1",
    "Q3"
  )
)


xcell_desc_reg <- compare_to_reference(
  
  current =
    xcell_summary,
  
  historical =
    ref06$xCell_descriptive,
  
  key_columns = c(
    "Feature",
    "ImmunePhenotype"
  ),
  
  numeric_columns = c(
    "N_nonmissing",
    "Zero_n",
    "Zero_pct",
    "Median",
    "Q1",
    "Q3"
  )
)


############################################################
## 6.30 — KW HISTORICAL REGRESSION
############################################################

estimate_kw_reg <- compare_to_reference(
  
  current =
    estimate_kw,
  
  historical =
    ref06$ESTIMATE_KW,
  
  key_columns =
    "Feature",
  
  numeric_columns = c(
    "KW_chisq",
    "df",
    "Pvalue",
    "BH_FDR"
  )
)


mcp_kw_reg <- compare_to_reference(
  
  current =
    mcp_kw,
  
  historical =
    ref06$MCP_KW,
  
  key_columns =
    "Feature",
  
  numeric_columns = c(
    "KW_chisq",
    "df",
    "Pvalue",
    "BH_FDR"
  )
)


xcell_kw_reg <- compare_to_reference(
  
  current =
    xcell_kw,
  
  historical =
    ref06$xCell_KW,
  
  key_columns =
    "Feature",
  
  numeric_columns = c(
    "KW_chisq",
    "df",
    "Pvalue",
    "BH_FDR",
    "Zero_n",
    "Zero_pct"
  )
)


############################################################
## 6.31 — DUNN HISTORICAL REGRESSION
############################################################

estimate_dunn_reg <- compare_to_reference(
  
  current =
    estimate_dunn,
  
  historical =
    ref06$ESTIMATE_Dunn,
  
  key_columns = c(
    "Feature",
    "Comparison"
  ),
  
  numeric_columns = c(
    "Z",
    "P.unadj",
    "P.adj"
  )
)


mcp_dunn_reg <- compare_to_reference(
  
  current =
    mcp_dunn,
  
  historical =
    ref06$MCP_Dunn,
  
  key_columns = c(
    "Feature",
    "Comparison"
  ),
  
  numeric_columns = c(
    "Z",
    "P.unadj",
    "P.adj"
  )
)


xcell_dunn_reg <- compare_to_reference(
  
  current =
    xcell_dunn,
  
  historical =
    ref06$xCell_Dunn,
  
  key_columns = c(
    "Feature",
    "Comparison"
  ),
  
  numeric_columns = c(
    "Z",
    "P.unadj",
    "P.adj"
  )
)


############################################################
## 6.32 — REGRESSION TOLERANCES
############################################################

DESCRIPTIVE_TOLERANCE <- 1e-10

KW_TOLERANCE <- 1e-10

DUNN_TOLERANCE <- 1e-10


descriptive_regression_pass <- all(
  
  estimate_desc_reg$max_numeric_difference <
    DESCRIPTIVE_TOLERANCE,
  
  mcp_desc_reg$max_numeric_difference <
    DESCRIPTIVE_TOLERANCE,
  
  xcell_desc_reg$max_numeric_difference <
    DESCRIPTIVE_TOLERANCE,
  
  estimate_desc_reg$exact_nonnumeric_match,
  
  mcp_desc_reg$exact_nonnumeric_match,
  
  xcell_desc_reg$exact_nonnumeric_match
)


kw_regression_pass <- all(
  
  estimate_kw_reg$max_numeric_difference <
    KW_TOLERANCE,
  
  mcp_kw_reg$max_numeric_difference <
    KW_TOLERANCE,
  
  xcell_kw_reg$max_numeric_difference <
    KW_TOLERANCE,
  
  estimate_kw_reg$exact_nonnumeric_match,
  
  mcp_kw_reg$exact_nonnumeric_match,
  
  xcell_kw_reg$exact_nonnumeric_match
)


dunn_regression_pass <- all(
  
  estimate_dunn_reg$max_numeric_difference <
    DUNN_TOLERANCE,
  
  mcp_dunn_reg$max_numeric_difference <
    DUNN_TOLERANCE,
  
  xcell_dunn_reg$max_numeric_difference <
    DUNN_TOLERANCE,
  
  estimate_dunn_reg$exact_nonnumeric_match,
  
  mcp_dunn_reg$exact_nonnumeric_match,
  
  xcell_dunn_reg$exact_nonnumeric_match
)


assert_true(
  descriptive_regression_pass,
  "Historical descriptive-table regression failed."
)


assert_true(
  kw_regression_pass,
  "Historical Kruskal-Wallis regression failed."
)


assert_true(
  dunn_regression_pass,
  "Historical Dunn-BH regression failed."
)


############################################################
## 6.33 — KNOWN FROZEN SCIENTIFIC CHECKPOINTS
############################################################

assert_true(
  sum(
    estimate_kw$BH_FDR < 0.001
  ) == 4L,
  "Expected all four ESTIMATE omnibus tests at FDR < 0.001."
)


assert_true(
  sum(
    estimate_dunn$P.adj < 0.05
  ) == 12L,
  "Expected all 12 ESTIMATE pairwise comparisons significant."
)


assert_true(
  sum(
    mcp_kw$BH_FDR < 0.001
  ) == 10L,
  "Expected all ten MCP-counter omnibus tests at FDR < 0.001."
)


assert_true(
  sum(
    xcell_kw$BH_FDR < 0.05
  ) == 25L,
  "Expected 25 significant xCell outputs at FDR < 0.05."
)


assert_true(
  sum(
    xcell_kw$BH_FDR < 0.01
  ) == 24L,
  "Expected 24 significant xCell outputs at FDR < 0.01."
)


assert_true(
  sum(
    xcell_kw$BH_FDR < 0.001
  ) == 21L,
  "Expected 21 significant xCell outputs at FDR < 0.001."
)


############################################################
## 6.34 — KNOWN FIBROBLAST CHECKPOINT
############################################################

fibroblast_mid_high <- mcp_dunn[
  
  mcp_dunn$Feature ==
    "Fibroblasts_MCPcounter" &
    
    grepl(
      "Immune_Mid",
      mcp_dunn$Comparison
    ) &
    
    grepl(
      "Immune_High",
      mcp_dunn$Comparison
    ),
  
  ,
  drop = FALSE
]


assert_true(
  nrow(
    fibroblast_mid_high
  ) == 1L,
  "Could not uniquely identify fibroblast Mid-vs-High comparison."
)


assert_true(
  abs(
    fibroblast_mid_high$P.adj -
      0.16709669369441832
  ) <
    1e-12,
  "Fibroblast Mid-vs-High adjusted P differs from frozen value."
)


############################################################
## 6.35 — SAVE PRODUCTION DATA BUNDLE
############################################################

phenotype_characterization <- list(
  
  phenotype_counts =
    phenotype_counts,
  
  method_missingness =
    method_missingness,
  
  ESTIMATE = list(
    
    descriptive =
      estimate_summary,
    
    Kruskal_Wallis =
      estimate_kw,
    
    Dunn_BH =
      estimate_dunn
  ),
  
  MCPcounter = list(
    
    descriptive =
      mcp_summary,
    
    Kruskal_Wallis =
      mcp_kw,
    
    Dunn_BH =
      mcp_dunn
  ),
  
  xCell = list(
    
    descriptive =
      xcell_summary,
    
    distribution_QC =
      xcell_distribution_QC,
    
    Kruskal_Wallis =
      xcell_kw,
    
    significant_features =
      xcell_sig_features,
    
    Dunn_BH =
      xcell_dunn,
    
    Low_vs_High =
      xcell_low_high,
    
    direction =
      xcell_direction
  )
)


save_rds_checked(
  phenotype_characterization,
  file.path(
    PATHS$processed,
    "06_tcga_phenotype_characterization.rds"
  )
)


############################################################
## 6.36 — SAVE INDIVIDUAL PRODUCTION TABLES
############################################################

write_csv_checked(
  method_missingness,
  file.path(
    table_output_dir,
    "06_TME_missingness.csv"
  )
)


write_csv_checked(
  estimate_summary,
  file.path(
    table_output_dir,
    "06_ESTIMATE_descriptive.csv"
  )
)


write_csv_checked(
  estimate_kw,
  file.path(
    table_output_dir,
    "06_ESTIMATE_KruskalWallis_BH.csv"
  )
)


write_csv_checked(
  estimate_dunn,
  file.path(
    table_output_dir,
    "06_ESTIMATE_Dunn_BH.csv"
  )
)


write_csv_checked(
  mcp_summary,
  file.path(
    table_output_dir,
    "06_MCPcounter_descriptive.csv"
  )
)


write_csv_checked(
  mcp_kw,
  file.path(
    table_output_dir,
    "06_MCPcounter_KruskalWallis_BH.csv"
  )
)


write_csv_checked(
  mcp_dunn,
  file.path(
    table_output_dir,
    "06_MCPcounter_Dunn_BH.csv"
  )
)


write_csv_checked(
  xcell_summary,
  file.path(
    table_output_dir,
    "06_xCell_descriptive.csv"
  )
)


write_csv_checked(
  xcell_distribution_QC,
  file.path(
    table_output_dir,
    "06_xCell_distribution_QC.csv"
  )
)


write_csv_checked(
  xcell_kw,
  file.path(
    table_output_dir,
    "06_xCell_KruskalWallis_BH.csv"
  )
)


write_csv_checked(
  xcell_dunn,
  file.path(
    table_output_dir,
    "06_xCell_Dunn_BH_significant_features.csv"
  )
)


write_csv_checked(
  xcell_low_high,
  file.path(
    table_output_dir,
    "06_xCell_Low_vs_High_Dunn_BH.csv"
  )
)


write_csv_checked(
  xcell_direction,
  file.path(
    table_output_dir,
    "06_xCell_direction.csv"
  )
)


############################################################
## 6.37 — REGRESSION QC TABLE
############################################################

module06_regression <- data.frame(
  
  Component = c(
    
    "ESTIMATE descriptive",
    "MCP-counter descriptive",
    "xCell descriptive",
    
    "ESTIMATE Kruskal-Wallis",
    "MCP-counter Kruskal-Wallis",
    "xCell Kruskal-Wallis",
    
    "ESTIMATE Dunn-BH",
    "MCP-counter Dunn-BH",
    "xCell Dunn-BH"
  ),
  
  Max_numeric_difference = c(
    
    estimate_desc_reg$max_numeric_difference,
    mcp_desc_reg$max_numeric_difference,
    xcell_desc_reg$max_numeric_difference,
    
    estimate_kw_reg$max_numeric_difference,
    mcp_kw_reg$max_numeric_difference,
    xcell_kw_reg$max_numeric_difference,
    
    estimate_dunn_reg$max_numeric_difference,
    mcp_dunn_reg$max_numeric_difference,
    xcell_dunn_reg$max_numeric_difference
  ),
  
  Exact_nonnumeric_match = c(
    
    estimate_desc_reg$exact_nonnumeric_match,
    mcp_desc_reg$exact_nonnumeric_match,
    xcell_desc_reg$exact_nonnumeric_match,
    
    estimate_kw_reg$exact_nonnumeric_match,
    mcp_kw_reg$exact_nonnumeric_match,
    xcell_kw_reg$exact_nonnumeric_match,
    
    estimate_dunn_reg$exact_nonnumeric_match,
    mcp_dunn_reg$exact_nonnumeric_match,
    xcell_dunn_reg$exact_nonnumeric_match
  ),
  
  stringsAsFactors = FALSE
)


write_csv_checked(
  module06_regression,
  file.path(
    PATHS$logs,
    "06_phenotype_characterization_historical_regression.csv"
  )
)


############################################################
## 6.38 — REFERENCE FILE HASH
############################################################

reference_md5 <- unname(
  tools::md5sum(
    reference_file
  )
)


############################################################
## 6.39 — MODULE QC SUMMARY
############################################################

module06_qc <- data.frame(
  
  Metric = c(
    
    "TCGA tumors",
    
    "Immune Low tumors",
    "Immune Mid tumors",
    "Immune High tumors",
    
    "ESTIMATE features",
    "MCP-counter features",
    "xCell features",
    
    "Total TME missing values",
    
    "ESTIMATE omnibus BH-FDR < 0.001",
    "ESTIMATE significant Dunn comparisons",
    
    "MCP omnibus BH-FDR < 0.001",
    
    "Fibroblast Mid-vs-High Dunn-BH P",
    
    "xCell omnibus BH-FDR < 0.05",
    "xCell omnibus BH-FDR < 0.01",
    "xCell omnibus BH-FDR < 0.001",
    
    "xCell Dunn-tested features",
    "xCell Dunn comparisons",
    
    "Descriptive historical regression",
    "Kruskal-Wallis historical regression",
    "Dunn-BH historical regression",
    
    "FSA version",
    
    "Frozen TME reference MD5",
    
    "Inferential scale",
    "Dunn implementation"
  ),
  
  Value = c(
    
    nrow(meta05),
    
    phenotype_counts["Immune_Low"],
    phenotype_counts["Immune_Mid"],
    phenotype_counts["Immune_High"],
    
    length(estimate_features),
    length(mcp_features),
    length(xcell_features),
    
    sum(
      method_missingness$Missing_n
    ),
    
    sum(
      estimate_kw$BH_FDR < 0.001
    ),
    
    sum(
      estimate_dunn$P.adj < 0.05
    ),
    
    sum(
      mcp_kw$BH_FDR < 0.001
    ),
    
    format(
      fibroblast_mid_high$P.adj,
      digits = 17
    ),
    
    sum(
      xcell_kw$BH_FDR < 0.05
    ),
    
    sum(
      xcell_kw$BH_FDR < 0.01
    ),
    
    sum(
      xcell_kw$BH_FDR < 0.001
    ),
    
    length(
      xcell_sig_features
    ),
    
    nrow(
      xcell_dunn
    ),
    
    descriptive_regression_pass,
    kw_regression_pass,
    dunn_regression_pass,
    
    as.character(
      packageVersion("FSA")
    ),
    
    reference_md5,
    
    "Original unstandardized ESTIMATE/MCP-counter/xCell scores",
    
    "FSA::dunnTest(method = 'bh')"
  ),
  
  stringsAsFactors = FALSE
)


write_csv_checked(
  module06_qc,
  file.path(
    PATHS$logs,
    "06_phenotype_characterization_QC.csv"
  )
)


############################################################
## 6.40 — SESSION INFORMATION
############################################################

capture.output(
  sessionInfo(),
  file = file.path(
    PATHS$logs,
    "06_phenotype_characterization_sessionInfo.txt"
  )
)


############################################################
## 6.41 — FINAL STATUS
############################################################

cat(
  "\n=============================================\n"
)


cat(
  "06_phenotype_characterization.R: PASS\n"
)


cat(
  "TCGA tumors: ",
  nrow(meta05),
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
  "ESTIMATE features: ",
  length(estimate_features),
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
  "Total TME missing values: ",
  sum(
    method_missingness$Missing_n
  ),
  "\n",
  sep = ""
)


cat(
  "ESTIMATE omnibus FDR < 0.001: ",
  sum(
    estimate_kw$BH_FDR < 0.001
  ),
  " / 4\n",
  sep = ""
)


cat(
  "MCP-counter omnibus FDR < 0.001: ",
  sum(
    mcp_kw$BH_FDR < 0.001
  ),
  " / 10\n",
  sep = ""
)


cat(
  "Fibroblast Mid-vs-High Dunn-BH P: ",
  format(
    fibroblast_mid_high$P.adj,
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "xCell omnibus FDR < 0.05: ",
  sum(
    xcell_kw$BH_FDR < 0.05
  ),
  " / 67\n",
  sep = ""
)


cat(
  "xCell omnibus FDR < 0.01: ",
  sum(
    xcell_kw$BH_FDR < 0.01
  ),
  " / 67\n",
  sep = ""
)


cat(
  "xCell omnibus FDR < 0.001: ",
  sum(
    xcell_kw$BH_FDR < 0.001
  ),
  " / 67\n",
  sep = ""
)


cat(
  "xCell Dunn-tested features: ",
  length(
    xcell_sig_features
  ),
  "\n",
  sep = ""
)


cat(
  "xCell Dunn comparisons: ",
  nrow(
    xcell_dunn
  ),
  "\n",
  sep = ""
)


cat(
  "Descriptive historical regression: VERIFIED\n"
)


cat(
  "Kruskal-Wallis historical regression: VERIFIED\n"
)


cat(
  "Dunn-BH historical regression: VERIFIED\n"
)


cat(
  "ESTIMATE/MCP-counter/xCell recomputation: NONE\n"
)


cat(
  "ssGSEA/PCA/phenotype recomputation: NONE\n"
)


cat(
  "Inferential statistics used original unstandardized TME scores\n"
)


cat(
  "=============================================\n"
)