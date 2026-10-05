############################################################
## 08_functional_enrichment.R
##
## PURPOSE
##
##   Reproducible production module for the validated
##   functional-enrichment analyses downstream of Module 07.
##
## ANALYSES PRESERVED
##
##   A. GO Biological Process ORA
##   B. KEGG ORA
##   C. Hallmark GSEA
##   D. Five-Hallmark leading-edge union for PPI analysis
##
##
## IMPORTANT REFACTORING RULE
##
##   This module does NOT redesign or update the previously
##   validated biological analyses.
##
##   It preserves the frozen corrected analysis while removing
##   dependence on:
##
##     - hidden R-session state
##     - changing online KEGG annotation
##     - changing MSigDB content
##     - stochastic GSEA reruns
##
##
## GO-BP
##
##   GO-BP is deterministically recomputed because the exact
##   corrected analysis has already been shown to reproduce
##   the frozen result with numerical difference = 0.
##
##
## KEGG
##
##   The corrected historical analysis used:
##
##       enrichKEGG(..., use_internal_data = FALSE)
##
##   which accesses changing online KEGG annotation.
##
##   The current audit showed annotation drift:
##
##       frozen effective universe = 8637
##       current effective universe = 8640
##
##   Therefore the validated frozen corrected KEGG result
##   is used as the production result.
##
##   The analytical method is documented but NOT rerun against
##   a newer KEGG database.
##
##
## HALLMARK GSEA
##
##   The exact Module-07 ranked gene list is reconstructed and
##   validated against the frozen GSEA rank.
##
##   Frozen MSigDB version:
##
##       2026.1.Hs
##
##   Frozen Hallmark TERM2GENE:
##
##       50 gene sets
##       7,322 pathway-gene pairs
##
##   The previously validated GSEA output is preserved rather
##   than replaced by a new stochastic/backend-dependent rerun.
##
##
## PPI INPUT
##
##   Five selected Hallmark leading edges are reconstructed
##   from the frozen validated GSEA core-enrichment results.
##
##   Expected:
##
##       memberships = 622
##       unique genes = 455
##
############################################################


############################################################
## 8.1 — REQUIRE PROJECT SETUP
############################################################

if (!exists("PATHS")) {
  source("R/00_setup.R")
}


############################################################
## 8.2 — REQUIRED PACKAGES
############################################################

require_package(
  "clusterProfiler"
)


require_package(
  "org.Hs.eg.db"
)


############################################################
## 8.3 — INPUT FILES
############################################################

deg_complete_file <- file.path(
  PATHS$processed,
  "07_tcga_deg_complete.rds"
)


deg_significant_file <- file.path(
  PATHS$processed,
  "07_tcga_deg_significant.rds"
)


reference_dir <- file.path(
  "external",
  "functional_enrichment"
)


go_reference_file <- file.path(
  reference_dir,
  "TCGA_BRCA_frozen_GO_BP_reference.rds"
)


kegg_reference_file <- file.path(
  reference_dir,
  "TCGA_BRCA_frozen_KEGG_reference.rds"
)


hallmark_reference_file <- file.path(
  reference_dir,
  "TCGA_BRCA_frozen_Hallmark_GSEA_reference.rds"
)


assert_file(
  deg_complete_file
)


assert_file(
  deg_significant_file
)


assert_file(
  go_reference_file
)


assert_file(
  kegg_reference_file
)


assert_file(
  hallmark_reference_file
)


############################################################
## 8.4 — OUTPUT DIRECTORY
############################################################

enrichment_output_dir <- file.path(
  "results",
  "tables",
  "functional_enrichment"
)


dir.create(
  enrichment_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
## 8.5 — LOAD MODULE-07 DEG RESULTS
############################################################

deg_complete <- readRDS(
  deg_complete_file
)


deg_significant <- readRDS(
  deg_significant_file
)


############################################################
## 8.6 — LOAD FROZEN ENRICHMENT REFERENCES
############################################################

S5 <- readRDS(
  go_reference_file
)


S6 <- readRDS(
  kegg_reference_file
)


S7 <- readRDS(
  hallmark_reference_file
)


############################################################
## 8.7 — MODULE-07 INPUT QC
############################################################

assert_true(
  nrow(
    deg_complete
  ) == 35009L,
  "Module-07 complete DEG table must contain 35,009 genes."
)


assert_true(
  nrow(
    deg_significant
  ) == 2830L,
  "Module-07 significant DEG table must contain 2,830 genes."
)


assert_unique(
  deg_complete$Gene,
  "Module-07 tested genes"
)


assert_unique(
  deg_significant$Gene,
  "Module-07 significant DEG genes"
)


assert_true(
  sum(
    deg_significant$Direction ==
      "Higher_in_Immune_High"
  ) == 2463L,
  "Expected 2,463 significant DEGs higher in Immune-High."
)


assert_true(
  sum(
    deg_significant$Direction ==
      "Higher_in_Immune_Low"
  ) == 367L,
  "Expected 367 significant DEGs higher in Immune-Low."
)


############################################################
## 8.8 — FROZEN INPUT REGRESSION
############################################################

significant_symbol_set_pass <- setequal(
  deg_significant$Gene,
  S5$significant_DEG_symbols
)


background_symbol_set_pass <- setequal(
  deg_complete$Gene,
  S5$supplied_background_symbols
)


assert_true(
  significant_symbol_set_pass,
  paste(
    "Module-07 significant DEG symbol set differs",
    "from frozen enrichment input."
  )
)


assert_true(
  background_symbol_set_pass,
  paste(
    "Module-07 tested-gene background differs",
    "from frozen enrichment input."
  )
)


############################################################
## 8.9 — FROZEN SYMBOL-to-ENTREZ MAPPING QC
############################################################

assert_true(
  length(
    S5$significant_DEG_symbols
  ) == 2830L,
  "Frozen significant DEG symbol count must be 2,830."
)


assert_true(
  nrow(
    S5$mapped_significant_DEG
  ) == 2818L,
  "Frozen mapped significant DEG count must be 2,818."
)


assert_true(
  length(
    S5$unmapped_significant_DEG
  ) == 12L,
  "Frozen unmapped significant DEG count must be 12."
)


assert_true(
  length(
    S5$supplied_background_symbols
  ) == 35009L,
  "Frozen tested-gene symbol background must contain 35,009 genes."
)


assert_true(
  length(
    S5$supplied_background_Entrez
  ) == 34627L,
  "Frozen supplied Entrez background must contain 34,627 IDs."
)


assert_unique(
  S5$mapped_significant_DEG$ENTREZID,
  "mapped significant DEG Entrez IDs"
)


assert_unique(
  S5$supplied_background_Entrez,
  "supplied background Entrez IDs"
)


assert_true(
  all(
    as.character(
      S5$mapped_significant_DEG$ENTREZID
    ) %in%
      as.character(
        S5$supplied_background_Entrez
      )
  ),
  "Significant DEG Entrez IDs are not all present in background."
)


############################################################
## 8.10 — GO-BP ORA
##
## Exact corrected historical method is recomputed.
############################################################

go_bp_object <- clusterProfiler::enrichGO(
  
  gene =
    as.character(
      S5$mapped_significant_DEG$ENTREZID
    ),
  
  universe =
    as.character(
      S5$supplied_background_Entrez
    ),
  
  OrgDb =
    org.Hs.eg.db::org.Hs.eg.db,
  
  keyType =
    "ENTREZID",
  
  ont =
    "BP",
  
  pAdjustMethod =
    "BH",
  
  pvalueCutoff =
    0.05,
  
  qvalueCutoff =
    0.05,
  
  readable =
    TRUE
)


go_bp_complete <- as.data.frame(
  go_bp_object
)


############################################################
## 8.11 — GO-BP CORE QC
############################################################

assert_true(
  nrow(
    go_bp_complete
  ) == 1145L,
  paste0(
    "Expected 1,145 corrected GO-BP terms; observed ",
    nrow(
      go_bp_complete
    ),
    "."
  )
)


assert_true(
  length(
    go_bp_object@universe
  ) == 17987L,
  paste0(
    "Expected effective GO-BP universe of 17,987; observed ",
    length(
      go_bp_object@universe
    ),
    "."
  )
)


assert_true(
  setequal(
    as.character(
      go_bp_object@universe
    ),
    as.character(
      S5$effective_GO_BP_universe
    )
  ),
  "Effective GO-BP universe differs from frozen corrected analysis."
)


assert_true(
  setequal(
    go_bp_complete$ID,
    S5$complete_GO_BP_results$ID
  ),
  "GO-BP significant term set differs from frozen analysis."
)


############################################################
## 8.12 — GO-BP EXACT NUMERICAL REGRESSION
############################################################

go_frozen_index <- match(
  go_bp_complete$ID,
  S5$complete_GO_BP_results$ID
)


assert_true(
  !anyNA(
    go_frozen_index
  ),
  "Unable to align current GO-BP terms to frozen result."
)


go_frozen_aligned <- S5$complete_GO_BP_results[
  go_frozen_index,
  ,
  drop = FALSE
]


go_numeric_columns <- c(
  "RichFactor",
  "FoldEnrichment",
  "zScore",
  "pvalue",
  "p.adjust",
  "qvalue",
  "Count"
)


go_regression <- data.frame(
  
  Metric =
    go_numeric_columns,
  
  Max_abs_difference =
    NA_real_,
  
  Mean_abs_difference =
    NA_real_,
  
  stringsAsFactors = FALSE
)


for (
  i in seq_along(
    go_numeric_columns
  )
) {
  
  metric <- go_numeric_columns[i]
  
  
  current_value <- as.numeric(
    go_bp_complete[[metric]]
  )
  
  
  frozen_value <- as.numeric(
    go_frozen_aligned[[metric]]
  )
  
  
  go_regression$Max_abs_difference[i] <-
    max(
      abs(
        current_value -
          frozen_value
      ),
      na.rm = TRUE
    )
  
  
  go_regression$Mean_abs_difference[i] <-
    mean(
      abs(
        current_value -
          frozen_value
      ),
      na.rm = TRUE
    )
}


go_gene_ratio_pass <- identical(
  
  as.character(
    go_bp_complete$GeneRatio
  ),
  
  as.character(
    go_frozen_aligned$GeneRatio
  )
)


go_bg_ratio_pass <- identical(
  
  as.character(
    go_bp_complete$BgRatio
  ),
  
  as.character(
    go_frozen_aligned$BgRatio
  )
)


GO_NUMERIC_TOLERANCE <- 1e-10


go_numerical_regression_pass <-
  max(
    go_regression$Max_abs_difference
  ) <
  GO_NUMERIC_TOLERANCE


assert_true(
  go_gene_ratio_pass,
  "GO-BP GeneRatio differs from frozen corrected result."
)


assert_true(
  go_bg_ratio_pass,
  "GO-BP BgRatio differs from frozen corrected result."
)


assert_true(
  go_numerical_regression_pass,
  paste0(
    "GO-BP numerical regression failed. Maximum difference = ",
    format(
      max(
        go_regression$Max_abs_difference
      ),
      scientific = TRUE,
      digits = 17
    ),
    "."
  )
)


############################################################
## 8.13 — GO-BP TOP-20
############################################################

go_bp_top20 <- S5$top20_GO_BP


assert_true(
  nrow(
    go_bp_top20
  ) == 20L,
  "Frozen GO-BP top-20 table must contain 20 terms."
)


assert_true(
  all(
    go_bp_top20$ID %in%
      go_bp_complete$ID
  ),
  "Frozen GO-BP top-20 contains terms absent from production result."
)


############################################################
## 8.14 — KEGG
##
## IMPORTANT:
##
## The validated corrected KEGG result is deliberately
## preserved from the frozen reference.
##
## It is NOT rerun against current online KEGG annotation.
############################################################

kegg_complete <- S6$complete_results


kegg_top20 <- S6$top20


kegg_qc <- S6$QC


kegg_historical_comparison <-
  S6$historical_comparison


kegg_changed_pathways <-
  S6$changed_pathways


############################################################
## 8.15 — KEGG FROZEN QC
############################################################

assert_true(
  identical(
    S6$analysis,
    "KEGG over-representation analysis"
  ),
  "Unexpected frozen KEGG analysis label."
)


assert_true(
  identical(
    S6$directionality,
    "Non-directional; combined significant DEG set"
  ),
  "Unexpected KEGG directionality."
)


assert_true(
  length(
    S6$significant_DEG_Entrez
  ) == 2818L,
  "Frozen KEGG significant Entrez input must contain 2,818 IDs."
)


assert_true(
  length(
    S6$supplied_background_Entrez
  ) == 34627L,
  "Frozen KEGG supplied background must contain 34,627 Entrez IDs."
)


assert_true(
  as.integer(
    S6$effective_KEGG_DEG_denominator
  ) == 1172L,
  "Frozen KEGG effective DEG denominator must be 1,172."
)


assert_true(
  length(
    S6$effective_KEGG_universe
  ) == 8637L,
  "Frozen corrected KEGG effective universe must contain 8,637 IDs."
)


assert_true(
  nrow(
    kegg_complete
  ) == 90L,
  "Frozen corrected KEGG result must contain 90 significant pathways."
)


assert_true(
  nrow(
    kegg_top20
  ) == 20L,
  "Frozen KEGG top-20 table must contain 20 pathways."
)


assert_true(
  setequal(
    as.character(
      S6$significant_DEG_Entrez
    ),
    as.character(
      S5$mapped_significant_DEG$ENTREZID
    )
  ),
  paste(
    "GO and KEGG must use the same frozen mapped",
    "significant DEG Entrez set."
  )
)


assert_true(
  setequal(
    as.character(
      S6$supplied_background_Entrez
    ),
    as.character(
      S5$supplied_background_Entrez
    )
  ),
  paste(
    "GO and KEGG must use the same frozen",
    "tested-gene Entrez background."
  )
)


############################################################
## 8.16 — HALLMARK RANK VECTOR FROM MODULE 07
############################################################

hallmark_rank <- deg_complete$logFC


names(
  hallmark_rank
) <- deg_complete$Gene


hallmark_rank <- sort(
  hallmark_rank,
  decreasing = TRUE
)


assert_true(
  length(
    hallmark_rank
  ) == 35009L,
  "Hallmark rank vector must contain 35,009 genes."
)


assert_unique(
  names(
    hallmark_rank
  ),
  "Hallmark ranked genes"
)


assert_true(
  sum(
    is.na(
      hallmark_rank
    )
  ) == 0L,
  "Missing values detected in Hallmark ranking vector."
)


assert_true(
  all(
    diff(
      hallmark_rank
    ) <= 0
  ),
  "Hallmark rank vector is not sorted in descending order."
)


############################################################
## 8.17 — EXACT FROZEN RANK REGRESSION
############################################################

hallmark_rank_gene_order_pass <- identical(
  
  names(
    hallmark_rank
  ),
  
  names(
    S7$ranked_gene_list
  )
)


hallmark_rank_max_difference <- max(
  abs(
    hallmark_rank -
      S7$ranked_gene_list
  )
)


hallmark_rank_value_pass <-
  hallmark_rank_max_difference <
  1e-12


assert_true(
  hallmark_rank_gene_order_pass,
  "Hallmark rank gene order differs from frozen GSEA rank."
)


assert_true(
  hallmark_rank_value_pass,
  paste0(
    "Hallmark rank values differ from frozen GSEA rank. Max difference = ",
    format(
      hallmark_rank_max_difference,
      scientific = TRUE,
      digits = 17
    ),
    "."
  )
)


############################################################
## 8.18 — FROZEN HALLMARK TERM2GENE
############################################################

hallmark_term2gene <- as.data.frame(
  S7$Hallmark_TERM2GENE
)


assert_true(
  nrow(
    hallmark_term2gene
  ) == 7322L,
  "Frozen Hallmark TERM2GENE must contain 7,322 pairs."
)


assert_true(
  ncol(
    hallmark_term2gene
  ) == 2L,
  "Frozen Hallmark TERM2GENE must contain two columns."
)


assert_true(
  length(
    unique(
      hallmark_term2gene[[1]]
    )
  ) == 50L,
  "Frozen Hallmark collection must contain 50 gene sets."
)


assert_true(
  sum(
    duplicated(
      hallmark_term2gene
    )
  ) == 0L,
  "Duplicated Hallmark pathway-gene pairs detected."
)


assert_true(
  identical(
    S7$MSigDB_version,
    "2026.1.Hs"
  ),
  "Frozen Hallmark MSigDB version is not 2026.1.Hs."
)


############################################################
## 8.19 — PRESERVE FROZEN HALLMARK GSEA RESULTS
##
## We DO NOT rerun GSEA here.
############################################################

hallmark_returned <-
  S7$returned_results


hallmark_significant <-
  S7$BH_significant_results


hallmark_leading_edge_summary <-
  S7$leading_edge_summary


hallmark_manifest <-
  S7$Hallmark_manifest


############################################################
## 8.20 — HALLMARK RESULT QC
############################################################

assert_true(
  nrow(
    hallmark_returned
  ) == 27L,
  "Frozen Hallmark returned-result table must contain 27 rows."
)


assert_true(
  nrow(
    hallmark_significant
  ) == 19L,
  "Frozen Hallmark significant-result table must contain 19 pathways."
)


assert_true(
  sum(
    hallmark_significant$p.adjust <
      0.05
  ) == 19L,
  "Expected exactly 19 Hallmark pathways with BH-FDR < 0.05."
)


assert_true(
  sum(
    hallmark_significant$NES >
      0
  ) == 19L,
  "All 19 significant Hallmark pathways must have positive NES."
)


assert_true(
  sum(
    hallmark_significant$NES <
      0
  ) == 0L,
  "No significant Hallmark pathway should have negative NES."
)


assert_true(
  all(
    hallmark_significant$Direction ==
      "Enriched_toward_Immune_High"
  ),
  paste(
    "All significant Hallmark pathways should be",
    "enriched toward Immune-High."
  )
)


assert_true(
  nrow(
    hallmark_leading_edge_summary
  ) == 19L,
  "Hallmark leading-edge summary must contain 19 pathways."
)


assert_true(
  nrow(
    hallmark_manifest
  ) == 50L,
  "Hallmark manifest must contain 50 gene sets."
)


############################################################
## 8.21 — DOCUMENT HISTORICAL GSEA PARAMETERS
############################################################

gsea_parameters <- data.frame(
  
  Parameter = c(
    "Comparison",
    "Ranking_metric",
    "Rank_direction",
    "Ranked_genes",
    "MSigDB_version",
    "Collection",
    "Hallmark_gene_sets",
    "TERM2GENE_pairs",
    "minGSSize",
    "maxGSSize",
    "pvalueCutoff",
    "pAdjustMethod",
    "eps",
    "nPermSimple",
    "seed",
    "exponent"
  ),
  
  Value = c(
    "Immune_High - Immune_Low",
    "limma logFC",
    "Descending",
    length(
      hallmark_rank
    ),
    S7$MSigDB_version,
    "H / Hallmark",
    length(
      unique(
        hallmark_term2gene[[1]]
      )
    ),
    nrow(
      hallmark_term2gene
    ),
    S7$GSEA_source_call_parameters$minGSSize,
    S7$GSEA_source_call_parameters$maxGSSize,
    S7$GSEA_source_call_parameters$pvalueCutoff,
    S7$GSEA_source_call_parameters$pAdjustMethod,
    S7$GSEA_source_call_parameters$eps,
    S7$GSEA_source_call_parameters$nPermSimple,
    S7$GSEA_source_call_parameters$seed,
    S7$stored_gseaResult_params$exponent
  ),
  
  stringsAsFactors = FALSE
)


############################################################
## 8.22 — FIVE SELECTED HALLMARKS FOR PPI
############################################################

ppi_selected_hallmarks <-
  S7$PPI_selected_Hallmarks


assert_true(
  nrow(
    ppi_selected_hallmarks
  ) == 5L,
  "Exactly five Hallmark pathways must be selected for PPI."
)


expected_ppi_hallmarks <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)


assert_true(
  setequal(
    ppi_selected_hallmarks$ID,
    expected_ppi_hallmarks
  ),
  "Five PPI-selected Hallmark pathways differ from frozen analysis."
)


############################################################
## 8.23 — EXPECTED FIVE-HALLMARK LEADING-EDGE SIZES
############################################################

expected_leading_edge_counts <- c(
  
  HALLMARK_INTERFERON_GAMMA_RESPONSE =
    151L,
  
  HALLMARK_ALLOGRAFT_REJECTION =
    131L,
  
  HALLMARK_INFLAMMATORY_RESPONSE =
    145L,
  
  HALLMARK_IL6_JAK_STAT3_SIGNALING =
    54L,
  
  HALLMARK_TNFA_SIGNALING_VIA_NFKB =
    141L
)


for (
  pathway in names(
    expected_leading_edge_counts
  )
) {
  
  observed <- ppi_selected_hallmarks$LeadingEdge_N[
    ppi_selected_hallmarks$ID ==
      pathway
  ]
  
  
  expected <-
    expected_leading_edge_counts[
      pathway
    ]
  
  
  assert_true(
    length(
      observed
    ) == 1L &&
      as.integer(
        observed
      ) ==
      as.integer(
        expected
      ),
    paste0(
      "Unexpected leading-edge size for ",
      pathway,
      "."
    )
  )
}


############################################################
## 8.24 — RECONSTRUCT FIVE-HALLMARK MEMBERSHIP
##
## Uses the frozen validated core_enrichment strings.
############################################################

selected_gsea_rows <- hallmark_returned[
  hallmark_returned$ID %in%
    expected_ppi_hallmarks,
  ,
  drop = FALSE
]


assert_true(
  nrow(
    selected_gsea_rows
  ) == 5L,
  paste(
    "Unable to identify all five selected Hallmarks",
    "in frozen GSEA result."
  )
)


membership_list <- lapply(
  
  seq_len(
    nrow(
      selected_gsea_rows
    )
  ),
  
  function(i) {
    
    pathway <-
      selected_gsea_rows$ID[i]
    
    
    genes <- strsplit(
      
      as.character(
        selected_gsea_rows$core_enrichment[i]
      ),
      
      "/",
      
      fixed = TRUE
      
    )[[1]]
    
    
    data.frame(
      
      Pathway =
        pathway,
      
      Gene =
        genes,
      
      stringsAsFactors = FALSE
    )
  }
)


ppi_membership_reconstructed <- do.call(
  rbind,
  membership_list
)


rownames(
  ppi_membership_reconstructed
) <- NULL


############################################################
## 8.25 — MEMBERSHIP REGRESSION
############################################################

assert_true(
  nrow(
    ppi_membership_reconstructed
  ) == 622L,
  paste0(
    "Expected 622 five-Hallmark pathway-gene memberships; observed ",
    nrow(
      ppi_membership_reconstructed
    ),
    "."
  )
)


membership_current_key <- paste(
  ppi_membership_reconstructed$Pathway,
  ppi_membership_reconstructed$Gene,
  sep = "|||"
)


membership_frozen_key <- paste(
  S7$PPI_leading_edge_membership$Pathway,
  S7$PPI_leading_edge_membership$Gene,
  sep = "|||"
)


ppi_membership_pass <- setequal(
  membership_current_key,
  membership_frozen_key
)


assert_true(
  ppi_membership_pass,
  "Reconstructed PPI leading-edge membership differs from frozen analysis."
)


############################################################
## 8.26 — PER-PATHWAY LEADING-EDGE COUNT REGRESSION
############################################################

ppi_pathway_counts <- table(
  ppi_membership_reconstructed$Pathway
)


for (
  pathway in names(
    expected_leading_edge_counts
  )
) {
  
  observed <- as.integer(
    ppi_pathway_counts[
      pathway
    ]
  )
  
  
  expected <- as.integer(
    expected_leading_edge_counts[
      pathway
    ]
  )
  
  
  assert_true(
    observed == expected,
    paste0(
      "Leading-edge membership count mismatch for ",
      pathway,
      "."
    )
  )
}


############################################################
## 8.27 — RECONSTRUCT 455-GENE UNION
############################################################

ppi_union_genes <- unique(
  ppi_membership_reconstructed$Gene
)


assert_true(
  length(
    ppi_union_genes
  ) == 455L,
  paste0(
    "Expected 455 unique PPI leading-edge genes; observed ",
    length(
      ppi_union_genes
    ),
    "."
  )
)


ppi_union_gene_set_pass <- setequal(
  ppi_union_genes,
  S7$PPI_leading_edge_union$Gene
)


assert_true(
  ppi_union_gene_set_pass,
  "Reconstructed 455-gene PPI union differs from frozen analysis."
)


############################################################
## 8.28 — RECONSTRUCT UNION ANNOTATION
############################################################

ppi_union_reconstructed <- lapply(
  
  sort(
    ppi_union_genes
  ),
  
  function(gene) {
    
    pathways <- unique(
      ppi_membership_reconstructed$Pathway[
        ppi_membership_reconstructed$Gene ==
          gene
      ]
    )
    
    
    pathways <- sort(
      pathways
    )
    
    
    data.frame(
      
      Gene =
        gene,
      
      N_Selected_Hallmarks =
        length(
          pathways
        ),
      
      Selected_Hallmarks =
        paste(
          pathways,
          collapse = ";"
        ),
      
      stringsAsFactors = FALSE
    )
  }
)


ppi_union_reconstructed <- do.call(
  rbind,
  ppi_union_reconstructed
)


rownames(
  ppi_union_reconstructed
) <- NULL


############################################################
## 8.29 — VALIDATE UNION COUNTS AGAINST FROZEN TABLE
############################################################

frozen_union_index <- match(
  ppi_union_reconstructed$Gene,
  S7$PPI_leading_edge_union$Gene
)


assert_true(
  !anyNA(
    frozen_union_index
  ),
  "Unable to align reconstructed PPI union to frozen table."
)


frozen_union_aligned <- S7$PPI_leading_edge_union[
  frozen_union_index,
  ,
  drop = FALSE
]


ppi_union_count_pass <- identical(
  
  as.integer(
    ppi_union_reconstructed$N_Selected_Hallmarks
  ),
  
  as.integer(
    frozen_union_aligned$N_Selected_Hallmarks
  )
)


assert_true(
  ppi_union_count_pass,
  "PPI union Hallmark-membership counts differ from frozen analysis."
)


############################################################
## 8.30 — PRESERVE AUTHORITATIVE FROZEN PPI UNION
##
## The frozen table preserves the original pathway text
## ordering as used in the validated supplement.
############################################################

ppi_leading_edge_membership <-
  S7$PPI_leading_edge_membership


ppi_leading_edge_union <-
  S7$PPI_leading_edge_union


############################################################
## 8.31 — ANALYSIS-PROVENANCE TABLE
############################################################

analysis_provenance <- data.frame(
  
  Analysis = c(
    "GO_BP_ORA",
    "KEGG_ORA",
    "Hallmark_GSEA",
    "PPI_leading_edge_input"
  ),
  
  Production_behavior = c(
    "Recomputed",
    "Frozen corrected result preserved",
    "Frozen validated result preserved",
    "Reconstructed and frozen-result verified"
  ),
  
  Reason = c(
    paste(
      "Exact corrected GO-BP computation is deterministic",
      "and reproduces frozen result exactly."
    ),
    
    paste(
      "Historical corrected KEGG used live online annotation;",
      "current KEGG annotation has drifted, so frozen corrected",
      "result is preserved."
    ),
    
    paste(
      "Exact Module-07 rank and frozen MSigDB TERM2GENE are",
      "verified; validated historical GSEA result is preserved",
      "rather than replaced by a backend/stochastic rerun."
    ),
    
    paste(
      "Five selected Hallmark leading edges are reconstructed",
      "from frozen validated GSEA core-enrichment genes."
    )
  ),
  
  stringsAsFactors = FALSE
)


############################################################
## 8.32 — GO SUMMARY
############################################################

go_summary <- data.frame(
  
  Metric = c(
    "Significant DEG symbols supplied",
    "Mapped significant DEG Entrez IDs",
    "Unmapped significant DEG symbols",
    "Tested/background symbols supplied",
    "Supplied background Entrez IDs",
    "Effective GO-BP DEG denominator",
    "Effective GO-BP universe",
    "Significant GO-BP terms",
    "Ontology",
    "Directionality",
    "Multiple-testing method",
    "P-value cutoff",
    "Q-value cutoff",
    "Historical numerical regression max difference"
  ),
  
  Value = c(
    2830,
    2818,
    12,
    35009,
    34627,
    2292,
    17987,
    1145,
    "Biological Process",
    "Combined significant DEG set; non-directional ORA",
    "Benjamini-Hochberg",
    0.05,
    0.05,
    format(
      max(
        go_regression$Max_abs_difference
      ),
      scientific = TRUE,
      digits = 17
    )
  ),
  
  stringsAsFactors = FALSE
)


############################################################
## 8.33 — KEGG SUMMARY
############################################################

kegg_summary <- data.frame(
  
  Metric = c(
    "Significant DEG Entrez IDs supplied",
    "Effective KEGG DEG denominator",
    "Tested/background Entrez IDs supplied",
    "Effective frozen KEGG universe",
    "Significant corrected KEGG pathways",
    "Organism",
    "Key type",
    "Multiple-testing method",
    "P-value cutoff",
    "Q-value cutoff",
    "Directionality",
    "Annotation mode",
    "Production rerun"
  ),
  
  Value = c(
    2818,
    1172,
    34627,
    8637,
    90,
    "hsa",
    "ncbi-geneid",
    "Benjamini-Hochberg",
    0.05,
    0.05,
    "Combined significant DEG set; non-directional ORA",
    "Historical corrected analysis used online KEGG annotation",
    "FALSE"
  ),
  
  stringsAsFactors = FALSE
)


############################################################
## 8.34 — HALLMARK SUMMARY
############################################################

hallmark_summary <- data.frame(
  
  Metric = c(
    "Comparison",
    "Ranking metric",
    "Ranked genes",
    "Rank direction",
    "Frozen MSigDB version",
    "Hallmark gene sets",
    "TERM2GENE pairs",
    "Returned GSEA rows",
    "BH-FDR significant pathways",
    "Significant positive NES pathways",
    "Significant negative NES pathways",
    "PPI-selected Hallmarks",
    "PPI pathway-gene memberships",
    "Unique PPI leading-edge genes",
    "Production GSEA rerun"
  ),
  
  Value = c(
    "Immune_High - Immune_Low",
    "limma logFC",
    35009,
    "Descending",
    "2026.1.Hs",
    50,
    7322,
    27,
    19,
    19,
    0,
    5,
    622,
    455,
    "FALSE"
  ),
  
  stringsAsFactors = FALSE
)


############################################################
## 8.35 — SAVE GO-BP OBJECTS
############################################################

save_rds_checked(
  go_bp_object,
  file.path(
    PATHS$processed,
    "08_go_bp_enrichResult.rds"
  )
)


save_rds_checked(
  go_bp_complete,
  file.path(
    PATHS$processed,
    "08_go_bp_complete.rds"
  )
)


write_csv_checked(
  go_bp_complete,
  file.path(
    enrichment_output_dir,
    "08_GO_BP_complete.csv"
  )
)


write_csv_checked(
  go_bp_top20,
  file.path(
    enrichment_output_dir,
    "08_GO_BP_top20.csv"
  )
)


write_csv_checked(
  go_summary,
  file.path(
    enrichment_output_dir,
    "08_GO_BP_summary.csv"
  )
)


############################################################
## 8.36 — SAVE KEGG FROZEN PRODUCTION OBJECTS
############################################################

save_rds_checked(
  kegg_complete,
  file.path(
    PATHS$processed,
    "08_kegg_complete_frozen_corrected.rds"
  )
)


write_csv_checked(
  kegg_complete,
  file.path(
    enrichment_output_dir,
    "08_KEGG_complete_corrected.csv"
  )
)


write_csv_checked(
  kegg_top20,
  file.path(
    enrichment_output_dir,
    "08_KEGG_top20_corrected.csv"
  )
)


write_csv_checked(
  kegg_summary,
  file.path(
    enrichment_output_dir,
    "08_KEGG_summary.csv"
  )
)


write_csv_checked(
  kegg_changed_pathways,
  file.path(
    enrichment_output_dir,
    "08_KEGG_historical_vs_corrected_changed_pathways.csv"
  )
)


############################################################
## 8.37 — SAVE HALLMARK RANK + TERM2GENE
############################################################

save_rds_checked(
  hallmark_rank,
  file.path(
    PATHS$processed,
    "08_hallmark_ranked_gene_list.rds"
  )
)


save_rds_checked(
  hallmark_term2gene,
  file.path(
    PATHS$processed,
    "08_hallmark_TERM2GENE_2026_1_Hs.rds"
  )
)


write_csv_checked(
  hallmark_term2gene,
  file.path(
    enrichment_output_dir,
    "08_Hallmark_TERM2GENE_2026_1_Hs.csv"
  )
)


############################################################
## 8.38 — SAVE FROZEN HALLMARK GSEA RESULTS
############################################################

save_rds_checked(
  hallmark_returned,
  file.path(
    PATHS$processed,
    "08_hallmark_gsea_returned_frozen.rds"
  )
)


save_rds_checked(
  hallmark_significant,
  file.path(
    PATHS$processed,
    "08_hallmark_gsea_significant_frozen.rds"
  )
)


write_csv_checked(
  hallmark_returned,
  file.path(
    enrichment_output_dir,
    "08_Hallmark_GSEA_returned.csv"
  )
)


write_csv_checked(
  hallmark_significant,
  file.path(
    enrichment_output_dir,
    "08_Hallmark_GSEA_BH_significant.csv"
  )
)


write_csv_checked(
  hallmark_leading_edge_summary,
  file.path(
    enrichment_output_dir,
    "08_Hallmark_GSEA_leading_edge_summary.csv"
  )
)


write_csv_checked(
  hallmark_manifest,
  file.path(
    enrichment_output_dir,
    "08_Hallmark_gene_set_manifest.csv"
  )
)


write_csv_checked(
  gsea_parameters,
  file.path(
    enrichment_output_dir,
    "08_Hallmark_GSEA_parameters.csv"
  )
)


write_csv_checked(
  hallmark_summary,
  file.path(
    enrichment_output_dir,
    "08_Hallmark_GSEA_summary.csv"
  )
)


############################################################
## 8.39 — SAVE PPI INPUTS
############################################################

save_rds_checked(
  ppi_selected_hallmarks,
  file.path(
    PATHS$processed,
    "08_ppi_selected_hallmarks.rds"
  )
)


save_rds_checked(
  ppi_leading_edge_membership,
  file.path(
    PATHS$processed,
    "08_ppi_leading_edge_membership.rds"
  )
)


save_rds_checked(
  ppi_leading_edge_union,
  file.path(
    PATHS$processed,
    "08_ppi_leading_edge_union_455.rds"
  )
)


write_csv_checked(
  ppi_selected_hallmarks,
  file.path(
    enrichment_output_dir,
    "08_PPI_selected_Hallmarks.csv"
  )
)


write_csv_checked(
  ppi_leading_edge_membership,
  file.path(
    enrichment_output_dir,
    "08_PPI_leading_edge_membership.csv"
  )
)


write_csv_checked(
  ppi_leading_edge_union,
  file.path(
    enrichment_output_dir,
    "08_PPI_leading_edge_union_455.csv"
  )
)


############################################################
## 8.40 — SAVE ANALYSIS PROVENANCE
############################################################

write_csv_checked(
  analysis_provenance,
  file.path(
    enrichment_output_dir,
    "08_functional_enrichment_provenance.csv"
  )
)


############################################################
## 8.41 — SAVE GO REGRESSION LOG
############################################################

write_csv_checked(
  go_regression,
  file.path(
    PATHS$logs,
    "08_GO_BP_frozen_regression.csv"
  )
)


############################################################
## 8.42 — REFERENCE HASHES
############################################################

reference_hashes <- tools::md5sum(
  c(
    go_reference_file,
    kegg_reference_file,
    hallmark_reference_file
  )
)


reference_hash_table <- data.frame(
  
  File =
    names(
      reference_hashes
    ),
  
  MD5 =
    unname(
      reference_hashes
    ),
  
  stringsAsFactors = FALSE
)


write_csv_checked(
  reference_hash_table,
  file.path(
    PATHS$logs,
    "08_functional_enrichment_reference_MD5.csv"
  )
)


############################################################
## 8.43 — MODULE QC
############################################################

module08_qc <- data.frame(
  
  Metric = c(
    "Module07 tested genes",
    "Module07 significant DEGs",
    "Significant DEG input regression",
    "Tested-gene background regression",
    "Mapped significant DEG Entrez IDs",
    "Mapped background Entrez IDs",
    "GO-BP effective universe",
    "GO-BP significant terms",
    "GO-BP numerical max difference",
    "GO-BP numerical regression",
    "KEGG effective frozen universe",
    "KEGG effective DEG denominator",
    "KEGG significant corrected pathways",
    "KEGG production mode",
    "Hallmark rank genes",
    "Hallmark rank maximum difference",
    "Hallmark rank regression",
    "Hallmark MSigDB version",
    "Hallmark gene sets",
    "Hallmark TERM2GENE pairs",
    "Hallmark returned rows",
    "Hallmark BH-significant pathways",
    "Hallmark positive significant NES",
    "Hallmark negative significant NES",
    "Hallmark production mode",
    "PPI-selected Hallmarks",
    "PPI pathway-gene memberships",
    "PPI unique leading-edge union genes",
    "PPI membership regression",
    "PPI union regression",
    "clusterProfiler version",
    "org.Hs.eg.db version"
  ),
  
  Value = c(
    nrow(
      deg_complete
    ),
    
    nrow(
      deg_significant
    ),
    
    significant_symbol_set_pass,
    
    background_symbol_set_pass,
    
    nrow(
      S5$mapped_significant_DEG
    ),
    
    length(
      S5$supplied_background_Entrez
    ),
    
    length(
      go_bp_object@universe
    ),
    
    nrow(
      go_bp_complete
    ),
    
    format(
      max(
        go_regression$Max_abs_difference
      ),
      scientific = TRUE,
      digits = 17
    ),
    
    go_numerical_regression_pass,
    
    length(
      S6$effective_KEGG_universe
    ),
    
    as.integer(
      S6$effective_KEGG_DEG_denominator
    ),
    
    nrow(
      kegg_complete
    ),
    
    "Frozen corrected result; no live KEGG rerun",
    
    length(
      hallmark_rank
    ),
    
    format(
      hallmark_rank_max_difference,
      scientific = TRUE,
      digits = 17
    ),
    
    all(
      hallmark_rank_gene_order_pass,
      hallmark_rank_value_pass
    ),
    
    S7$MSigDB_version,
    
    length(
      unique(
        hallmark_term2gene[[1]]
      )
    ),
    
    nrow(
      hallmark_term2gene
    ),
    
    nrow(
      hallmark_returned
    ),
    
    nrow(
      hallmark_significant
    ),
    
    sum(
      hallmark_significant$NES >
        0
    ),
    
    sum(
      hallmark_significant$NES <
        0
    ),
    
    "Frozen validated result; no production GSEA rerun",
    
    nrow(
      ppi_selected_hallmarks
    ),
    
    nrow(
      ppi_leading_edge_membership
    ),
    
    nrow(
      ppi_leading_edge_union
    ),
    
    ppi_membership_pass,
    
    ppi_union_gene_set_pass,
    
    as.character(
      packageVersion(
        "clusterProfiler"
      )
    ),
    
    as.character(
      packageVersion(
        "org.Hs.eg.db"
      )
    )
  ),
  
  stringsAsFactors = FALSE
)


write_csv_checked(
  module08_qc,
  file.path(
    PATHS$logs,
    "08_functional_enrichment_QC.csv"
  )
)


############################################################
## 8.44 — SESSION INFORMATION
############################################################

capture.output(
  sessionInfo(),
  file = file.path(
    PATHS$logs,
    "08_functional_enrichment_sessionInfo.txt"
  )
)


############################################################
## 8.45 — FINAL STATUS
############################################################

cat(
  "\n=============================================\n"
)


cat(
  "08_functional_enrichment.R: PASS\n"
)


cat(
  "\nGO-BP ORA\n"
)


cat(
  "Significant DEG symbols: ",
  length(
    S5$significant_DEG_symbols
  ),
  "\n",
  sep = ""
)


cat(
  "Mapped significant DEG Entrez IDs: ",
  nrow(
    S5$mapped_significant_DEG
  ),
  "\n",
  sep = ""
)


cat(
  "Effective GO-BP universe: ",
  length(
    go_bp_object@universe
  ),
  "\n",
  sep = ""
)


cat(
  "Significant GO-BP terms: ",
  nrow(
    go_bp_complete
  ),
  "\n",
  sep = ""
)


cat(
  "GO-BP frozen max difference: ",
  format(
    max(
      go_regression$Max_abs_difference
    ),
    scientific = TRUE,
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "GO-BP historical regression: VERIFIED\n"
)


cat(
  "\nKEGG ORA\n"
)


cat(
  "Effective frozen KEGG universe: ",
  length(
    S6$effective_KEGG_universe
  ),
  "\n",
  sep = ""
)


cat(
  "Effective KEGG DEG denominator: ",
  as.integer(
    S6$effective_KEGG_DEG_denominator
  ),
  "\n",
  sep = ""
)


cat(
  "Significant corrected KEGG pathways: ",
  nrow(
    kegg_complete
  ),
  "\n",
  sep = ""
)


cat(
  "KEGG production result: FROZEN CORRECTED RESULT PRESERVED\n"
)


cat(
  "Live KEGG database queried in production: FALSE\n"
)


cat(
  "\nHallmark GSEA\n"
)


cat(
  "Ranked genes: ",
  length(
    hallmark_rank
  ),
  "\n",
  sep = ""
)


cat(
  "Rank maximum difference: ",
  format(
    hallmark_rank_max_difference,
    scientific = TRUE,
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "Frozen MSigDB version: ",
  S7$MSigDB_version,
  "\n",
  sep = ""
)


cat(
  "Hallmark gene sets: ",
  length(
    unique(
      hallmark_term2gene[[1]]
    )
  ),
  "\n",
  sep = ""
)


cat(
  "TERM2GENE pairs: ",
  nrow(
    hallmark_term2gene
  ),
  "\n",
  sep = ""
)


cat(
  "BH-significant Hallmarks: ",
  nrow(
    hallmark_significant
  ),
  "\n",
  sep = ""
)


cat(
  "Positive significant NES: ",
  sum(
    hallmark_significant$NES >
      0
  ),
  "\n",
  sep = ""
)


cat(
  "Negative significant NES: ",
  sum(
    hallmark_significant$NES <
      0
  ),
  "\n",
  sep = ""
)


cat(
  "Hallmark rank regression: VERIFIED\n"
)


cat(
  "Production GSEA rerun: FALSE\n"
)


cat(
  "\nPPI input\n"
)


cat(
  "Selected Hallmarks: ",
  nrow(
    ppi_selected_hallmarks
  ),
  "\n",
  sep = ""
)


cat(
  "Leading-edge memberships: ",
  nrow(
    ppi_leading_edge_membership
  ),
  "\n",
  sep = ""
)


cat(
  "Unique leading-edge union genes: ",
  nrow(
    ppi_leading_edge_union
  ),
  "\n",
  sep = ""
)


cat(
  "PPI leading-edge reconstruction: VERIFIED\n"
)


cat(
  "\nGO/KEGG ORA directionality: NON-DIRECTIONAL\n"
)


cat(
  "Hallmark GSEA direction: HIGH-minus-LOW logFC\n"
)


cat(
  "Upstream DEG analysis recomputed: FALSE\n"
)


cat(
  "Immune phenotype recomputed: FALSE\n"
)


cat(
  "=============================================\n"
)