############################################################
## 09_ppi_hubs.R
##
## PURPOSE
##
## Production/refactoring module for the validated
## STRING -> Cytoscape -> CytoHubba PPI analysis.
##
## IMPORTANT
##
## This module DOES NOT:
##
##   - query STRING again
##   - rebuild the PPI network
##   - rerun Cytoscape
##   - rerun CytoHubba
##   - recalculate MCC scores
##   - change the historical network threshold
##
## The validated frozen network artifact is preserved and
## checked against the exact 455-gene input produced by
## Module 08.
##
##
## VALIDATED HISTORICAL WORKFLOW
##
## Input:
##   Union of leading-edge genes from five selected
##   Hallmark GSEA pathways
##
## Input genes:
##   455
##
## STRING representation:
##   Exact symbol matches              388
##   Alias remapping                     1
##       RIGI -> DDX58
##   Represented after reconciliation  389
##   Absent from exported network       66
##
## Preserved network:
##   Nodes                              389
##   Original reciprocal rows        8,088
##   Unique undirected edges         4,044
##   combined_score range        0.700-0.999
##
## Hub prioritization:
##   Cytoscape / CytoHubba
##   MCC = Maximal Clique Centrality
##   Top 20 genes retained
##
##
## INTERPRETATION GUARDRAIL
##
## Hub genes are highly connected immune-associated genes
## within this curated immune-focused PPI network.
##
## They must NOT be described as causal master regulators.
############################################################


############################################################
## 9.1 — REQUIRE PROJECT SETUP
############################################################

if (!exists("PATHS")) {
  source("R/00_setup.R")
}


############################################################
## 9.2 — INPUT FILES
############################################################

ppi_union_file <- file.path(
  PATHS$processed,
  "08_ppi_leading_edge_union_455.rds"
)


deg_file <- file.path(
  PATHS$processed,
  "07_tcga_deg_complete.rds"
)


ppi_reference_file <- file.path(
  "external",
  "ppi",
  "TCGA_BRCA_frozen_STRING_Cytoscape_CytoHubba_reference.rds"
)


assert_file(
  ppi_union_file
)


assert_file(
  deg_file
)


assert_file(
  ppi_reference_file
)


############################################################
## 9.3 — OUTPUT DIRECTORY
############################################################

ppi_output_dir <- file.path(
  "results",
  "tables",
  "ppi_hubs"
)


dir.create(
  ppi_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
## 9.4 — LOAD INPUTS
############################################################

ppi_union <- readRDS(
  ppi_union_file
)


deg07 <- readRDS(
  deg_file
)


S8 <- readRDS(
  ppi_reference_file
)


############################################################
## 9.5 — MODULE-08 PPI INPUT QC
############################################################

assert_true(
  nrow(
    ppi_union
  ) == 455L,
  "Module-08 PPI union must contain 455 rows."
)


assert_true(
  "Gene" %in%
    colnames(
      ppi_union
    ),
  "Module-08 PPI union must contain a Gene column."
)


assert_unique(
  ppi_union$Gene,
  "Module-08 PPI union genes"
)


assert_true(
  length(
    unique(
      ppi_union$Gene
    )
  ) == 455L,
  "Module-08 PPI union must contain 455 unique genes."
)


############################################################
## 9.6 — FROZEN S8 STRUCTURE QC
############################################################

required_S8_objects <- c(
  "input_genes",
  "node_mapping",
  "network_nodes",
  "undirected_edges",
  "top20_MCC_hubs",
  "network_QC",
  "artifact_provenance",
  "method_provenance",
  "identifier_reconciliation",
  "CytoHubba_algorithm",
  "interpretation"
)


assert_true(
  all(
    required_S8_objects %in%
      names(
        S8
      )
  ),
  "Frozen S8 reference is missing required objects."
)


############################################################
## 9.7 — 455-GENE INPUT REGRESSION
############################################################

assert_true(
  length(
    S8$input_genes
  ) == 455L,
  "Frozen S8 input must contain 455 genes."
)


assert_unique(
  S8$input_genes,
  "Frozen S8 PPI input genes"
)


ppi_input_pass <- setequal(
  ppi_union$Gene,
  S8$input_genes
)


assert_true(
  ppi_input_pass,
  paste(
    "Module-08 455-gene PPI input differs",
    "from the frozen STRING/Cytoscape analysis."
  )
)


############################################################
## 9.8 — NODE-MAPPING STRUCTURE
############################################################

required_mapping_columns <- c(
  "Input_Gene",
  "Network_Label",
  "Mapping_Status",
  "Represented_in_network"
)


assert_true(
  all(
    required_mapping_columns %in%
      colnames(
        S8$node_mapping
      )
  ),
  "Frozen PPI node-mapping table is missing required columns."
)


assert_true(
  nrow(
    S8$node_mapping
  ) == 455L,
  "PPI node-mapping table must contain 455 input genes."
)


assert_unique(
  S8$node_mapping$Input_Gene,
  "PPI node-mapping input genes"
)


assert_true(
  setequal(
    S8$node_mapping$Input_Gene,
    S8$input_genes
  ),
  "Node-mapping input genes differ from frozen PPI input."
)


############################################################
## 9.9 — NODE-MAPPING STATUS QC
############################################################

mapping_status_counts <- table(
  S8$node_mapping$Mapping_Status
)


assert_true(
  "Exact_symbol_match" %in%
    names(
      mapping_status_counts
    ),
  "Exact_symbol_match category absent from node mapping."
)


assert_true(
  "Alias_remapped_RIGI_to_DDX58" %in%
    names(
      mapping_status_counts
    ),
  "Expected RIGI-to-DDX58 alias category absent."
)


assert_true(
  "Absent_from_exported_network" %in%
    names(
      mapping_status_counts
    ),
  "Expected absent-network category missing."
)


assert_true(
  as.integer(
    mapping_status_counts[
      "Exact_symbol_match"
    ]
  ) == 388L,
  "Expected 388 exact input/network symbol matches."
)


assert_true(
  as.integer(
    mapping_status_counts[
      "Alias_remapped_RIGI_to_DDX58"
    ]
  ) == 1L,
  "Expected exactly one alias reconciliation."
)


assert_true(
  as.integer(
    mapping_status_counts[
      "Absent_from_exported_network"
    ]
  ) == 66L,
  "Expected 66 input genes absent from exported network."
)


############################################################
## 9.10 — REPRESENTED NODE QC
############################################################

represented_count <- sum(
  S8$node_mapping$Represented_in_network
)


absent_count <- sum(
  !S8$node_mapping$Represented_in_network
)


assert_true(
  represented_count == 389L,
  "Expected 389 PPI input genes represented after reconciliation."
)


assert_true(
  absent_count == 66L,
  "Expected 66 genes absent from the preserved network."
)


############################################################
## 9.11 — IDENTIFIER RECONCILIATION
############################################################

assert_true(
  length(
    S8$identifier_reconciliation
  ) == 1L,
  "Expected exactly one identifier reconciliation."
)


assert_true(
  identical(
    names(
      S8$identifier_reconciliation
    ),
    "RIGI"
  ),
  "Expected reconciled Hallmark identifier RIGI."
)


assert_true(
  identical(
    unname(
      S8$identifier_reconciliation
    ),
    "DDX58"
  ),
  "Expected RIGI to map to DDX58."
)


alias_row <- S8$node_mapping[
  S8$node_mapping$Input_Gene ==
    "RIGI",
  ,
  drop = FALSE
]


assert_true(
  nrow(
    alias_row
  ) == 1L,
  "Expected exactly one RIGI node-mapping row."
)


assert_true(
  identical(
    as.character(
      alias_row$Network_Label
    ),
    "DDX58"
  ),
  "RIGI node mapping does not point to DDX58."
)


assert_true(
  isTRUE(
    alias_row$Represented_in_network
  ),
  "RIGI/DDX58 reconciled node must be represented in network."
)


############################################################
## 9.12 — NETWORK NODE QC
############################################################

assert_true(
  length(
    S8$network_nodes
  ) == 389L,
  "Preserved STRING network must contain 389 nodes."
)


assert_unique(
  S8$network_nodes,
  "STRING network nodes"
)


represented_labels <- unique(
  S8$node_mapping$Network_Label[
    S8$node_mapping$Represented_in_network
  ]
)


assert_true(
  setequal(
    represented_labels,
    S8$network_nodes
  ),
  paste(
    "Represented node-mapping labels differ",
    "from preserved STRING network nodes."
  )
)


############################################################
## 9.13 — NETWORK EDGE STRUCTURE QC
############################################################

required_edge_columns <- c(
  "node1",
  "node2",
  "node1_string_id",
  "node2_string_id",
  "combined_score"
)


assert_true(
  all(
    required_edge_columns %in%
      colnames(
        S8$undirected_edges
      )
  ),
  "Frozen STRING edge table is missing required columns."
)


assert_true(
  nrow(
    S8$undirected_edges
  ) == 4044L,
  "Preserved network must contain 4,044 unique undirected edges."
)


assert_true(
  all(
    S8$undirected_edges$node1 %in%
      S8$network_nodes
  ),
  "Some edge node1 labels are absent from network node set."
)


assert_true(
  all(
    S8$undirected_edges$node2 %in%
      S8$network_nodes
  ),
  "Some edge node2 labels are absent from network node set."
)


assert_true(
  !any(
    S8$undirected_edges$node1 ==
      S8$undirected_edges$node2
  ),
  "Self-edges detected in preserved undirected network."
)


############################################################
## 9.14 — VERIFY UNDIRECTED-EDGE UNIQUENESS
############################################################

edge_key <- vapply(
  
  seq_len(
    nrow(
      S8$undirected_edges
    )
  ),
  
  function(i) {
    
    paste(
      sort(
        c(
          S8$undirected_edges$node1[i],
          S8$undirected_edges$node2[i]
        )
      ),
      collapse = "|||"
    )
  },
  
  character(1)
)


assert_true(
  !anyDuplicated(
    edge_key
  ),
  "Duplicated undirected STRING interactions remain."
)


############################################################
## 9.15 — STRING SCORE QC
############################################################

minimum_score <- min(
  S8$undirected_edges$combined_score,
  na.rm = TRUE
)


maximum_score <- max(
  S8$undirected_edges$combined_score,
  na.rm = TRUE
)


edges_below_0700 <- sum(
  S8$undirected_edges$combined_score <
    0.700,
  na.rm = TRUE
)


edges_exact_0700 <- sum(
  S8$undirected_edges$combined_score ==
    0.700,
  na.rm = TRUE
)


assert_true(
  isTRUE(
    all.equal(
      minimum_score,
      0.700,
      tolerance = 1e-12
    )
  ),
  "Expected minimum preserved STRING combined_score = 0.700."
)


assert_true(
  isTRUE(
    all.equal(
      maximum_score,
      0.999,
      tolerance = 1e-12
    )
  ),
  "Expected maximum preserved STRING combined_score = 0.999."
)


assert_true(
  edges_below_0700 == 0L,
  "Preserved STRING network contains scores below 0.700."
)


assert_true(
  edges_exact_0700 == 22L,
  "Expected 22 preserved edges with combined_score exactly 0.700."
)


############################################################
## 9.16 — CYTOHUBBA METHOD QC
############################################################

assert_true(
  identical(
    S8$CytoHubba_algorithm,
    "MCC"
  ),
  "Frozen CytoHubba algorithm must be MCC."
)


assert_true(
  nrow(
    S8$top20_MCC_hubs
  ) == 20L,
  "Frozen CytoHubba result must contain 20 hubs."
)


############################################################
## 9.17 — EXPECTED TOP-20 HUB ORDER
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


hub_order_pass <- identical(
  as.character(
    S8$top20_MCC_hubs$Gene
  ),
  expected_hubs
)


assert_true(
  hub_order_pass,
  "Frozen top-20 MCC hub order differs from validated analysis."
)


assert_true(
  identical(
    as.integer(
      S8$top20_MCC_hubs$MCC_Rank
    ),
    1:20
  ),
  "Frozen MCC ranks must be exactly 1 through 20."
)


assert_true(
  all(
    expected_hubs %in%
      S8$network_nodes
  ),
  "One or more final hub genes are absent from preserved network."
)


assert_true(
  all(
    expected_hubs %in%
      ppi_union$Gene
  ),
  "One or more final hub genes are absent from Module-08 PPI input."
)


############################################################
## 9.18 — MODULE-07 DEG INPUT QC
############################################################

assert_true(
  nrow(
    deg07
  ) == 35009L,
  "Module-07 DEG table must contain 35,009 genes."
)


assert_unique(
  deg07$Gene,
  "Module-07 DEG genes"
)


hub_deg_index <- match(
  expected_hubs,
  deg07$Gene
)


assert_true(
  !anyNA(
    hub_deg_index
  ),
  "One or more hub genes are absent from Module-07 DEG table."
)


hub_deg <- deg07[
  hub_deg_index,
  ,
  drop = FALSE
]


############################################################
## 9.19 — REATTACH CURRENT MODULE-07 DEG STATISTICS
##
## MCC ranking and score remain frozen.
## DEG statistics are linked directly to Module 07 so that
## production outputs do not depend on stale duplicated values.
############################################################

hub_table <- data.frame(
  
  MCC_Rank =
    as.integer(
      S8$top20_MCC_hubs$MCC_Rank
    ),
  
  Gene =
    as.character(
      S8$top20_MCC_hubs$Gene
    ),
  
  MCC_Score =
    S8$top20_MCC_hubs$MCC_Score,
  
  logFC_High_vs_Low =
    hub_deg$logFC,
  
  DEG_FDR_BH =
    hub_deg$FDR_BH,
  
  AveExpr =
    hub_deg$AveExpr,
  
  limma_t =
    hub_deg$t,
  
  stringsAsFactors = FALSE
)


############################################################
## 9.20 — REGRESSION AGAINST FROZEN HUB TABLE
############################################################

assert_true(
  identical(
    hub_table$Gene,
    S8$top20_MCC_hubs$Gene
  ),
  "Production hub gene order differs from frozen S8."
)


assert_true(
  identical(
    hub_table$MCC_Rank,
    as.integer(
      S8$top20_MCC_hubs$MCC_Rank
    )
  ),
  "Production MCC ranks differ from frozen S8."
)


############################################################
## MCC scores are preserved directly rather than recalculated.
############################################################

mcc_score_pass <- identical(
  hub_table$MCC_Score,
  S8$top20_MCC_hubs$MCC_Score
)


assert_true(
  mcc_score_pass,
  "Production MCC scores differ from frozen S8 values."
)


############################################################
## Numerical DEG regression
############################################################

hub_logfc_max_diff <- max(
  abs(
    hub_table$logFC_High_vs_Low -
      S8$top20_MCC_hubs$logFC_High_vs_Low
  )
)


hub_fdr_max_diff <- max(
  abs(
    hub_table$DEG_FDR_BH -
      S8$top20_MCC_hubs$DEG_FDR_BH
  )
)


hub_aveexpr_max_diff <- max(
  abs(
    hub_table$AveExpr -
      S8$top20_MCC_hubs$AveExpr
  )
)


hub_t_max_diff <- max(
  abs(
    hub_table$limma_t -
      S8$top20_MCC_hubs$limma_t
  )
)


HUB_DEG_TOLERANCE <- 1e-10


hub_deg_regression_pass <- all(
  hub_logfc_max_diff <
    HUB_DEG_TOLERANCE,
  
  hub_fdr_max_diff <
    HUB_DEG_TOLERANCE,
  
  hub_aveexpr_max_diff <
    HUB_DEG_TOLERANCE,
  
  hub_t_max_diff <
    HUB_DEG_TOLERANCE
)


assert_true(
  hub_deg_regression_pass,
  paste0(
    "Hub DEG-statistic regression failed. ",
    "Maximum differences: logFC=",
    hub_logfc_max_diff,
    ", FDR=",
    hub_fdr_max_diff,
    ", AveExpr=",
    hub_aveexpr_max_diff,
    ", t=",
    hub_t_max_diff,
    "."
  )
)


############################################################
## 9.21 — PRESERVE FROZEN NETWORK COMPONENTS
############################################################

ppi_node_mapping <-
  S8$node_mapping


ppi_network_nodes <-
  S8$network_nodes


ppi_edges <-
  S8$undirected_edges


ppi_network_qc <-
  S8$network_QC


ppi_method_provenance <-
  S8$method_provenance


ppi_artifact_provenance <-
  S8$artifact_provenance


ppi_identifier_reconciliation <-
  S8$identifier_reconciliation


ppi_interpretation <-
  S8$interpretation


############################################################
## 9.22 — PRODUCTION METHOD SUMMARY
############################################################

ppi_method_summary <- data.frame(
  
  Item = c(
    "PPI_input_derivation",
    "PPI_input_genes",
    "Exact_symbol_matches",
    "Alias_remappings",
    "Alias_reconciliation",
    "Genes_represented_after_reconciliation",
    "Genes_absent_from_exported_network",
    "Network_nodes",
    "Original_reciprocal_interaction_rows",
    "Unique_undirected_edges",
    "Minimum_preserved_combined_score",
    "Maximum_preserved_combined_score",
    "Edges_below_combined_score_0.700",
    "Edges_exactly_combined_score_0.700",
    "CytoHubba_algorithm",
    "Selected_hubs",
    "STRING_network_recomputed",
    "Cytoscape_recomputed",
    "CytoHubba_MCC_recomputed",
    "Interpretation"
  ),
  
  Value = c(
    paste(
      "Union of leading-edge genes from five",
      "selected Hallmark GSEA pathways"
    ),
    455,
    388,
    1,
    "RIGI -> DDX58",
    389,
    66,
    389,
    8088,
    4044,
    minimum_score,
    maximum_score,
    edges_below_0700,
    edges_exact_0700,
    "MCC (Maximal Clique Centrality)",
    20,
    "FALSE",
    "FALSE",
    "FALSE",
    paste(
      "Connectivity-based prioritization within",
      "an immune-focused PPI network; non-causal"
    )
  ),
  
  stringsAsFactors = FALSE
)


############################################################
## 9.23 — PRODUCTION PROVENANCE NOTE
##
## Deliberately conservative wording.
##
## We know the preserved network contains only interactions
## with combined_score >= 0.700.
##
## We do NOT claim an exact historical STRING web-interface
## preset because that interface setting was not independently
## recovered.
############################################################

production_provenance <- data.frame(
  
  Component = c(
    "Input_gene_set",
    "STRING_network",
    "Identifier_reconciliation",
    "Cytoscape_session",
    "CytoHubba_ranking",
    "DEG_statistics",
    "Interpretation"
  ),
  
  Production_behavior = c(
    "Reconstructed from Module 08 and verified against frozen S8",
    "Frozen validated exported network preserved",
    "Frozen reconciliation preserved: RIGI -> DDX58",
    "Historical Cytoscape analysis preserved; not rerun",
    "Historical MCC ranking preserved; not recalculated",
    "Reattached directly from validated Module 07 output",
    "Connectivity-based prioritization only; non-causal"
  ),
  
  stringsAsFactors = FALSE
)


############################################################
## 9.24 — SAVE PROCESSED RDS OBJECTS
############################################################

save_rds_checked(
  hub_table,
  file.path(
    PATHS$processed,
    "09_top20_MCC_hub_genes.rds"
  )
)


save_rds_checked(
  ppi_node_mapping,
  file.path(
    PATHS$processed,
    "09_ppi_node_mapping.rds"
  )
)


save_rds_checked(
  ppi_network_nodes,
  file.path(
    PATHS$processed,
    "09_ppi_network_nodes.rds"
  )
)


save_rds_checked(
  ppi_edges,
  file.path(
    PATHS$processed,
    "09_ppi_undirected_edges.rds"
  )
)


############################################################
## 9.25 — EXPORT TOP-20 HUB TABLE
############################################################

write_csv_checked(
  hub_table,
  file.path(
    ppi_output_dir,
    "09_MCC_top20_hub_genes.csv"
  )
)


############################################################
## 9.26 — EXPORT NODE MAPPING
############################################################

write_csv_checked(
  ppi_node_mapping,
  file.path(
    ppi_output_dir,
    "09_STRING_input_node_mapping.csv"
  )
)


############################################################
## 9.27 — EXPORT NETWORK NODE LIST
############################################################

network_node_table <- data.frame(
  Gene =
    ppi_network_nodes,
  stringsAsFactors = FALSE
)


write_csv_checked(
  network_node_table,
  file.path(
    ppi_output_dir,
    "09_STRING_network_nodes.csv"
  )
)


############################################################
## 9.28 — EXPORT UNDIRECTED EDGE TABLE
############################################################

write_csv_checked(
  ppi_edges,
  file.path(
    ppi_output_dir,
    "09_STRING_undirected_edges.csv"
  )
)


############################################################
## 9.29 — EXPORT NETWORK QC / PROVENANCE
############################################################

write_csv_checked(
  ppi_network_qc,
  file.path(
    ppi_output_dir,
    "09_STRING_network_QC.csv"
  )
)


write_csv_checked(
  ppi_method_provenance,
  file.path(
    ppi_output_dir,
    "09_PPI_method_provenance_frozen.csv"
  )
)


write_csv_checked(
  ppi_artifact_provenance,
  file.path(
    ppi_output_dir,
    "09_PPI_historical_artifact_provenance.csv"
  )
)


write_csv_checked(
  ppi_method_summary,
  file.path(
    ppi_output_dir,
    "09_PPI_method_summary.csv"
  )
)


write_csv_checked(
  production_provenance,
  file.path(
    ppi_output_dir,
    "09_PPI_production_provenance.csv"
  )
)


############################################################
## 9.30 — IDENTIFIER RECONCILIATION TABLE
############################################################

identifier_reconciliation_table <- data.frame(
  
  Input_Gene =
    names(
      ppi_identifier_reconciliation
    ),
  
  Network_Label =
    unname(
      ppi_identifier_reconciliation
    ),
  
  stringsAsFactors = FALSE
)


write_csv_checked(
  identifier_reconciliation_table,
  file.path(
    ppi_output_dir,
    "09_PPI_identifier_reconciliation.csv"
  )
)


############################################################
## 9.31 — HUB DEG REGRESSION TABLE
############################################################

hub_deg_regression <- data.frame(
  
  Metric = c(
    "logFC_High_vs_Low",
    "DEG_FDR_BH",
    "AveExpr",
    "limma_t"
  ),
  
  Max_abs_difference = c(
    hub_logfc_max_diff,
    hub_fdr_max_diff,
    hub_aveexpr_max_diff,
    hub_t_max_diff
  ),
  
  stringsAsFactors = FALSE
)


write_csv_checked(
  hub_deg_regression,
  file.path(
    PATHS$logs,
    "09_hub_DEG_regression.csv"
  )
)


############################################################
## 9.32 — REFERENCE MD5
############################################################

ppi_reference_md5 <- tools::md5sum(
  ppi_reference_file
)


ppi_reference_hash <- data.frame(
  
  File =
    names(
      ppi_reference_md5
    ),
  
  MD5 =
    unname(
      ppi_reference_md5
    ),
  
  stringsAsFactors = FALSE
)


write_csv_checked(
  ppi_reference_hash,
  file.path(
    PATHS$logs,
    "09_PPI_reference_MD5.csv"
  )
)


############################################################
## 9.33 — MODULE QC TABLE
############################################################

module09_qc <- data.frame(
  
  Metric = c(
    "Module08 PPI input genes",
    "PPI input regression",
    "Exact symbol matches",
    "Alias remappings",
    "Alias mapping",
    "Represented input genes",
    "Absent input genes",
    "Network nodes",
    "Unique undirected edges",
    "Minimum combined score",
    "Maximum combined score",
    "Edges below 0.700",
    "Edges exactly 0.700",
    "CytoHubba algorithm",
    "Selected hub genes",
    "Exact hub order regression",
    "MCC score preservation",
    "Hub DEG regression",
    "Hub logFC max difference",
    "Hub FDR max difference",
    "Hub AveExpr max difference",
    "Hub limma-t max difference",
    "STRING network recomputed",
    "Cytoscape recomputed",
    "CytoHubba MCC recomputed"
  ),
  
  Value = c(
    455,
    ppi_input_pass,
    388,
    1,
    "RIGI -> DDX58",
    represented_count,
    absent_count,
    length(
      ppi_network_nodes
    ),
    nrow(
      ppi_edges
    ),
    minimum_score,
    maximum_score,
    edges_below_0700,
    edges_exact_0700,
    S8$CytoHubba_algorithm,
    nrow(
      hub_table
    ),
    hub_order_pass,
    mcc_score_pass,
    hub_deg_regression_pass,
    format(
      hub_logfc_max_diff,
      scientific = TRUE,
      digits = 17
    ),
    format(
      hub_fdr_max_diff,
      scientific = TRUE,
      digits = 17
    ),
    format(
      hub_aveexpr_max_diff,
      scientific = TRUE,
      digits = 17
    ),
    format(
      hub_t_max_diff,
      scientific = TRUE,
      digits = 17
    ),
    FALSE,
    FALSE,
    FALSE
  ),
  
  stringsAsFactors = FALSE
)


write_csv_checked(
  module09_qc,
  file.path(
    PATHS$logs,
    "09_ppi_hubs_QC.csv"
  )
)


############################################################
## 9.34 — SESSION INFORMATION
############################################################

capture.output(
  sessionInfo(),
  file = file.path(
    PATHS$logs,
    "09_ppi_hubs_sessionInfo.txt"
  )
)


############################################################
## 9.35 — FINAL STATUS
############################################################

cat(
  "\n=============================================\n"
)


cat(
  "09_ppi_hubs.R: PASS\n"
)


cat(
  "\nPPI input\n"
)


cat(
  "Module-08 leading-edge union genes: ",
  nrow(
    ppi_union
  ),
  "\n",
  sep = ""
)


cat(
  "455-gene input regression: VERIFIED\n"
)


cat(
  "\nSTRING network\n"
)


cat(
  "Exact symbol matches: 388\n"
)


cat(
  "Alias remappings: 1\n"
)


cat(
  "Identifier reconciliation: RIGI -> DDX58\n"
)


cat(
  "Represented after reconciliation: ",
  represented_count,
  "\n",
  sep = ""
)


cat(
  "Absent from exported network: ",
  absent_count,
  "\n",
  sep = ""
)


cat(
  "Unique network nodes: ",
  length(
    ppi_network_nodes
  ),
  "\n",
  sep = ""
)


cat(
  "Unique undirected edges: ",
  nrow(
    ppi_edges
  ),
  "\n",
  sep = ""
)


cat(
  "Preserved combined-score range: ",
  minimum_score,
  " - ",
  maximum_score,
  "\n",
  sep = ""
)


cat(
  "Edges below combined_score 0.700: ",
  edges_below_0700,
  "\n",
  sep = ""
)


cat(
  "Exact historical STRING web-interface preset recovered: FALSE\n"
)


cat(
  "\nCytoHubba\n"
)


cat(
  "Ranking algorithm: MCC (Maximal Clique Centrality)\n"
)


cat(
  "Selected hubs: ",
  nrow(
    hub_table
  ),
  "\n",
  sep = ""
)


cat(
  "Exact top-20 hub order: VERIFIED\n"
)


cat(
  "MCC scores: FROZEN / PRESERVED\n"
)


cat(
  "Hub DEG statistics regression: VERIFIED\n"
)


cat(
  "Maximum hub DEG-statistic difference: ",
  format(
    max(
      hub_deg_regression$Max_abs_difference
    ),
    scientific = TRUE,
    digits = 17
  ),
  "\n",
  sep = ""
)


cat(
  "\nProduction behavior\n"
)


cat(
  "STRING network recomputed: FALSE\n"
)


cat(
  "Cytoscape analysis recomputed: FALSE\n"
)


cat(
  "CytoHubba MCC recomputed: FALSE\n"
)


cat(
  "Module-07 DEG statistics reattached: TRUE\n"
)


cat(
  "Interpretation: highly connected immune-associated genes ",
  "within the curated immune-focused PPI network; non-causal\n",
  sep = ""
)


cat(
  "=============================================\n"
)