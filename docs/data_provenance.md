# Data Provenance

## Overview

This repository contains analysis scripts, selected processed outputs,
validation workflows, and reproducibility documentation for the
TCGA-BRCA immune-state analysis project.

The repository is designed as a computational research record. Raw
sequencing datasets and restricted external reference objects are not
redistributed when they are controlled by their original repositories or
subject to access limitations.

------------------------------------------------------------------------

# TCGA-BRCA Dataset

## Dataset

The primary analysis uses The Cancer Genome Atlas Breast Cancer
(TCGA-BRCA) cohort.

## Data Type

The analysis operates on processed gene-expression matrices and
associated sample metadata.

Raw sequencing processing steps (such as FASTQ alignment and transcript
quantification) are outside the scope of this repository.

------------------------------------------------------------------------

# Expression Data Processing

The downstream analysis begins from processed expression data generated
before the analytical modules.

The workflow uses:

-   processed gene-expression matrix
-   sample metadata
-   immune phenotype annotations

The differential-expression analysis does not recompute upstream
expression preprocessing steps.

------------------------------------------------------------------------

# Immune Phenotype Definition

Samples are categorized according to the project immune-state framework.

For differential-expression analysis:

-   Immune_High samples are compared against Immune_Low samples.
-   Immune_Mid samples are excluded from the DEG comparison.

The primary contrast is:

Immune_High - Immune_Low

------------------------------------------------------------------------

# Differential Expression Reference Objects

The repository includes reference DEG files used for regression and
reproducibility checks.

These files are used to verify regenerated results against validated
outputs.

Reference files include:

-   historical DEG result reference
-   complete DEG reference output

------------------------------------------------------------------------

# External Validation Dataset

External validation uses the SCAN-B breast cancer cohort.

The validation workflow uses:

-   SCAN-B expression data
-   clinical metadata
-   frozen validation reference objects for reproducibility checking

------------------------------------------------------------------------

# Reproducibility Notes

The repository provides:

-   analysis scripts
-   processed intermediate objects where applicable
-   validation scripts
-   quality-control outputs

Users wishing to reproduce the complete analysis should obtain original
datasets from their respective repositories and follow the documented
workflow.
