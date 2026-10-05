# Run Instructions

## Overview

This repository contains the computational workflow for TCGA-BRCA
immune-state analysis.

The scripts are organized as sequential modules. Each module performs a
defined analytical step and generates intermediate objects,
quality-control files, or final outputs.

------------------------------------------------------------------------

# Requirements

## Software

Recommended:

-   R version 4.x
-   RStudio (recommended)

## Required Packages

The workflow requires R packages used by individual modules, including:

-   limma
-   GSVA
-   survival
-   survminer
-   tidyverse and related data-processing packages

Package requirements are checked during execution.

------------------------------------------------------------------------

# Repository Structure

    R/
    ├── 00_setup.R
    ├── 01_*.R
    ├── 02_*.R
    ...
    ├── 13_scanb_validation.R

    data_raw/
    data_processed/
    external/
    results/
    logs/
    docs/

------------------------------------------------------------------------

# Execution Workflow

Run scripts sequentially according to the numbered module order:

    00_setup.R

    01_import_tcga.R

    02_expression_processing.R

    03_*.R

    04_*.R

    05_*.R

    06_*.R

    07_differential_expression.R

    08_*.R

    09_*.R

    10_*.R

    11_*.R

    12_*.R

    13_scanb_validation.R

The exact script names should match the files present in the repository.

------------------------------------------------------------------------

# Differential Expression Module

The differential-expression workflow:

1.  Loads finalized expression data.
2.  Loads immune phenotype metadata.
3.  Selects Immune_High and Immune_Low groups.
4.  Performs limma differential expression analysis.
5.  Generates complete and significant DEG tables.
6.  Performs regression checks against reference DEG outputs.

The primary contrast is:

    Immune_High - Immune_Low

------------------------------------------------------------------------

# SCAN-B Validation Module

The SCAN-B validation module:

1.  Loads external validation cohort files.
2.  Applies the previously defined immune-state framework.
3.  Performs validation analyses.
4.  Compares outputs against frozen validation references where
    available.

------------------------------------------------------------------------

# Outputs

Generated files are saved into:

    results/
    logs/
    data_processed/

Outputs include:

-   differential-expression tables
-   quality-control reports
-   validation summaries
-   processed analysis objects

------------------------------------------------------------------------

# Notes on Reproducibility

This repository is intended to provide transparent documentation of the
computational workflow.

Raw datasets are not included when they are available through original
repositories.

Users should obtain original datasets from the appropriate sources and
then follow the documented module sequence.
