# clim_risk_phylosdm

> **Status:** ![v1](https://img.shields.io/badge/version-v1-blue) ![Unpublished](https://img.shields.io/badge/status-unpublished-orange) ![Work in Progress](https://img.shields.io/badge/work--in--progress-yellow)
>
> This is pre-publication research code. Results and methods are subject to change. Please do not cite or redistribute without permission.

## Overview

This repository contains the analysis pipeline for estimating **climate change risk in amphibians** using phylogenetically-informed species distribution models (SDMs). The framework integrates occurrence data, expert range maps, and phylogenetic trees to model habitat suitability, then decomposes climate risk into three components:

- **Exposure** — magnitude of climate change within a species' predicted range
- **Sensitivity** — partial derivative of habitat suitability with respect to climate variables (from model coefficients)
- **Vulnerability** — combined measure capturing projected habitat suitability losses (sensitivity × exposure)

Models are fit using either **Stan** (LGCP) or **TMB**, and the pipeline is designed to run at scale on an HPC cluster (Yale McCleary / Slurm).

## Repository Structure

```
scripts/
├── 00-harmonize.R          # Harmonize GBIF occurrences, expert ranges, and phylogeny
├── 00-makeClusters.R       # Build species clusters for parallel jobs
├── 01-config.yaml          # Experiment configuration (model type, HPC specs, etc.)
├── 02-build_jobs.R         # Generate DSQ job arrays
├── 03-*.R                  # Data preparation and soft-clipping
├── 04-*.R                  # Package data for model input
├── 05-run_model.R          # Fit SDM (Stan or TMB)
├── 05-stan_model.R         # Stan model specification (LGCP with background)
├── 10-15-*.R               # Prediction: conditional, continuous, thresholding, binary
├── 20-exposure.R           # Compute climate exposure (SSP3-7.0, 2041–2060)
├── 21-sensitivity.R        # Compute climate sensitivity from model coefficients
├── 22-vulnerability.R      # Compute vulnerability = sensitivity × exposure
├── figures/                # Figure-generating scripts
│   ├── data_assembly/      # Compile master dataset, species lists, quantiles
│   └── *.R                 # Richness maps, phylo plots, case studies, etc.
├── handlers/               # Shell scripts for fetching data on HPC
├── sim/                    # Simulation scripts
└── 999-model_compare.R     # TMB vs Stan comparison diagnostics
```

## Data

Input data (not tracked in this repository) includes:

- Amphibian occurrence records (GBIF, harmonized)
- Expert range maps (IUCN / MOL)
- Amphibian phylogeny (MCC tree, 100 trees)
- Current and future climate layers (SSP3-7.0, 2041–2060)

Data paths are configured for the Yale HPC environment (`/vast/palmer/pi/jetz/ss4224/`). To run locally, update root paths in each script.

## Dependencies

R packages: `rstan`, `TMB`, `terra`, `ape`, `yaml`

## Running the Pipeline

1. Edit `scripts/01-config.yaml` with experiment settings
2. Run `scripts/02-build_jobs.R` to generate job arrays
3. Submit via DSQ on Slurm, or run scripts sequentially in order (00 → 22)

## Contact

Shubhi Sharma — Yale University ([shubhi.sharma@yale.edu](mailto:shubhi.sharma@yale.edu))