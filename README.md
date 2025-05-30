# Sola et al. 2025 — *Counteracting Cascades Challenge the Heterogeneity–Stability Relationship*

This repository provides the data and R code used in the study:

**Sola, J., Fairchild, T.P., Perkins, M.J., Bull, J.C., & Griffin, J.N. (2025).** *Counteracting cascades challenge the heterogeneity–stability relationship.* *Ecology Letters.*

The study is based on a three-year field experiment on a temperate rocky shore, testing the widely held ecological assumption that spatial heterogeneity promotes community temporal stability. Contrary to theoretical expectations, we found that heterogeneity had no net stabilising effect, due to opposing cascades involving species richness, population variability, and suppression of dominant species. This repository enables full reproducibility of all analyses, figures, and models presented in the paper.

## Repository Structure

The repository is organised into four main folders:

- **`Data/`**  
  Contains the processed datasets used in the analyses, as well as scripts for calculating temporal stability and its components (e.g., species synchrony, variability, richness).

- **`Models/`**  
  Includes the full suite of statistical models used in the study, including mixed-effects models and structural equation models (SEMs) examining richness, composition, and the roles of dominant taxa.

- **`Plots/`**  
  Contains scripts to generate all main-text figures, based on model outputs.

- **`SupplementaryMaterials/`**  
  Provides additional code for reproducing supplementary figures, running sensitivity analyses, and exploring further model outputs.

## Reproducibility Instructions

1. Begin by running the scripts in the **`Models/`** folder. These rely on data files located in **`Data/`**.
2. Once the models are executed, use the scripts in **`Plots/`** to generate the figures included in the main text.
3. Supplementary analyses and figures can be reproduced using the scripts in **`SupplementaryMaterials/`**.

## Notes on Data

The datasets in the folder **`Datasets/`** provided here are derived from a subset of the full raw dataset. However, the full processed dataset used for the published analyses is included to ensure full reproducibility.

## Package Version Compatibility

Please note that all analyses were originally conducted using R package versions from 2020. Running the code with newer package versions may lead to errors or incompatibilities, particularly with packages such as nlme and piecewiseSEM (see for example: GitHub Issue #763).

To ensure full reproducibility, we recommend managing package versions using tools such as remotes, renv, or checkpoint. These allow users to recreate the original analysis environment with the appropriate package versions.

We plan to include a future update with code to help users automatically replicate the exact setup used in the original analyses.
