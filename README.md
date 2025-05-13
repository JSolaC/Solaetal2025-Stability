# Sola et al. 2025 — Counteracting Cascades Challenge the Heterogeneity–Stability Relationship

This repository contains the data and R code used for the study:
Sola, J., Fairchild, T.P., Perkins, M.J., Bull, J.C., & Griffin, J.N. (2025). Counteracting cascades challenge the heterogeneity–stability relationship.

The study presents a three-year field experiment conducted on a rocky shore to test the ecological assumption that spatial heterogeneity enhances community temporal stability. Contrary to theoretical expectations, we found no net effect of heterogeneity on stability, due to counteracting cascades involving species richness, population dynamics, and dominant species suppression. This repository enables full reproducibility of the analyses, figures, and models presented in the paper.

The repository is organised into four main folders. The Data folder contains the processed data used in the analyses, along with scripts for calculating temporal stability metrics and related components. The Models folder includes mixed-effects models and Structural Equation Models used to assess species richness, community composition, and dominant taxa. The Plots folder provides code to generate all figures from the main text. The SupplementaryMaterials folder contains additional scripts used to reproduce supplementary figures, conduct sensitivity analyses, and perform further assessments.

To reproduce the results, begin by running the scripts in the Models folder, which rely on the data files in the Data folder. Once the models have been executed, use the code in the Plots folder to visualise the results. Supplementary analyses and additional figures can be reproduced using the scripts in the SupplementaryMaterials folder.
