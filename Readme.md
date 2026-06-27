# Robustness of ANOVA, Welch's ANOVA, and Kruskal–Wallis in Small Samples

## Overview
This repository contains all simulation code, results, and figures for the manuscript:

**"Robustness of ANOVA, Welch's ANOVA, and Kruskal–Wallis in Small Samples Under Distributional and Variance Violations"**

Submitted to: *Communications in Statistics – Simulation and Computation*

---

## Repository Structure
├── code/

│   ├── h0_hepsi_son.R              # Type I error (H0) simulations

│   ├── Final h1.R                  # Power (H1) simulations - main

│   ├── tekrar.R                    # Negbin heterogeneous rerun (seed fix)

│   ├── eksikler.R                  # Negbin n=10 k=8 extended rerun

│   ├── Birlestirme.R               # Final dataset assembly

│   ├── Grafik.R                    # Figures 2–4 (Empirical Power)

│   └── Type I Error Grafik.R       # Figure 1 (Type I Error)

│

├── results/

│   ├── H1_FINAL.csv                # Final H1 results (240 conditions)

│   ├── H0_FINAL_RESULTS_ALL_SCENARIOS.xlsx  # Final H0 results

│   └── intermediate/

│       ├── negbin_het_rerun_fixed.csv

│       ├── negbin_het_incomplete_fixed.csv

│       └── poisson_negbin_rerun_clean.csv

│

└── figures/

├── Figure1_Type1Error.png

├── Figure2_ANOVA_Power.png

├── Figure3_Welch_Power.png

└── Figure4_KW_Power.png

---

## Simulation Design

| Factor | Levels |
|---|---|
| Tests | Classical ANOVA, Welch's ANOVA, Kruskal–Wallis |
| Distributions | Normal, Exponential, Lognormal, Poisson, Negative Binomial |
| Variance structure | Homogeneous, Heterogeneous (SD ratio 1:2, variance ratio 4:1) |
| Sample size per group | n ∈ {5, 6, 8, 10} |
| Number of groups | k ∈ {3, 5, 8} |
| Design balance | Balanced, Unbalanced |
| Replications | 10,000 per condition |
| Target effect size | Cohen's f = 0.69 |
| Nominal alpha | 0.05 |

---

## Reproducibility Notes

- Fixed seeds were used for all simulations
- The negbin heterogeneous conditions were rerun after correcting a seed generation error in the original script that failed to incorporate `var_type` in the seed formula, causing homogeneous and heterogeneous conditions to receive identical seeds
- Four conditions (negbin, n = 10, k = 8, balanced and unbalanced × homogeneous and heterogeneous) required extended attempts (MAX_ATTEMPTS = 15,000,000) due to low acceptance rates (~0.17%); all four reached the target of 10,000 accepted replications
- The final dataset (`H1_FINAL.csv`) contains 240 rows with no missing values for any condition

---

## Requirements

- R version ≥ 4.0  
- Packages: `dplyr`, `ggplot2`, `tidyr`, `openxlsx`, `car`

---

## Citation

If you use this code or data, please cite the associated manuscript (citation will be updated upon publication).
