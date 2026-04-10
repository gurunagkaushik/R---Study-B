# Study B: Hospital Throughput Analysis

**Jain et al. 2025 — Stage I Endometrial Cancer | RAS vs LPS Surgery**

## Overview

This R script models the operational benefits of switching from laparoscopic (LPS)
to robotic-assisted (RAS) surgery for Stage I endometrial cancer. Using aggregate
data from Jain et al. (2025), it answers two questions:

1. How many extra surgeries can a hospital fit per year? (OR Throughput)
2. How many bed-days are freed per year? (Bed Liberation)

## Data Source

Jain et al., Cureus 2025. DOI: 10.7759/cureus.90782
Only published means and SDs are used — no patient-level data required.

## Requirements

R 4.0+ and the following packages:

install.packages(c("tidyverse", "MASS", "fitdistrplus", "ggplot2", "scales"))

## How to Run

1. Open study_B_throughput.R in RStudio or any R environment
2. Install packages above if not already installed
3. Run the script top to bottom (Ctrl+Alt+R / Cmd+Alt+R)

All outputs are saved to your working directory automatically.

## Outputs

| File                            | Description                              |
|---------------------------------|------------------------------------------|
| study_B_results.csv             | Summary table of all scenarios           |
| or_time_distribution.png        | Distribution of OR time saved per case   |
| additional_cases_by_volume.png  | Extra surgeries per year by volume       |
| beddays_freed_by_volume.png     | Bed-days freed per year by volume        |

## Methods Summary

- Gamma distributions fitted to OR time and LOS via Method of Moments
- Monte Carlo simulation: 10,000 iterations
- Three hospital volume scenarios: 100, 300, and 500 cases/year
- Sensitivity analysis on OR day length: 420 / 450 / 480 usable minutes
- Random seed: 2026

## Key Results (Base Case: 450 OR min/day)

| Hospital Volume  | Extra Surgeries/Year | Bed-Days Freed/Year |
|------------------|----------------------|----------------------|
| 100 cases/year   | ~20                  | ~106                 |
| 300 cases/year   | ~60                  | ~318                 |
| 500 cases/year   | ~100                 | ~530                 |

Note: Values above are approximate. Run the script for exact means and 95% CIs.

## Reference

Jain et al. Robotic-Assisted vs Laparoscopic Surgery for Stage I Endometrial Cancer.
Cureus 2025. DOI: 10.7759/cureus.90782
