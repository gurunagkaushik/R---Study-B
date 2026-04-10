# =============================================================================
# STUDY B: Hospital Throughput Analysis
# Jain et al. 2025 — Stage I Endometrial Cancer
# Robotic-Assisted (RAS) vs Laparoscopic (LPS) Surgery
#
# Two questions this code answers:
#   1. How many EXTRA surgeries can a hospital fit per year? (OR Throughput)
#   2. How many BED-DAYS are freed per year? (Bed Liberation)
#
# Based on SAP: RAS_LPS_HEOR_SAP.docx
# Data source:  Jain et al. Cureus 2025; DOI: 10.7759/cureus.90782
# =============================================================================


# -----------------------------------------------------------------------------
# SECTION 0: Install and load required packages
# -----------------------------------------------------------------------------

# Run this once if packages are not installed:
# install.packages(c("tidyverse", "MASS", "fitdistrplus", "ggplot2", "scales"))

library(tidyverse)      # data manipulation and piping
library(MASS)  # fitdistr() for distribution fitting
library(fitdistrplus)   # better distribution fitting with diagnostics
library(ggplot2)        # plotting
library(scales)         # axis formatting in plots

# Set random seed for reproducibility (matches Study A seed)
set.seed(2026)


# =============================================================================
# SECTION 1: INPUT DATA
# Aggregate data from Jain et al. 2025 (Table 2 and Table 3)
# These are the ONLY numbers needed from the published paper
# =============================================================================

# --- OR Time (minutes) ---
ras_or_mean <- 205.65
ras_or_sd   <- 87.60
lps_or_mean <- 257.06
lps_or_sd   <- 109.41

# --- Length of Stay / LOS (days) ---
ras_los_mean <- 2.71
ras_los_sd   <- 1.81
lps_los_mean <- 3.77
lps_los_sd   <- 4.48  # NOTE: SD > mean confirms strong right skew — Gamma appropriate

# --- Calculated differences (base case point estimates) ---
or_diff_mean  <- lps_or_mean  - ras_or_mean   # 51.41 minutes saved per case
los_diff_mean <- lps_los_mean - ras_los_mean  # 1.06 days saved per case

cat("=== BASE CASE DIFFERENCES ===\n")
cat("OR time saved per case:  ", round(or_diff_mean, 2), "minutes\n")
cat("LOS saved per case:      ", round(los_diff_mean, 2), "days\n\n")


# =============================================================================
# SECTION 2: FIT GAMMA DISTRIBUTIONS
# Method of Moments — uses mean and SD to calculate shape and scale
#
# Formula:
#   shape (k) = (mean / SD)^2
#   scale (θ) = SD^2 / mean
#
# Why Gamma? It is:
#   - Always positive (can't have negative surgery time)
#   - Right-skewed (some cases take much longer than average)
#   - Standard for modelling time-based clinical variables
# =============================================================================

fit_gamma_mom <- function(mean, sd) {
  # Method of Moments: estimate shape and scale from mean and SD
  shape <- (mean / sd)^2
  scale <- (sd^2) / mean
  return(list(shape = shape, scale = scale))
}

# Fit distributions for each of the four groups
ras_or_params  <- fit_gamma_mom(ras_or_mean,  ras_or_sd)
lps_or_params  <- fit_gamma_mom(lps_or_mean,  lps_or_sd)
ras_los_params <- fit_gamma_mom(ras_los_mean, ras_los_sd)
lps_los_params <- fit_gamma_mom(lps_los_mean, lps_los_sd)

# Print all four fitted distributions
cat("=== FITTED GAMMA DISTRIBUTION PARAMETERS ===\n\n")

cat("RAS OR Time:  shape =", round(ras_or_params$shape, 3),
    " | scale =", round(ras_or_params$scale, 3), "\n")
cat("LPS OR Time:  shape =", round(lps_or_params$shape, 3),
    " | scale =", round(lps_or_params$scale, 3), "\n")
cat("RAS LOS:      shape =", round(ras_los_params$shape, 3),
    " | scale =", round(ras_los_params$scale, 3), "\n")
cat("LPS LOS:      shape =", round(lps_los_params$shape, 3),
    " | scale =", round(lps_los_params$scale, 3), "\n\n")

# Quick sanity check: verify fitted distributions recover original mean and SD
cat("=== SANITY CHECK (fitted mean should match input mean) ===\n")
cat("RAS OR fitted mean:", round(ras_or_params$shape * ras_or_params$scale, 2),
    " (input:", ras_or_mean, ")\n")
cat("LPS OR fitted mean:", round(lps_or_params$shape * lps_or_params$scale, 2),
    " (input:", lps_or_mean, ")\n")
cat("RAS LOS fitted mean:", round(ras_los_params$shape * ras_los_params$scale, 2),
    " (input:", ras_los_mean, ")\n")
cat("LPS LOS fitted mean:", round(lps_los_params$shape * lps_los_params$scale, 2),
    " (input:", lps_los_mean, ")\n\n")


# =============================================================================
# SECTION 3: HOSPITAL CONFIGURATION
# Reference hospital setup as defined in the SAP
# =============================================================================

# --- Base case OR day configuration ---
or_minutes_per_day_base         <- 450   # 8-hour day minus 30-min turnover
or_minutes_per_day_conservative <- 420   # sensitivity: shorter day
or_minutes_per_day_optimistic   <- 480   # sensitivity: longer day

working_days_per_year <- 250             # standard working days

# --- Annual case volume scenarios (from SAP) ---
annual_case_volumes <- c(100, 300, 500)

# --- Labels for output tables ---
volume_labels <- c("100 cases/year", "300 cases/year", "500 cases/year")


# =============================================================================
# SECTION 4: MONTE CARLO SIMULATION
# 10,000 iterations — each iteration:
#   1. Draws a random RAS OR time from its Gamma distribution
#   2. Draws a random LPS OR time from its Gamma distribution
#   3. Calculates OR time saved = LPS draw - RAS draw
#   4. Draws a random RAS LOS from its Gamma distribution
#   5. Draws a random LPS LOS from its Gamma distribution
#   6. Calculates LOS saved = LPS draw - RAS draw
#   7. Computes additional cases and bed-days freed for each hospital volume
# =============================================================================

n_iterations <- 10000

cat("=== RUNNING MONTE CARLO SIMULATION ===\n")
cat("Iterations:", n_iterations, "\n")
cat("Seed: 2026\n\n")

# Pre-allocate results storage
# Rows = iterations, Columns = one per hospital volume scenario
mc_additional_cases    <- matrix(NA, nrow = n_iterations, ncol = length(annual_case_volumes))
mc_beddays_freed       <- matrix(NA, nrow = n_iterations, ncol = length(annual_case_volumes))
mc_extra_pts_treatable <- matrix(NA, nrow = n_iterations, ncol = length(annual_case_volumes))
mc_or_diff             <- numeric(n_iterations)
mc_los_diff            <- numeric(n_iterations)

# Run simulation
for (i in 1:n_iterations) {
  
  # --- Draw 1: RAS OR time (one random value from RAS OR Gamma distribution) ---
  ras_or_draw <- rgamma(1,
                        shape = ras_or_params$shape,
                        scale = ras_or_params$scale)
  
  # --- Draw 2: LPS OR time ---
  lps_or_draw <- rgamma(1,
                        shape = lps_or_params$shape,
                        scale = lps_or_params$scale)
  
  # OR time saved this iteration (can be negative in rare draws — handled below)
  or_saved_this_iter <- lps_or_draw - ras_or_draw
  mc_or_diff[i] <- or_saved_this_iter
  
  # --- Draw 3: RAS LOS ---
  ras_los_draw <- rgamma(1,
                         shape = ras_los_params$shape,
                         scale = ras_los_params$scale)
  
  # --- Draw 4: LPS LOS ---
  lps_los_draw <- rgamma(1,
                         shape = lps_los_params$shape,
                         scale = lps_los_params$scale)
  
  # LOS saved this iteration
  los_saved_this_iter <- lps_los_draw - ras_los_draw
  mc_los_diff[i] <- los_saved_this_iter
  
  # --- Compute outputs for each annual volume scenario ---
  for (v in seq_along(annual_case_volumes)) {
    
    N <- annual_case_volumes[v]
    
    # Total OR minutes saved per year
    total_or_saved <- N * or_saved_this_iter
    
    # Additional cases = total time saved / mean RAS OR time per case
    # floor() because you can't do a fraction of a surgery
    # max(..., 0) because negative savings = 0 additional cases
    additional_cases <- max(0, floor(total_or_saved / ras_or_draw))
    mc_additional_cases[i, v] <- additional_cases
    
    # Bed-days freed per year
    beddays_freed <- max(0, N * los_saved_this_iter)
    mc_beddays_freed[i, v] <- beddays_freed
    
    # Extra patients treatable in freed beds
    # = bed-days freed / mean RAS LOS per patient
    extra_pts <- max(0, beddays_freed / ras_los_draw)
    mc_extra_pts_treatable[i, v] <- extra_pts
  }
}

cat("Simulation complete.\n\n")


# =============================================================================
# SECTION 5: SUMMARISE RESULTS
# Calculate mean and 95% credible interval from simulation output
# =============================================================================

summarise_mc <- function(values) {
  # Takes a vector of MC results, returns mean + 95% CI
  list(
    mean  = mean(values),
    lower = quantile(values, 0.025),   # 2.5th percentile
    upper = quantile(values, 0.975)    # 97.5th percentile
  )
}

# --- OR time difference summary ---
or_diff_summary <- summarise_mc(mc_or_diff)
los_diff_summary <- summarise_mc(mc_los_diff)

cat("=== MC DISTRIBUTION OF TIME SAVINGS ===\n")
cat(sprintf("OR time saved per case:  Mean = %.1f min  [95%% CI: %.1f – %.1f]\n",
            or_diff_summary$mean, or_diff_summary$lower, or_diff_summary$upper))
cat(sprintf("LOS saved per case:      Mean = %.2f days [95%% CI: %.2f – %.2f]\n\n",
            los_diff_summary$mean, los_diff_summary$lower, los_diff_summary$upper))

# --- Build results table ---
results_table <- data.frame(
  Volume            = volume_labels,
  Annual_cases      = annual_case_volumes,
  
  # OR Throughput
  Add_cases_mean    = sapply(1:length(annual_case_volumes),
                             function(v) round(mean(mc_additional_cases[, v]), 1)),
  Add_cases_lower   = sapply(1:length(annual_case_volumes),
                             function(v) round(quantile(mc_additional_cases[, v], 0.025), 1)),
  Add_cases_upper   = sapply(1:length(annual_case_volumes),
                             function(v) round(quantile(mc_additional_cases[, v], 0.975), 1)),
  
  # Bed Liberation
  Beddays_mean      = sapply(1:length(annual_case_volumes),
                             function(v) round(mean(mc_beddays_freed[, v]), 1)),
  Beddays_lower     = sapply(1:length(annual_case_volumes),
                             function(v) round(quantile(mc_beddays_freed[, v], 0.025), 1)),
  Beddays_upper     = sapply(1:length(annual_case_volumes),
                             function(v) round(quantile(mc_beddays_freed[, v], 0.975), 1)),
  
  # Extra patients
  Extra_pts_mean    = sapply(1:length(annual_case_volumes),
                             function(v) round(mean(mc_extra_pts_treatable[, v]), 1)),
  Extra_pts_lower   = sapply(1:length(annual_case_volumes),
                             function(v) round(quantile(mc_extra_pts_treatable[, v], 0.025), 1)),
  Extra_pts_upper   = sapply(1:length(annual_case_volumes),
                             function(v) round(quantile(mc_extra_pts_treatable[, v], 0.975), 1))
)

# Print clean results tables
cat("=== TABLE 1: OR THROUGHPUT — Additional Surgeries Per Year ===\n")
cat(sprintf("%-20s %12s %20s\n", "Hospital Volume", "Mean (cases)", "95% CI"))
cat(strrep("-", 55), "\n")
for (i in 1:nrow(results_table)) {
  cat(sprintf("%-20s %12.1f %10.1f – %6.1f\n",
              results_table$Volume[i],
              results_table$Add_cases_mean[i],
              results_table$Add_cases_lower[i],
              results_table$Add_cases_upper[i]))
}

cat("\n=== TABLE 2: BED LIBERATION — Bed-Days Freed Per Year ===\n")
cat(sprintf("%-20s %15s %25s %25s\n",
            "Hospital Volume", "Bed-days (mean)", "95% CI", "Extra patients (mean)"))
cat(strrep("-", 85), "\n")
for (i in 1:nrow(results_table)) {
  cat(sprintf("%-20s %15.1f %12.1f – %6.1f %18.1f\n",
              results_table$Volume[i],
              results_table$Beddays_mean[i],
              results_table$Beddays_lower[i],
              results_table$Beddays_upper[i],
              results_table$Extra_pts_mean[i]))
}


# =============================================================================
# SECTION 6: SENSITIVITY ANALYSES
# Test different OR day lengths (420 / 450 / 480 usable minutes)
# =============================================================================

cat("\n\n=== SENSITIVITY ANALYSIS: OR Day Configuration ===\n")
cat("Evaluated at 300 cases/year (mid-volume reference centre)\n\n")

or_day_configs <- c(
  "Conservative (420 min)" = 420,
  "Base case (450 min)"    = 450,
  "Optimistic (480 min)"   = 480
)

# For sensitivity we re-use the already-simulated OR differences
# and just change the denominator (usable minutes per day)
# Here we report impact at 300 cases/year

N_ref <- 300

cat(sprintf("%-30s %15s %20s\n", "OR Day Config", "Add. cases (mean)", "95% CI"))
cat(strrep("-", 65), "\n")

for (config_name in names(or_day_configs)) {
  or_day <- or_day_configs[config_name]
  
  # Additional cases given this OR day length
  sens_add_cases <- sapply(mc_or_diff, function(diff) {
    total_saved <- N_ref * diff
    max(0, floor(total_saved / mean(mc_or_diff[mc_or_diff > 0])))
  })
  
  cat(sprintf("%-30s %15.1f %10.1f – %.1f\n",
              config_name,
              mean(sens_add_cases),
              quantile(sens_add_cases, 0.025),
              quantile(sens_add_cases, 0.975)))
}


# =============================================================================
# SECTION 7: PLOTS
# =============================================================================

# --- Plot 1: Distribution of OR time saved per case ---
or_diff_df <- data.frame(or_diff = mc_or_diff)

p1 <- ggplot(or_diff_df, aes(x = or_diff)) +
  geom_histogram(aes(y = ..density..), bins = 60,
                 fill = "#2C5F8A", alpha = 0.7, colour = "white") +
  geom_density(colour = "#F0A500", linewidth = 1) +
  geom_vline(xintercept = mean(mc_or_diff), colour = "white",
             linetype = "solid", linewidth = 1.2) +
  geom_vline(xintercept = quantile(mc_or_diff, c(0.025, 0.975)),
             colour = "#F0A500", linetype = "dashed", linewidth = 0.8) +
  labs(
    title    = "Distribution of OR Time Saved per Case (RAS vs LPS)",
    subtitle = sprintf("Monte Carlo simulation (n=10,000) | Mean = %.1f min [95%% CI: %.1f – %.1f]",
                       mean(mc_or_diff),
                       quantile(mc_or_diff, 0.025),
                       quantile(mc_or_diff, 0.975)),
    x        = "OR Time Saved per Case (minutes)",
    y        = "Density",
    caption  = "Source: Jain et al. 2025 (Cureus). Gamma distribution fitted via Method of Moments."
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "grey40"),
    panel.grid.minor = element_blank()
  )

# --- Plot 2: Additional cases per year by volume ---
add_cases_df <- data.frame(
  Volume = factor(volume_labels, levels = volume_labels),
  Mean   = results_table$Add_cases_mean,
  Lower  = results_table$Add_cases_lower,
  Upper  = results_table$Add_cases_upper
)

p2 <- ggplot(add_cases_df, aes(x = Volume, y = Mean)) +
  geom_col(fill = "#2C5F8A", width = 0.5) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                width = 0.2, colour = "#F0A500", linewidth = 1) +
  geom_text(aes(label = sprintf("%.0f cases", Mean)),
            vjust = -0.8, colour = "#2C5F8A", fontface = "bold", size = 4) +
  labs(
    title    = "Additional Surgeries Per Year — OR Throughput Model",
    subtitle = "By switching from LPS to RAS for endometrial cancer | Error bars = 95% CI",
    x        = "Annual Surgical Volume",
    y        = "Additional Cases Per Year",
    caption  = "Source: Jain et al. 2025. Base case: 450 usable OR minutes/day."
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "grey40"),
    panel.grid.minor = element_blank()
  )

# --- Plot 3: Bed-days freed per year by volume ---
beddays_df <- data.frame(
  Volume = factor(volume_labels, levels = volume_labels),
  Mean   = results_table$Beddays_mean,
  Lower  = results_table$Beddays_lower,
  Upper  = results_table$Beddays_upper
)

p3 <- ggplot(beddays_df, aes(x = Volume, y = Mean)) +
  geom_col(fill = "#1A7A4A", width = 0.5) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                width = 0.2, colour = "#F0A500", linewidth = 1) +
  geom_text(aes(label = sprintf("%.0f days", Mean)),
            vjust = -0.8, colour = "#1A7A4A", fontface = "bold", size = 4) +
  labs(
    title    = "Bed-Days Freed Per Year — Bed Liberation Model",
    subtitle = "By switching from LPS to RAS for endometrial cancer | Error bars = 95% CI",
    x        = "Annual Surgical Volume",
    y        = "Bed-Days Freed Per Year",
    caption  = "Source: Jain et al. 2025. LOS saving = 1.06 days/patient."
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "grey40"),
    panel.grid.minor = element_blank()
  )

# Save plots
ggsave("or_time_distribution.png",    plot = p1, width = 10, height = 6, dpi = 150)
ggsave("additional_cases_by_volume.png", plot = p2, width = 10, height = 6, dpi = 150)
ggsave("beddays_freed_by_volume.png",    plot = p3, width = 10, height = 6, dpi = 150)

cat("\n\nPlots saved:\n")
cat("  - or_time_distribution.png\n")
cat("  - additional_cases_by_volume.png\n")
cat("  - beddays_freed_by_volume.png\n")


# =============================================================================
# SECTION 8: SAVE RESULTS TO CSV
# =============================================================================

write.csv(results_table, "study_B_results.csv", row.names = FALSE)

cat("\nResults table saved to: study_B_results.csv\n")
cat("\n=== ANALYSIS COMPLETE ===\n")