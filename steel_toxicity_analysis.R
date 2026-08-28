# Steel's (Dunnett's) Test Analysis for Drug Toxicity Data - MCF7A1H1
# Analysis of cell viability across drug concentrations with repeated measures (occasions)
# 3 technical replicates per concentration per occasion
# 4 occasions (different days, same cell line)

# Clear workspace
rm(list = ls())

# ============================================================================
# 1. LOAD REQUIRED LIBRARIES
# ============================================================================

packages <- c("tidyverse", "multcomp", "car", "lme4", "emmeans", "ggplot2", "gridExtra")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# ============================================================================
# 2. LOAD YOUR DATA
# ============================================================================

# Load MCF7A1H1 data from CSV file
data_toxicity <- read.csv("MCF7A1H1.csv")

# Ensure columns are correctly formatted
data_toxicity$occasion <- factor(data_toxicity$occasion)
data_toxicity$replicate <- factor(data_toxicity$replicate)
data_toxicity$concentration <- factor(
  data_toxicity$concentration,
  levels = c("Control", "50", "100", "200", "500")
)

# Display first rows
cat("\n========== DATA LOADED ==========\n")
cat("Cell line: MCF7A1H1\n")
cat("File: MCF7A1H1.csv\n")
cat("Data dimensions:", nrow(data_toxicity), "rows,", ncol(data_toxicity), "columns\n\n")

head(data_toxicity, 15)
print("Data structure:")
str(data_toxicity)

# ============================================================================
# 3. DESCRIPTIVE STATISTICS
# ============================================================================

cat("\n========== DESCRIPTIVE STATISTICS ==========\n")

# Summary by concentration
summary_by_conc <- data_toxicity %>%
  group_by(concentration) %>%
  summarise(
    N = n(),
    Mean = mean(viability, na.rm = TRUE),
    SD = sd(viability, na.rm = TRUE),
    Min = min(viability, na.rm = TRUE),
    Max = max(viability, na.rm = TRUE),
    SE = SD / sqrt(N),
    .groups = "drop"
  )

print("Summary by Concentration:")
print(summary_by_conc)

# Summary by occasion
summary_by_occasion <- data_toxicity %>%
  group_by(occasion) %>%
  summarise(
    N = n(),
    Mean = mean(viability, na.rm = TRUE),
    SD = sd(viability, na.rm = TRUE),
    .groups = "drop"
  )

print("\nSummary by Occasion:")
print(summary_by_occasion)

# Summary by concentration and occasion
summary_by_both <- data_toxicity %>%
  group_by(concentration, occasion) %>%
  summarise(
    N = n(),
    Mean = mean(viability, na.rm = TRUE),
    SD = sd(viability, na.rm = TRUE),
    .groups = "drop"
  )

print("\nSummary by Concentration and Occasion:")
print(summary_by_both)

# ============================================================================
# 4. VISUALIZATION
# ============================================================================

cat("\n========== CREATING VISUALIZATIONS ==========\n")

# Plot 1: Cell viability by concentration and occasion (single plot with occasion color-coded)
plot1 <- ggplot(data_toxicity, aes(x = concentration, y = viability, fill = occasion)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.75)) +
  geom_jitter(aes(color = occasion), position = position_dodge(width = 0.75), alpha = 0.4, size = 2) +
  labs(
    title = "MCF7A1H1 - Cell Viability by Concentration and Occasion",
    x = "Drug Concentration (µg/ml)",
    y = "Viability",
    fill = "Occasion",
    color = "Occasion"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

# Plot 2: Line plot showing dose-response across occasions
plot2 <- data_toxicity %>%
  group_by(concentration, occasion) %>%
  summarise(Mean = mean(viability), SE = sd(viability) / sqrt(n()), .groups = "drop") %>%
  ggplot(aes(x = concentration, y = Mean, color = occasion, group = occasion)) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  labs(
    title = "MCF7A1H1 - Mean Viability by Concentration (with SE)",
    x = "Drug Concentration (µg/ml)",
    y = "Mean Viability",
    color = "Occasion"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot 3: Dose-response curve (collapsed across occasions)
plot3 <- data_toxicity %>%
  group_by(concentration) %>%
  summarise(Mean = mean(viability), SE = sd(viability) / sqrt(n()), .groups = "drop") %>%
  ggplot(aes(x = concentration, y = Mean)) +
  geom_point(size = 4, color = "steelblue") +
  geom_line(aes(group = 1), color = "steelblue", size = 1) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2, color = "steelblue") +
  labs(
    title = "MCF7A1H1 - Dose-Response Curve (All Occasions Combined)",
    x = "Drug Concentration (µg/ml)",
    y = "Mean Viability"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine and display plots
combined_plots <- gridExtra::grid.arrange(plot1, plot2, plot3, nrow = 2)
print(combined_plots)

# Save plots
ggsave("MCF7A1H1_toxicity_analysis_plots.pdf", combined_plots, width = 14, height = 10)
cat("Plots saved to: MCF7A1H1_toxicity_analysis_plots.pdf\n")

# ============================================================================
# 5. LINEAR MIXED EFFECTS MODEL (for repeated measures)
# ============================================================================

cat("\n========== LINEAR MIXED EFFECTS MODEL ==========\n")

# Fit mixed model with random intercept for occasion (repeated measures)
# Fixed effects: concentration
# Random effects: occasion (accounts for between-day variation)

model_lme <- lme4::lmer(
  viability ~ concentration + (1 | occasion),
  data = data_toxicity,
  REML = TRUE
)

print(summary(model_lme))

# Model diagnostics
cat("\n--- Model Diagnostics ---\n")
print("Residual plot:")
plot(model_lme)

# Check assumptions
residuals_lme <- residuals(model_lme)
shapiro_test <- shapiro.test(residuals_lme)
cat("\nShapiro-Wilk test for normality of residuals:\n")
print(shapiro_test)

# Variance homogeneity test (Levene's test)
data_toxicity$fitted_values <- predict(model_lme)
levene_test <- car::leveneTest(residuals_lme ~ data_toxicity$concentration)
cat("\nLevene's test for homogeneity of variance:\n")
print(levene_test)

# ============================================================================
# 6. STEEL'S TEST (Dunnett's Test) - Compare Each Concentration to Control
# ============================================================================

cat("\n========== STEEL'S TEST (DUNNETT'S TEST) ==========\n")
cat("Testing whether each concentration differs significantly from Control\n")

# Using emmeans for contrasts
emmeans_fit <- emmeans::emmeans(model_lme, ~concentration)

# Dunnett test: compare all treatments to control
dunnett_contrasts <- emmeans::contrast(emmeans_fit, method = "dunnett", ref = 1)
print(summary(dunnett_contrasts, adjust = "dunnett"))

# Convert to a more readable format
dunnett_results <- summary(dunnett_contrasts, adjust = "dunnett") %>%
  as.data.frame() %>%
  mutate(
    Significant = ifelse(p.value < 0.05, "***", ""),
    p.value = round(p.value, 4)
  )

cat("\n--- Dunnett's Test Results ---\n")
cat("Null hypothesis: Each concentration has the same viability as Control\n")
cat("Adjusted p-values (Dunnett correction):\n\n")
print(dunnett_results)

# ============================================================================
# 7. TREND TEST (Linear Contrast for Dose-Response)
# ============================================================================

cat("\n========== TREND TEST (Linear Dose-Response) ==========\n")

# Create numeric dose variable (assuming ordered concentrations)
data_toxicity_dose <- data_toxicity %>%
  mutate(dose = case_when(
    concentration == "Control" ~ 0,
    concentration == "50" ~ 1,
    concentration == "100" ~ 2,
    concentration == "200" ~ 3,
    concentration == "500" ~ 4
  ))

# Fit model with dose as continuous variable
model_dose <- lme4::lmer(
  viability ~ dose + (1 | occasion),
  data = data_toxicity_dose,
  REML = TRUE
)

print(summary(model_dose))

# Extract dose effect
dose_coef <- fixef(model_dose)["dose"]
cat(sprintf("\nDose coefficient: %.4f\n", dose_coef))
cat("Interpretation: For each unit increase in concentration,\n")
cat(sprintf("viability changes by approximately %.2f units\n", dose_coef))

# Test if trend is significant
cat("\nLinear trend test (from dose model):\n")
trend_summary <- summary(model_dose)
print(trend_summary)

# ============================================================================
# 8. COMPREHENSIVE RESULTS TABLE
# ============================================================================

cat("\n========== COMPREHENSIVE RESULTS TABLE ==========\n")

results_table <- as.data.frame(emmeans_fit) %>%
  mutate(
    Mean_Viability = round(emmean, 2),
    SE_Value = round(SE, 2),
    CI_Lower = round(lower.CL, 2),
    CI_Upper = round(upper.CL, 2)
  ) %>%
  dplyr::select(concentration, Mean_Viability, SE_Value, CI_Lower, CI_Upper)

print(results_table)

# ============================================================================
# 9. SAVE RESULTS
# ============================================================================

cat("\n========== SAVING RESULTS ==========\n")

# Save results table
write.csv(results_table, "MCF7A1H1_toxicity_results_table.csv", row.names = FALSE)
write.csv(dunnett_results, "MCF7A1H1_dunnett_test_results.csv", row.names = FALSE)
write.csv(summary_by_conc, "MCF7A1H1_summary_statistics.csv", row.names = FALSE)

cat("Results saved:\n")
cat("  - MCF7A1H1_toxicity_results_table.csv\n")
cat("  - MCF7A1H1_dunnett_test_results.csv\n")
cat("  - MCF7A1H1_summary_statistics.csv\n")

# ============================================================================
# 10. SUMMARY AND INTERPRETATION
# ============================================================================

cat("\n========== SUMMARY AND INTERPRETATION (MCF7A1H1) ==========\n")
cat("\n1. DOSE-RESPONSE RELATIONSHIP:\n")
if (abs(dose_coef) > 0.1 && summary(model_dose)$coefficients["dose", "Pr(>|t|)"] < 0.05) {
  cat("   ✓ Significant linear trend detected\n")
  cat("   ✓ Cell viability significantly changes with drug concentration\n")
} else {
  cat("   ✗ No significant linear dose-response trend detected\n")
}

cat("\n2. CONCENTRATION-SPECIFIC EFFECTS (vs. Control):\n")
sig_contrasts <- dunnett_results %>% filter(p.value < 0.05)
if (nrow(sig_contrasts) > 0) {
  cat(sprintf("   ✓ %d concentration(s) show significant difference from control:\n", nrow(sig_contrasts)))
  for (i in 1:nrow(sig_contrasts)) {
    cat(sprintf("     - %s (p = %.4f)\n", sig_contrasts$contrast[i], sig_contrasts$p.value[i]))
  }
} else {
  cat("   ✗ No concentrations show significant difference compared to control\n")
}

cat("\n3. REPEATED MEASURES EFFECT (Occasion variability):\n")
random_effect_var <- as.numeric(VarCorr(model_lme)$occasion)
residual_var <- attr(VarCorr(model_lme), "sc")^2
cat(sprintf("   - Between-occasion variance: %.2f\n", random_effect_var))
cat(sprintf("   - Residual variance: %.2f\n", residual_var))
cat(sprintf("   - Ratio: %.2f%%\n", 100 * random_effect_var / (random_effect_var + residual_var)))

cat("\n========== ANALYSIS COMPLETE ==========\n")
cat("Cell line: MCF7A1H1\n")
cat("Analysis date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
