# Analysis of cell viability across drug concentrations with repeated measures (occasions)
# 3 technical replicates per concentration per occasion
# 4 occasions (different days, same cell line)

# Clear workspace
rm(list = ls())

# ============================================================================
# 1. LOAD REQUIRED LIBRARIES
# ============================================================================

packages <- c("tidyverse", "multcomp", "car", "lme4", "emmeans", "ggplot2", "gridExtra", "clinfun")

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

print(head(data_toxicity, 15))
cat("\nData structure:\n")
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
    SE = sd(viability, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

print(summary_by_conc)

# ============================================================================
# 4. STATISTICAL TESTS
# ============================================================================

cat("\n========== STATISTICAL TESTS ==========\n")

# --- DUNNETT'S TEST ---
cat("\n--- Dunnett's Test: Compare Each Concentration to Control ---\n")
model_mixed <- lme4::lmer(viability ~ concentration + (1 | occasion), data = data_toxicity)
dunnett_results <- emmeans::emmeans(model_mixed, dunnett ~ concentration)
p_val_data_dunnett <- as.data.frame(dunnett_results$contrast)
print(p_val_data_dunnett)

# --- JONCKHEERE TEST FOR DOSE DEPENDENCY ---
cat("\n--- Jonckheere Test: Dose-Dependency Trend Analysis ---\n")

# Create numeric dose levels (ordered)
dose_numeric <- as.numeric(data_toxicity$concentration) - 1  # 0, 1, 2, 3, 4

# Perform Jonckheere-Terpstra test
jonckheere_test <- clinfun::jonckheere.test(
  x = data_toxicity$viability,
  g = dose_numeric,
  alternative = "two.sided"
)

cat("Jonckheere-Terpstra Test Results:\n")
cat("Statistic:", jonckheere_test$statistic, "\n")
cat("P-value:", jonckheere_test$p.value, "\n")

# Create a summary dataframe for visualization
jonckheere_summary <- data.frame(
  Test = "Jonckheere-Terpstra",
  Statistic = jonckheere_test$statistic,
  P_Value = jonckheere_test$p.value,
  Significant = jonckheere_test$p.value < 0.05
)

print(jonckheere_summary)

# ============================================================================
# 5. VISUALIZATION
# ============================================================================

cat("\n========== CREATING VISUALIZATIONS ==========\n")

n_levels <- length(levels(data_toxicity$concentration))

# Fixed y-axis range for all raw value plots
y_axis_limits <- c(50, 180)

# ------- Plot 1: Cell viability by concentration and occasion (Boxplot) -------
plot1 <- ggplot(data_toxicity, aes(x = as.numeric(concentration), y = viability, group = concentration)) +
  geom_vline(xintercept = seq(1.5, n_levels - 0.5, by = 1), color = "grey80", linetype = "dashed") +
  geom_boxplot(aes(fill = occasion, group = interaction(concentration, occasion)), 
               alpha = 0.7, position = position_dodge(width = 0.75)) +
  geom_jitter(aes(color = occasion, group = interaction(concentration, occasion)), 
              position = position_dodge(width = 0.75), alpha = 0.4, size = 2) +
  scale_x_continuous(
    breaks = 1:n_levels,
    labels = levels(data_toxicity$concentration)
  ) +
  scale_y_continuous(limits = y_axis_limits) +
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
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

print(plot1)
ggsave("MCF7A1H1_plot1_viability_by_concentration.pdf", plot1, width = 10, height = 7)
cat("Plot 1 saved: MCF7A1H1_plot1_viability_by_concentration.pdf\n")

# ------- Plot 2: Line plot showing dose-response across occasions -------
plot2_data <- data_toxicity %>%
  group_by(concentration, occasion) %>%
  summarise(Mean = mean(viability), SE = sd(viability) / sqrt(n()), .groups = "drop")

plot2 <- ggplot(plot2_data, aes(x = as.numeric(concentration), y = Mean, color = occasion, group = occasion)) +
  geom_vline(xintercept = seq(1.5, n_levels - 0.5, by = 1), color = "grey80", linetype = "dashed") +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  scale_x_continuous(
    breaks = 1:n_levels,
    labels = levels(data_toxicity$concentration)
  ) +
  scale_y_continuous(limits = y_axis_limits) +
  labs(
    title = "MCF7A1H1 - Mean Viability by Concentration (with SE)",
    x = "Drug Concentration (µg/ml)",
    y = "Mean Viability",
    color = "Occasion"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

print(plot2)
ggsave("MCF7A1H1_plot2_mean_viability.pdf", plot2, width = 10, height = 7)
cat("Plot 2 saved: MCF7A1H1_plot2_mean_viability.pdf\n")

# ------- Plot 3: Dose-response curve (collapsed across occasions) -------
plot3_data <- data_toxicity %>%
  group_by(concentration) %>%
  summarise(Mean = mean(viability), SE = sd(viability) / sqrt(n()), .groups = "drop")

plot3 <- ggplot(plot3_data, aes(x = as.numeric(concentration), y = Mean)) +
  geom_vline(xintercept = seq(1.5, n_levels - 0.5, by = 1), color = "grey80", linetype = "dashed") +
  geom_point(size = 4, color = "steelblue") +
  geom_line(aes(group = 1), color = "steelblue", linewidth = 1) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2, color = "steelblue") +
  scale_x_continuous(
    breaks = 1:n_levels,
    labels = levels(data_toxicity$concentration)
  ) +
  scale_y_continuous(limits = y_axis_limits) +
  labs(
    title = "MCF7A1H1 - Dose-Response Curve (All Occasions Combined)",
    x = "Drug Concentration (µg/ml)",
    y = "Mean Viability"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

print(plot3)
ggsave("MCF7A1H1_plot3_dose_response_curve.pdf", plot3, width = 10, height = 7)
cat("Plot 3 saved: MCF7A1H1_plot3_dose_response_curve.pdf\n")

# ------- Plot 4: P-value Significance Visualization (Dunnett's Test) -------
plot4 <- ggplot(p_val_data_dunnett, aes(x = contrast, y = p.value)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0, ymax = 0.05), 
            fill = "green", alpha = 0.01) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_segment(aes(x = contrast, xend = contrast, y = 0, yend = p.value), color = "grey50") +
  geom_point(aes(color = p.value < 0.05), size = 4) +
  scale_color_manual(values = c("TRUE" = "darkgreen", "FALSE" = "firebrick")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(
    title = "Dunnett's Test: P-Values vs. Control Group",
    x = "Comparison",
    y = "Adjusted P-Value",
    color = "Significant (p < 0.05)"
  ) +
  annotate("text", x = 1, y = 0.08, label = "Significance Threshold (p = 0.05)", 
           color = "red", fontface = "italic", size = 3, adj = 0) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

print(plot4)
ggsave("MCF7A1H1_plot4_dunnett_pvalues.pdf", plot4, width = 10, height = 7)
cat("Plot 4 saved: MCF7A1H1_plot4_dunnett_pvalues.pdf\n")

# ------- Plot 5: Jonckheere Test P-Value Visualization -------
plot5 <- ggplot(jonckheere_summary, aes(x = Test, y = P_Value)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0, ymax = 0.05), 
            fill = "green", alpha = 0.01) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_segment(aes(x = Test, xend = Test, y = 0, yend = P_Value), color = "grey50") +
  geom_point(aes(color = P_Value < 0.05), size = 6) +
  scale_color_manual(values = c("TRUE" = "darkgreen", "FALSE" = "firebrick")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(
    title = "Jonckheere-Terpstra Test: Dose-Dependency Analysis",
    x = "Test",
    y = "P-Value",
    color = "Significant (p < 0.05)"
  ) +
  annotate("text", x = 1.15, y = 0.08, label = "Significance Threshold (p = 0.05)", 
           color = "red", fontface = "italic", size = 3, adj = 0) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

print(plot5)
ggsave("MCF7A1H1_plot5_jonckheere_pvalues.pdf", plot5, width = 10, height = 7)
cat("Plot 5 saved: MCF7A1H1_plot5_jonckheere_pvalues.pdf\n")

# ============================================================================
# 6. SUMMARY REPORT
# ============================================================================

cat("\n========== ANALYSIS SUMMARY ==========\n")
cat("\nCell line: MCF7A1H1\n")
cat("Analysis date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\nDunnett's Test (vs. Control):\n")
print(p_val_data_dunnett)
cat("\nJonckheere-Terpstra Test (Dose-Dependency):\n")
print(jonckheere_summary)
cat("\n========== ANALYSIS COMPLETE ==========\n")
