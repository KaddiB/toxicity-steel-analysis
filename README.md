# Steel's Test Analysis for Drug Toxicity Data

This repository contains an R script for analyzing drug toxicity data using **Steel's test (Dunnett's test)** with **repeated measures** design.

## Study Design

- **4 drug concentrations** tested against a **control**
- **3 technical replicates** per concentration per test occasion
- **4 test occasions** (different days, same cell line → quasi-biological replicates)
- **Outcome**: Cell viability (%)
- **Research question**: Does viability decrease with increasing drug concentration?

## Analysis Overview

The script performs a complete statistical analysis including:

### 1. **Descriptive Statistics**
   - Summary statistics by concentration and occasion
   - Grouped statistics showing variability across replicates and occasions

### 2. **Visualization**
   - Box plots with individual data points (by concentration and occasion)
   - Line plots showing dose-response trends across occasions
   - Combined dose-response curve with error bars
   - Automatically saved as PDF

### 3. **Linear Mixed Effects Model**
   - Accounts for **repeated measures** (occasion as random effect)
   - Concentration as fixed effect
   - Model diagnostics: residual plots, normality tests, homogeneity of variance

### 4. **Steel's Test (Dunnett's Test)**
   - Compares each concentration directly to **Control**
   - Adjusted p-values (Dunnett correction for multiple comparisons)
   - Tests: Is viability significantly reduced at each concentration vs. control?

### 5. **Trend Test (Linear Dose-Response)**
   - Tests for a significant **linear relationship** between concentration and viability
   - Quantifies the effect size: how much viability decreases per concentration unit
   - Answers: Is there a statistically significant dose-dependent effect?

### 6. **Results Summary**
   - Comprehensive results table with means, standard errors, and confidence intervals
   - Automated interpretation of findings

## How to Use

### Step 1: Prepare Your Data

Replace the sample data in **Section 2** with your actual data. Your data should have columns:
- `occasion`: Factor (e.g., "Day_1", "Day_2", etc.)
- `replicate`: Factor (e.g., "Rep_1", "Rep_2", "Rep_3")
- `concentration`: Factor with levels: "Control", "Conc_1", "Conc_2", "Conc_3", "Conc_4"
- `viability`: Numeric (% viability, 0-100)

**Example data structure:**
```
  occasion replicate concentration viability
1    Day_1      Rep_1       Control      96.2
2    Day_1      Rep_2       Control      97.8
3    Day_1      Rep_3       Control      98.1
4    Day_1      Rep_1       Conc_1      92.5
...
```

**To import from CSV:**
```r
data_toxicity <- read.csv("your_data.csv")
```

### Step 2: Run the Script

```r
source("steel_toxicity_analysis.R")
```

Or run it line-by-line in RStudio.

### Step 3: Interpret Results

The script outputs:

1. **Steel's Test (Dunnett's Test) Results:**
   - Compares each concentration to Control
   - `p.value < 0.05` → Significant difference from control
   - `estimate` → Mean difference in viability (%)

2. **Linear Trend Test Results:**
   - `dose coefficient` → Change in viability per concentration unit
   - `p.value < 0.05` → Significant dose-response relationship

3. **Summary Interpretation Section:**
   - ✓ or ✗ indicating whether dose-response is significant
   - Which concentrations differ significantly from control
   - Repeated measures effect size

## Key Output Files

- **`toxicity_analysis_plots.pdf`** — Visualization of all results
- **`toxicity_results_table.csv`** — Mean viability, SE, confidence intervals by concentration
- **`dunnett_test_results.csv`** — Steel's test results with p-values
- **`summary_statistics.csv`** — Descriptive statistics by concentration

## Statistical Details

### Why Linear Mixed Effects Model?

Your data has a **repeated measures structure**:
- Same cell line tested on 4 different days
- This introduces correlation within days (not independent observations)
- A standard ANOVA would violate independence assumptions
- Mixed effects model accounts for this with a random intercept by occasion

### Why Steel's Test?

Steel's test (Dunnett's test) is ideal for:
- **Multiple comparisons to a control group** (not pairwise comparisons)
- **Correction for multiple testing** (maintains family-wise error rate at α = 0.05)
- **Planned comparisons** (you're interested in vs. control, not all pairwise)

### Why Trend Test?

Tests the specific hypothesis: **"Viability decreases linearly with dose"**
- More powerful than testing individual concentrations separately
- Quantifies dose-response relationship
- Appropriate for ordered categorical doses

## Required R Packages

The script automatically installs (if needed) and loads:
- `tidyverse` — Data manipulation and visualization
- `multcomp` — Multiple comparison tests
- `car` — Statistical tests (Levene's test)
- `lme4` — Linear mixed effects models
- `emmeans` — Estimated marginal means and contrasts
- `ggplot2` — Advanced visualization
- `gridExtra` — Combine plots

## Customization

### Change Control Group
In Section 6, modify the `ref` argument:
```r
dunnett_contrasts <- contrast(emmeans_fit, method = "dunnett", ref = 1)
# ref = 1 means first level is control (default)
```

### Adjust Significance Level
Change `0.05` to your desired α-level throughout:
```r
if (p.value < 0.01) { # Use 0.01 for α = 0.01
```

### Add Covariates
Extend the model with additional fixed or random effects:
```r
# Example: if you have a "cell_batch" covariate
model_lme <- lmer(viability ~ concentration + cell_batch + (1 | occasion), data = data)
```

### Transformation
If your data violates normality, apply a transformation:
```r
data_toxicity$viability_log <- log(data_toxicity$viability)
# Then use viability_log in the model
```

## Troubleshooting

**Issue: "Error in lme4::lmer() ... REML=TRUE"**
- Update R and all packages: `update.packages()`

**Issue: "No contrasts found"**
- Ensure `concentration` is a factor with the correct levels
- Check: `levels(data_toxicity$concentration)`

**Issue: "Singular fit"**
- Random effects variance is very small; try removing random intercept:
```r
model_lme <- lm(viability ~ concentration, data = data_toxicity)
```

**Issue: Plots not saving**
- Check your working directory: `getwd()`
- Ensure write permissions in that directory

## Contact & Questions

For questions about the analysis or adjustments needed for your specific data, feel free to open an issue in this repository.

---

**References:**
- Dunnett, C. W. (1955). A multiple comparison procedure for comparing several treatments with a control. Journal of the American Statistical Association, 50(271), 1096–1121.
- Pinheiro, J., & Bates, D. (2000). Mixed-effects models in S and S-PLUS. Springer.
- R Core Team (2023). R: A language and environment for statistical computing.
