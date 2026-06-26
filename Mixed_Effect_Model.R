############################################
########    MIXED EFFECT MODEL    ##########
############################################


## Step 1: Load Required Libraries
library(dplyr)       # Data manipulation
library(reshape2)    # Wide to long transformation
library(lme4)        # Linear mixed model (lmer) — used to confirm model justification
library(lmerTest)    # F-tests with Satterthwaite degrees of freedom for lmer
library(nlme)        # Linear mixed model (lme) — used for covariance structure selection
library(emmeans)     # Estimated marginal means and post-hoc comparisons
library(performance) # ICC calculation
library(ggplot2)     # Data visualization


## Step 2: Import Data
setwd("/Users/sushmavuppala/Documents/RCT_Project")

data <- read.csv("RCT_Data.csv")

head(data)


## Step 3: Data Inspection and Cleaning
# Review structure and distributions
str(data)
summary(data)

# Check for missing values and duplicate records
sapply(data, function(x) sum(is.na(x)))
sum(duplicated(data))

# Convert participant identifier and treatment group to categorical factors
data$id    <- factor(data$id)
data$group <- factor(data$group,
                     levels = c("1", "2"),
                     labels = c("FK506 (Treatment)", "CsA (Control)"))

str(data)


## Step 4: Reshape to Long Format
# Convert from wide format (one row per subject) to long format (one row per observation)
# This is required for longitudinal mixed-effects modeling
long_data <- melt(
  data,
  id.vars       = c("id", "group"),
  measure.vars  = c("chol0", "chol1", "chol3", "chol6"),
  variable.name = "time",
  value.name    = "cholesterol"
)

head(long_data)


## Step 5: Data Transformation 
# Relabel time points with clinically meaningful labels
long_data$time <- factor(
  long_data$time,
  levels = c("chol0", "chol1", "chol3", "chol6"),
  labels = c("Baseline", "Month 1", "Month 3", "Month 6")
)

# Create numeric time variable using true month values (0, 1, 3, 6)
# This correctly reflects the unequal spacing between visits
# Used by SP(EXP) and AR(1) covariance structures in Step 7
long_data$time_numeric <- ifelse(long_data$time == "Baseline", 0,
                                 ifelse(long_data$time == "Month 1",  1,
                                        ifelse(long_data$time == "Month 3",  3, 6)))

# Create consecutive integer time variable (1, 2, 3, 4)
# Required specifically by corSymm (Unstructured) which needs consecutive integers
long_data$time_int <- as.integer(long_data$time)



## Step 6: Data Validation
# Verify final structure
str(long_data)

# Confirm no missing values after transformation
sum(is.na(long_data$cholesterol))

# Confirm balanced design — equal observations per group and time point
table(long_data$group, long_data$time)


## Step 7: Exploratory Data Analysis
# Summary statistics with 95% confidence intervals for each group × time combination
eda_summary <- long_data %>%
  group_by(group, time) %>%
  summarise(
    n        = n(),
    mean     = mean(cholesterol, na.rm = TRUE),
    sd       = sd(cholesterol, na.rm = TRUE),
    se       = sd / sqrt(n),
    ci_lower = mean - 1.96 * se,
    ci_upper = mean + 1.96 * se,
    .groups  = "drop"
  )

eda_summary

# Visualize mean cholesterol trajectories over time by treatment group
ggplot(long_data, aes(x = time, y = cholesterol, color = group, group = group)) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2) +
  scale_color_manual(values = c("FK506 (Treatment)" = "#1B6CA8",
                                "CsA (Control)"     = "#C43C3C")) +
  labs(
    title = "Mean Total Cholesterol Over Time by Treatment Group",
    x     = "Time Point",
    y     = "Total Cholesterol (mg/dL)",
    color = "Treatment Group"
  ) +
  theme_minimal() +
  theme(
    plot.title      = element_text(hjust = 0.5),
    legend.position = "bottom"
  )


## Step 8: Linear Mixed Model — Confirming Model Justification (lmer)
model_lmer <- lmer(cholesterol ~ group * time + (1 | id), data = long_data)

summary(model_lmer)


# Type III F-tests for fixed effects using Satterthwaite degrees of freedom
anova(model_lmer)

# Intraclass Correlation Coefficient (ICC)
# Quantifies the proportion of total variance attributable to between-subject differences
# ICC > 0.1 generally justifies use of a mixed model
icc(model_lmer)

# Likelihood ratio test for the random intercept
# Compares model with vs. without random intercept to confirm it is needed
rand(model_lmer)

## Step 9: Covariance Structure Selection (lme)

# VC: Variance Components
# Assumes complete independence between time points — no within-subject correlation
model_VC <- lme(cholesterol ~ group * time,
                random = ~1 | id,
                data = long_data)

# CS: Compound Symmetry
# Assumes all pairs of time points share the same correlation
# Equivalent to the lmer assumption in Step 7
model_CS <- lme(cholesterol ~ group * time,
                random = ~1 | id,
                correlation = corCompSymm(form = ~time_numeric | id),
                data = long_data)

# AR(1): First-Order Autoregressive
# Correlation decays exponentially with time — closer visits more correlated
# Assumes equal spacing between visits, which does not hold here (0,1,3,6 months)
model_AR <- lme(cholesterol ~ group * time,
                random = ~1 | id,
                correlation = corAR1(form = ~time_numeric | id),
                data = long_data)

# SP(EXP): Spatial Exponential
# Correlation decays with time and correctly accounts for unequal visit spacing
# Most appropriate structure given visits at months 0, 1, 3, and 6
model_SP <- lme(cholesterol ~ group * time,
                random = ~1 | id,
                correlation = corSpatial(form = ~time_numeric | id, type = "exp"),
                data = long_data)

# UN: Unstructured
# Each pair of time points has its own freely estimated correlation
# Most flexible but most parameter-intensive
# Note: corSymm requires consecutive integers (1,2,3,4) — uses time_int not time_numeric
model_UN <- lme(cholesterol ~ group * time,
                random = ~1 | id,
                correlation = corSymm(form = ~time_int | id),
                weights = varIdent(form = ~1 | time),
                data = long_data)

# Compare all five structures using AIC and BIC
# Lower values indicate better fit, penalized for model complexity
model_comparison <- data.frame(
  Model       = c("VC", "CS", "AR(1)", "SP(EXP)", "UN"),
  Description = c("Independence", "Equal correlations",
                  "Decaying correlation (equal spacing)",
                  "Decaying correlation (unequal spacing)",
                  "Freely estimated correlations"),
  AIC = c(AIC(model_VC), AIC(model_CS), AIC(model_AR),
          AIC(model_SP), AIC(model_UN)),
  BIC = c(BIC(model_VC), BIC(model_CS), BIC(model_AR),
          BIC(model_SP), BIC(model_UN))
)

model_comparison

# Final model selected based on AIC/BIC and clinical appropriateness
final_model <- model_SP

summary(final_model)
anova(final_model)

## Step 10: Estimated Marginal Means and Post-Hoc Comparisons
# Extract LS Means from final model
emm <- emmeans(final_model, ~ group * time)
emm

# Tukey-adjusted pairwise comparisons across all group × time combinations
contrast(emm, "pairwise", adjust = "tukey")

# 95% confidence intervals for all pairwise contrasts
confint(contrast(emm, "pairwise", adjust = "tukey"))


## Step 11: Visualization of Estimated Marginal Means
emm_df <- as.data.frame(emm)

ggplot(emm_df, aes(x = time, y = emmean, color = group, group = group)) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  scale_color_manual(values = c("FK506 (Treatment)" = "brown",
                                "CsA (Control)"     = "darkgreen")) +
  labs(
    title = "Model-Adjusted Mean Total Cholesterol by Group and Time",
    x     = "Time Point",
    y     = "LS Mean Total Cholesterol (mg/dL)",
    color = "Treatment Group"
  ) +
  theme_minimal() +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )


## Step 12: Model Diagnostics

residuals_final <- residuals(final_model)
fitted_final    <- fitted(final_model)

# Residuals vs Fitted Values
# Checks the assumption of constant variance (homoscedasticity)
# Residuals should be randomly scattered around zero with no pattern
plot(fitted_final, residuals_final,
     main = "Residuals vs Fitted Values",
     xlab = "Fitted Values",
     ylab = "Residuals",
     pch  = 16, col = "steelblue")
abline(h = 0, col = "red", lwd = 2)
grid(col = "gray", lty = "dotted")

# Q-Q Plot of Residuals
# Checks the normality assumption — points should follow the diagonal reference line
qqnorm(residuals_final,
       main = "Q-Q Plot of Residuals",
       pch  = 16, col = "black", cex = 0.7)
qqline(residuals_final, col = "red", lwd = 2)
grid(col = "gray", lty = "dotted")

# Shapiro-Wilk Test
# Formal test of normality — p > 0.05 indicates no significant departure from normality
shapiro.test(residuals_final)