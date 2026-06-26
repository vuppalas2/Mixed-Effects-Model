# Mixed Effects Model — Longitudinal Cholesterol Analysis in Renal Transplant Recipients

## Background

Cyclosporine A (CsA) is an immunosuppressive agent widely used in renal transplant 
recipients but is associated with hyperlipidemia, significantly raising total 
cholesterol, LDL cholesterol, and apoprotein B levels. In 1994, Tacrolimus (FK506) 
received FDA approval as an alternative immunosuppressant shown to have less adverse 
effect on post-transplant lipid profiles, with studies in liver and kidney transplant 
recipients demonstrating significantly lower triglyceride, cholesterol, and LDL 
cholesterol levels in FK506-treated patients compared to CsA.

This analysis uses data from a randomized controlled trial (RCT) investigating:

1. Whether it is safe to convert stable renal transplant recipients (at least 12 months 
post-transplant, total cholesterol > 240 mg/dL) from CsA to FK506
2. Whether conversion from CsA to FK506 results in a clinically significant reduction 
in total cholesterol, defined as a drop of at least 40 mg/dL

---

## Research Questions

**a.** Do mean total cholesterol levels differ between the FK506 (Treatment) and CsA 
(Control) groups at baseline?

**b.** Do mean total cholesterol levels differ between the FK506 (Treatment) and CsA 
(Control) groups at six months?

**c.** Does the mean change in total cholesterol from baseline to six months differ 
between the FK506 and CsA groups?

---

## Data

- **File:** `Data.csv`
- **Sample size:** 54 renal transplant recipients
- **Groups:** FK506 (Treatment, n = 27) vs CsA (Control, n = 27)
- **Outcome:** Total cholesterol (mg/dL)
- **Time points:** Baseline, Month 1, Month 3, Month 6
- **Design:** Fully balanced, no missing data

---

## Analysis Workflow

| Step | Description |
|------|-------------|
| 1 | Load required R libraries |
| 2 | Import dataset |
| 3 | Data inspection and cleaning — structure, missing values, duplicates, factor conversion |
| 4 | Reshape wide to long format for longitudinal modeling |
| 5 | Data transformation — relabel time points, create numeric time variables |
| 6 | Data validation — confirm balanced design and no missing values |
| 7 | Exploratory data analysis — summary statistics and visualization |
| 8 | Linear mixed model with lmer — confirm model justification via ICC and random effects test |
| 9 | Covariance structure selection with lme — compare VC, CS, AR(1), SP(EXP), UN via AIC/BIC |
| 10 | Estimated marginal means and Tukey-adjusted pairwise comparisons |
| 11 | Visualization of model-adjusted means with 95% confidence intervals |
| 12 | Model diagnostics — residuals vs fitted, Q-Q plot, Shapiro-Wilk test |

---

## Methods

A two-stage linear mixed-effects modeling approach was used:

**Stage 1 — Model Justification (lmer)**  
An initial linear mixed model was fit with treatment group, time, and their interaction 
as fixed effects, and a subject-specific random intercept to account for repeated 
measurements. The ICC and likelihood ratio test confirmed the mixed-effects approach 
was statistically justified (ICC = 0.506, LRT p = 3.3e-15).

**Stage 2 — Covariance Structure Selection (lme)**  
Five candidate within-subject covariance structures were compared using AIC and BIC:

| Structure | Description |
|-----------|-------------|
| VC | Independence — no correlation between time points |
| CS | Compound Symmetry — equal correlation between all time points |
| AR(1) | First-order autoregressive — correlation decays with time (assumes equal spacing) |
| SP(EXP) | Spatial Exponential — correlation decays with time, handles unequal spacing |
| UN | Unstructured — freely estimated correlations and variances |

The spatial exponential structure SP(EXP) was selected as most appropriate given the 
unequal visit spacing (0, 1, 3, 6 months). Estimated marginal means were extracted 
using the emmeans package with Tukey adjustment for multiple comparisons.

---

## Key Findings

- At baseline, cholesterol did not differ significantly between groups (p = 0.934), confirming successful randomization
- FK506 group showed a mean reduction of 71.6 mg/dL from Baseline to Month 1 (p < 0.001), exceeding the pre-specified 40 mg/dL clinical threshold
- Between-group difference at Month 6 was 58.6 mg/dL (p = 0.0001), favoring FK506
- The group x time interaction was highly significant (F = 11.55, p < 0.001)
- Model residuals were normally distributed (Shapiro-Wilk W = 0.996, p = 0.866)

---

## R Packages

| Package | Purpose |
|---------|---------|
| lme4 + lmerTest | Linear mixed model and Satterthwaite F-tests |
| nlme | Flexible covariance structure modeling |
| emmeans | Estimated marginal means and pairwise contrasts |
| performance | ICC calculation |
| ggplot2 | Data visualization |
| reshape2 | Wide to long data transformation |
| dplyr | Data manipulation |
