# Conformal Sensitivity Analysis for Individual Treatment Effects

This project applies a **Conformal Sensitivity Analysis (CSA)** framework to construct valid prediction intervals for Individual Treatment Effects (ITE) in clinical observational datasets. It compares three outcome models within the CSA framework and evaluates their performance via simulation and real-data experiments.

## Problem Setting

In observational studies, standard conformal inference for ITE requires the unconfoundedness assumption. This project relaxes that via the **Marginal Sensitivity Model** (Rosenbaum): the constructed intervals remain valid as long as the odds ratio of unmeasured selection bias is bounded by Γ. Larger Γ means more conservative (wider) intervals.

The key question for each test unit: does the ITE interval contain zero — i.e., is the treatment **effective** under the assumed confounding level Γ?

## Methods Compared

| Method | Description |
|---|---|
| **CSA-Mean** | Nested CSA with mean regression (Random Forest) — best overall interval precision in simulation |
| **CSA-Quantile** | Nested CSA with conformalized quantile regression (CQR, Quantile RF) |
| **CSA-Huber** | Nested CSA with gradient boosting under Huber loss — robust to heavy-tailed errors |
| **ITE-NUC** | Inexact nested approach from `cfcausal` assuming no unmeasured confounding (Γ=1 baseline) |
| **ITE-Bonf.** | Bonferroni/naive: union of two single-arm sensitivity intervals |

All nested CSA methods follow a two-group pipeline: fit outcome models on group 1, calibrate sensitivity-adjusted conformal scores on group 2, then smooth the resulting ITE bounds via quantile regression on the test set.

## Experiments

### Simulation Study (`run_tests/syn/`, `results/synthetic/`)

Real datasets (VK, VD, Glutathione) provide the **covariate distribution** only. Synthetic outcomes are generated from known `taufun` / `taufun0` with three error distributions (normal, logistic, polluted normal) across Γ ∈ {1, 2, 3}. This gives a ground-truth ITE for evaluating:
- **Coverage**: fraction of test units whose true ITE falls in the predicted interval
- **Interval length**: efficiency of the intervals across methods and Γ values

Simulation result: among methods with satisfactory predictive accuracy, **CSA-Mean achieves the highest interval precision**.

### Real Data Application (`run_tests/`, `results/ITE/`)

Applies CSA-Mean to three clinical datasets. Produces sensitivity intervals as Γ varies, and classifies test units as "effectively treated" if 0 is excluded from the interval.

| Dataset | Clinical Context | Outcome | Finding |
|---|---|---|---|
| **VK** | Vitamin K2 supplementation in hemodialysis patients (Elshinnawy 2023) | MGP (vascular calcification marker) | K2 significantly alleviates vascular calcification; robust across Γ |
| **VD** | Vitamin D supplementation in chronic hepatitis C after DAA treatment (Sriphoosanaphan 2021) | TGF-β1 (liver fibrosis marker) | No significant effect; robust across Γ |
| **Glutathione** | Sublingual glutathione in COPD patients (Farag 2023) | SOD3 (antioxidant marker) | Positive effects only for a subset of patients; high individual variability, poor robustness |

## Datasets

| Dataset | File | n | Usage |
|---|---|---|---|
| VK2 | `data/VK2.csv` | 120 | Simulation (covariate distribution) + real data application |
| VD | `data/VD.csv` | 75 | Simulation (covariate distribution) + real data application |
| Glutathione | `data/data1.xlsx` | 50 | Simulation (covariate distribution) + real data application |

In simulation, covariates are fitted with a multivariate normal and used to generate synthetic samples; real outcome data is not used.

## Repository Structure

```
R/                          # Core algorithm implementations
│   ├── nested_conformalSA.R    # Two-group nested CSA fit
│   ├── conformal_SA.R          # Single-arm conformal with sensitivity weighting
│   ├── cutoff_SA.R             # Sensitivity-adjusted conformal cutoff under Γ
│   ├── fit_and_predict_band.R  # Smooth & extrapolate ITE bounds to test set
│   ├── predict.nested.R        # Get ITE bands from a nested_conformalSA object
│   ├── conformal_learners.R    # RF, Quantile RF, Huber Boosting wrappers
│   ├── conformalIte.R          # cfcausal ITE interface (nest/naive/counterfactual)
│   └── util_SA.R               # Simulation utilities (samplecf, summary_CI, ...)
data/                       # Observational datasets (CSV/xlsx)
run_tests/                  # Experiment entry scripts
│   ├── VD.r / vk.R / data1.R  # Real data experiments (per dataset)
│   └── syn/                    # Simulation experiments
│       ├── VD-syn.R / VD-huber.R
│       ├── VK-syn.R / VK-huber.R
│       └── data1-syn.R / data1-huber.R
results/                    # Output CSVs (gitignored)
│   ├── ITE/                # Per-dataset real-data ITE results
│   └── synthetic/          # Simulation coverage & length results
figures/                    # Output figures (gitignored)
plot_figures/               # Plotting scripts
ref/                        # Reference experiment scripts (exp-cf, exp-ite, exp-fish)
```

## Quickstart

All scripts run from the **repository root** (the R package root).

**Install dependencies:**
```r
install.packages(c("devtools", "randomForest", "grf", "gbm", "h2o",
                   "MASS", "argparse", "dplyr", "ggplot2", "readxl"))
devtools::install_github("lihualei71/cfcausal")
```

**Run simulation experiment (VK covariate distribution):**
```bash
Rscript run_tests/syn/VK-syn.R \
  --alpha 0.2 \
  --cftype 2 \
  --ntrial 50 \
  --ntrain 2000 \
  --ntest 5000 \
  --seed 1
```

**Run real data experiment:**
```bash
Rscript run_tests/vk.R \
  --data_name VK \
  --method mean \
  --gmm_star 3 \
  --alpha 0.2 \
  --ntrial 10 \
  --seed1 123 --seed2 45
```

**Key CLI arguments:**

| Argument | Description |
|---|---|
| `--gmm_star` | Sensitivity parameter Γ (≥ 1; Γ=1 means no unmeasured confounding assumed) |
| `--method` | Inference type: `mean` or `cqr` |
| `--alpha` | Miscoverage level (produces 1−α prediction intervals) |
| `--cftype` | Confounding type for simulation (1–4, see `util_SA.R`) |
| `--ntrain` / `--ntest` | Training / testing sample sizes for simulation |
| `--huber_alpha` | Huber loss quantile parameter ∈ [0,1] |

Results are saved to `results/synthetic/<dataset>/<errdist>/gmm_<Γ>/` for simulation, and `results/ITE/<data_name>-<method>/alpha_<α>_gmm_<Γ>/` for real data.

## References

- Lei, L. & Candès, E. (2021). Conformal inference of counterfactuals and individual treatment effects. *JRSS-B*.
- Yin, M. et al. (2024). Conformal sensitivity analysis for individual treatment effects. *JASA*, 119(545): 122–135.
- Jin, Y. et al. (2023). Sensitivity analysis of individual treatment effects: A robust conformal inference approach. *PNAS*, 120(6).
- Elshinnawy, H. A. et al. (2023). Different vitamin K forms in hemodialysis patients. *Egyptian Journal of Internal Medicine*, 35(1): 4.
- Sriphoosanaphan, S. et al. (2021). Effect of vitamin D supplementation in patients with chronic hepatitis C. *PeerJ*, 9: e10709.
- Farag, A. et al. (2023). Evaluation of the antioxidant and anti-inflammatory effect of sublingual glutathione on COPD patients. *Journal of Medicine and Life*, 16(12): 1796–1801.
