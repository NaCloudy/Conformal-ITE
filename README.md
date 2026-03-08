# Conformal Sensitivity Analysis for Individual Treatment Effects

This project proposes and evaluates **Conformal Sensitivity Analysis (CSA)** methods for constructing valid prediction intervals of Individual Treatment Effects (ITE) under unmeasured confounding. The key contribution is a **CSA-Huber** variant that uses Huber regression for robustness, compared against several baselines.

## Problem Setting

In observational studies, standard conformal inference for ITE requires unconfoundedness. This project relaxes that assumption via a **Γ-sensitivity model**: the constructed intervals remain valid as long as the odds ratio of unmeasured selection bias is bounded by Γ. Larger Γ means more conservative (wider) intervals.

The key question for each test unit: does the ITE interval contain zero — i.e., is the treatment **effective** under the assumed confounding level Γ?

## Methods Compared

| Method | Description |
|---|---|
| **CSA-Huber** | Nested CSA with Huber regression outcome model (robust to heavy tails/outliers) — *main contribution* |
| **CSA-M** | Nested CSA with conditional mean regression (Random Forest) |
| **CSA-Q** | Nested CSA with conformalized quantile regression (CQR) |
| **ITE-NUC** | Inexact nested approach from cfcausal assuming no unmeasured confounding |
| **CSA-B** | Bonferroni/naive: union of two single-arm sensitivity intervals |

All nested CSA methods follow a two-group pipeline: fit outcome models on group 1, calibrate sensitivity-adjusted conformal scores on group 2, then smooth the resulting ITE bounds via quantile regression.

## Experiments

### Simulation Study (`run_tests/syn/`, `results/synthetic/`)

Real datasets (VD, VK) provide the **covariate distribution** only. Synthetic outcomes are generated from known `taufun` / `taufun0` with various error distributions (normal, logistic, Student-t) and confounding types (4 cases). This gives a ground-truth ITE for evaluating:
- **Coverage**: fraction of test units whose true ITE falls in the predicted interval
- **Interval length**: efficiency of the intervals across methods and Γ values

### Real Data Application (`run_tests/`, `results/ITE/`)

Applies the methods to actual observational datasets with no known ground truth ITE. Produces sensitivity intervals as Γ varies, and classifies test units as "effectively treated" if 0 is excluded from the interval.

## Datasets

| Dataset | Description |
|---|---|
| VD | Clinical trial data with liver fibrosis markers; outcome: TGF-β at 6 months |
| VK2 | Clinical dataset (kidney-related) |
| bweight | Birth weight observational study |
| drugged_AS / drugged_TS | Drug study, two arms |
| data1 / data3 / data9 / data30 | Additional observational datasets |
| fish | Environmental study (from Zenodo) |

In simulation, VD/VK covariates are fitted with a multivariate normal and used to generate synthetic samples; the real outcome data is not used.

## Repository Structure

```
code/
├── R/                      # Core algorithm implementations
│   ├── nested_conformalSA.R    # Two-group nested CSA fit
│   ├── conformal_SA.R          # Single-arm conformal with sensitivity weighting
│   ├── cutoff_SA.R             # Sensitivity-adjusted conformal cutoff under Γ
│   ├── fit_and_predict_band.R  # Smooth & extrapolate ITE bounds to test set
│   ├── predict.nested.R        # Get ITE bands from a nested_conformalSA object
│   ├── conformal_learners.R    # Huber, RF, quantRF, Boosting wrappers
│   ├── conformalIte.R          # cfcausal ITE interface (nest/naive/counterfactual)
│   └── util_SA.R               # Simulation utilities (samplecf, summary_CI, ...)
├── data/                   # Observational datasets (CSV/xlsx)
├── run_tests/              # Experiment entry scripts
│   ├── VD.R / VK.R / bweight.R / ...   # Real data experiments
│   └── syn/                             # Simulation experiments using real covariate distributions
│       ├── VD-syn.R / VD-huber.R
│       └── VK-syn.R / VK-huber.R / ...
├── results/                # Output CSVs (gitignored)
│   ├── ITE/                # Per-dataset real-data ITE results
│   └── synthetic/          # Simulation coverage & length results
├── figures/                # Output figures (gitignored)
├── plot_figures/           # Plotting scripts
│   └── or_plot/            # Outcome regression diagnostics & final figures
├── exp-cf/                 # Counterfactual conformal experiments (cfcausal baseline)
├── exp-ite/                # ITE simulation experiments
└── exp-fish/               # Fish dataset experiment
```

## Quickstart

All scripts run from the `code/` directory (the R package root).

**Install dependencies:**
```r
install.packages(c("devtools", "randomForest", "grf", "gbm", "h2o",
                   "MASS", "argparse", "dplyr", "ggplot2", "readxl"))
devtools::install_github("lihualei71/cfcausal")
```

**Run simulation experiment (VD covariate distribution):**
```bash
cd code
Rscript run_tests/syn/VD-huber.R \
  --alpha 0.2 \
  --cftype 2 \
  --ntrial 50 \
  --ntrain 2000 \
  --ntest 5000 \
  --seed 1
```

**Run real data experiment:**
```bash
cd code
Rscript run_tests/VD.R \
  --data_name VD \
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
| `--huber_alpha` | Huber loss quantile parameter α ∈ [0,1] |

Results are saved to `results/synthetic/<dataset>/<errdist>/gmm_<Γ>/` for simulation, and `results/ITE/<data_name>-<method>/alpha_<α>_gmm_<Γ>/` for real data.

## References

- Lei, L. & Candès, E. (2021). Conformal inference of counterfactuals and individual treatment effects. *JRSS-B*. [GitHub](https://github.com/lihualei71/cfcausal)
- Yin, M. et al. Conformal sensitivity analysis for individual treatment effects. [GitHub](https://github.com/mingzhang-yin/conformal-sensitivity-analysis)
- Fish dataset: [Zenodo 5794567](https://zenodo.org/records/5794567)
- Additional data: [Zenodo 4404285](https://zenodo.org/records/4404285)
