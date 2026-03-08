# Conformal Sensitivity Analysis for Individual Treatment Effects

This project applies **Conformal Sensitivity Analysis (CSA)** to estimate prediction intervals for Individual Treatment Effects (ITE) on observational datasets, with validity guarantees under unmeasured confounding parameterized by a sensitivity parameter Γ.

## Background

Standard conformal inference for ITE assumes unconfoundedness (no unmeasured confounders). This project extends the nested conformal approach from [Lei & Candes (2021)](https://github.com/lihualei71/cfcausal) with a Γ-sensitivity model: the intervals remain valid as long as the true odds ratio of selection bias does not exceed Γ. The key question asked for each test unit is whether its ITE interval contains zero — i.e., whether the treatment is **effective** under the assumed confounding level.

## Methods

The core pipeline (`code/R/`) implements:

| Function | Description |
|---|---|
| `nested_conformalSA` | Two-group nested approach: fits outcome models on group 1, calibrates sensitivity-adjusted conformal scores on group 2 |
| `predict.nested` | Produces pointwise ITE bands from a fitted `nested_conformalSA` object |
| `fit_and_predict_band` | Fits a regression model to smooth the ITE bounds and extrapolates to test set |
| `conformal_SA` | Single-arm conformal inference with sensitivity weighting |
| `cutoff_SA` | Computes the sensitivity-adjusted conformal cutoff under Γ |

Outcome models supported: Random Forest (`RF`), Quantile RF (`quantRF`), Boosting, Quantile Boosting, BART, Quantile BART, or any user-supplied function.

Two inference types are supported:
- **mean**: standard conformal inference based on conditional mean regression
- **CQR**: conformalized quantile regression

## Repository Structure

```
code/
├── R/                      # Core algorithm implementations
│   ├── nested_conformalSA.R
│   ├── conformal_SA.R
│   ├── cutoff_SA.R
│   ├── fit_and_predict_band.R
│   ├── predict.nested.R
│   ├── conformalIte.R      # Original cfcausal ITE interface
│   ├── conformalCf.R
│   ├── util_SA.R           # Sensitivity simulation utilities
│   └── ...
├── data/                   # Datasets
│   ├── VD.csv              # Clinical dataset (liver fibrosis, outcome: TGF6)
│   ├── VK2.csv             # Clinical dataset (kidney-related)
│   ├── bweight.csv         # Birth weight observational study
│   ├── drugged_AS.csv      # Drug study (AS arm)
│   ├── drugged_TS.csv      # Drug study (TS arm)
│   ├── data1/3/9/30.csv    # Additional datasets
│   └── fish.csv            # Fish dataset (from Zenodo)
├── run_tests/              # Experiment entry scripts (one per dataset)
│   ├── VD.R
│   ├── VK.R
│   ├── bweight.R
│   ├── bweight_Gamma.R     # Sensitivity parameter sweep for bweight
│   ├── drug_AS.R / drug_TS.R
│   ├── data1/3/9/30.R
│   └── syn/                # Synthetic data experiments (VD-syn, VK-syn, huber variants)
├── results/                # Output CSV files
│   ├── ITE/                # Per-dataset ITE interval results
│   └── synthetic/          # Simulation results
├── figures/                # Output figures
│   ├── ITE/
│   └── cov_len_shrinkage/  # Coverage & interval length vs Γ
├── plot_figures/           # Plotting scripts
│   └── or_plot/            # Outcome regression diagnostics & final plots
├── exp-cf/                 # Counterfactual conformal experiments (simulation)
├── exp-ite/                # ITE conformal experiments (simulation)
└── exp-fish/               # Fish dataset experiment
```

## Quickstart

All scripts are run from the `code/` directory (the R package root).

**Install dependencies:**
```r
install.packages(c("devtools", "randomForest", "grf", "gbm", "argparse",
                   "dplyr", "ggplot2", "readxl"))
devtools::install_github("lihualei71/cfcausal")
```

**Run an experiment (e.g., VD dataset):**
```bash
cd code
Rscript run_tests/VD.R \
  --data_name VD \
  --method mean \
  --gmm_star 3 \
  --alpha 0.2 \
  --ntrial 10 \
  --seed1 123 \
  --seed2 45
```

Results are saved to `results/ITE/<data_name>-<method>/alpha_<alpha>_gmm_<gmm_star>/`.

**Key CLI arguments:**

| Argument | Default | Description |
|---|---|---|
| `--data_name` | `VD` | Dataset name (matches filename in `data/`) |
| `--method` | `mean` | Inference type: `mean` or `cqr` |
| `--gmm_star` | `3` | Sensitivity parameter Γ (≥ 1; Γ=1 means no unmeasured confounding) |
| `--alpha` | `0.2` | Miscoverage level (produces 1−α prediction intervals) |
| `--ntrial` | `10` | Number of random trials |
| `--seed1` | `123` | Data split seed |
| `--seed2` | `45` | Model fitting seed |

## Simulation Experiments

Synthetic experiments under `run_tests/syn/` and `exp-ite/` evaluate coverage and interval length across varying Γ values for homoscedastic and heteroscedastic error settings.

The `exp-cf/` folder contains counterfactual conformal experiments from the original cfcausal paper, included for comparison.

## References

- Lei, L. & Candès, E. (2021). Conformal inference of counterfactuals and individual treatment effects. *JRSS-B*. [GitHub](https://github.com/lihualei71/cfcausal)
- Yin, M. et al. Conformal sensitivity analysis for individual treatment effects. [GitHub](https://github.com/mingzhang-yin/conformal-sensitivity-analysis)
- Fish dataset: [Zenodo 5794567](https://zenodo.org/records/5794567)
- Additional data: [Zenodo 4404285](https://zenodo.org/records/4404285)
