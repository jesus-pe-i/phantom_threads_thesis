# Phantom Threads

Replication code for ***Phantom Threads: Structured Priors for Network Recovery on Block VAR Models***.

This repository contains the estimation, simulation, benchmarking, validation, and empirical-analysis code developed for the thesis. The project studies high-dimensional VARs in which variables are grouped into multivariate units, so that the VAR coefficient matrices define a **directed, lag-specific network between units**.

The main contribution is **M3**, a hierarchical Bayesian shrinkage prior that aligns shrinkage with the unit-to-unit-lag blocks defining the network. The thesis also adapts the **GIGG** prior to the Block VAR setting and compares both methods against Bayesian, penalised, and structurally restricted alternatives.

## Network target

The VAR contains `K` multivariate units, each observed through `m` series. Each lag coefficient matrix is therefore divided into `K × K` blocks.

A block `A[i, j, l]` contains the lag-`l` predictive effects from unit `j` to unit `i`:

**sender `j` → receiver `i`**

Network strength is measured from the Frobenius norm of each coefficient block, scaled by the number of variables per unit.

The benchmark treats these blocks as the main inferential target and evaluates separately:

* network topology recovery;
* lag attribution;
* active interaction-strength recovery; and
* leakage into inactive relationships.

## Models

Seven estimators are implemented behind a common interface:

* **M3** — proposed block-aligned Bayesian shrinkage prior
* **GIGG** — Block VAR adaptation of the group inverse-gamma gamma prior
* **Half-t** — coefficient-wise Bayesian global-local baseline
* **gEN** — blockwise group elastic net via `adelie`
* **HLAG** — hierarchical lag regularisation via `BigVAR`
* **MAR** — Matrix Autoregression
* **NIRVAR** — adapted Network-Informed Restricted VAR

The Bayesian models use purpose-built C++ Gibbs samplers through `RcppArmadillo`.

## Repository structure

```text
phantom_threads_thesis/
├── R/                          # Model, simulation and benchmark infrastructure
├── src/                        # C++ Bayesian sampler backends
├── data/
│   └── dgp_bank/               # Frozen simulation DGP bank
├── scripts/
│   ├── analysis_mexican_banks.R
│   ├── run_bayesian_benchmarks.R
│   ├── run_bayesian_full_noncomplex.R
│   ├── run_non_bayesian_benchmarks.R
│   ├── run_non_bayesian_full_noncomplex.R
│   ├── simulation/
│   │   └── review_dgp_bank.R
│   └── validation/
│       ├── validate_var_core.R
│       ├── validate_bayesian_models.R
│       └── validate_benchmarks.R
└── README.md
```

The main Bayesian model implementations follow the pattern

```text
R/m3_structure.R
R/m3_sampler.R
R/m3_fit.R

R/gigg_structure.R
R/gigg_sampler.R
R/gigg_fit.R

R/half_t_structure.R
R/half_t_sampler.R
R/half_t_fit.R
```

while `R/benchmarks.R` provides the common seven-model fitting interface.

## Requirements

The repository is research code rather than an R package and should be run from the project root.

Main dependencies include:

```r
install.packages(
  c(
    "Rcpp",
    "RcppArmadillo",
    "data.table",
    "urca",
    "adelie",
    "BigVAR",
    "mclust",
    "posterior"
  )
)
```

A working C++ toolchain is also required for the Bayesian samplers.

## Basic usage

Load the benchmark infrastructure:

```r
source("R/benchmarks.R")
source_benchmark_models()
```

Input data are supplied as a list of unit-level matrices:

```r
Y_list <- list(
  unit_1,
  unit_2,
  unit_3
)
```

where each unit contains the same number of observations and variables.

Fit M3 directly:

```r
fit <- fit_m3(
  Y_list = Y_list,
  p_lags = 2L,
  chains = 4L,
  burnin = 500L,
  draws = 1000L,
  seed = 991L
)
```

or fit several estimators through the common interface:

```r
fits <- fit_benchmarks(
  models = c(
    "m3",
    "gigg",
    "half_t",
    "adelie_gen",
    "bigvar_hlag",
    "mar",
    "nirvar"
  ),
  Y_list = Y_list,
  p_lags = 2L
)
```

Common outputs include estimated VAR coefficients and lag-specific and unit-level block-network strengths.

## Simulation benchmark

The frozen DGP bank used by the public simulation infrastructure is included under:

```text
data/dgp_bank/
```

The DGP definitions can be reviewed with:

```r
source("scripts/simulation/review_dgp_bank.R")
```

The public benchmark campaign runners are:

```text
scripts/run_bayesian_benchmarks.R
scripts/run_non_bayesian_benchmarks.R
scripts/run_bayesian_full_noncomplex.R
scripts/run_non_bayesian_full_noncomplex.R
```

These scripts use the shared benchmark infrastructure in `R/` and the frozen DGP bank in `data/dgp_bank/`.

## Mexican banking application

The thesis applies the framework to a monthly panel of 14 Mexican banks observed through four balance-sheet and credit-risk variables, producing a 56-dimensional Block VAR system.

The public empirical analysis is:

```text
scripts/analysis_mexican_banks.R
```

It reconstructs the final banking panel, performs deterministic preprocessing, compares lag orders using information criteria and expanding-window forecasts, fits gEN, HLAG, GIGG and M3, and reports headline cross-bank network summaries.

### Mexican banking data

The frozen CNBV dataset used by the empirical application is distributed separately from the Git repository because of its size.

Download the dataset from the GitHub Release:

[data-mex-v1](https://github.com/jesus-pe-i/phantom_threads_thesis/releases/tag/data-mex-v1)

Extract the release asset so that the required files are located at:

```text
data/mexican_banks/
├── sh_datos_40.csv
└── cat_instituciones_40.csv
```

The empirical application can then be run from the repository root with:

```r
source("scripts/analysis_mexican_banks.R")
```

The script reproduces the public empirical pipeline directly from the frozen data and does not write fitted objects or result files to disk.

## Validation

The repository includes permanent validation scripts:

```bash
Rscript scripts/validation/validate_var_core.R
Rscript scripts/validation/validate_bayesian_models.R
Rscript scripts/validation/validate_benchmarks.R
```

These cover the common VAR representation, Bayesian sampler infrastructure, fixed-seed reproducibility, model loading and dispatch, coefficient and network reconstruction, nested benchmark records, benchmark metrics, Bayesian settings, runtimes, and direct-versus-dispatch equivalence.

## Reproducibility note

This repository contains the reusable model, simulation, benchmark, validation, and empirical-analysis infrastructure used in the thesis.

The frozen simulation DGP bank is included directly in the repository. The Mexican banking data are provided separately through the `data-mex-v1` GitHub Release because of their size.

Full simulation campaign outputs and thesis figure/table artifacts are not included. The repository is therefore intended to reproduce the modelling and analysis pipelines rather than serve as an archive of every intermediate or generated thesis artifact.

## Citation

If you use this code, please cite:

> Pinera Esquivel, Jesus Antonio (2026).  
> ***Phantom Threads: Structured Priors for Network Recovery on Block VAR Models.***
