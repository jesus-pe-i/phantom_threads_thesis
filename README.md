# Phantom Threads

Replication code for ***Phantom Threads: Structured Priors for Network Recovery on Block VAR Models***.

This repository contains the estimation, simulation, benchmarking, and validation code developed for the thesis. The project studies high-dimensional VARs in which variables are grouped into multivariate units, so that the VAR coefficient matrices define a **directed, lag-specific network between units**.

The main contribution is **M3**, a hierarchical Bayesian shrinkage prior that aligns shrinkage with the unit-to-unit-lag blocks defining the network. The thesis also adapts the **GIGG** prior to the Block VAR setting and compares both methods against Bayesian, penalised, and structurally restricted alternatives.

## Network target

For a VAR

[
y_t = \sum_{\ell=1}^{p} A_\ell y_{t-\ell} + u_t,
]

each lag matrix is partitioned into blocks

[
A_{ij,\ell} \in \mathbb{R}^{m\times m},
]

where (A_{ij,\ell}) represents the lag-(\ell) predictive relationship

[
j \rightarrow i.
]

Interaction strength is measured as

[
s_{ij,\ell}
===========

\frac{\lVert A_{ij,\ell}\rVert_F}{m}.
]

The benchmark therefore evaluates network topology, lag attribution, active-signal recovery, and inactive leakage separately.

## Models

Seven estimators are implemented behind a common interface:

* **M3** — proposed block-aligned Bayesian shrinkage prior
* **GIGG** — Block VAR adaptation of the group inverse-gamma gamma prior
* **Half-t** — coefficient-wise Bayesian global-local baseline
* **gEN** — blockwise group elastic net via `adelie`
* **HLAG** — hierarchical lag regularisation via `BigVAR`
* **MAR** — Matrix Autoregression
* **NIRVAR** — Network-Informed Restricted VAR

The Bayesian models use purpose-built C++ Gibbs samplers through `RcppArmadillo`.

## Repository structure

```text
phantom_threads_thesis/
├── R/                  # R model, simulation and benchmark infrastructure
├── src/                # C++ Bayesian sampler backends
├── scripts/
│   ├── simulation/     # DGP review utilities
│   └── validation/     # Integration and reproducibility checks
└── README.md
```

The main model implementations follow the pattern

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
    "adelie",
    "BigVAR",
    "mclust"
  )
)
```

A working C++ toolchain is also required.

## Basic usage

Load the benchmark infrastructure:

```r
source("R/benchmarks.R")
source_benchmark_models()
```

Compile the Bayesian engines:

```r
load_half_t_cpp()
load_gigg_cpp()
load_m3_cpp()
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

## Validation

The repository includes permanent validation scripts:

```bash
Rscript scripts/validation/validate_var_core.R
Rscript scripts/validation/validate_bayesian_models.R
Rscript scripts/validation/validate_benchmarks.R
```

These check the common VAR representation, the three Bayesian engines, fixed-seed reproducibility, and the unified benchmark interface.

## Thesis results

The simulation study considers block-diagonal, community, core-periphery, hub, and larger complex networks.

Overall, no estimator dominates every inferential target. M3 provides particularly strong control of inactive network leakage and strong topology recovery, while GIGG allows greater flexibility in recovering active magnitudes. The group elastic net often provides a strong compromise, while HLAG performs well for active magnitudes but can be weaker for lag attribution.

The thesis also applies the framework to a 56-dimensional panel of 14 Mexican banks, finding dominant own-bank dynamics alongside a structured cross-bank predictive network.

## Reproducibility note

This repository contains the reusable model, simulation, benchmark, and validation infrastructure.

The full frozen DGP bank, simulation campaign outputs, Mexican banking data, and thesis figure/table artifacts are not currently included, so the repository alone does not reproduce every numerical result in the thesis.

## Citation

If you use this code, please cite:

> Piñera Esquivel, Jesus Antonio (2026).
> ***Phantom Threads: Structured Priors for Network Recovery on Block VAR Models.***

## Author

**Jesus Antonio Piñera Esquivel**
