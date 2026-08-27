# Frozen DGP bank

Each directory contains one frozen Block VAR data-generating process.

Directory names follow the pattern `<family>_n<n_units>_m<m>_p<p_lags>`.

`A_list.rds` contains the lag coefficient matrices used by the public simulation and benchmark infrastructure. `meta.rds` stores construction metadata, while `recipe.R` records the original DGP construction recipe.

The bank can be reviewed with `scripts/simulation/review_dgp_bank.R`.
