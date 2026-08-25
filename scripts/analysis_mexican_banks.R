# Mexican banking empirical application -----
#
# Reconstructs the final CNBV banking panel, selects the VAR lag
# order, fits the four empirical estimators and reports headline
# network summaries.


# Setup -----

## Project sources -----

project_files <- c(
  "R/var_data.R", "R/var_coefficients.R", "R/network_blocks.R",
  "R/m3_structure.R", "R/m3_sampler.R", "R/m3_fit.R",
  "R/gigg_structure.R", "R/gigg_sampler.R", "R/gigg_fit.R",
  "R/fit_adelie.R", "R/fit_bigvar.R"
)

invisible(lapply(project_files, source))
load_m3_cpp(rebuild = FALSE)
load_gigg_cpp(rebuild = FALSE)

## Packages and settings -----

library(data.table)

data_path <- "data/mexican_banks"
data_file <- file.path(data_path, "sh_datos_40.csv")
bank_file <- file.path(data_path, "cat_instituciones_40.csv")

start_date <- as.IDate("2003-01-01")
end_date <- as.IDate("2026-05-01")

model_vars <- c("imor", "eprc", "liquidity", "capital")

concept_ids <- c(
  imor = 40200017L,
  eprc = 40200118L,
  cash = 40100045L,
  assets = 40100001L,
  equity = 40100002L
)


# Data preparation -----

## Load CNBV data -----

bank_cat <- fread(bank_file, encoding = "Latin-1")
bank_cat[, entidad := as.integer(entidad)]

bank_map <- unique(
  bank_cat[entidad >= 40000L, .(entidad, bank_name = nombre_entidad)],
  by = "entidad"
)

raw <- fread(
  data_file,
  select = c("sector", "idconcepto", "entidad", "periodo", "saldo", "valor"),
  na.strings = c("", "NA", "N/A")
)

raw[, `:=`(
  sector = as.integer(sector),
  idconcepto = as.integer(idconcepto),
  entidad = as.integer(entidad),
  periodo = as.integer(periodo),
  saldo = as.integer(saldo),
  valor = as.numeric(valor)
)]

raw <- raw[
  sector == 40L &
    entidad >= 40000L &
    saldo == 133L &
    idconcepto %in% concept_ids
]

raw[, date := as.IDate(paste0(periodo, "01"), format = "%Y%m%d")]
raw <- raw[date >= start_date & date <= end_date]
raw[, item := names(concept_ids)[match(idconcepto, concept_ids)]]


## Construct variables and bank panel -----

wide <- dcast(raw, entidad + date ~ item, value.var = "valor")

wide[, `:=`(
  imor = 100 * imor,
  eprc = 100 * eprc,
  liquidity = 100 * cash / assets,
  capital = 100 * equity / assets
)]

dates <- seq(start_date, end_date, by = "month")

availability <- wide[, .(
  months = uniqueN(date),
  complete = sum(complete.cases(.SD))
), by = entidad, .SDcols = model_vars]

keep <- availability[
  months == length(dates) &
    complete == length(dates),
  entidad
]

keep <- setdiff(keep, c(40108L, 40112L))

sample <- wide[
  entidad %in% keep,
  .(mean_assets = mean(assets)),
  by = entidad
]

sample <- merge(sample, bank_map, by = "entidad", all.x = TRUE)
sample[, asset_share := 100 * mean_assets / sum(mean_assets)]
setorder(sample, -asset_share)
sample[, unit_index := .I]

panel <- merge(
  wide[entidad %in% keep],
  sample[, .(entidad, unit_index, bank_name, mean_assets, asset_share)],
  by = "entidad"
)

panel <- melt(
  panel,
  id.vars = c("date", "unit_index", "entidad", "bank_name", "mean_assets", "asset_share"),
  measure.vars = model_vars,
  variable.name = "variable",
  value.name = "value",
  variable.factor = FALSE
)

panel[, var_order := match(variable, model_vars)]
setorder(panel, unit_index, var_order, date)
panel[, var_order := NULL]


## Detrend, deseasonalise and standardise -----

trend_grid <- data.table(
  trend = c("none", "linear", "log", "quadratic"),
  complexity = c(0L, 1L, 1L, 2L)
)

fit_deterministic <- function(y, trend) {
  tt <- seq_along(y)
  
  switch(
    trend,
    none = lm(y ~ 1),
    linear = lm(y ~ tt),
    log = lm(y ~ log(tt)),
    quadratic = lm(y ~ tt + I(tt^2))
  )
}

get_10 <- function(x) {
  if (is.matrix(x) && "10pct" %in% colnames(x)) return(as.numeric(x[1L, "10pct"]))
  if (is.matrix(x) && "10pct" %in% rownames(x)) return(as.numeric(x["10pct", 1L]))
  if ("10pct" %in% names(x)) return(as.numeric(x["10pct"]))
  as.numeric(x[1L])
}

stationarity_pass <- function(x) {
  adf <- urca::ur.df(x, type = "none", lags = 12L, selectlags = "BIC")
  kpss <- urca::ur.kpss(x, type = "mu", lags = "short")
  
  c(
    adf = as.numeric(adf@teststat[1L]) < get_10(adf@cval),
    kpss = as.numeric(kpss@teststat[1L]) < get_10(kpss@cval)
  )
}

choose_trend <- function(y) {
  candidates <- rbindlist(lapply(seq_len(nrow(trend_grid)), function(i) {
    fit <- fit_deterministic(y, trend_grid$trend[i])
    pass <- stationarity_pass(residuals(fit))
    
    data.table(
      trend = trend_grid$trend[i],
      complexity = trend_grid$complexity[i],
      bic = BIC(fit),
      adf_pass = pass["adf"],
      kpss_pass = pass["kpss"]
    )
  }))
  
  admissible <- candidates[adf_pass == TRUE & kpss_pass == TRUE]
  basis <- "ADF + KPSS"
  
  if (!nrow(admissible)) {
    admissible <- candidates[adf_pass == TRUE]
    basis <- "ADF only"
  }
  
  if (!nrow(admissible)) return(list(trend = NA_character_, basis = NA_character_))
  
  setorder(admissible, complexity, bic)
  list(trend = admissible$trend[1L], basis = basis)
}

trend_choice <- panel[
  ,
  choose_trend(value),
  by = .(unit_index, variable)
]

# Only Invex EPRC reaches this fallback in the thesis sample.
trend_choice[is.na(trend), `:=`(
  trend = "linear",
  basis = "manual fallback"
)]

panel <- merge(
  panel,
  trend_choice,
  by = c("unit_index", "variable"),
  all.x = TRUE,
  sort = FALSE
)

panel[, resid := residuals(fit_deterministic(value, first(trend))), by = .(unit_index, variable)]
panel[, month := as.integer(format(date, "%m"))]
panel[, seasonal := median(resid), by = .(unit_index, variable, month)]
panel[, deseason := resid - seasonal]
panel[, z := (deseason - mean(deseason)) / sd(deseason), by = .(unit_index, variable)]

series_map <- unique(
  panel[, .(
    unit_index,
    entidad,
    bank_name,
    mean_assets,
    asset_share,
    variable,
    trend,
    basis
  )]
)

series_map[, var_order := match(variable, model_vars)]
setorder(series_map, unit_index, var_order)
series_map[, series_name := sprintf("b%02d_%s", unit_index, variable)]

panel[
  series_map[, .(unit_index, variable, series_name)],
  on = .(unit_index, variable),
  series_name := i.series_name
]

wide_empirical <- dcast(panel, date ~ series_name, value.var = "z")
series_columns <- series_map$series_name

Z_empirical <- as.matrix(wide_empirical[, ..series_columns])
dates_empirical <- wide_empirical$date
colnames(Z_empirical) <- series_columns

n_units_empirical <- nrow(sample)
m_empirical <- length(model_vars)
N_empirical <- ncol(Z_empirical)

data_summary <- data.table(
  start = min(dates_empirical),
  end = max(dates_empirical),
  observations = nrow(Z_empirical),
  banks = n_units_empirical,
  variables_per_bank = m_empirical,
  system_dimension = N_empirical
)

data_summary
sample

# Lag selection -----

## Helpers -----

lag_grid <- 1:3
ridge_grid <- 10^seq(-4, 1, length.out = 30)
cv_months <- 60L

make_var_data <- function(Z, p) {
  tt <- (p + 1L):nrow(Z)
  X <- do.call(cbind, lapply(seq_len(p), function(l) Z[tt - l, , drop = FALSE]))
  list(X = X, Y = Z[tt, , drop = FALSE])
}

fit_ridge_var <- function(Z, p) {
  dat <- make_var_data(Z, p)
  Xc <- sweep(dat$X, 2L, colMeans(dat$X))
  Yc <- sweep(dat$Y, 2L, colMeans(dat$Y))
  XtX <- crossprod(Xc)
  XtY <- crossprod(Xc, Yc)
  scale <- mean(diag(XtX))
  lambdas <- ridge_grid * scale
  d <- svd(Xc, nu = 0, nv = 0)$d
  
  gcv <- vapply(lambdas, function(lambda) {
    B <- solve(XtX + diag(lambda, ncol(Xc)), XtY)
    U <- Yc - Xc %*% B
    df <- sum(d^2 / (d^2 + lambda))
    sum(U^2) / (nrow(Xc) - df)^2
  }, numeric(1L))
  
  best <- which.min(gcv)
  B <- solve(XtX + diag(lambdas[best], ncol(Xc)), XtY)
  
  list(B = B, x_mean = colMeans(dat$X), y_mean = colMeans(dat$Y), lambda_ratio = ridge_grid[best])
}


## Information criteria -----

information_criteria <- rbindlist(lapply(lag_grid, function(p) {
  dat <- make_var_data(Z_empirical, p)
  
  B <- solve(crossprod(dat$X), crossprod(dat$X, dat$Y))
  U <- dat$Y - dat$X %*% B
  
  n <- nrow(U)
  k <- ncol(dat$X) * ncol(dat$Y)
  logdet <- as.numeric(determinant(crossprod(U) / n, logarithm = TRUE)$modulus)
  
  data.table(
    p = p,
    AIC = logdet + 2 * k / n,
    BIC = logdet + log(n) * k / n
  )
}))


## Expanding-window forecasts -----

forecast_rows <- seq.int(nrow(Z_empirical) - cv_months + 1L, nrow(Z_empirical))
cv_rows <- list()

for (p in lag_grid) {
  for (t in forecast_rows) {
    train <- Z_empirical[seq_len(t - 1L), , drop = FALSE]
    fit <- fit_ridge_var(train, p)
    x_new <- unlist(lapply(seq_len(p), function(l) Z_empirical[t - l, ]))
    forecast <- fit$y_mean + as.numeric((x_new - fit$x_mean) %*% fit$B)
    error <- Z_empirical[t, ] - forecast
    
    cv_rows[[length(cv_rows) + 1L]] <- data.table(
      p = p, date = dates_empirical[t],
      MSE = mean(error^2), MAE = mean(abs(error))
    )
  }
}

cv_results <- rbindlist(cv_rows)

forecast_summary <- cv_results[, .(
  MSE = mean(MSE),
  SE = sd(MSE) / sqrt(.N),
  MAE = mean(MAE)
), by = p]


## Select lag order -----

best_forecast <- forecast_summary[which.min(MSE)]
one_se_cutoff <- best_forecast$MSE + best_forecast$SE
p_empirical <- min(forecast_summary[MSE <= one_se_cutoff, p])

lag_selection <- merge(information_criteria, forecast_summary, by = "p")
lag_selection[, selected := p == p_empirical]

# Fit models -----

## Common model input -----

Y_list_empirical <- lapply(seq_len(n_units_empirical), function(i) {
  cols <- ((i - 1L) * m_empirical + 1L):(i * m_empirical)
  out <- Z_empirical[, cols, drop = FALSE]
  colnames(out) <- model_vars
  out
})

names(Y_list_empirical) <- sample[order(unit_index)]$bank_name


## Group elastic net -----

fit_gen <- fit_adelie(
  Y_list = Y_list_empirical, p_lags = p_empirical,
  alpha = 0.5, lambda_rule = "lambda.min", n_folds = 5L,
  standardize = TRUE, intercept = FALSE, min_ratio = 0.01,
  lmda_path_size = 100L, n_threads = 1L, seed = 991L,
  keep_fits = FALSE
)


## HLAG -----

fit_hlag <- fit_bigvar(
  Y_list = Y_list_empirical, p_lags = p_empirical,
  struct = "HLAGOO", gran = c(200, 20), h = 1L,
  cv = "Rolling", IC = FALSE, T1 = NULL, T2 = NULL,
  tol = 1e-4, verbose = FALSE, selected_eps = 1e-12,
  keep_cv_fit = FALSE
)


## GIGG -----

fit_gigg_empirical <- fit_gigg(
  Y_list = Y_list_empirical, p_lags = p_empirical,
  chains = 4L, burnin = 1000L, draws = 10000L, thin = 1L,
  seed = 991L, beta_algorithm = "auto",
  use_asis = TRUE, asis_every = 1L,
  keep_chain_results = FALSE
)


## M3 -----

fit_m3_empirical <- fit_m3(
  Y_list = Y_list_empirical, p_lags = p_empirical,
  chains = 4L, burnin = 1000L, draws = 10000L, thin = 1L,
  seed = 991L, beta_algorithm = "auto",
  q_update = "gibbs_end",
  use_c_asis = TRUE, c_asis_every = 1L,
  keep_chain_results = FALSE, keep_all = FALSE
)

## Final fits object and checks

fits <- list(
  gEN = fit_gen,
  HLAG = fit_hlag,
  GIGG = fit_gigg_empirical,
  M3 = fit_m3_empirical
)

lag_selection
p_empirical
names(fits)


# Network summaries -----

## Block strengths -----

bank_names <- sample[order(unit_index), bank_name]

foreign_edges <- rbindlist(lapply(names(fits), function(model) {
  S <- fits[[model]]$s_hat_lag[, , 1L]
  
  grid <- CJ(receiver_unit = seq_len(n_units_empirical),
             sender_unit = seq_len(n_units_empirical))
  
  grid[, strength := S[cbind(receiver_unit, sender_unit)]]
  
  grid[, `:=`(
    model = model,
    receiver_bank = bank_names[receiver_unit],
    sender_bank = bank_names[sender_unit]
  )]
  
  grid[receiver_unit != sender_unit]
  
}))

model_summary <- rbindlist(lapply(names(fits), function(model) {
  
  S <- fits[[model]]$s_hat_lag[, , 1L]
  foreign <- S[row(S) != col(S)]
  self <- diag(S)
  
  data.table(
    model = model,
    median_cross = median(foreign),
    max_cross = max(foreign),
    median_self = median(self),
    self_cross_ratio = median(self) / median(foreign)
  )
}))


## Sender and receiver rankings -----

bank_roles <- rbindlist(lapply(names(fits), function(model) {
  
  S <- fits[[model]]$s_hat_lag[, , 1L]
  diag(S) <- 0
  
  data.table(
    model = model,
    unit_index = seq_len(n_units_empirical),
    bank_name = bank_names,
    sender_strength = colSums(S),
    receiver_strength = rowSums(S)
  )
}))

bank_roles[, sender_rank := frank(-sender_strength, ties.method = "average"), by = model]
bank_roles[, receiver_rank := frank(-receiver_strength, ties.method = "average"), by = model]

model_summary[
  bank_roles[sender_rank == 1, .(top_sender = paste(bank_name, collapse = ", ")), by = model],
  on = "model",
  top_sender := i.top_sender
]

model_summary[
  bank_roles[receiver_rank == 1, .(top_receiver = paste(bank_name, collapse = ", ")), by = model],
  on = "model",
  top_receiver := i.top_receiver
]


## Cross-model agreement -----

edge_wide <- dcast(
  foreign_edges,
  sender_unit + receiver_unit + sender_bank + receiver_bank ~ model,
  value.var = "strength"
)

bayesian_edge_spearman <- cor(edge_wide$GIGG, edge_wide$M3, method = "spearman")


## Selected relationships -----

selected_pairs <- data.table(
  sender_bank = c("Banamex", "Banco del Bajío", "Inbursa", "Santander", "Banregio"),
  receiver_bank = c("Invex", "Afirme", "HSBC", "HSBC", "Banca Mifel")
)

selected_relationships <- merge(
  foreign_edges,
  selected_pairs,
  by = c("sender_bank", "receiver_bank")
)

selected_relationships[, relationship := paste(sender_bank, "->", receiver_bank)]

selected_relationships <- dcast(
  selected_relationships,
  relationship ~ model,
  value.var = "strength"
)

setcolorder(selected_relationships, c("relationship", "gEN", "HLAG", "GIGG", "M3"))


# Headline results -----

cat("\nMEXICAN BANKING EMPIRICAL APPLICATION\n")
cat("=====================================\n\n")

cat("DATA\n")
print(data_summary)

cat("\nLAG SELECTION\n")
print(lag_selection, digits = 5)
cat("\nSelected lag order:", p_empirical, "\n")

cat("\nNETWORK SUMMARY\n")
print(model_summary, digits = 4)

cat("\nSELECTED CROSS-BANK RELATIONSHIPS\n")
print(selected_relationships, digits = 4)

cat(
  "\nGIGG-M3 cross-bank strength rank correlation:",
  round(bayesian_edge_spearman, 3),
  "\n"
)

