# Runs selected Bayesian benchmark estimators over the ordinary DGP bank.
#
# Models are run sequentially and records in each model are fitted in parallel.
# Chain draws are retained only long enough to compute compact diagnostics.


rm(
  list = ls()
)


# Settings -----

models <- c(
  "m3"
)


dgp_ids <- c(
  
  
  # Block diagonal -----
  
  "blockdiag_n6_m3_p2",
  
  "blockdiag_n6_m3_p1",
  "blockdiag_n6_m3_p5",
  
  "blockdiag_n6_m2_p2",
  "blockdiag_n6_m6_p2",
  
  "blockdiag_n4_m3_p2",
  "blockdiag_n10_m3_p2",
  
  
  # Communities -----
  
  "communities_n6_m3_p2",
  
  "communities_n6_m3_p1",
  "communities_n6_m3_p5",
  
  "communities_n6_m2_p2",
  "communities_n6_m6_p2",
  
  "communities_n4_m3_p2",
  "communities_n10_m3_p2",
  
  
  # Core-periphery -----
  
  "coreperiphery_n6_m3_p2",
  
  "coreperiphery_n6_m3_p1",
  "coreperiphery_n6_m3_p5",
  
  "coreperiphery_n6_m2_p2",
  "coreperiphery_n6_m6_p2",
  
  "coreperiphery_n4_m3_p2",
  "coreperiphery_n10_m3_p2",
  
  
  # Hubs -----
  
  "hubs_n6_m3_p2",
  
  "hubs_n6_m3_p1",
  "hubs_n6_m3_p5",
  
  "hubs_n6_m2_p2",
  "hubs_n6_m6_p2",
  
  "hubs_n4_m3_p2",
  "hubs_n10_m3_p2"
)


T_grid <- c(
  50L,
  150L,
  400L
)


seeds <- c(
  101L,
  202L,
  303L,
  404L,
  505L,
  606L,
  707L,
  808L,
  909L,
  1000L,
  1101L,
  1202L,
  1303L,
  1404L,
  1505L,
  1606L,
  1707L,
  1808L,
  1909L,
  2000L
)


workers <- 14L

simulation_burn_in <- 300L
sigma <- 1

fit_seed <- 991L

mixing_diagnostics_seed <- 2718L
max_monitor_blocks <- 100L
extra_monitor_beta <- 200L


model_args <- list(
  
  m3 = list(
    chains = 4L,
    burnin = 1000L,
    draws = 8000L,
    thin = 1L,
    beta_algorithm = "auto",
    q_update = "gibbs_end",
    use_c_asis = TRUE,
    c_asis_every = 1L
  ),
  
  half_t = list(
    global_grouping = "self_diagonal",
    chains = 4L,
    burnin = 1000L,
    draws = 8000L,
    thin = 1L,
    beta_algorithm = "auto",
    use_asis = TRUE,
    asis_every = 1L
  ),
  
  gigg = list(
    chains = 4L,
    burnin = 1000L,
    draws = 8000L,
    thin = 1L,
    beta_algorithm = "auto",
    use_asis = TRUE,
    asis_every = 1L
  )
)


output_root <- file.path(
  "outputs",
  "benchmark_runs"
)


# Checks -----

valid_models <- names(
  model_args
)

unknown_models <- setdiff(
  models,
  valid_models
)

if (length(unknown_models) > 0L) {
  stop(
    "Unknown Bayesian benchmark models: ",
    paste(
      unknown_models,
      collapse = ", "
    )
  )
}


missing_dgps <- dgp_ids[
  !dir.exists(
    file.path(
      "data",
      "dgp_bank",
      dgp_ids
    )
  )
]

if (length(missing_dgps) > 0L) {
  stop(
    "Missing DGPs from data/dgp_bank: ",
    paste(
      missing_dgps,
      collapse = ", "
    )
  )
}


# Parallel settings -----

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)


# Sources -----

source(
  "R/benchmarks.R"
)

source(
  "R/var_data.R"
)

source(
  "R/var_coefficients.R"
)

source(
  "R/network_blocks.R"
)

source(
  "R/simulation.R"
)

source(
  "R/benchmark_records.R"
)

source(
  "R/benchmark_metrics.R"
)

source(
  "R/benchmark_bayesian_diagnostics.R"
)


# Packages -----

required_packages <- c(
  "Rcpp",
  "RcppArmadillo",
  "posterior"
)

package_available <- vapply(
  required_packages,
  requireNamespace,
  logical(1L),
  quietly = TRUE
)

if (!all(package_available)) {
  stop(
    "Missing packages: ",
    paste(
      required_packages[!package_available],
      collapse = ", "
    )
  )
}


# Model backends -----

model_worker_files <- list(
  
  m3 = c(
    "R/m3_structure.R",
    "R/m3_sampler.R",
    "R/m3_fit.R"
  ),
  
  half_t = c(
    "R/half_t_structure.R",
    "R/half_t_sampler.R",
    "R/half_t_fit.R"
  ),
  
  gigg = c(
    "R/gigg_structure.R",
    "R/gigg_sampler.R",
    "R/gigg_fit.R"
  )
)


native_loaders <- c(
  m3 = "load_m3_cpp",
  half_t = "load_half_t_cpp",
  gigg = "load_gigg_cpp"
)


# Worker -----

fit_record_bayesian <- function(
    record,
    model,
    fit_seed,
    model_args,
    mixing_diagnostics_seed,
    max_monitor_blocks,
    extra_monitor_beta) {
  
  plan <- make_bayesian_monitor_plan(
    record = record,
    seed = mixing_diagnostics_seed,
    max_blocks = max_monitor_blocks,
    extra_beta = extra_monitor_beta
  )
  
  monitor_args <- make_bayesian_monitor_args(
    model = model,
    plan = plan
  )
  
  args <- utils::modifyList(
    model_args,
    monitor_args
  )
  
  args$keep_chain_results <- TRUE
  
  fit <- fit_benchmark(
    model = model,
    Y_list = record$data$Y_list,
    p_lags = record$dimensions$p_lags,
    seed = fit_seed,
    args = args
  )
  
  metrics <- evaluate_benchmark_fit(
    fit = fit,
    truth = record$truth
  )
  
  diagnostics <- summarise_bayesian_diagnostics(
    fit
  )
  
  diagnostics <- cbind(
    data.frame(
      record_id = record$id,
      dgp_id = record$design$dgp_id,
      seed = record$design$seed,
      T_obs = record$design$T_obs,
      model = model,
      stringsAsFactors = FALSE
    ),
    diagnostics
  )
  
  result <- cbind(
    data.frame(
      record_id = record$id,
      dgp_id = record$design$dgp_id,
      seed = record$design$seed,
      T_obs = record$design$T_obs,
      model = model,
      backend = fit$backend,
      n_units = record$dimensions$n_units,
      m = record$dimensions$m,
      N = record$dimensions$N,
      p_lags = record$dimensions$p_lags,
      k = (
        record$dimensions$N *
          record$dimensions$p_lags
      ),
      runtime_sec = fit$runtime_seconds,
      fit_seed = fit$fit_seed,
      stringsAsFactors = FALSE
    ),
    metrics
  )
  
  strengths <- fit$s_hat_lag
  
  rm(
    fit
  )
  
  gc()
  
  list(
    result = result,
    strengths = strengths,
    diagnostics = diagnostics
  )
}


# Campaign -----

root <- normalizePath(
  ".",
  winslash = "/",
  mustWork = TRUE
)

common_worker_files <- c(
  "R/benchmarks.R",
  "R/var_data.R",
  "R/var_coefficients.R",
  "R/network_blocks.R",
  "R/benchmark_metrics.R",
  "R/benchmark_bayesian_diagnostics.R"
)

campaign_start <- Sys.time()


for (model in models) {
  
  model_start <- Sys.time()
  
  current_model_args <- model_args[[model]]
  
  worker_files <- c(
    common_worker_files,
    model_worker_files[[model]]
  )
  
  native_loader <- native_loaders[[model]]
  
  
  # Cluster -----
  
  cluster <- parallel::makePSOCKcluster(
    workers
  )
  
  parallel::clusterExport(
    cluster,
    c(
      "root",
      "worker_files",
      "native_loader"
    )
  )
  
  invisible(
    parallel::clusterEvalQ(
      cluster,
      {
        
        Sys.setenv(
          OMP_NUM_THREADS = "1",
          OPENBLAS_NUM_THREADS = "1",
          MKL_NUM_THREADS = "1",
          VECLIB_MAXIMUM_THREADS = "1",
          NUMEXPR_NUM_THREADS = "1"
        )
        
        setwd(
          root
        )
        
        for (file in worker_files) {
          source(
            file.path(
              root,
              file
            )
          )
        }
        
        get(
          native_loader,
          mode = "function"
        )(
          rebuild = FALSE
        )
        
        if (
          !requireNamespace(
            "posterior",
            quietly = TRUE
          )
        ) {
          stop(
            "Package 'posterior' is required on workers"
          )
        }
        
        NULL
      }
    )
  )
  
  
  # Run -----
  
  tryCatch(
    {
      
      for (dgp_id in dgp_ids) {
        
        dgp_start <- Sys.time()
        
        cat(
          "\n",
          model,
          " | ",
          dgp_id,
          "\n",
          sep = ""
        )
        
        records <- make_benchmark_records(
          dgp_id = dgp_id,
          T_grid = T_grid,
          seeds = seeds,
          burn_in = simulation_burn_in,
          sigma = sigma
        )
        
        cat(
          length(records),
          " records\n",
          sep = ""
        )
        
        output <- parallel::parLapplyLB(
          cluster,
          records,
          fit_record_bayesian,
          model = model,
          fit_seed = fit_seed,
          model_args = current_model_args,
          mixing_diagnostics_seed =
            mixing_diagnostics_seed,
          max_monitor_blocks =
            max_monitor_blocks,
          extra_monitor_beta =
            extra_monitor_beta
        )
        
        
        # Save -----
        
        results <- do.call(
          rbind,
          lapply(
            output,
            function(record_output) {
              record_output$result
            }
          )
        )
        
        rownames(results) <- NULL
        
        strengths <- lapply(
          output,
          function(record_output) {
            record_output$strengths
          }
        )
        
        names(strengths) <- names(records)
        
        diagnostics <- do.call(
          rbind,
          lapply(
            output,
            function(record_output) {
              record_output$diagnostics
            }
          )
        )
        
        rownames(diagnostics) <- NULL
        
        benchmark <- list(
          config = list(
            model = model,
            dgp_id = dgp_id,
            T_grid = T_grid,
            seeds = seeds,
            simulation_burn_in =
              simulation_burn_in,
            sigma = sigma,
            fit_seed = fit_seed,
            workers = workers,
            model_args =
              current_model_args,
            mixing_diagnostics_seed =
              mixing_diagnostics_seed,
            max_monitor_blocks =
              max_monitor_blocks,
            extra_monitor_beta =
              extra_monitor_beta
          ),
          results = results,
          strengths = strengths,
          diagnostics = diagnostics
        )
        
        output_dir <- file.path(
          output_root,
          model,
          dgp_id
        )
        
        dir.create(
          output_dir,
          recursive = TRUE,
          showWarnings = FALSE
        )
        
        saveRDS(
          benchmark,
          file.path(
            output_dir,
            "benchmark.rds"
          )
        )
        
        cat(
          "Saved ",
          nrow(results),
          " fits and ",
          nrow(diagnostics),
          " diagnostic rows in ",
          round(
            as.numeric(
              difftime(
                Sys.time(),
                dgp_start,
                units = "mins"
              )
            ),
            2L
          ),
          " minutes.\n",
          sep = ""
        )
        
        rm(
          records,
          output,
          results,
          strengths,
          diagnostics,
          benchmark
        )
        
        gc()
      }
      
    },
    finally = {
      
      parallel::stopCluster(
        cluster
      )
    }
  )
  
  
  # Model runtime -----
  
  cat(
    "\nReduced ",
    model,
    " campaign wall time: ",
    round(
      as.numeric(
        difftime(
          Sys.time(),
          model_start,
          units = "hours"
        )
      ),
      3L
    ),
    " hours.\n",
    sep = ""
  )
}


# Runtime -----

cat(
  "\nTotal Bayesian campaign wall time: ",
  round(
    as.numeric(
      difftime(
        Sys.time(),
        campaign_start,
        units = "hours"
      )
    ),
    3L
  ),
  " hours.\n",
  sep = ""
)