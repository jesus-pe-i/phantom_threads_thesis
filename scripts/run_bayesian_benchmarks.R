# Runs one Bayesian benchmark estimator over frozen simulation records.
#
# Records are fitted in parallel. Each worker retains chain draws only long
# enough to compute mixing diagnostics, then returns compact benchmark output.


rm(
  list = ls()
)


# Settings -----

dgp_ids <- c(
  "blockdiag_n6_m3_p2",
  "communities_n6_m3_p2",
  "blockdiag_n6_m4_p2",
  "communities_n6_m4_p2"
)

model <- "m3"

T_grid <- c(
  50L,
  150L,
  400L
)

seeds <- c(
  101L,
  202L,
  303L,
  404L
)


dgp_ids <- "communities_n6_m3_p2"

T_grid <- c(
  50L,
  150L
)

seeds <- c(
  101L,
  202L
)

workers <- 12L

simulation_burn_in <- 300L
sigma <- 1

fit_seed <- 991L

mixing_diagnostics_seed <- 2718L
max_monitor_blocks <- 100L
extra_monitor_beta <- 200L

model_args <- switch(
  model,
  
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
  ),
  
  stop(
    "Unknown Bayesian benchmark model: ",
    model
  )
)

output_root <- file.path(
  "outputs",
  "benchmark_runs"
)


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


# Cluster -----

root <- normalizePath(
  ".",
  winslash = "/",
  mustWork = TRUE
)

common_worker_files <- c(
  "R/benchmarks.R",
  "R/var_data.R",
  "R/var_coefficients.R",
  "R/network_blocks.R"
)

model_worker_files <- switch(
  model,
  
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

native_loader <- switch(
  model,
  m3 = "load_m3_cpp",
  half_t = "load_half_t_cpp",
  gigg = "load_gigg_cpp"
)

worker_files <- c(
  common_worker_files,
  model_worker_files,
  "R/benchmark_metrics.R",
  "R/benchmark_bayesian_diagnostics.R"
)

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
        model_args = model_args,
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
          model_args = model_args,
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
        " diagnostic rows.\n",
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