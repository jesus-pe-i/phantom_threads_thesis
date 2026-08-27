# Runs the non-Bayesian benchmark estimators over the full ordinary DGP bank.
#
# Records are fitted in parallel. Each worker fits all requested estimators
# sequentially, retains benchmark metrics and lag-specific network strengths,
# and returns only the compact benchmark output.


rm(
  list = ls()
)


# Settings -----

dgp_ids <- c(
  
  
  # Block diagonal -----
  
  "blockdiag_n6_m3_p2",
  
  "blockdiag_n6_m3_p1",
  "blockdiag_n6_m3_p3",
  "blockdiag_n6_m3_p4",
  "blockdiag_n6_m3_p5",
  
  "blockdiag_n6_m2_p2",
  "blockdiag_n6_m4_p2",
  "blockdiag_n6_m5_p2",
  "blockdiag_n6_m6_p2",
  
  "blockdiag_n4_m3_p2",
  "blockdiag_n8_m3_p2",
  "blockdiag_n10_m3_p2",
  
  
  # Communities -----
  
  "communities_n6_m3_p2",
  
  "communities_n6_m3_p1",
  "communities_n6_m3_p3",
  "communities_n6_m3_p4",
  "communities_n6_m3_p5",
  
  "communities_n6_m2_p2",
  "communities_n6_m4_p2",
  "communities_n6_m5_p2",
  "communities_n6_m6_p2",
  
  "communities_n4_m3_p2",
  "communities_n8_m3_p2",
  "communities_n10_m3_p2",
  
  
  # Core-periphery -----
  
  "coreperiphery_n6_m3_p2",
  
  "coreperiphery_n6_m3_p1",
  "coreperiphery_n6_m3_p3",
  "coreperiphery_n6_m3_p4",
  "coreperiphery_n6_m3_p5",
  
  "coreperiphery_n6_m2_p2",
  "coreperiphery_n6_m4_p2",
  "coreperiphery_n6_m5_p2",
  "coreperiphery_n6_m6_p2",
  
  "coreperiphery_n4_m3_p2",
  "coreperiphery_n8_m3_p2",
  "coreperiphery_n10_m3_p2",
  
  
  # Hubs -----
  
  "hubs_n6_m3_p2",
  
  "hubs_n6_m3_p1",
  "hubs_n6_m3_p3",
  "hubs_n6_m3_p4",
  "hubs_n6_m3_p5",
  
  "hubs_n6_m2_p2",
  "hubs_n6_m4_p2",
  "hubs_n6_m5_p2",
  "hubs_n6_m6_p2",
  
  "hubs_n4_m3_p2",
  "hubs_n8_m3_p2",
  "hubs_n10_m3_p2"
)

models <- c(
  "adelie_gen",
  "bigvar_hlag",
  "mar",
  "nirvar"
)

T_grid <- c(
  50L,
  100L,
  150L,
  200L,
  250L,
  300L,
  350L,
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

workers <- 12L

simulation_burn_in <- 300L
sigma <- 1

fit_seed <- 991L

model_args <- list(
  adelie_gen = list(
    n_threads = 1L
  ),
  bigvar_hlag = list(),
  mar = list(
    n_starts = 3L,
    max_iter = 100L
  ),
  nirvar = list(
    k = 2L
  )
)

output_root <- file.path(
  "outputs",
  "benchmark_runs"
)


# DGP bank -----

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


# Packages -----

required_packages <- c(
  "adelie",
  "BigVAR",
  "mclust"
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

fit_record_non_bayesian <- function(
    record,
    models,
    fit_seed,
    model_args) {
  
  output <- setNames(
    vector(
      "list",
      length(models)
    ),
    models
  )
  
  for (model in models) {
    
    args <- model_args[[model]]
    
    if (is.null(args)) {
      args <- list()
    }
    
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
    
    output[[model]] <- list(
      result = result,
      strengths = fit$s_hat_lag
    )
  }
  
  output
}


# Cluster -----

root <- normalizePath(
  ".",
  winslash = "/",
  mustWork = TRUE
)

worker_files <- c(
  "R/benchmarks.R",
  "R/var_data.R",
  "R/var_coefficients.R",
  "R/network_blocks.R",
  "R/fit_adelie.R",
  "R/fit_bigvar.R",
  "R/fit_mar.R",
  "R/fit_nirvar.R",
  "R/benchmark_metrics.R"
)

campaign_start <- Sys.time()

cluster <- parallel::makePSOCKcluster(
  workers
)

parallel::clusterExport(
  cluster,
  c(
    "root",
    "worker_files"
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
      
      for (file in worker_files) {
        
        source(
          file.path(
            root,
            file
          )
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
        " records x ",
        length(models),
        " models\n",
        sep = ""
      )
      
      output <- parallel::parLapplyLB(
        cluster,
        records,
        fit_record_non_bayesian,
        models = models,
        fit_seed = fit_seed,
        model_args = model_args
      )
      
      
      # Save -----
      
      for (model in models) {
        
        results <- do.call(
          rbind,
          lapply(
            output,
            function(record_output) {
              record_output[[model]]$result
            }
          )
        )
        
        rownames(results) <- NULL
        
        strengths <- lapply(
          output,
          function(record_output) {
            record_output[[model]]$strengths
          }
        )
        
        names(strengths) <- names(records)
        
        args <- model_args[[model]]
        
        if (is.null(args)) {
          args <- list()
        }
        
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
            model_args = args
          ),
          results = results,
          strengths = strengths
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
      }
      
      cat(
        "Saved in ",
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
        output
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


# Runtime -----

cat(
  "\nFull campaign wall time: ",
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