# NIRVAR benchmark for block VAR network recovery.
#
# Builds feature-wise correlation networks, estimates an unfolded spectral
# embedding, clusters feature-unit nodes, constructs NIRVAR restrictions,
# and fits the resulting restricted VAR by equation-wise least squares.


# Data construction -----

make_nirvar_array <- function(
    Y,
    n_units,
    m) {
  
  T_obs <- nrow(Y)
  
  Y_array <- array(
    NA_real_,
    dim = c(
      T_obs,
      n_units,
      m
    )
  )
  
  for (unit in seq_len(n_units)) {
    
    columns <- (
      (unit - 1L) * m + 1L
    ):(
      unit * m
    )
    
    Y_array[, unit, ] <- Y[
      ,
      columns,
      drop = FALSE
    ]
  }
  
  Y_array
}


# Correlation networks -----

nirvar_pearson_correlations <- function(
    Y_array) {
  
  n_units <- dim(Y_array)[2L]
  m <- dim(Y_array)[3L]
  
  correlation <- array(
    NA_real_,
    dim = c(
      m,
      n_units,
      n_units
    )
  )
  
  for (feature in seq_len(m)) {
    
    Y_feature <- Y_array[
      ,
      ,
      feature,
      drop = FALSE
    ][, , 1L]
    
    C <- stats::cor(
      Y_feature,
      use = "everything",
      method = "pearson"
    )
    
    if (any(!is.finite(C))) {
      stop(
        "NIRVAR produced non-finite correlations"
      )
    }
    
    C <- (
      C +
        t(C)
    ) / 2
    
    correlation[feature, , ] <- C
  }
  
  correlation
}


# Embedding dimension -----

nirvar_mp_dimension <- function(
    correlation,
    T_obs,
    d = NULL) {
  
  m <- dim(correlation)[1L]
  n_units <- dim(correlation)[2L]
  
  mp_cutoff <- (
    1 +
      sqrt(
        n_units /
          T_obs
      )
  )^2
  
  eigenvalues_by_feature <- matrix(
    NA_real_,
    nrow = m,
    ncol = n_units
  )
  
  for (feature in seq_len(m)) {
    
    C <- correlation[
      feature,
      ,
      ,
      drop = FALSE
    ][1L, , ]
    
    C <- (
      C +
        t(C)
    ) / 2
    
    eigenvalues_by_feature[feature, ] <-
      eigen(
        C,
        symmetric = TRUE,
        only.values = TRUE
      )$values
  }
  
  d_mp_by_feature <- as.integer(
    rowSums(
      eigenvalues_by_feature >
        mp_cutoff
    )
  )
  
  d_cap <- max(
    1L,
    n_units - 1L
  )
  
  d_mp <- max(
    d_mp_by_feature
  )
  
  if (is.null(d)) {
    
    d_used <- max(
      1L,
      min(
        d_cap,
        d_mp
      )
    )
    
    dimension_rule <-
      "featurewise_max_mp"
    
    minimum_dimension_fallback <-
      all(
        d_mp_by_feature == 0L
      )
    
  } else {
    
    d_used <- as.integer(d)
    
    if (
      d_used < 1L ||
      d_used > d_cap
    ) {
      stop(
        "d must be between 1 and n_units - 1"
      )
    }
    
    dimension_rule <- "fixed"
    minimum_dimension_fallback <- FALSE
  }
  
  list(
    d = as.integer(d_used),
    d_mp = as.integer(d_mp),
    d_mp_by_feature = d_mp_by_feature,
    mp_cutoff = as.numeric(mp_cutoff),
    eigenvalues_by_feature = eigenvalues_by_feature,
    signal_features = which(
      d_mp_by_feature > 0L
    ),
    dimension_rule = dimension_rule,
    minimum_dimension_fallback =
      minimum_dimension_fallback,
    d_cap = as.integer(d_cap)
  )
}


# UASE embedding -----

nirvar_uase_embedding <- function(
    correlation,
    d) {
  
  m <- dim(correlation)[1L]
  n_units <- dim(correlation)[2L]
  
  flat_correlation <- do.call(
    rbind,
    lapply(
      seq_len(m),
      function(feature) {
        correlation[feature, , ]
      }
    )
  )
  
  svd_fit <- base::svd(
    flat_correlation,
    nu = d,
    nv = 0L
  )
  
  embedding_matrix <- sweep(
    svd_fit$u[
      ,
      seq_len(d),
      drop = FALSE
    ],
    2L,
    svd_fit$d[
      seq_len(d)
    ],
    `*`
  )
  
  embedding <- array(
    NA_real_,
    dim = c(
      m,
      n_units,
      d
    )
  )
  
  for (feature in seq_len(m)) {
    
    rows <- (
      (feature - 1L) *
        n_units +
        1L
    ):(
      feature *
        n_units
    )
    
    embedding[feature, , ] <-
      embedding_matrix[
        rows,
        ,
        drop = FALSE
      ]
  }
  
  list(
    embedding = embedding,
    embedding_matrix = embedding_matrix,
    flat_correlation = flat_correlation,
    singular_values = as.numeric(
      svd_fit$d
    )
  )
}


make_nirvar_embedding <- function(
    Y_array,
    d = NULL) {
  
  correlation <-
    nirvar_pearson_correlations(
      Y_array
    )
  
  dimension_fit <-
    nirvar_mp_dimension(
      correlation = correlation,
      T_obs = dim(Y_array)[1L],
      d = d
    )
  
  embedding_fit <-
    nirvar_uase_embedding(
      correlation = correlation,
      d = dimension_fit$d
    )
  
  list(
    correlation = correlation,
    
    embedding =
      embedding_fit$embedding,
    
    embedding_matrix =
      embedding_fit$embedding_matrix,
    
    flat_correlation =
      embedding_fit$flat_correlation,
    
    singular_values =
      embedding_fit$singular_values,
    
    d =
      dimension_fit$d,
    
    d_mp =
      dimension_fit$d_mp,
    
    d_mp_by_feature =
      dimension_fit$d_mp_by_feature,
    
    mp_cutoff =
      dimension_fit$mp_cutoff,
    
    eigenvalues_by_feature =
      dimension_fit$eigenvalues_by_feature,
    
    signal_features =
      dimension_fit$signal_features,
    
    dimension_rule =
      dimension_fit$dimension_rule,
    
    minimum_dimension_fallback =
      dimension_fit$minimum_dimension_fallback,
    
    d_cap =
      dimension_fit$d_cap
  )
}


# Clustering -----

nirvar_canonicalize_labels <- function(
    labels,
    embedding_matrix) {
  
  labels <- as.integer(
    labels
  )
  
  cluster_ids <- sort(
    unique(
      labels
    )
  )
  
  centers <- do.call(
    rbind,
    lapply(
      cluster_ids,
      function(cluster) {
        
        colMeans(
          embedding_matrix[
            labels == cluster,
            ,
            drop = FALSE
          ]
        )
      }
    )
  )
  
  centers <- as.matrix(
    centers
  )
  
  order_arguments <- lapply(
    seq_len(
      ncol(centers)
    ),
    function(column) {
      centers[, column]
    }
  )
  
  ordering <- do.call(
    order,
    c(
      order_arguments,
      list(
        cluster_ids
      )
    )
  )
  
  ordered_cluster_ids <-
    cluster_ids[
      ordering
    ]
  
  label_map <- setNames(
    seq_along(
      ordered_cluster_ids
    ),
    ordered_cluster_ids
  )
  
  as.integer(
    label_map[
      as.character(
        labels
      )
    ]
  )
}


nirvar_cluster_embedding <- function(
    embedding_matrix,
    m,
    n_units,
    k,
    seed = 991L) {
  
  if (!requireNamespace(
    "mclust",
    quietly = TRUE
  )) {
    stop(
      "Package 'mclust' is required for NIRVAR"
    )
  }
  
  embedding_matrix <- as.matrix(
    embedding_matrix
  )
  
  k <- as.integer(k)
  seed <- as.integer(seed)
  
  n_rows <- m *
    n_units
  
  distinct_rows <- nrow(
    unique(
      as.data.frame(
        embedding_matrix
      )
    )
  )
  
  k_attempted <- min(
    k,
    n_rows,
    distinct_rows
  )
  
  fallback_reason <- NULL
  bic <- NA_real_
  loglik <- NA_real_
  
  uncertainty <- rep(
    NA_real_,
    n_rows
  )
  
  
  # Degenerate single-cluster case -----
  
  if (k_attempted <= 1L) {
    
    labels <- rep(
      1L,
      n_rows
    )
    
    method <- if (
      k == 1L
    ) {
      "single_cluster"
    } else {
      "single_cluster_distinct_rows"
    }
    
    model_name <- NA_character_
    
    if (k > 1L) {
      fallback_reason <- paste0(
        "Only ",
        distinct_rows,
        " distinct embedding row(s)"
      )
    }
    
    
    # Gaussian mixture -----
    
  } else {
    
    d <- ncol(
      embedding_matrix
    )
    
    model_name <- if (
      d == 1L
    ) {
      "V"
    } else {
      "VVV"
    }
    
    gmm_data <- if (
      d == 1L
    ) {
      as.numeric(
        embedding_matrix[, 1L]
      )
    } else {
      embedding_matrix
    }
    
    set.seed(
      seed
    )
    
    gmm_result <- tryCatch(
      {
        
        fit <- suppressWarnings(
          mclust::Mclust(
            data = gmm_data,
            G = k_attempted,
            modelNames = model_name,
            verbose = FALSE
          )
        )
        
        if (
          is.null(fit) ||
          is.null(
            fit$classification
          )
        ) {
          stop(
            "mclust returned no classification"
          )
        }
        
        labels_candidate <- as.integer(
          fit$classification
        )
        
        if (
          length(labels_candidate) !=
          n_rows ||
          anyNA(labels_candidate) ||
          length(
            unique(
              labels_candidate
            )
          ) != k_attempted
        ) {
          stop(
            "mclust returned an invalid classification"
          )
        }
        
        list(
          success = TRUE,
          fit = fit,
          labels = labels_candidate,
          error = NULL
        )
      },
      error = function(error) {
        
        list(
          success = FALSE,
          fit = NULL,
          labels = NULL,
          error = conditionMessage(
            error
          )
        )
      }
    )
    
    
    # GMM result -----
    
    if (isTRUE(
      gmm_result$success
    )) {
      
      labels <-
        gmm_result$labels
      
      method <- "gmm"
      
      bic_values <- as.numeric(
        gmm_result$fit$bic
      )
      
      if (
        any(
          is.finite(
            bic_values
          )
        )
      ) {
        
        bic <- max(
          bic_values[
            is.finite(
              bic_values
            )
          ]
        )
      }
      
      if (!is.null(
        gmm_result$fit$loglik
      )) {
        
        loglik <- as.numeric(
          gmm_result$fit$loglik
        )
      }
      
      if (!is.null(
        gmm_result$fit$uncertainty
      )) {
        
        uncertainty <- as.numeric(
          gmm_result$fit$uncertainty
        )
      }
      
      
      # K-means fallback -----
      
    } else {
      
      fallback_reason <-
        gmm_result$error
      
      set.seed(
        seed
      )
      
      kmeans_result <- tryCatch(
        {
          
          fit <- suppressWarnings(
            stats::kmeans(
              x = embedding_matrix,
              centers = k_attempted,
              iter.max = 100L,
              nstart = 1L,
              algorithm = "Hartigan-Wong"
            )
          )
          
          labels_candidate <- as.integer(
            fit$cluster
          )
          
          if (
            length(labels_candidate) !=
            n_rows ||
            anyNA(labels_candidate) ||
            length(
              unique(
                labels_candidate
              )
            ) != k_attempted
          ) {
            stop(
              "kmeans returned an invalid classification"
            )
          }
          
          list(
            success = TRUE,
            labels = labels_candidate,
            error = NULL
          )
        },
        error = function(error) {
          
          list(
            success = FALSE,
            labels = NULL,
            error = conditionMessage(
              error
            )
          )
        }
      )
      
      if (isTRUE(
        kmeans_result$success
      )) {
        
        labels <-
          kmeans_result$labels
        
        method <-
          "kmeans_fallback"
        
      } else {
        
        labels <- rep(
          1L,
          n_rows
        )
        
        method <-
          "single_cluster_fallback"
        
        fallback_reason <- paste(
          "GMM:",
          fallback_reason,
          "| kmeans:",
          kmeans_result$error
        )
      }
    }
  }
  
  
  # Canonical labels -----
  
  labels <- nirvar_canonicalize_labels(
    labels = labels,
    embedding_matrix =
      embedding_matrix
  )
  
  k_used <- length(
    unique(
      labels
    )
  )
  
  labels_matrix <- matrix(
    labels,
    nrow = m,
    ncol = n_units,
    byrow = TRUE
  )
  
  cluster_sizes <- as.integer(
    table(
      factor(
        labels,
        levels = seq_len(
          k_used
        )
      )
    )
  )
  
  if (
    length(uncertainty) !=
    n_rows
  ) {
    
    uncertainty <- rep(
      NA_real_,
      n_rows
    )
  }
  
  list(
    labels = labels_matrix,
    labels_vector = labels,
    
    k_requested =
      k,
    
    k_attempted =
      as.integer(
        k_attempted
      ),
    
    k_used =
      as.integer(
        k_used
      ),
    
    cluster_sizes =
      cluster_sizes,
    
    method =
      method,
    
    model_name =
      model_name,
    
    seed =
      seed,
    
    bic =
      bic,
    
    loglik =
      loglik,
    
    uncertainty =
      uncertainty,
    
    mean_uncertainty = if (
      any(
        is.finite(
          uncertainty
        )
      )
    ) {
      
      mean(
        uncertainty[
          is.finite(
            uncertainty
          )
        ]
      )
      
    } else {
      
      NA_real_
    },
    
    max_uncertainty = if (
      any(
        is.finite(
          uncertainty
        )
      )
    ) {
      
      max(
        uncertainty[
          is.finite(
            uncertainty
          )
        ]
      )
      
    } else {
      
      NA_real_
    },
    
    fallback_reason =
      fallback_reason,
    
    distinct_embedding_rows =
      as.integer(
        distinct_rows
      )
  )
}


# Restrictions -----

make_nirvar_restriction <- function(
    labels) {
  
  labels <- as.matrix(
    labels
  )
  
  m <- nrow(labels)
  n_units <- ncol(labels)
  
  restriction <- array(
    0L,
    dim = c(
      m,
      m,
      n_units,
      n_units
    ),
    dimnames = list(
      target_feature =
        seq_len(m),
      
      predictor_feature =
        seq_len(m),
      
      receiver =
        seq_len(n_units),
      
      sender =
        seq_len(n_units)
    )
  )
  
  for (target_feature in seq_len(m)) {
    
    for (predictor_feature in seq_len(m)) {
      
      restriction[
        target_feature,
        predictor_feature,
        ,
      ] <- outer(
        labels[
          target_feature,
        ],
        labels[
          predictor_feature,
        ],
        FUN = "=="
      )
    }
  }
  
  storage.mode(
    restriction
  ) <- "integer"
  
  restriction
}


make_nirvar_restrictions <- function(
    Y_array,
    k,
    d = NULL,
    seed = 991L) {
  
  embedding_fit <-
    make_nirvar_embedding(
      Y_array = Y_array,
      d = d
    )
  
  cluster_fit <-
    nirvar_cluster_embedding(
      embedding_matrix =
        embedding_fit$embedding_matrix,
      m = dim(Y_array)[3L],
      n_units = dim(Y_array)[2L],
      k = k,
      seed = seed
    )
  
  restriction <-
    make_nirvar_restriction(
      labels =
        cluster_fit$labels
    )
  
  list(
    embedding = embedding_fit,
    clustering = cluster_fit,
    restriction = restriction
  )
}


# Restricted least squares -----

nirvar_ols_no_intercept <- function(
    X,
    y,
    tol = 1e-10) {
  
  X <- as.matrix(
    X
  )
  
  y <- as.numeric(
    y
  )
  
  if (ncol(X) == 0L) {
    
    fitted <- rep(
      0,
      length(y)
    )
    
    return(
      list(
        coefficients =
          numeric(0L),
        
        fitted =
          fitted,
        
        residuals =
          y,
        
        rank =
          0L,
        
        singular_values =
          numeric(0L),
        
        condition_number =
          NA_real_,
        
        cutoff =
          NA_real_,
        
        solver =
          "no_predictors"
      )
    )
  }
  
  svd_fit <- base::svd(
    X,
    nu = min(
      dim(X)
    ),
    nv = min(
      dim(X)
    )
  )
  
  largest_singular_value <- max(
    svd_fit$d
  )
  
  cutoff <-
    tol *
    max(
      dim(X)
    ) *
    largest_singular_value
  
  retained <- which(
    svd_fit$d >
      cutoff
  )
  
  rank <- length(
    retained
  )
  
  if (rank == 0L) {
    
    coefficients <- rep(
      0,
      ncol(X)
    )
    
  } else {
    
    projected_response <- crossprod(
      svd_fit$u[
        ,
        retained,
        drop = FALSE
      ],
      y
    )
    
    scaled_projection <-
      projected_response /
      svd_fit$d[
        retained
      ]
    
    coefficients <- drop(
      svd_fit$v[
        ,
        retained,
        drop = FALSE
      ] %*%
        scaled_projection
    )
  }
  
  fitted <- drop(
    X %*%
      coefficients
  )
  
  residuals <-
    y -
    fitted
  
  condition_number <- if (
    rank >= 1L
  ) {
    
    largest_singular_value /
      min(
        svd_fit$d[
          retained
        ]
      )
    
  } else {
    
    NA_real_
  }
  
  list(
    coefficients =
      coefficients,
    
    fitted =
      fitted,
    
    residuals =
      residuals,
    
    rank =
      as.integer(rank),
    
    singular_values =
      as.numeric(
        svd_fit$d
      ),
    
    condition_number =
      as.numeric(
        condition_number
      ),
    
    cutoff =
      as.numeric(
        cutoff
      ),
    
    solver =
      "svd_minimum_norm"
  )
}


fit_nirvar_restricted_var <- function(
    prepared,
    restriction,
    tol = 1e-10) {
  
  m <- prepared$m
  n_units <- prepared$n_units
  N <- prepared$n_series
  p_lags <- prepared$p_lags
  T_eff <- prepared$T_eff
  
  restriction <- as.array(
    restriction
  )
  
  restriction <-
    restriction != 0
  
  A_list <- lapply(
    seq_len(p_lags),
    function(lag) {
      
      matrix(
        0,
        nrow = N,
        ncol = N
      )
    }
  )
  
  coefficient_array <- array(
    0,
    dim = c(
      p_lags,
      m,
      m,
      n_units,
      n_units
    ),
    dimnames = list(
      lag =
        seq_len(p_lags),
      
      target_feature =
        seq_len(m),
      
      predictor_feature =
        seq_len(m),
      
      receiver =
        seq_len(n_units),
      
      sender =
        seq_len(n_units)
    )
  )
  
  fitted_values <- matrix(
    0,
    nrow = T_eff,
    ncol = N
  )
  
  residuals <- matrix(
    0,
    nrow = T_eff,
    ncol = N
  )
  
  equation_n_predictors <-
    integer(N)
  
  equation_rank <-
    integer(N)
  
  equation_condition_number <- rep(
    NA_real_,
    N
  )
  
  equation_solver <-
    character(N)
  
  
  # Equation-wise restricted OLS -----
  
  for (receiver in seq_len(n_units)) {
    
    for (target_feature in seq_len(m)) {
      
      equation <- (
        receiver - 1L
      ) * m +
        target_feature
      
      active_base <- logical(
        N
      )
      
      for (sender in seq_len(n_units)) {
        
        for (predictor_feature in seq_len(m)) {
          
          predictor <- (
            sender - 1L
          ) * m +
            predictor_feature
          
          active_base[predictor] <-
            restriction[
              target_feature,
              predictor_feature,
              receiver,
              sender
            ]
        }
      }
      
      active_base_indices <- which(
        active_base
      )
      
      active_columns <- unlist(
        lapply(
          seq_len(p_lags),
          function(lag) {
            
            (
              lag - 1L
            ) * N +
              active_base_indices
          }
        ),
        use.names = FALSE
      )
      
      X_equation <- prepared$X[
        ,
        active_columns,
        drop = FALSE
      ]
      
      y_equation <- prepared$Y[
        ,
        equation
      ]
      
      equation_fit <-
        nirvar_ols_no_intercept(
          X = X_equation,
          y = y_equation,
          tol = tol
        )
      
      fitted_values[, equation] <-
        equation_fit$fitted
      
      residuals[, equation] <-
        equation_fit$residuals
      
      equation_n_predictors[equation] <-
        length(
          active_columns
        )
      
      equation_rank[equation] <-
        equation_fit$rank
      
      equation_condition_number[equation] <-
        equation_fit$condition_number
      
      equation_solver[equation] <-
        equation_fit$solver
      
      for (index in seq_along(active_columns)) {
        
        design_column <-
          active_columns[index]
        
        lag <- (
          (
            design_column -
              1L
          ) %/%
            N
        ) + 1L
        
        predictor <- (
          (
            design_column -
              1L
          ) %%
            N
        ) + 1L
        
        sender <- (
          (
            predictor -
              1L
          ) %/%
            m
        ) + 1L
        
        predictor_feature <- (
          (
            predictor -
              1L
          ) %%
            m
        ) + 1L
        
        coefficient <-
          equation_fit$coefficients[
            index
          ]
        
        A_list[[lag]][
          equation,
          predictor
        ] <- coefficient
        
        coefficient_array[
          lag,
          target_feature,
          predictor_feature,
          receiver,
          sender
        ] <- coefficient
      }
    }
  }
  
  Sigma_hat <- crossprod(
    residuals
  ) /
    T_eff
  
  equation_diagnostics <- data.frame(
    equation =
      seq_len(N),
    
    receiver =
      rep(
        seq_len(n_units),
        each = m
      ),
    
    target_feature =
      rep(
        seq_len(m),
        times = n_units
      ),
    
    n_predictors =
      equation_n_predictors,
    
    rank =
      equation_rank,
    
    rank_deficient =
      equation_rank <
      equation_n_predictors,
    
    condition_number =
      equation_condition_number,
    
    solver =
      equation_solver,
    
    stringsAsFactors =
      FALSE
  )
  
  list(
    A_list =
      A_list,
    
    coefficient_array =
      coefficient_array,
    
    fitted_values =
      fitted_values,
    
    residuals =
      residuals,
    
    Sigma_hat =
      Sigma_hat,
    
    equation_diagnostics =
      equation_diagnostics
  )
}


# Model fit -----

fit_nirvar <- function(
    Y_list,
    p_lags,
    k = 2L,
    d = NULL,
    seed = 991L,
    ols_tol = 1e-10,
    selected_eps = 1e-12,
    keep_details = FALSE) {
  
  
  # Data -----
  
  prepared <- prepare_var_data(
    Y_list = Y_list,
    p_lags = p_lags,
    standardize = FALSE
  )
  
  Y_array <- make_nirvar_array(
    Y = prepared$Y_full,
    n_units =
      prepared$n_units,
    m =
      prepared$m
  )
  
  
  # Network restrictions -----
  
  restriction_fit <-
    make_nirvar_restrictions(
      Y_array = Y_array,
      k = k,
      d = d,
      seed = seed
    )
  
  
  # Restricted VAR -----
  
  regression_fit <-
    fit_nirvar_restricted_var(
      prepared = prepared,
      restriction =
        restriction_fit$restriction,
      tol = ols_tol
    )
  
  A_hat_lag <-
    regression_fit$A_list
  
  beta_hat <- A_list_to_beta(
    A_hat_lag
  )
  
  
  # Network estimate -----
  
  s_hat_lag <- A_list_to_s_lag(
    A_list =
      A_hat_lag,
    n_units =
      prepared$n_units,
    m =
      prepared$m
  )
  
  s_hat_unit_max <- s_lag_to_unit(
    s_lag =
      s_hat_lag,
    method =
      "max"
  )
  
  s_hat_unit_rms <- s_lag_to_unit(
    s_lag =
      s_hat_lag,
    method =
      "rms"
  )
  
  
  # Output -----
  
  out <- list(
    model =
      "nirvar",
    
    backend =
      "nirvar",
    
    beta_hat =
      beta_hat,
    
    A_hat_lag =
      A_hat_lag,
    
    s_hat_lag =
      s_hat_lag,
    
    s_hat_unit_max =
      s_hat_unit_max,
    
    s_hat_unit_rms =
      s_hat_unit_rms,
    
    score_lag =
      s_hat_lag,
    
    selected_lag =
      s_hat_lag >
      selected_eps,
    
    coefficient_array =
      regression_fit$coefficient_array,
    
    labels =
      restriction_fit$clustering$labels,
    
    restriction =
      restriction_fit$restriction,
    
    Sigma_hat =
      regression_fit$Sigma_hat,
    
    equation_diagnostics =
      regression_fit$equation_diagnostics,
    
    d =
      restriction_fit$embedding$d,
    
    d_mp =
      restriction_fit$embedding$d_mp,
    
    d_mp_by_feature =
      restriction_fit$embedding$d_mp_by_feature,
    
    mp_cutoff =
      restriction_fit$embedding$mp_cutoff,
    
    dimension_rule =
      restriction_fit$embedding$dimension_rule,
    
    minimum_dimension_fallback =
      restriction_fit$embedding$minimum_dimension_fallback,
    
    k_requested =
      restriction_fit$clustering$k_requested,
    
    k_attempted =
      restriction_fit$clustering$k_attempted,
    
    k_used =
      restriction_fit$clustering$k_used,
    
    cluster_sizes =
      restriction_fit$clustering$cluster_sizes,
    
    cluster_method =
      restriction_fit$clustering$method,
    
    cluster_model =
      restriction_fit$clustering$model_name,
    
    cluster_bic =
      restriction_fit$clustering$bic,
    
    cluster_loglik =
      restriction_fit$clustering$loglik,
    
    mean_cluster_uncertainty =
      restriction_fit$clustering$mean_uncertainty,
    
    max_cluster_uncertainty =
      restriction_fit$clustering$max_uncertainty,
    
    cluster_fallback_reason =
      restriction_fit$clustering$fallback_reason,
    
    fit_seed =
      as.integer(seed),
    
    dimensions = list(
      n_units =
        prepared$n_units,
      
      m =
        prepared$m,
      
      N =
        prepared$n_series,
      
      p_lags =
        prepared$p_lags,
      
      k =
        prepared$n_coef,
      
      T_full =
        prepared$T_full,
      
      T_p =
        prepared$T_eff
    ),
    
    preprocessing = list(
      standardize =
        FALSE,
      
      centering =
        FALSE,
      
      scale =
        prepared$scale,
      
      series_names =
        prepared$series_names
    ),
    
    control = list(
      k =
        as.integer(k),
      
      d =
        d,
      
      ols_tol =
        as.numeric(
          ols_tol
        ),
      
      selected_eps =
        as.numeric(
          selected_eps
        )
    )
  )
  
  if (keep_details) {
    
    out$embedding <-
      restriction_fit$embedding
    
    out$clustering <-
      restriction_fit$clustering
    
    out$fitted_values <-
      regression_fit$fitted_values
    
    out$residuals <-
      regression_fit$residuals
  }
  
  out
}