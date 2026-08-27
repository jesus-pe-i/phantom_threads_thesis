# Reviews the frozen DGP bank and writes one descriptive row per DGP.
#
# Reads each A_list.rds directly, recovers dimensions from the DGP ID,
# checks coefficient dimensions and stability, and summarizes network
# density, block strengths, lag structure, sender/receiver activity, and
# weakly connected component structure.


rm(
  list = ls()
)


# Sources -----

source(
  "R/network_blocks.R"
)

source(
  "R/simulation.R"
)


# Paths -----

dgp_root <- file.path(
  "data",
  "dgp_bank"
)

output_dir <- file.path(
  "data",
  "outputs"
)

output_path <- file.path(
  output_dir,
  "dgp_bank_review.csv"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Settings -----

truth_eps <- 1e-12


# Helpers -----

parse_dgp_id <- function(
    dgp_id) {
  
  pattern <- paste0(
    "^(.*)_n",
    "([0-9]+)",
    "_m",
    "([0-9]+)",
    "_p",
    "([0-9]+)$"
  )
  
  match <- regexec(
    pattern,
    dgp_id
  )
  
  parts <- regmatches(
    dgp_id,
    match
  )[[1L]]
  
  if (length(parts) != 5L) {
    stop(
      "Unable to parse DGP ID: ",
      dgp_id
    )
  }
  
  list(
    family =
      parts[2L],
    
    n_units =
      as.integer(
        parts[3L]
      ),
    
    m =
      as.integer(
        parts[4L]
      ),
    
    p_lags =
      as.integer(
        parts[5L]
      )
  )
}


summarize_dgp <- function(
    dgp_dir) {
  
  dgp_id <- basename(
    dgp_dir
  )
  
  parsed <- parse_dgp_id(
    dgp_id
  )
  
  A_path <- file.path(
    dgp_dir,
    "A_list.rds"
  )
  
  if (!file.exists(
    A_path
  )) {
    stop(
      "Missing A_list.rds for ",
      dgp_id
    )
  }
  
  A_list <- readRDS(
    A_path
  )
  
  validate_A_list(
    A_list
  )
  
  n_units <- parsed$n_units
  m <- parsed$m
  p_lags <- parsed$p_lags
  
  N <-
    n_units *
    m
  
  
  ## Dimensions -----
  
  actual_p <- length(
    A_list
  )
  
  actual_N <- nrow(
    A_list[[1L]]
  )
  
  dimensions_match <-
    actual_p == p_lags &&
    actual_N == N &&
    all(
      vapply(
        A_list,
        function(A) {
          
          identical(
            dim(A),
            as.integer(
              c(
                N,
                N
              )
            )
          )
        },
        logical(1L)
      )
    )
  
  
  ## Stability -----
  
  stability <- check_A_radius(
    A_list
  )
  
  
  ## Network truth -----
  
  s_lag <- A_list_to_s_lag(
    A_list =
      A_list,
    n_units =
      n_units,
    m =
      m
  )
  
  G_lag <-
    s_lag >
    truth_eps
  
  G_unit <- G_lag_to_G_unit(
    G_lag
  )
  
  s_unit_max <- s_lag_to_unit(
    s_lag =
      s_lag,
    method =
      "max"
  )
  
  s_unit_rms <- s_lag_to_unit(
    s_lag =
      s_lag,
    method =
      "rms"
  )
  
  
  ## Block masks -----
  
  self_mask <- row(
    G_unit
  ) ==
    col(
      G_unit
    )
  
  foreign_mask <- !self_mask
  
  active_self_unit <-
    G_unit &
    self_mask
  
  active_foreign_unit <-
    G_unit &
    foreign_mask
  
  active_lag <- G_lag
  
  lag_self_mask <- array(
    self_mask,
    dim =
      c(
        n_units,
        n_units,
        p_lags
      )
  )
  
  lag_foreign_mask <- !lag_self_mask
  
  active_self_lag <-
    active_lag &
    lag_self_mask
  
  active_foreign_lag <-
    active_lag &
    lag_foreign_mask
  
  
  ## Strength summaries -----
  
  active_strengths <- s_lag[
    active_lag
  ]
  
  self_strengths <- s_lag[
    active_self_lag
  ]
  
  foreign_strengths <- s_lag[
    active_foreign_lag
  ]
  
  mean_active_strength <- if (
    length(
      active_strengths
    ) > 0L
  ) {
    mean(
      active_strengths
    )
  } else {
    0
  }
  
  median_active_strength <- if (
    length(
      active_strengths
    ) > 0L
  ) {
    median(
      active_strengths
    )
  } else {
    0
  }
  
  min_active_strength <- if (
    length(
      active_strengths
    ) > 0L
  ) {
    min(
      active_strengths
    )
  } else {
    0
  }
  
  max_active_strength <- if (
    length(
      active_strengths
    ) > 0L
  ) {
    max(
      active_strengths
    )
  } else {
    0
  }
  
  mean_self_strength <- if (
    length(
      self_strengths
    ) > 0L
  ) {
    mean(
      self_strengths
    )
  } else {
    0
  }
  
  mean_foreign_strength <- if (
    length(
      foreign_strengths
    ) > 0L
  ) {
    mean(
      foreign_strengths
    )
  } else {
    0
  }
  
  
  ## Degree structure -----
  
  foreign_unit <- G_unit
  
  diag(
    foreign_unit
  ) <- FALSE
  
  receiver_degree <- rowSums(
    foreign_unit
  )
  
  sender_degree <- colSums(
    foreign_unit
  )
  
  zero_foreign_degree <- (
    receiver_degree +
      sender_degree
  ) == 0L
  
  undirected_foreign <-
    foreign_unit |
    t(
      foreign_unit
    )
  
  visited <- rep(
    FALSE,
    n_units
  )
  
  weak_component_sizes <- integer(0L)
  
  for (start in seq_len(
    n_units
  )) {
    
    if (visited[start]) {
      next
    }
    
    queue <- start
    visited[start] <- TRUE
    component_size <- 0L
    
    while (length(
      queue
    ) > 0L) {
      
      node <- queue[1L]
      queue <- queue[-1L]
      component_size <-
        component_size +
        1L
      
      neighbours <- which(
        undirected_foreign[
          node,
          
        ]
      )
      
      new_neighbours <- neighbours[
        !visited[neighbours]
      ]
      
      if (length(
        new_neighbours
      ) > 0L) {
        
        visited[new_neighbours] <- TRUE
        queue <- c(
          queue,
          new_neighbours
        )
      }
    }
    
    weak_component_sizes <- c(
      weak_component_sizes,
      component_size
    )
  }
  
  weak_component_count <- length(
    weak_component_sizes
  )
  
  largest_weak_component_size <- max(
    weak_component_sizes
  )
  
  largest_weak_component_share <-
    largest_weak_component_size /
    n_units
  
  
  ## Main row -----
  
  output <- data.frame(
    dgp_id =
      dgp_id,
    
    family =
      parsed$family,
    
    n_units =
      n_units,
    
    m =
      m,
    
    p_lags =
      p_lags,
    
    N =
      N,
    
    n_coef =
      N *
      N *
      p_lags,
    
    dimensions_match =
      dimensions_match,
    
    spectral_radius =
      stability$radius,
    
    stable =
      stability$radius <
      1,
    
    active_lag_blocks =
      sum(
        active_lag
      ),
    
    total_lag_blocks =
      n_units *
      n_units *
      p_lags,
    
    lag_density =
      mean(
        active_lag
      ),
    
    active_unit_blocks =
      sum(
        G_unit
      ),
    
    total_unit_blocks =
      n_units *
      n_units,
    
    unit_density =
      mean(
        G_unit
      ),
    
    active_self_lag_blocks =
      sum(
        active_self_lag
      ),
    
    active_foreign_lag_blocks =
      sum(
        active_foreign_lag
      ),
    
    active_self_unit_blocks =
      sum(
        active_self_unit
      ),
    
    active_foreign_unit_blocks =
      sum(
        active_foreign_unit
      ),
    
    foreign_unit_density =
      sum(
        active_foreign_unit
      ) /
      (
        n_units *
          (
            n_units -
              1L
          )
      ),
    
    min_active_strength =
      min_active_strength,
    
    mean_active_strength =
      mean_active_strength,
    
    median_active_strength =
      median_active_strength,
    
    max_active_strength =
      max_active_strength,
    
    mean_self_strength =
      mean_self_strength,
    
    mean_foreign_strength =
      mean_foreign_strength,
    
    mean_unit_rms_strength =
      mean(
        s_unit_rms[
          G_unit
        ]
      ),
    
    max_unit_rms_strength =
      max(
        s_unit_rms
      ),
    
    max_unit_strength =
      max(
        s_unit_max
      ),
    
    mean_receiver_degree =
      mean(
        receiver_degree
      ),
    
    max_receiver_degree =
      max(
        receiver_degree
      ),
    
    mean_sender_degree =
      mean(
        sender_degree
      ),
    
    max_sender_degree =
      max(
        sender_degree
      ),
    
    active_receivers =
      sum(
        receiver_degree >
          0L
      ),
    
    active_senders =
      sum(
        sender_degree >
          0L
      ),
    
    zero_foreign_degree_units =
      sum(
        zero_foreign_degree
      ),
    
    weak_component_count =
      weak_component_count,
    
    largest_weak_component_size =
      largest_weak_component_size,
    
    largest_weak_component_share =
      largest_weak_component_share,
    
    stringsAsFactors =
      FALSE
  )
  
  
  ## Lag-specific summaries -----
  
  for (lag in seq_len(
    p_lags
  )) {
    
    lag_active <-
      G_lag[
        ,
        ,
        lag
      ]
    
    lag_strength <-
      s_lag[
        ,
        ,
        lag
      ]
    
    lag_foreign <-
      lag_active &
      foreign_mask
    
    active_values <-
      lag_strength[
        lag_active
      ]
    
    foreign_values <-
      lag_strength[
        lag_foreign
      ]
    
    output[[
      paste0(
        "lag",
        lag,
        "_active_blocks"
      )
    ]] <- sum(
      lag_active
    )
    
    output[[
      paste0(
        "lag",
        lag,
        "_foreign_blocks"
      )
    ]] <- sum(
      lag_foreign
    )
    
    output[[
      paste0(
        "lag",
        lag,
        "_density"
      )
    ]] <- mean(
      lag_active
    )
    
    output[[
      paste0(
        "lag",
        lag,
        "_mean_active_strength"
      )
    ]] <- if (
      length(
        active_values
      ) > 0L
    ) {
      mean(
        active_values
      )
    } else {
      0
    }
    
    output[[
      paste0(
        "lag",
        lag,
        "_mean_foreign_strength"
      )
    ]] <- if (
      length(
        foreign_values
      ) > 0L
    ) {
      mean(
        foreign_values
      )
    } else {
      0
    }
  }
  
  output
}


# DGP bank -----

dgp_dirs <- list.dirs(
  dgp_root,
  recursive = FALSE,
  full.names = TRUE
)

dgp_dirs <- dgp_dirs[
  file.exists(
    file.path(
      dgp_dirs,
      "A_list.rds"
    )
  )
]

if (length(
  dgp_dirs
) == 0L) {
  stop(
    "No frozen DGPs found in ",
    dgp_root
  )
}

dgp_review <- lapply(
  dgp_dirs,
  summarize_dgp
)


# Combine lag columns -----

all_columns <- unique(
  unlist(
    lapply(
      dgp_review,
      names
    )
  )
)

dgp_review <- lapply(
  dgp_review,
  function(row) {
    
    missing_columns <- setdiff(
      all_columns,
      names(
        row
      )
    )
    
    for (column in missing_columns) {
      row[[column]] <- NA
    }
    
    row[
      all_columns
    ]
  }
)

dgp_review <- do.call(
  rbind,
  dgp_review
)

rownames(
  dgp_review
) <- NULL


# Ordering -----

dgp_review <- dgp_review[
  order(
    dgp_review$family,
    dgp_review$n_units,
    dgp_review$m,
    dgp_review$p_lags,
    dgp_review$dgp_id
  ),
  ,
  drop = FALSE
]


# Save -----

utils::write.csv(
  dgp_review,
  output_path,
  row.names = FALSE
)


# Console report -----

separator <- paste0(
  rep(
    "=",
    76L
  ),
  collapse = ""
)

cat(
  "\n",
  separator,
  "\nDGP BANK REVIEW\n",
  separator,
  "\n\n",
  sep = ""
)

cat(
  sprintf(
    "DGPs reviewed ......................... %d\n",
    nrow(
      dgp_review
    )
  ),
  sprintf(
    "Families .............................. %d\n",
    length(
      unique(
        dgp_review$family
      )
    )
  ),
  sprintf(
    "n_units range ......................... %d -- %d\n",
    min(
      dgp_review$n_units
    ),
    max(
      dgp_review$n_units
    )
  ),
  sprintf(
    "m range ............................... %d -- %d\n",
    min(
      dgp_review$m
    ),
    max(
      dgp_review$m
    )
  ),
  sprintf(
    "p range ............................... %d -- %d\n",
    min(
      dgp_review$p_lags
    ),
    max(
      dgp_review$p_lags
    )
  ),
  sprintf(
    "N range ............................... %d -- %d\n",
    min(
      dgp_review$N
    ),
    max(
      dgp_review$N
    )
  ),
  sprintf(
    "Spectral radius range ................. %.6f -- %.6f\n",
    min(
      dgp_review$spectral_radius
    ),
    max(
      dgp_review$spectral_radius
    )
  ),
  "\n",
  sep = ""
)


# Family summary -----

cat(
  "Families\n"
)

family_table <- table(
  dgp_review$family
)

for (family in names(
  family_table
)) {
  
  cat(
    sprintf(
      "  %-24s %d\n",
      family,
      family_table[[family]]
    )
  )
}


# Potential issues -----

dimension_failures <- dgp_review[
  !dgp_review$dimensions_match,
  ,
  drop = FALSE
]

unstable <- dgp_review[
  !dgp_review$stable,
  ,
  drop = FALSE
]

near_boundary <- dgp_review[
  dgp_review$spectral_radius >= 0.98,
  ,
  drop = FALSE
]

very_low_radius <- dgp_review[
  dgp_review$spectral_radius < 0.50,
  ,
  drop = FALSE
]

cat(
  "\nPotential issues\n",
  sprintf(
    "  Dimension mismatches ................ %d\n",
    nrow(
      dimension_failures
    )
  ),
  sprintf(
    "  Unstable DGPs ....................... %d\n",
    nrow(
      unstable
    )
  ),
  sprintf(
    "  Radius >= 0.98 ...................... %d\n",
    nrow(
      near_boundary
    )
  ),
  sprintf(
    "  Radius < 0.50 ....................... %d\n",
    nrow(
      very_low_radius
    )
  ),
  "\n",
  sep = ""
)

if (nrow(
  dimension_failures
) > 0L) {
  
  cat(
    "Dimension mismatches:\n"
  )
  
  print(
    dimension_failures[
      ,
      c(
        "dgp_id",
        "n_units",
        "m",
        "p_lags",
        "N"
      ),
      drop = FALSE
    ],
    row.names = FALSE
  )
  
  cat(
    "\n"
  )
}

if (nrow(
  unstable
) > 0L) {
  
  cat(
    "Unstable DGPs:\n"
  )
  
  print(
    unstable[
      ,
      c(
        "dgp_id",
        "spectral_radius"
      ),
      drop = FALSE
    ],
    row.names = FALSE
  )
  
  cat(
    "\n"
  )
}

if (nrow(
  near_boundary
) > 0L) {
  
  cat(
    "Near-boundary DGPs:\n"
  )
  
  print(
    near_boundary[
      ,
      c(
        "dgp_id",
        "spectral_radius"
      ),
      drop = FALSE
    ],
    row.names = FALSE
  )
  
  cat(
    "\n"
  )
}


# Compact inventory -----

cat(
  "Compact inventory\n\n"
)

inventory <- dgp_review[
  ,
  c(
    "dgp_id",
    "family",
    "n_units",
    "m",
    "p_lags",
    "N",
    "spectral_radius",
    "active_lag_blocks",
    "active_foreign_lag_blocks",
    "lag_density",
    "mean_active_strength",
    "mean_foreign_strength",
    "zero_foreign_degree_units",
    "weak_component_count",
    "largest_weak_component_size",
    "largest_weak_component_share"
  ),
  drop = FALSE
]

inventory$spectral_radius <- round(
  inventory$spectral_radius,
  4L
)

inventory$lag_density <- round(
  inventory$lag_density,
  4L
)

inventory$mean_active_strength <- round(
  inventory$mean_active_strength,
  5L
)

inventory$mean_foreign_strength <- round(
  inventory$mean_foreign_strength,
  5L
)

inventory$largest_weak_component_share <- round(
  inventory$largest_weak_component_share,
  4L
)

print(
  inventory,
  row.names = FALSE
)

cat(
  "\n",
  separator,
  "\nDGP BANK REVIEW COMPLETE\n",
  separator,
  "\n",
  sprintf(
    "CSV: %s\n",
    output_path
  ),
  sep = ""
)