# This script prepares the data for the exploratory analysis.
# It does four things:
# 1. Find the reference draws file path for each posterior.
# 2. Compute geometry features from the reference draws.
# 3. Read sampler runs results information.
# 4. Compute RMSE, ESS/GE, and posterior-level performance.





# ============================================================
# Load in packages and data directory
# ============================================================

# Load packages
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(tibble)
library(posteriordb)
library(PDMPSamplersR)
library(posterior)


# Directory with posteriorDB unconstrained reference draws
pdb_ref_dir <- "/Users/zhangsr/Desktop/UvA/Thesis/programming/Pipeline/unconstrained reference draws"

# Directory with reference draws that I sampled
own_ref_dir <- "/Users/zhangsr/Desktop/UvA/Thesis/programming/Pipeline/try_to_find_reference"

# Directory with sampler fitting results
results_dir <- "/Users/zhangsr/Desktop/UvA/Thesis/programming/Pipeline/benchmark_runs_final"

# Directory for saving output tables
output_dir <- "/Users/zhangsr/Desktop/UvA/Thesis/programming/Pipeline"





# ============================================================
# Define own reference posteriors
# ============================================================

# These are the posteriors for which my own MCMC reference draws
# passed the reference-draw quality check. 
# (i.e., requirement from PosterioDB paper):
# 1. 10 000 draws per parameter in the model (or more)
# 2. all parameters: abs(mean autocorrelation at lag 1) < 0.05
# 3. all parameters: R hat < 1.01
# 4. all expected fraction of missing information (E-FMI) > 0.2 
# (the paper said below 0.2, but the reference posterior from 
# posteriorDB also has E-FMI > 0.2)
# 5. no divergent transitions
posterior_groups_ref <- list(
  simple = c(
    "Mth_data-Mth_model",
    "bones_data-bones_model",
    "dogs-dogs_hierarchical",
    "dogs-dogs",
    "dugongs_data-dugongs_model",
    "GLM_Binomial_data-GLM_Binomial_model",
    "GLM_Poisson_Data-GLM_Poisson_model",
    "GLMM_data-GLMM1_model",
    "M0_data-M0_model",
    "Mb_data-Mb_model",
    "Mt_data-Mt_model",
    "Mtbh_data-Mtbh_model",
    "nes_logit_data-nes_logit_model",
    "radon_all-radon_pooled",
    "radon_mn-radon_pooled",
    "Rate_1_data-Rate_1_model",
    "Rate_2_data-Rate_2_model",
    "Rate_3_data-Rate_3_model",
    "Rate_4_data-Rate_4_model",
    "Rate_5_data-Rate_5_model",
    "sesame_data-sesame_one_pred_a",
    "wells_data-wells_daae_c_model",
    "wells_data-wells_dae_c_model",
    "wells_data-wells_dae_inter_model",
    "wells_data-wells_dae_model",
    "wells_data-wells_dist100_model",
    "wells_data-wells_dist100ars_model",
    "wells_data-wells_interaction_c_model",
    "wells_data-wells_interaction_model"
  )
)

# Flatten the list into one vector of posterior names
own_reference_posteriors <- unlist(posterior_groups_ref)





# ============================================================
# Build reference table
# ============================================================

# Read posteriorDB reference files.
pdb_ref_files <- tibble(
  ref_file_path = list.files(
    pdb_ref_dir,
    pattern = "_unconstrained_ref\\.rds$",
    full.names = TRUE
  )
) %>%
  mutate(
    file_name = basename(ref_file_path),
    posterior = str_remove(file_name, "_unconstrained_ref\\.rds$"),
    reference = 1L,
    reference_label = "posteriordb"
  ) %>%
  dplyr::select(
    posterior,
    reference,
    reference_label,
    ref_file_path
  )


# Read my own reference files.
# Only keep the posteriors that passed the reference-draw quality check.
own_ref_files <- tibble(
  ref_file_path = list.files(
    own_ref_dir,
    pattern = "_MCMC_seed4711_unconstrained_matrix\\.rds$",
    full.names = TRUE
  )
) %>%
  mutate(
    file_name = basename(ref_file_path),
    posterior = str_remove(file_name, "_MCMC_seed4711_unconstrained_matrix\\.rds$"),
    reference = 2L,
    reference_label = "own_reference"
  ) %>%
  filter(
    posterior %in% own_reference_posteriors
  ) %>%
  dplyr::select(
    posterior,
    reference,
    reference_label,
    ref_file_path
  )


# Check whether all expected own reference files exist.
own_reference_check <- tibble(
  posterior = own_reference_posteriors
) %>%
  mutate(
    expected_file_path = file.path(
      own_ref_dir,
      paste0(posterior, "_MCMC_seed4711_unconstrained_matrix.rds")
    ),
    file_exists = file.exists(expected_file_path)
  )


# Combine posteriorDB references and own references.
reference_table <- bind_rows(
  pdb_ref_files,
  own_ref_files
) %>%
  mutate(
    priority = case_when(
      reference == 2L ~ 1L,
      reference == 1L ~ 2L,
      TRUE ~ 99L
    )
  ) %>%
  arrange(posterior, priority) %>%
  group_by(posterior) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    file_exists = file.exists(ref_file_path)
  ) %>%
  dplyr::select(
    posterior,
    reference,
    reference_label,
    ref_file_path,
    file_exists
  ) %>%
  arrange(reference, posterior)





# ============================================================
# Functions for geometry features
# ============================================================

# Compute sample skewness for one parameter.
# This is used to measure marginal asymmetry.
compute_skewness <- function(x) {
  
  # Keep only finite values
  x <- x[is.finite(x)]
  n <- length(x)
  
  # Skewness needs at least 3 values
  if (n < 3) return(NA_real_)
  
  m <- mean(x)
  s <- sd(x)
  
  # If the parameter is constant, return 0 skewness
  if (is.na(s) || s == 0) return(0)
  
  mean(((x - m) / s)^3)
}


# Compute condition number and dimension from posterior draws.
#
# The condition number is computed from the covariance matrix of the 
# unconstrained reference draws, which is defined as:
# largest covariance eigenvalue / smallest covariance eigenvalue
compute_cond_dim <- function(draws) {
  
  # Convert posterior draws to matrix format
  draws_mat <- as_draws_matrix(draws)
  draws_mat <- as.matrix(draws_mat)
  
  # Remove log-density columns if they exist
  drop_cols <- grepl("^lp__", colnames(draws_mat)) |
    grepl("^lp_approx__", colnames(draws_mat))
  
  draws_mat <- draws_mat[, !drop_cols, drop = FALSE]
  
  # Remove constant parameters
  keep <- apply(draws_mat, 2, sd, na.rm = TRUE) > 0
  draws_mat <- draws_mat[, keep, drop = FALSE]
  
  # Number of unconstrained parameters
  dim_unconstrained <- ncol(draws_mat)
  
  # If no parameter remains, return NA
  if (dim_unconstrained < 1) {
    return(list(cond = NA_real_, dim = 0))
  }
  
  # Compute covariance matrix
  cov_mat <- tryCatch(
    cov(draws_mat),
    error = function(e) NULL
  )
  
  # Return NA if covariance matrix is invalid
  if (is.null(cov_mat) || any(!is.finite(cov_mat))) {
    return(list(cond = NA_real_, dim = dim_unconstrained))
  }
  
  # Compute covariance eigenvalues
  eigvals <- tryCatch(
    eigen(cov_mat, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NULL
  )
  
  if (is.null(eigvals)) {
    return(list(cond = NA_real_, dim = dim_unconstrained))
  }
  
  # Remove non-finite and near-zero eigenvalues
  eigvals <- eigvals[is.finite(eigvals) & eigvals > 1e-12]
  
  # Compute condition number
  if (length(eigvals) == 1) {
    cond <- 1
  } else if (length(eigvals) < 1) {
    cond <- NA_real_
  } else {
    cond <- max(eigvals) / min(eigvals)
  }
  
  list(cond = cond, dim = dim_unconstrained)
}





# ============================================================
# Compute geometry features from one reference file
# ============================================================

# This function computes geometry features from the unconstrained reference draws.
compute_geometry_features <- function(ref_file_path) {
  
  # Read reference draws
  X <- readRDS(ref_file_path)
  X <- as.matrix(X)
  storage.mode(X) <- "numeric"
  
  n_draws <- nrow(X)
  d <- ncol(X)
  
  # Return NA if the reference draw matrix is too small
  if (n_draws < 10 || d < 1) {
    return(tibble::tibble(
      n_ref_draws = n_draws,
      dim = d,
      mean_abs_skewness = NA_real_,
      tail_ratio = NA_real_,
      mean_abs_correlation = NA_real_,
      condition_number = NA_real_,
      log_condition_number = NA_real_
    ))
  }
  
  # 1: mean absolute skewness
  skew_values <- apply(X, 2, compute_skewness)
  mean_abs_skewness <- mean(abs(skew_values), na.rm = TRUE)
  
  # 2: tail ratio
  # Mean ratio of the 1%-99% quantile range to the IQR.
  # For each parameter: tail_ratio = (q99 - q01) / IQR
  # Then average across parameters.
  q99 <- apply(X, 2, quantile, probs = 0.99, na.rm = TRUE)
  q01 <- apply(X, 2, quantile, probs = 0.01, na.rm = TRUE)
  q75 <- apply(X, 2, quantile, probs = 0.75, na.rm = TRUE)
  q25 <- apply(X, 2, quantile, probs = 0.25, na.rm = TRUE)
  
  iqr <- q75 - q25
  tail_ratio_each <- (q99 - q01) / (iqr + 1e-12)
  tail_ratio <- mean(tail_ratio_each, na.rm = TRUE)
  
  # 3: mean absolute correlation
  # This is the mean of the absolute off-diagonal entries
  # of the posterior correlation matrix.
  # Larger values mean stronger average posterior dependence.
  if (d > 1) {
    R <- suppressWarnings(stats::cor(X))
    off_diag <- R[upper.tri(R)]
    mean_abs_correlation <- mean(abs(off_diag), na.rm = TRUE)
  } else {
    mean_abs_correlation <- 0
  }
  
  # 4: condition number and dimension
  cond_res <- compute_cond_dim(X)
  
  condition_number <- cond_res$cond
  dim_unconstrained <- cond_res$dim
  
  # Use log scale for condition number
  # Because condition numbers can be very large.
  log_condition_number <- ifelse(
    is.na(condition_number) || condition_number <= 0,
    NA_real_,
    log(condition_number)
  )
  
  tibble::tibble(
    n_ref_draws = n_draws,
    dim = dim_unconstrained,
    mean_abs_skewness = mean_abs_skewness,
    tail_ratio = tail_ratio,
    mean_abs_correlation = mean_abs_correlation,
    condition_number = condition_number,
    log_condition_number = log_condition_number
  )
}





# ============================================================
# Build geometry table for all reference posteriors
# ============================================================

# Compute geometry features for every posterior with a reference file.
geometry_table <- reference_table %>%
  dplyr::filter(file_exists) %>%
  dplyr::mutate(
    geometry = purrr::map(
      ref_file_path,
      ~ tryCatch(
        compute_geometry_features(.x),
        error = function(e) {
          tibble::tibble(
            n_ref_draws = NA_integer_,
            dim = NA_integer_,
            mean_abs_skewness = NA_real_,
            tail_ratio = NA_real_,
            mean_abs_correlation = NA_real_,
            condition_number = NA_real_,
            log_condition_number = NA_real_,
            geometry_error = conditionMessage(e)
          )
        }
      )
    )
  ) %>%
  tidyr::unnest(geometry)


# Save geometry table
saveRDS(
  geometry_table,
  file.path(output_dir, "geometry_reference_table.rds")
)





# ============================================================
# Read sampler run information
# ============================================================

# This function extracts metadata from a result file name.
get_runs_info <- function(file, results_dir) {
  
  x <- basename(file)
  
  # Algorithms used in the benchmark
  methods <- c("MCMC", "ZigZag", "BPS", "Boomerang")
  method_pattern <- paste(methods, collapse = "|")
  
  # Extract algorithm
  algorithm <- str_match(
    x,
    paste0("_(", method_pattern, ")_seed[0-9]+\\.rds$")
  )[, 2]
  
  # Extract seed
  seed <- as.integer(
    str_match(x, "_seed([0-9]+)\\.rds$")[, 2]
  )
  
  # Remove algorithm and seed from file name
  prefix <- str_remove(
    x,
    paste0("_(", method_pattern, ")_seed[0-9]+\\.rds$")
  )
  
  # Extract group
  group <- str_extract(prefix, "^[^_]+")
  
  # Extract posterior name
  posterior <- prefix %>%
    str_remove("^(simple|moderate|complex)_")
  
  # For MCMC, store the path to the unconstrained draw matrix.
  # This file is needed to compute RMSE and ESS for MCMC.
  mcmc_draws_file <- NA_character_
  
  if (!is.na(algorithm) && algorithm == "MCMC") {
    draws_candidate <- file.path(
      results_dir,
      paste0(
        group,
        "_",
        posterior,
        "_MCMC_seed",
        seed,
        "_unconstrained_matrix.rds"
      )
    )
    
    if (file.exists(draws_candidate)) {
      mcmc_draws_file <- draws_candidate
    }
  }
  
  tibble(
    posterior = posterior,
    group = group,
    algorithm = algorithm,
    seed = seed,
    fit_file = file,
    mcmc_draws_file = mcmc_draws_file
  )
}


# List all RDS files in the result directory
all_rds_files <- list.files(
  results_dir,
  pattern = "\\.rds$",
  full.names = TRUE
)

# Keep only sampler result files.
# Remove x0 files and MCMC unconstrained draw matrices.
runs <- all_rds_files[
  !grepl("_x0\\.rds$", all_rds_files) &
    !grepl("_unconstrained_matrix\\.rds$", all_rds_files)
]

# Extract metadata for all runs
runs_info <- bind_rows(
  lapply(runs, get_runs_info, results_dir = results_dir)
)





# ============================================================
# Filter runs for the exploratory analysis
# ============================================================

# This table contains the number of successful seeds for each posterior.
posterior_method_seed_counts <- readRDS(
  file.path(output_dir, "posterior_method_seed_counts.rds")
)


# Keep only: posteriors with reference draws 
# and posteriors with more than 5 seeds across all methods.
runs_info_filtered <- runs_info %>%
  dplyr::semi_join(
    reference_table %>%
      dplyr::distinct(posterior),
    by = "posterior"
  ) %>%
  dplyr::anti_join(
    posterior_method_seed_counts %>%
      dplyr::filter(min_seeds_across_methods <= 5) %>%
      dplyr::distinct(posterior),
    by = "posterior"
  )


# Check how many runs and posteriors remain
tibble::tibble(
  step = c("original runs_info", "after filtering"),
  n_rows = c(nrow(runs_info), nrow(runs_info_filtered)),
  n_posteriors = c(
    dplyr::n_distinct(runs_info$posterior),
    dplyr::n_distinct(runs_info_filtered$posterior)
  )
)





# ============================================================
# Functions for performance measures
# ============================================================

# Compute global RMSE of posterior mean estimate.
#
# For MCMC, the mean estimate is computed from the MCMC unconstrained draws.
# For PDMP, the mean estimate is computed from the PDMP fit object.
# The reference mean is computed from the reference draws.
# The global RMSE is: sqrt(mean((estimated_mean - reference_mean)^2))
compute_pee <- function(draws = NULL,
                        ref_draws,
                        algorithm = c("MCMC", "ZigZag", "BPS", "Boomerang"),
                        pdmp_fit = NULL,
                        type = c("global", "per_param", "both")) {
  
  type <- match.arg(type)
  algorithm <- match.arg(algorithm)
  
  # Convert reference draws to matrix
  if (inherits(ref_draws, "draws_matrix") ||
      inherits(ref_draws, "draws") ||
      inherits(ref_draws, "array") ||
      is.list(ref_draws)) {
    ref_draws <- as.matrix(ref_draws)
  }
  
  # Reference posterior mean
  ref_mean <- colMeans(ref_draws)
  
  # Estimated posterior mean
  if (algorithm == "MCMC") {
    
    if (is.null(draws)) {
      stop("For MCMC, 'draws' must be provided.")
    }
    
    # Convert MCMC draws to matrix
    if (inherits(draws, "draws_matrix") ||
        inherits(draws, "draws") ||
        inherits(draws, "array") ||
        is.list(draws)) {
      draws <- as.matrix(draws)
    }
    
    if (ncol(draws) != ncol(ref_draws)) {
      stop("draws and ref_draws must have the same number of parameters")
    }
    
    draws_mean <- colMeans(draws)
    
  } else {
    
    if (is.null(pdmp_fit)) {
      stop("For PDMP algorithms, 'pdmp_fit' must be provided.")
    }
    
    # Get posterior mean estimate from the PDMP result
    draws_mean <- PDMPSamplersR:::mean.pdmp_result(pdmp_fit)
    draws_mean <- as.numeric(draws_mean)
    
    if (length(draws_mean) != length(ref_mean)) {
      stop("PDMP mean and ref_draws must have the same number of parameters")
    }
    
    names(draws_mean) <- colnames(ref_draws)
  }
  
  # Squared error for each parameter
  per_param_pee <- (draws_mean - ref_mean)^2
  
  if (type == "per_param") {
    return(per_param_pee)
  } else if (type == "global") {
    global_pee <- sqrt(mean(per_param_pee))
    return(global_pee)
  } else if (type == "both") {
    global_pee <- sqrt(mean(per_param_pee))
    return(list(global = global_pee, per_param = per_param_pee))
  }
}


# Compute ESS for one sampler run.
#
# For MCMC, ESS is computed from the unconstrained MCMC draw matrix.
# For PDMP, ESS is computed from the PDMP fit object.
# The returned summary uses the minimum ESS across parameters.
compute_ess <- function(fit_file, algorithm, mcmc_draws_file = NULL) {
  
  if (algorithm == "MCMC") {
    
    if (is.null(mcmc_draws_file) || !file.exists(mcmc_draws_file)) {
      return(NA_real_)
    }
    
    draws <- readRDS(mcmc_draws_file)
    ess_per_param <- apply(as.matrix(draws), 2, posterior::ess_bulk)
    
    return(list(
      per_param = ess_per_param,
      min = min(ess_per_param, na.rm = TRUE)
    ))
    
  } else {
    
    if (!file.exists(fit_file)) {
      return(NA_real_)
    }
    
    fit <- readRDS(fit_file)
    ess_per_param <- ess(fit)
    
    return(list(
      per_param = ess_per_param,
      min = min(ess_per_param, na.rm = TRUE)
    ))
  }
}


# Compute gradient evaluations for one sampler run.
#
# For MCMC, this is the total number of leapfrog steps after warmup.
# For PDMP, this is the number of gradient calls stored in the fit object.
compute_ge <- function(fit_file, algorithm) {
  
  if (!file.exists(fit_file)) return(NA_real_)
  
  fit <- readRDS(fit_file)
  
  if (algorithm == "MCMC") {
    
    diag <- fit$sampler_diagnostics(inc_warmup = FALSE)
    nlf <- diag[, , "n_leapfrog__"]
    ge <- sum(nlf, na.rm = TRUE)
    
  } else {
    
    ge <- fit$stats$gradient_calls
  }
  
  ge
}


# Add reference file information to each run
runs_info_filtered_with_ref <- runs_info_filtered %>%
  dplyr::left_join(
    geometry_table %>%
      dplyr::select(
        posterior,
        reference,
        reference_label,
        ref_file_path
      ) %>%
      dplyr::distinct(),
    by = "posterior"
  )





# ============================================================
# Compute run-level performance
# ============================================================

# Compute RMSE, minimum ESS, and gradient evaluations for each run.
# (It will take about 4 mins)
runs_performance <- runs_info_filtered_with_ref %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    
    # Global RMSE of posterior mean estimate
    rmse_global = tryCatch({
      
      ref_draws <- readRDS(ref_file_path)
      
      if (algorithm == "MCMC") {
        
        draws <- readRDS(mcmc_draws_file)
        
        compute_pee(
          draws = draws,
          ref_draws = ref_draws,
          algorithm = algorithm,
          type = "global"
        )
        
      } else {
        
        pdmp_fit <- readRDS(fit_file)
        
        compute_pee(
          ref_draws = ref_draws,
          algorithm = algorithm,
          pdmp_fit = pdmp_fit,
          type = "global"
        )
      }
      
    }, error = function(e) NA_real_),
    
    
    # Minimum ESS across parameters
    ess_min = tryCatch({
      
      ess_res <- compute_ess(
        fit_file = fit_file,
        algorithm = algorithm,
        mcmc_draws_file = mcmc_draws_file
      )
      
      if (is.list(ess_res)) {
        ess_res$min
      } else {
        NA_real_
      }
      
    }, error = function(e) NA_real_),
    
    
    # Number of gradient evaluations
    gradient_evaluations = tryCatch({
      
      compute_ge(
        fit_file = fit_file,
        algorithm = algorithm
      )
      
    }, error = function(e) NA_real_)
  ) %>%
  dplyr::ungroup()


# Compute ESS per ge.
# Larger values mean better sampling efficiency.
runs_performance <- runs_performance %>%
  dplyr::mutate(
    ess_per_ge = ess_min / gradient_evaluations
  )


# Save run-level performance table
saveRDS(
  runs_performance,
  file.path(output_dir, "runs_performance_with_rmse_ess.rds")
)





# ============================================================
# Aggregate performance across seeds
# ============================================================

# Aggregate each posterior and algorithm across seeds.
# The median is used because it is more robust to bad seeds or unstable runs than the mean.
posterior_performance <- runs_performance %>%
  dplyr::group_by(group, posterior, algorithm) %>%
  dplyr::summarise(
    
    # Number of attempted runs
    n_attempted = dplyr::n(),
    
    # Number of successful RMSE values
    n_success_rmse = sum(!is.na(rmse_global)),
    
    # Number of successful ESS/GE values
    n_success_ess_ge = sum(!is.na(ess_per_ge) & is.finite(ess_per_ge)),
    
    # Median performance across successful seeds
    median_rmse = median(rmse_global, na.rm = TRUE),
    median_ess_per_ge = median(ess_per_ge, na.rm = TRUE),
    
    # Bad performance outcomes used in the random forest analysis
    # Larger bad_rmse means larger posterior mean error.
    # Larger bad_efficiency means lower ESS per gradient evaluation.
    bad_rmse = log(median_rmse),
    bad_efficiency = -log(median_ess_per_ge),
    
    .groups = "drop"
  )


# Add reference-based geometry features.
posterior_performance <- posterior_performance %>%
  dplyr::left_join(
    geometry_table %>%
      dplyr::select(
        posterior,
        reference,
        reference_label,
        n_ref_draws,
        dim,
        mean_abs_skewness,
        tail_ratio,
        mean_abs_correlation,
        condition_number,
        log_condition_number
      ),
    by = "posterior"
  )


# Save final table for the random forest analysis
saveRDS(
  posterior_performance,
  file.path(output_dir, "posterior_algorithm_performance_with_geometry.rds")
)



