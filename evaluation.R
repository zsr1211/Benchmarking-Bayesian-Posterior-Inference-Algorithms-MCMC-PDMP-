# This script is about defining functions to compute evaluation measures, 
# and biuld the pipeline to get the evaluation measures.


# load in packages
library(posteriordb)
library(posterior)
library(stringr)
library(purrr)
library(dplyr)
library(tibble)
library(cmdstanr)
library(Matrix)
library(PDMPSamplersR)
library(transport)
library(ggplot2)
library(kldest)
library(RANN)


##########################################################################################
# The evaluation functions about accuracy-------------------------------------------------
##########################################################################################

# get fitting result file info
# load in posteriordb and the generated fitting result files
pdb <- pdb_local("/Users/zhangsr/Desktop/UvA/Thesis/posteriorDB/posteriordb")
results_dir <- "/Users/zhangsr/Desktop/UvA/Thesis/programming/Pipeline/benchmark_runs_final"
getwd()
setwd(results_dir)
getwd()



# for (f in files) {
#   algo <- str_extract(basename(f), "(MCMC|ZigZag|BPS|Boomerang)")
#   if (is.na(algo)) {
#     message("File has unknown algorithm: ", f)
#   }
# }
# function to extract file information
get_file_info <- function(file, results_dir, ref_dir) {
  
  # the file name
  x <- basename(file)
  
  # posterior name
  posterior <- x %>%
    str_remove("^(simple|moderate|complex)_") %>%
    str_remove("_(MCMC|ZigZag|BPS|Boomerang)_seed\\d+\\.rds$")
  
  # the posterior comlexity group
  group <- str_extract(x, "^[^_]+")
  
  # algorithm
  algorithm <- str_extract(x, "(MCMC|ZigZag|BPS|Boomerang)")
  
  # seed
  seed <- as.numeric(str_extract(x, "(?<=seed)\\d+"))
  
  # unconstrained reference draws file path
  ref_file_candidate <- file.path(ref_dir, paste0(posterior, "_unconstrained_ref.rds"))
  # if no reference file, marked as NA
  ref_file <- if (file.exists(ref_file_candidate)) ref_file_candidate else NA_character_
  
  # MCMC unconstrained draws file path
  mcmc_draws_file <- NA_character_
  if (algorithm == "MCMC") {
    # get the path of MCMC draws file
    draws_candidate <- file.path(results_dir, paste0(group, "_", posterior, "_MCMC_seed", seed, "_unconstrained_matrix.rds"))
    if (file.exists(draws_candidate)) {
      mcmc_draws_file <- draws_candidate
    }
  }
  
  # combine all the information together
  tibble(
    posterior = posterior,
    group = group,
    algorithm = algorithm,
    seed = seed,
    fit_file = file,               
    ref_file = ref_file,       
    has_ref = !is.na(ref_file),
    mcmc_draws_file = mcmc_draws_file
  )
}



# take one posterior as an example to show the following functions about computing measures
draws <- readRDS("/Users/zhangsr/Desktop/UvA/Thesis/programming/Pipeline/benchmark_runs_final/simple_eight_schools-eight_schools_noncentered_MCMC_seed8016_unconstrained_matrix.rds")
ref_draws <- readRDS("/Users/zhangsr/Desktop/UvA/Thesis/programming/Pipeline/unconstrained reference draws/eight_schools-eight_schools_noncentered_unconstrained_ref.rds")
file.exists("/Users/zhangsr/Desktop/UvA/Thesis/programming/Pipeline/benchmark_runs_final/simple_eight_schools-eight_schools_noncentered_MCMC_seed8016_unconstrained_matrix.rds")
draws_pdmp <- readRDS("/Users/zhangsr/Desktop/UvA/Thesis/programming/Pipeline/benchmark_runs_final/simple_eight_schools-eight_schools_noncentered_BPS_seed8016.rds")


# RMSE
# compute_pee <- function(draws, ref_draws, type = c("global", "per_param", "both")) {
#   type <- match.arg(type)
#   
#   # transform to matrix
#   if (inherits(draws, "draws_matrix") || inherits(draws, "draws") || inherits(draws, "array")) draws <- as.matrix(draws)
#   if (inherits(ref_draws, "draws_matrix") || inherits(ref_draws, "draws") || inherits(ref_draws, "array")) ref_draws <- as.matrix(ref_draws)
#   
#   # check if column numbers of draws and reference draws are the same 
#   if (ncol(draws) != ncol(ref_draws)) stop("draws and ref_draws must have the same number of parameters")
#   
#   # compute mean of each parameters
#   ref_mean <- colMeans(ref_draws)
#   draws_mean <- colMeans(draws)
#   
#   # get PEE for each parameters
#   per_param_pee <- (draws_mean - ref_mean)^2
#   
#   # return based on argument "type" 
#   if (type == "per_param") {
#     return(per_param_pee)
#   } else if (type == "global") {
#     global_pee <- sqrt(mean(per_param_pee))
#     return(global_pee)
#   } else if (type == "both") {
#     global_pee <- sqrt(mean(per_param_pee))
#     return(list(global = global_pee, per_param = per_param_pee))
#   }
# }


compute_pee <- function(draws = NULL,
                        ref_draws,
                        algorithm = c("MCMC", "ZigZag", "BPS", "Boomerang"),
                        pdmp_fit = NULL,
                        type = c("global", "per_param", "both")) {
  type <- match.arg(type)
  algorithm <- match.arg(algorithm)
  
  # transform reference draws to matrix
  if (inherits(ref_draws, "draws_matrix") ||
      inherits(ref_draws, "draws") ||
      inherits(ref_draws, "array") ||
      is.list(ref_draws)) {
    ref_draws <- as.matrix(ref_draws)
  }
  
  # reference mean
  ref_mean <- colMeans(ref_draws)
  
  # estimated mean from MCMC or PDMP
  if (algorithm == "MCMC") {
    
    if (is.null(draws)) {
      stop("For MCMC, 'draws' must be provided.")
    }
    
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
    
    draws_mean <- PDMPSamplersR:::mean.pdmp_result(pdmp_fit)
    
    # make sure it is numeric vector
    draws_mean <- as.numeric(draws_mean)
    
    if (length(draws_mean) != length(ref_mean)) {
      stop("PDMP mean and ref_draws must have the same number of parameters")
    }
    
    # keep parameter names if possible
    names(draws_mean) <- colnames(ref_draws)
  }
  
  # squared error for each parameter
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

compute_pee(ref_draws = ref_draws, algorithm = "BPS", pdmp_fit = draws_pdmp, type = "per_param")
compute_pee(ref_draws = ref_draws, algorithm = "BPS", pdmp_fit = draws_pdmp)
compute_pee(draws = draws, ref_draws = ref_draws, type = "per_param")
compute_pee(draws = draws, ref_draws = ref_draws)
ref_draws <- as.matrix(ref_draws)
colMeans(ref_draws)

# only return global
compute_pee(draws, ref_draws, type = "global")
# only per_param
compute_pee(draws, ref_draws, type = "per_param")
# both
compute_pee(draws, ref_draws, type = "both")
# will return list: $global and $per_param



# W1 and W2
# compute Wasserstein distance（1-Wasserstein or 2-Wasserstein）
compute_wasserstein <- function(draws, ref_draws, p = 1, type = c("global", "per_param", "both")) {
  type <- match.arg(type)
  
  # transform to matrix
  if (inherits(draws, "draws_matrix") || inherits(draws, "draws") || inherits(draws, "array")) draws <- as.matrix(draws)
  if (inherits(ref_draws, "draws_matrix") || inherits(ref_draws, "draws") || inherits(ref_draws, "array")) ref_draws <- as.matrix(ref_draws)
  
  # check if column numbers of draws and reference draws are the same 
  if (ncol(draws) != ncol(ref_draws)) stop("draws and ref_draws must have the same number of parameters")
  
  n_params <- ncol(draws)
  vals <- numeric(n_params)
  
  # for each parameter, compute the W1 or W2 distance
  for (j in seq_len(n_params)) {
    vals[j] <- wasserstein1d(draws[, j], ref_draws[, j], p = p)
  }
  
  # return based on the argument "type"
  if (type == "per_param") {
    return(vals)
  } else if (type == "global") {
    return(mean(vals))
  } else if (type == "both") {
    return(list(global = mean(vals), per_param = vals))
  }
}

# 1-Wasserstein
draws_pdmp <- discretize(draws_pdmp, dt = 0.1)
compute_wasserstein(draws, ref_draws, p = 1, type = "per_param")
compute_wasserstein(draws, ref_draws, p = 1, type = "global")
compute_wasserstein(draws_pdmp, ref_draws, p = 1, type = "both")
# 2-Wasserstein
compute_wasserstein(draws_pdmp, ref_draws, p = 2, type = "global")



# KL divergence
compute_kl <- function(X,
                       Y,
                       jitter = TRUE,
                       jitter_eps = 1e-8,
                       seed = 1234,
                       vartype = NULL) {
  if (!requireNamespace("kldest", quietly = TRUE)) {
    stop("Package 'kldest' is required.")
  }
  
  # Convert to ordinary matrices
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  # Basic checks
  if (ncol(X) != ncol(Y)) {
    stop("X and Y must have the same number of columns.")
  }
  
  
  if (!all(is.finite(X))) {
    stop("X contains NA, NaN, or Inf.")
  }
  
  if (!all(is.finite(Y))) {
    stop("Y contains NA, NaN, or Inf.")
  }
  
  # If vartype is not provided, treat all columns as continuous
  if (is.null(vartype)) {
    vartype <- rep("c", ncol(X))
  }
  
  if (length(vartype) != ncol(X)) {
    stop("vartype must have length equal to the number of columns.")
  }
  
  # Joint standardization
  Z <- rbind(X, Y)
  mu <- colMeans(Z)
  sdv <- apply(Z, 2, sd)
  
  if (any(sdv == 0)) {
    stop("At least one column has zero standard deviation after combining X and Y.")
  }
  
  Xs <- scale(X, center = mu, scale = sdv)
  Ys <- scale(Y, center = mu, scale = sdv)
  
  # Add tiny jitter to avoid zero nearest-neighbour distances from duplicates
  if (jitter) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) .Random.seed else NULL
    set.seed(seed)
    
    Xs <- Xs + matrix(rnorm(length(Xs), sd = jitter_eps),
                      nrow = nrow(Xs), ncol = ncol(Xs))
    Ys <- Ys + matrix(rnorm(length(Ys), sd = jitter_eps),
                      nrow = nrow(Ys), ncol = ncol(Ys))
    
    if (!is.null(old_seed)) {
      .Random.seed <<- old_seed
    }
  }
  
  # KL(X || Y)
  kl <- kldest::kld_est(
    X = Xs,
    Y = Ys,
    vartype = vartype
  )
  
  return(kl)
}
draws1 <- readRDS("simple_dogs-dogs_MCMC_seed8016_unconstrained_matrix.rds")
draws2 <- readRDS("simple_dogs-dogs_BPS_seed8016.rds")
draws2 <- discretize(draws2, dt = 0.01)
kl_value <- compute_kl(draws1, draws2)
print(kl_value)






# get accuracy measures -------------------------------

# load files



# get directory of the fitting results and reference files
results_dir <- "/Users/zhangsr/Desktop/UvA/Thesis/programming/Pipeline/benchmark_runs_final"
ref_dir <- "/Users/zhangsr/Desktop/UvA/Thesis/programming/Pipeline/unconstrained reference draws"


# for the fitting results folder, we include fitting results and the MCMC unconstrained draws file
# get list of all fitting .rds file
files <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)
# remove auxiliary files:
# x0 initialization files
# MCMC unconstrained draws files
files <- files[
  !grepl("_x0\\.rds$", files) &
    !grepl("_unconstrained_matrix\\.rds$", files)
]
# get metadata about the posterior info and file path
files_info <- bind_rows(
  lapply(files, get_file_info, results_dir = results_dir, ref_dir = ref_dir)
)
files_info
str(files_info)

# select posterior
target_posteriors <- c(
  # complex
  # ref
  # "kidiq-kidscore_momiq",
  # "earnings-logearn_height",
  "eight_schools-eight_schools_centered",
  "kidiq-kidscore_momhsiq",
  "mesquite-mesquite",
  "earnings-logearn_interaction",
  "kidiq-kidscore_interaction",
  
  # unref
  "pilots-pilots",
  "seeds_data-seeds_centered_model",
  # "lsat_data-lsat_model",
  
  
  # moderate
  # ref
  "gp_pois_regr-gp_pois_regr",
  "hmm_example-hmm_example",
  "nes2000-nes",
  # "kidiq-kidscore_momhs",
  "earnings-logearn_logheight_male",
  "mesquite-logmesquite_logvas",
  
  # unref
  # "low_dim_gauss_mix_collapse-low_dim_gauss_mix_collapse",
  "seeds_data-seeds_model",
  "surgical_data-surgical_model",
  
  
  # simple
  # ref
  "arK-arK",
  "bball_drive_event_0-hmm_drive_0",
  "eight_schools-eight_schools_noncentered",
  "garch-garch11",
  "gp_pois_regr-gp_regr",
  # "low_dim_gauss_mix-low_dim_gauss_mix",
  
  # unref
  "dogs-dogs_log",
  "normal_2-normal_mixture"
  # "Mh_data-Mh_model"
)

files_info <- files_info %>%
  dplyr::filter(posterior %in% target_posteriors)


# initialize the measure columns
files_info <- files_info %>%
  mutate(
    RMSE = NA_real_,
    W1 = NA_real_,
    W2 = NA_real_,
    KL = NA_real_
  )


batch_size <- 100
results_dir <- "/Users/zhangsr/Desktop/UvA/Thesis/programming/Pipeline/benchmark_runs_final"

# 新建保存中间结果的文件夹
save_dir <- "/Users/zhangsr/Desktop/UvA/Thesis/programming/Pipeline/measures_batchs"
dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)

pb <- txtProgressBar(min = 0, max = nrow(files_info), style = 3)

# for loops to get accuracy measures
for (i in seq_len(nrow(files_info))) {
  
  setTxtProgressBar(pb, i)
  
  # get the row number
  row <- files_info[i, ]
  
  # get information
  algorithm <- row$algorithm[[1]]
  has_ref <- row$has_ref[[1]]
  posterior_i <- row$posterior[[1]]
  seed_i <- row$seed[[1]]
  
  # get file path
  fit_file <- row$fit_file[[1]]
  mcmc_draws_file <- row$mcmc_draws_file[[1]]
  ref_file <- row$ref_file[[1]]
  
  # get draws or pdmp results
  fit <- NULL
  draws <- NULL
  
  if (algorithm == "MCMC") {
    if (!is.na(mcmc_draws_file)) {
      draws <- readRDS(mcmc_draws_file)
    } else next
  } else {
    if (!is.na(fit_file) && file.exists(fit_file)) {
      fit <- readRDS(fit_file)
      draws <- discretize(fit, dt = 0.01)
    } else next
  }
  
  # get reference draws
  ref_draws <- NULL
  if (has_ref && !is.na(ref_file) && file.exists(ref_file)) {
    ref_draws <- readRDS(ref_file)
    if (is.list(ref_draws)) ref_draws <- posterior::as_draws_matrix(ref_draws)
  }
  
  # accuracy measures
  if (!is.null(draws)) {
    
    # if has reference
    if (has_ref && !is.null(ref_draws)) {
      
      # RMSE
      files_info$RMSE[i] <- compute_pee(
        draws = draws,
        ref_draws = ref_draws,
        algorithm = algorithm,
        pdmp_fit = fit,
        type = "global"
      )
      
      # W1, W2
      files_info$W1[i] <- compute_wasserstein(draws, ref_draws, p = 1, type = "global")
      files_info$W2[i] <- compute_wasserstein(draws, ref_draws, p = 2, type = "global")
      files_info$KL[i] <- NA_real_
      
    } else {
      
      # if no reference
      files_info$RMSE[i] <- NA_real_
      files_info$W1[i] <- NA_real_
      files_info$W2[i] <- NA_real_
      
      # if no reference & is pdmp
      if (algorithm != "MCMC") {
        
        # get the mcmc draws with same posterior and seed
        mcmc_row <- files_info %>%
          filter(
            posterior == posterior_i,
            seed == seed_i,
            algorithm == "MCMC"
          )
        
        if (nrow(mcmc_row) == 1 && !is.na(mcmc_row$mcmc_draws_file[[1]])) {
          mcmc_draws <- readRDS(mcmc_row$mcmc_draws_file[[1]])
          files_info$KL[i] <- compute_kl(X = mcmc_draws, Y = draws)
        } else {
          files_info$KL[i] <- NA_real_
        }
      } else {
        files_info$KL[i] <- NA_real_     # if MCMC, KL is NA
      }
    }
  }
  
  # save intermediate results
  if (i %% batch_size == 0 || i == nrow(files_info)) {
    save_path <- file.path(save_dir, paste0("accuracy_progress_up_to_", i, ".rds"))
    saveRDS(files_info[1:i, ], save_path)
    cat("\nSaved intermediate results up to file ", i, " -> ", save_path, "\n", sep = "")
  }
  
}





###############
# The evaluation functions about efficiency-------------------------------------------------
###############

# ESS function
compute_ess <- function(fit_file, algorithm, mcmc_draws_file = NULL) {
  
  # Compute ESS for MCMC
  if (algorithm == "MCMC") {
    # we get ESS using the unconstrained MCMC draws
    if (is.null(mcmc_draws_file) || !file.exists(mcmc_draws_file)) return(NA_real_)
    draws <- readRDS(mcmc_draws_file)
    ess_per_param <- apply(as.matrix(draws), 2, posterior::ess_bulk)
    return(list(per_param = ess_per_param, 
                min = min(ess_per_param, na.rm = TRUE)))
  } else {
    
    # get ESS of PDMP
    if (!file.exists(fit_file)) return(NA_real_)
    fit <- readRDS(fit_file)
    ess_per_param <- ess(fit)
    return(list(per_param = ess_per_param, 
                min = min(ess_per_param, na.rm = TRUE)))
  }
}



# Gradient evaluation function (GE)
compute_ge <- function(fit_file, algorithm) {
  if (!file.exists(fit_file)) return(NA_real_)
  fit <- readRDS(fit_file)
  
  # for MCMC, get GE by sum "n_leapfrog__"
  if (algorithm == "MCMC") {
    diag <- fit$sampler_diagnostics(inc_warmup = FALSE)
    nlf <- diag[, , "n_leapfrog__"]  # iterations × chains
    ge <- sum(nlf, na.rm = TRUE)
  } else {
    
    # for PDMP
    ge <- fit$stats$gradient_calls
  }
  
  ge
}



# Sampling time function (include warm up time)
compute_time <- function(fit_file, algorithm) {
  if (!file.exists(fit_file)) return(NA_real_)
  fit <- readRDS(fit_file)
  
  # MCMC
  if (algorithm == "MCMC") {
    time_info <- fit$time()
    sampling_time <- time_info$total
  } else {
    
    # PDMP
    sampling_time <- fit$stats$elapsed_time
  }
  
  sampling_time
}



# compute efficiency measures ----------------------------
# initialize the efficiency measure columns
files_info <- files_info %>%
  mutate(
    ESS_min = NA_real_,
    GE = NA_real_,
    time = NA_real_
  )
# library(dplyr)
# files_info <- files_info %>% select(-ESS_median)

# compute efficiency measures
pb <- txtProgressBar(min = 0, max = nrow(files_info), style = 3)

for (i in seq_len(nrow(files_info))) {
  setTxtProgressBar(pb, i)
  row <- files_info[i, ]
  files_info$ESS_min[i] <- compute_ess(row$fit_file[[1]], row$algorithm[[1]], row$mcmc_draws_file[[1]])$min
  files_info$GE[i] <- compute_ge(row$fit_file[[1]], row$algorithm[[1]])
  files_info$time[i] <- compute_time(row$fit_file[[1]], row$algorithm[[1]])
}

# add the ESS/GE, ESS/time to files_info
files_info <- files_info %>%
  mutate(
    ESS_per_GE   = ifelse(!is.na(ESS_min) & !is.na(GE) & GE > 0,
                          ESS_min / GE,
                          NA_real_),
    ESS_per_time = ifelse(!is.na(ESS_min) & !is.na(time) & time > 0,
                          ESS_min / time,
                          NA_real_)
  )
saveRDS(
  files_info,
  file = "/Users/zhangsr/Desktop/UvA/Thesis/programming/Pipeline/files_info1.rds"
)
files_info <- readRDS("/Users/zhangsr/Desktop/UvA/Thesis/programming/Pipeline/files_info.rds")
str(files_info)
library(dplyr)

# 假设你的 tibble 是 files_info

# for each posterior, select the median w1 and w2 of mcmc as baseline
baseline <- files_info %>%
  filter(algorithm == "MCMC") %>%
  group_by(posterior) %>%
  summarize(
    W1_median = median(W1, na.rm = TRUE),
    W2_median = median(W2, na.rm = TRUE)
  )

# combine baseline to original tibble
files_info <- files_info %>%
  left_join(baseline, by = "posterior")

# relative Wasserstein distance
files_info <- files_info %>%
  mutate(
    W1_relative = W1 / W1_median,
    W2_relative = W2 / W2_median
  )
saveRDS(
  files_info,
  file = "/Users/zhangsr/Desktop/UvA/Thesis/programming/Pipeline/files_info2.rds"
)
files_info2 <- readRDS("/Users/zhangsr/Desktop/UvA/Thesis/programming/Pipeline/files_info2.rds")
str(files_info1)

# new relative distance colomns:
# W1_relative & W2_relative
# # conpoute ESS for one posterior
# row <- files_info[9, ]
# compute_ess(row$fit_file[[1]], row$algorithm[[1]], row$mcmc_draws_file[[1]])
# colnames(files_info)



# average measures by the seeds
files_summary <- files_info %>%
  group_by(posterior, algorithm, group) %>%  
  summarise(
    RMSE_median = median(RMSE, na.rm = TRUE),
    W1_relative_median = median(W1_relative, na.rm = TRUE),
    W2_relative_median = median(W2_relative, na.rm = TRUE),
    KL_median = median(KL, na.rm = TRUE),
    ESS_per_GE_median = median(ESS_per_GE, na.rm = TRUE),
    ESS_per_time_median = median(ESS_per_time, na.rm = TRUE),
    .groups = "drop"
  )

# check files_summary
files_summary













# visualization
measures <- c("RMSE_mean", "W1_mean", "W2_mean", "KL_mean", "ESS_per_GE_mean", "ESS_per_time_mean")

for (measure in measures) {
  
  p <- ggplot(files_summary, aes(x = group, y = .data[[measure]], color = algorithm, group = algorithm)) +
    geom_point(position = position_dodge(width = 0.5), size = 3, na.rm = TRUE) +
    geom_line(position = position_dodge(width = 0.5), na.rm = TRUE) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    ) +
    labs(
      x = "Posterior",
      y = measure,
      color = "Algorithm",
      title = paste0("Comparison of algorithms: ", measure)
    )
  
  # save to environment，each plot is named by measure
  assign(paste0("plot_", measure), p)
}


# see one plot
plot_RMSE_mean
plot_W1_mean
plot_W2_mean
plot_KL_mean
plot_ESS_per_GE_mean
plot_ESS_per_time_mean


files_summary <- files_info %>%
  group_by(posterior, algorithm, group) %>%  
  summarise(
    RMSE_median = median(RMSE, na.rm = TRUE),
    W1_relative_median = median(W1_relative, na.rm = TRUE),
    W2_relative_median = median(W2_relative, na.rm = TRUE),
    KL_median = median(KL, na.rm = TRUE),
    ESS_per_GE_median = median(ESS_per_GE, na.rm = TRUE),
    ESS_per_time_median = median(ESS_per_time, na.rm = TRUE),
    .groups = "drop"
  )

# check files_summary
files_summary



library(ggplot2)
library(dplyr)

# posterior_label
plot_df <- files_summary %>%
  mutate(
    posterior_label = interaction(group, posterior, sep = " / "),
    posterior_label = factor(
      posterior_label,
      levels = unique(posterior_label[order(group, posterior)])
    )
  )

# W1_relative_median
plot_W1 <- plot_df %>%
  filter(!is.na(W1_relative_median)) %>%
  ggplot(aes(x = posterior_label, y = W1_relative_median, color = algorithm)) +
  geom_point(position = position_dodge(width = 0.6), size = 2.5) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.6) +
  labs(x = "Posterior (group / name)", y = "Relative W1 distance", color = "Algorithm") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8))

# W2_relative_median
plot_W2 <- plot_df %>%
  filter(!is.na(W2_relative_median)) %>%
  ggplot(aes(x = posterior_label, y = W2_relative_median, color = algorithm)) +
  geom_point(position = position_dodge(width = 0.6), size = 2.5) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.6) +
  labs(x = "Posterior (group / name)", y = "Relative W2 distance", color = "Algorithm") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8))

# RMSE_median
plot_RMSE <- plot_df %>%
  filter(!is.na(RMSE_median)) %>%
  ggplot(aes(x = posterior_label, y = RMSE_median, color = algorithm)) +
  geom_point(position = position_dodge(width = 0.6), size = 2.5) +
  labs(x = "Posterior (group / name)", y = "Median RMSE", color = "Algorithm") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8))

# ESS_per_GE_median
plot_ESS_GE <- plot_df %>%
  filter(!is.na(ESS_per_GE_median)) %>%
  ggplot(aes(x = posterior_label, y = ESS_per_GE_median, color = algorithm)) +
  geom_point(position = position_dodge(width = 0.6), size = 2.5) +
  labs(x = "Posterior (group / name)", y = "Median ESS per GE", color = "Algorithm") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8))

# ESS_per_time_median
plot_ESS_time <- plot_df %>%
  filter(!is.na(ESS_per_time_median)) %>%
  ggplot(aes(x = posterior_label, y = ESS_per_time_median, color = algorithm)) +
  geom_point(position = position_dodge(width = 0.6), size = 2.5) +
  labs(x = "Posterior (group / name)", y = "Median ESS per time", color = "Algorithm") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8))

# KL_median
plot_KL <- plot_df %>%
  filter(!is.na(KL_median)) %>%
  ggplot(aes(x = posterior_label, y = KL_median, color = algorithm)) +
  geom_point(position = position_dodge(width = 0.6), size = 2.5) +
  labs(x = "Posterior (group / name)", y = "Median KL", color = "Algorithm") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8))



plot_W1
plot_W2
plot_RMSE
plot_KL
plot_ESS_GE
plot_ESS_time



ggplot(plot_df %>% filter(!is.na(KL_median)), aes(x = posterior_label, y = KL_median, color = algorithm)) +
  geom_boxplot(position = position_dodge(width = 0.7)) +
  labs(x = "Posterior (group / name)", y = "KL", color = "Algorithm") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))



# 假设你使用 plot_df，有每个 seed 的 W1_relative 或 KL_median
df <- plot_df %>%
  filter(!is.na(KL_median)) %>%
  group_by(algorithm, posterior_label) %>%
  summarise(
    median_KL = median(KL_median, na.rm = TRUE),
    q25 = quantile(KL_median, 0.25, na.rm = TRUE),
    q75 = quantile(KL_median, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(df, aes(x = posterior_label, y = median_KL, color = algorithm)) +
  geom_point(position = position_dodge(width = 0.6), size = 2.5) +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.3, position = position_dodge(width = 0.6)) +
  labs(x = "Posterior (group / name)", y = "KL (median ± IQR)", color = "Algorithm") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))



ggplot(plot_df %>% filter(!is.na(KL_median)), aes(x = posterior_label, y = KL_median, color = algorithm)) +
  geom_point(position = position_jitter(width = 0.15), alpha = 0.6) +   # 每个 seed
  stat_summary(fun = median, geom = "point", size = 5, shape = 18) +      # 中位数
  labs(x = "Posterior (group / name)", y = "KL", color = "Algorithm") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))


# log 






plot_df <- files_summary %>%
  mutate(
    posterior_label = interaction(group, posterior, sep = " / "),
    posterior_label = factor(
      posterior_label,
      levels = unique(posterior_label[order(group, posterior)])
    )
  )

# W1
plot_W1 <- plot_df %>%
  filter(!is.na(W1_relative_median)) %>%
  ggplot(aes(x = posterior_label, y = W1_relative_median, color = algorithm)) +
  geom_point(position = position_jitter(width = 0.15), alpha = 0.65) +
  stat_summary(fun = median, geom = "point", size = 3, shape = 18) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.6) +
  labs(
    x = "Posterior (group / name)",
    y = "Relative W1 distance",
    color = "Algorithm"
  ) +
  scale_y_continuous(transform = "log") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

# W2
plot_W2 <- plot_df %>%
  filter(!is.na(W2_relative_median)) %>%
  ggplot(aes(x = posterior_label, y = W2_relative_median, color = algorithm)) +
  geom_point(position = position_jitter(width = 0.15), alpha = 0.65) +
  stat_summary(fun = median, geom = "point", size = 3, shape = 18) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.6) +
  labs(
    x = "Posterior (group / name)",
    y = "Relative W2 distance",
    color = "Algorithm"
  ) +
  scale_y_continuous(transform = "log") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

# RMSE
plot_RMSE <- plot_df %>%
  filter(!is.na(RMSE_median)) %>%
  ggplot(aes(x = posterior_label, y = RMSE_median, color = algorithm)) +
  geom_point(position = position_jitter(width = 0.15), alpha = 0.65) +
  stat_summary(fun = median, geom = "point", size = 3, shape = 18) +
  labs(
    x = "Posterior (group / name)",
    y = "RMSE",
    color = "Algorithm"
  ) +
  scale_y_continuous(transform = "log") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

# KL
plot_KL <- plot_df %>%
  filter(!is.na(KL_median)) %>%
  ggplot(aes(x = posterior_label, y = KL_median, color = algorithm)) +
  geom_point(position = position_jitter(width = 0.15), alpha = 0.65) +
  stat_summary(fun = median, geom = "point", size = 3, shape = 18) +
  labs(
    x = "Posterior (group / name)",
    y = "KL",
    color = "Algorithm"
  ) +
  scale_y_continuous(transform = "log") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

# ESS per GE
plot_ESS_GE <- plot_df %>%
  filter(!is.na(ESS_per_GE_median)) %>%
  ggplot(aes(x = posterior_label, y = ESS_per_GE_median, color = algorithm)) +
  geom_point(position = position_jitter(width = 0.15), alpha = 0.65) +
  stat_summary(fun = median, geom = "point", size = 3, shape = 18) +
  labs(
    x = "Posterior (group / name)",
    y = "ESS per GE",
    color = "Algorithm"
  ) +
  scale_y_continuous(transform = "log") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

# ESS per time
plot_ESS_time <- plot_df %>%
  filter(!is.na(ESS_per_time_median)) %>%
  ggplot(aes(x = posterior_label, y = ESS_per_time_median, color = algorithm)) +
  geom_point(position = position_jitter(width = 0.15), alpha = 0.65) +
  stat_summary(fun = median, geom = "point", size = 3, shape = 18) +
  labs(
    x = "Posterior (group / name)",
    y = "ESS per time",
    color = "Algorithm"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

plot_W1
plot_W2
plot_RMSE
plot_KL
plot_ESS_GE
plot_ESS_time
