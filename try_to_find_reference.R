# This script runs only MCMC for posteriorDB posteriors with has_ref == 0,
# grouped by complexity_withNA, using seed 4711.
# It saves both the CmdStanMCMC fit object and unconstrained draws.
# try rto find reference draws by ourselves

library(posteriordb)
library(cmdstanr)
library(posterior)
library(stringr)
library(purrr)
library(dplyr)

# set posteriordb
pdb <- pdb_local("/Users/zhangsr/Desktop/UvA/Thesis/posteriorDB/posteriordb")

# keep only posteriors where has_ref == 0
summary_df_no_ref <- summary_df[summary_df$has_ref == 0, ]

# Posterior groups
posterior_groups <- list(
  simple = summary_df_no_ref$posterior[
    summary_df_no_ref$complexity_withNA == "simple"
  ],
  moderate = summary_df_no_ref$posterior[
    summary_df_no_ref$complexity_withNA == "moderate"
  ],
  complex = summary_df_no_ref$posterior[
    summary_df_no_ref$complexity_withNA == "complex"
  ]
)

# check
posterior_groups

# fixed seed
seed <- 4711

# output folder
results_dir <- "/Users/zhangsr/Desktop/UvA/Thesis/programming/Pipeline/try_to_find_reference"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# optional: keep csv files too
# csv_dir <- file.path(results_dir, "cmdstan_csv")
# dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)

# helper for safe filenames
safe_name <- function(x) {
  gsub("[^A-Za-z0-9_\\-]+", "_", x)
}

# main loops
for (group_name in names(posterior_groups)) {
  
  for (post_name in posterior_groups[[group_name]]) {
    
    cat("\n==========\n")
    cat("Group:", group_name, "Posterior:", post_name, "\n")
    
    post_safe <- safe_name(post_name)
    
    # output paths
    fit_file <- file.path(
      results_dir,
      paste0(group_name, "_", post_safe, "_MCMC_seed", seed, "_fit.rds")
    )
    
    unconstrained_file <- file.path(
      results_dir,
      paste0(group_name, "_", post_safe, "_MCMC_seed", seed, "_unconstrained_matrix.rds")
    )
    
    # skip if both files already exist
    if (file.exists(fit_file) && file.exists(unconstrained_file)) {
      cat("Fit and unconstrained draws already exist -> skip\n")
      next
    }
    
    tryCatch({
      
      # get posterior
      po <- posterior(post_name, pdb)
      
      # get Stan model and data path
      model_path <- stan_code_file_path(po)
      data_path <- data_file_path(po)
      
      cat("Model:", model_path, "\n")
      cat("Data:", data_path, "\n")
      
      # compile model
      mod <- cmdstan_model(model_path, force_recompile = TRUE)
      
      # run MCMC using the reference settings:
      # rstan iter = 20000, warmup = 10000, thin = 10
      # cmdstanr equivalent:
      # iter_warmup = 10000, iter_sampling = 10000, thin = 10
      fit <- mod$sample(
        data = data_path,
        chains = 10,
        iter_warmup = 10000,
        iter_sampling = 10000,
        thin = 10,
        seed = seed,
        adapt_delta = 0.8
      )
      
      # save fit object
      # Better than plain saveRDS(fit, ...) for CmdStanMCMC objects
      fit$save_object(file = fit_file)
      cat("Saved fit:", fit_file, "\n")
      
      # extract and save unconstrained draws
      unconstrained_matrix <- fit$unconstrain_draws(format = "matrix")
      saveRDS(unconstrained_matrix, unconstrained_file)
      cat("Saved unconstrained draws:", unconstrained_file, "\n")
      
    }, error = function(e) {
      cat("ERROR for posterior:", post_name, "\n")
      cat("Message:", e$message, "\n")
    })
  }
}

# try1 <- readRDS("/Users/zhangsr/Desktop/UvA/Thesis/programming/Pipeline/try_to_find_reference/simple_dogs-dogs_hierarchical_MCMC_seed4711_fit.rds")
# try1$diagnostic_summary()
# try1$summary()
# try1$ndraws(variables = variables)
# variables <- c("a", "b")
# summ <- try1$summary(variables = variables)
# try1$metadata()$model_params
# draws <- try1$draws(variables = variables)
# posterior::ndraws(draws)
# posterior::nchains(draws)
# summ$rhat
# summ$ess_bulk
# summ$ess_tail
# 
# draws_array <- posterior::as_draws_array(
#   try1$draws(variables = variables)
# )
# 
# mean_lag1_ac <- sapply(variables, function(v) {
#   
#   ac_by_chain <- sapply(seq_len(dim(draws_array)[2]), function(ch) {
#     
#     x <- draws_array[, ch, v]
#     
#     if (length(x) < 2 || sd(x) == 0) {
#       return(NA_real_)
#     }
#     
#     cor(x[-length(x)], x[-1])
#   })
#   
#   mean(ac_by_chain, na.rm = TRUE)
# })
# 
# mean_lag1_ac
# abs(mean_lag1_ac) < 0.05



a <- c(2.06787947520188, 1.98406323497227, 1.96730344820381, 1.9577588927253, 1.95885957253166, 1.93786829596686, 
        1.98889939415628, 2.00217411797618, 2.05892672314893, 2.09048669550594)
mean(a)

a1 <- c(1.99368134018704, 2.08543412544146, 1.93962475544812, 1.9660605916653, 1.98421406267258, 
        1.9971157673996, 2.04360655068323, 1.9977567315046, 2.05054124370642, 1.99763763366879)
mean(a1)




# check if the sampled results meet the standard of reference draws

# function: get parameter names from Stan code
get_stan_parameter_names <- function(posterior_name, pdb) {
  
  # get posterior object
  po <- posterior(posterior_name, pdb)
  
  # get Stan model code as a character
  stan_code_text <- stan_code(po)
  
  # remove comments in the Stan code
  code <- str_replace_all(stan_code_text, "//.*", "")
  code <- str_replace_all(code, "/\\*.*?\\*/", "")
  
  # extract parameters block
  m <- str_match(code, "parameters\\s*\\{([\\s\\S]*?)\\n\\}")
  
  if (is.na(m[1, 2])) {
    return(character(0))
  }
  
  block <- m[1, 2]
  
  # split declarations by ;
  declarations <- str_split(block, ";", simplify = FALSE)[[1]]
  declarations <- str_trim(declarations)
  declarations <- declarations[declarations != ""]
  
  # remove constraints, e.g. <lower=0>
  declarations <- str_replace_all(declarations, "<[^>]*>", "")
  declarations <- str_squish(declarations)
  
  # extract parameter names
  param_names <- map_chr(declarations, function(x) {
    tokens <- str_split(x, "\\s+", simplify = TRUE)
    last <- tokens[length(tokens)]
    str_replace(last, "\\[.*\\]$", "")
  })
  
  unique(param_names)
}


# function: expand parameter names
expand_parameter_names <- function(fit, param_names) {
  
  all_vars <- fit$metadata()$model_params
  
  all_vars <- all_vars[all_vars != "lp__"]
  
  variables <- unlist(lapply(param_names, function(p) {
    all_vars[all_vars == p | str_detect(all_vars, paste0("^", p, "\\["))]
  }))
  
  unique(variables)
}


# function: compute mean lag-1 autocorrelation
compute_mean_lag1_ac <- function(fit, variables) {
  
  draws_array <- posterior::as_draws_array(
    fit$draws(variables = variables)
  )
  
  mean_lag1_ac <- sapply(variables, function(v) {
    
    ac_by_chain <- sapply(seq_len(dim(draws_array)[2]), function(ch) {
      
      x <- draws_array[, ch, v]
      
      if (length(x) < 2 || sd(x) == 0) {
        return(NA_real_)
      }
      
      cor(x[-length(x)], x[-1])
    })
    
    mean(ac_by_chain, na.rm = TRUE)
  })
  
  names(mean_lag1_ac) <- variables
  
  mean_lag1_ac
}


# function: compute diagnostics for one fit
compute_checks_one_fit <- function(fit, variables) {
  
  # ESS and Rhat
  summ <- fit$summary(variables = variables)
  
  ess_bulk <- summ$ess_bulk
  ess_tail <- summ$ess_tail
  r_hat <- summ$rhat
  
  # divergence and E-BFMI
  diag_sum <- fit$diagnostic_summary()
  
  divergent_transitions <- diag_sum$num_divergent
  efmi <- diag_sum$ebfmi
  
  # lag-1 autocorrelation
  mean_lag1_ac <- compute_mean_lag1_ac(fit, variables)
  
  # draw info
  draws <- fit$draws(variables = variables)
  
  ndraws <- posterior::ndraws(draws)
  nchains <- posterior::nchains(draws)
  
  # checks
  checks <- list(
    ndraws_is_10k = ndraws >= 10000,
    nchains_is_gte_4 = nchains >= 4,
    ess_within_bounds = all(
      is.finite(ess_bulk),
      is.finite(ess_tail),
      ess_bulk > 0,
      ess_tail > 0,
      na.rm = TRUE
    ),
    r_hat_below_1_01 = all(r_hat < 1.01, na.rm = TRUE),
    efmi_above_0_2 = all(efmi > 0.2, na.rm = TRUE),
    abs_mean_lag1_ac_below_0_05 = all(abs(mean_lag1_ac) < 0.05, na.rm = TRUE),
    no_divergent_transitions = all(divergent_transitions == 0)
  )
  
  # full diagnostics, if you also want to save them
  diagnostics <- list(
    diagnostic_information = list(
      names = variables
    ),
    ndraws = ndraws,
    nchains = nchains,
    effective_sample_size_bulk = ess_bulk,
    effective_sample_size_tail = ess_tail,
    r_hat = r_hat,
    divergent_transitions = divergent_transitions,
    expected_fraction_of_missing_information = efmi,
    mean_lag1_ac = mean_lag1_ac,
    checks_made = checks
  )
  
  diagnostics
}


# main loop: read fit.rds and make table

# diagnostics_dir <- file.path(results_dir, "diagnostics")
# dir.create(diagnostics_dir, showWarnings = FALSE, recursive = TRUE)

# only use fit files that already exist
fit_files <- list.files(
  results_dir,
  pattern = "_MCMC_seed4711_fit\\.rds$",
  full.names = TRUE
)

summary_rows <- list()

for (fit_file in fit_files) {
  
  cat("\n==========\n")
  cat("Fit file:", fit_file, "\n")
  
  file_base <- basename(fit_file)
  
  # remove suffix
  name_part <- sub("_MCMC_seed4711_fit\\.rds$", "", file_base)
  
  # infer group from filename
  group_name <- sub("_.*$", "", name_part)
  
  # infer safe posterior name from filename
  post_safe <- sub(paste0("^", group_name, "_"), "", name_part)
  
  # find original posterior name by matching safe_name()
  all_posteriors <- unlist(posterior_groups, use.names = FALSE)
  matched <- all_posteriors[safe_name(all_posteriors) == post_safe]
  
  if (length(matched) == 0) {
    cat("Cannot match posterior name, skip\n")
    
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      group = group_name,
      posterior = NA_character_,
      fit_file = fit_file,
      fit_exists = TRUE,
      n_parameters = NA_integer_,
      ndraws = NA_integer_,
      nchains = NA_integer_,
      ndraws_is_10k = NA,
      nchains_is_gte_4 = NA,
      ess_within_bounds = NA,
      r_hat_below_1_01 = NA,
      efmi_above_0_2 = NA,
      abs_mean_lag1_ac_below_0_05 = NA,
      no_divergent_transitions = NA,
      all_checks_pass = FALSE,
      max_r_hat = NA_real_,
      max_abs_mean_lag1_ac = NA_real_,
      min_efmi = NA_real_,
      total_divergent = NA_integer_,
      error = "cannot match posterior name"
    )
    
    next
  }
  
  post_name <- matched[1]
  
  # diagnostics_file <- file.path(
  #   diagnostics_dir,
  #   paste0(name_part, "_diagnostics.rds")
  # )
  
  cat("Posterior:", post_name, "\n")
  
  tryCatch({
    
    fit <- readRDS(fit_file)
    
    # get parameter names from Stan code
    param_names <- get_stan_parameter_names(post_name, pdb)
    
    # expand vector/matrix parameters to beta[1], beta[2], etc.
    variables <- expand_parameter_names(fit, param_names)
    
    cat("Parameters:", paste(variables, collapse = ", "), "\n")
    
    diagnostics <- compute_checks_one_fit(fit, variables)
    
    # save full diagnostics for this posterior
    # saveRDS(diagnostics, diagnostics_file)
    
    checks <- diagnostics$checks_made
    
    all_checks_pass <- all(
      checks$ndraws_is_10k,
      checks$nchains_is_gte_4,
      checks$ess_within_bounds,
      checks$r_hat_below_1_01,
      checks$efmi_above_0_2,
      checks$abs_mean_lag1_ac_below_0_05,
      checks$no_divergent_transitions
    )
    
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      group = group_name,
      posterior = post_name,
      fit_file = fit_file,
      fit_exists = TRUE,
      n_parameters = length(variables),
      ndraws = diagnostics$ndraws,
      nchains = diagnostics$nchains,
      ndraws_is_10k = checks$ndraws_is_10k,
      nchains_is_gte_4 = checks$nchains_is_gte_4,
      ess_within_bounds = checks$ess_within_bounds,
      r_hat_below_1_01 = checks$r_hat_below_1_01,
      efmi_above_0_2 = checks$efmi_above_0_2,
      abs_mean_lag1_ac_below_0_05 = checks$abs_mean_lag1_ac_below_0_05,
      no_divergent_transitions = checks$no_divergent_transitions,
      all_checks_pass = all_checks_pass,
      max_r_hat = max(diagnostics$r_hat, na.rm = TRUE),
      max_abs_mean_lag1_ac = max(abs(diagnostics$mean_lag1_ac), na.rm = TRUE),
      min_efmi = min(diagnostics$expected_fraction_of_missing_information, na.rm = TRUE),
      total_divergent = sum(diagnostics$divergent_transitions),
      error = NA_character_
    )
    
  }, error = function(e) {
    
    cat("ERROR:", e$message, "\n")
    
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      group = group_name,
      posterior = post_name,
      fit_file = fit_file,
      fit_exists = TRUE,
      n_parameters = NA_integer_,
      ndraws = NA_integer_,
      nchains = NA_integer_,
      ndraws_is_10k = NA,
      nchains_is_gte_4 = NA,
      ess_within_bounds = NA,
      r_hat_below_1_01 = NA,
      efmi_above_0_2 = NA,
      abs_mean_lag1_ac_below_0_05 = NA,
      no_divergent_transitions = NA,
      all_checks_pass = FALSE,
      max_r_hat = NA_real_,
      max_abs_mean_lag1_ac = NA_real_,
      min_efmi = NA_real_,
      total_divergent = NA_integer_,
      error = e$message
    )
  })
}

diagnostics_summary_df <- dplyr::bind_rows(summary_rows)

summary_file <- file.path(results_dir, "reference_check_summary_existing_fits.csv")
write.csv(diagnostics_summary_df, summary_file, row.names = FALSE)
file <- read.csv("/Users/zhangsr/Desktop/UvA/Thesis/programming/Pipeline/try_to_find_reference/reference_check_summary_existing_fits.csv")
diagnostics_summary_df
