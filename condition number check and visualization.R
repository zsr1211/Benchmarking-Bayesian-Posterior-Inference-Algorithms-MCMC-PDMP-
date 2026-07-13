# This script is about classifing the posterior based on their complexity


# load in packages
library(cmdstanr)
library(posterior)
library(posteriordb)
library(R.utils)   # for withTimeout()

# set seed
set.seed(123)



# ----------------------------------------
# Load local PosteriorDB
pdb <- pdb_local("/Users/zhangsr/Desktop/UvA/Thesis/posteriorDB/posteriordb")
# Get all posterior IDs available in the database
posterior_ids <- posterior_names(pdb)
length(posterior_ids)
getwd()
setwd("/Users/zhangsr/Desktop/UvA/Thesis/programming/Pipeline")



# ----------------------------------------
# Set up result storage
results_file <- "results_condition_number.rds"
results <- readRDS(results_file)
# If results already exist, load them
if (file.exists(results_file)) {
  results <- readRDS(results_file)
  cat("Loaded existing results:", length(results), "posteriors\n")
} else {
  results <- list()
}




# ----------------------------------------
# Function to compute condition number
compute_cond_dim <- function(draws) {
  
  # Convert posterior draws to matrix format
  draws_mat <- as_draws_matrix(draws)
  draws_mat <- as.matrix(draws_mat)
  
  # Remove log-probability related columns
  drop_cols <- grepl("^lp__", colnames(draws_mat)) |
    grepl("^lp_approx__", colnames(draws_mat))
  draws_mat <- draws_mat[, !drop_cols, drop = FALSE]
  
  # Remove constant parameters
  keep <- apply(draws_mat, 2, sd, na.rm = TRUE) > 0
  draws_mat <- draws_mat[, keep, drop = FALSE]
  
  # Get parameter dimension
  dim_unconstrained <- ncol(draws_mat)
  
  # If no parameters remain, return NA
  if (dim_unconstrained < 1) {
    return(list(cond = NA_real_, dim = 0))
  }
  
  # Covariance matrix of posterior samples
  cov_mat <- tryCatch(cov(draws_mat), error = function(e) NULL)
  
  # Numerical stability check
  if (is.null(cov_mat) || any(!is.finite(cov_mat))) {
    return(list(cond = NA_real_, dim = dim_unconstrained))
  }
  
  # Eigenvalues
  eigvals <- tryCatch(
    eigen(cov_mat, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NULL
  )
  if (is.null(eigvals)) {
    return(list(cond = NA_real_, dim = dim_unconstrained))
  }
  
  # Remove non-finite or near-zero eigenvalues
  eigvals <- eigvals[is.finite(eigvals) & eigvals > 1e-12]
  
  # condition number
  if (length(eigvals) == 1) {
    cond <- 1
  } else if (length(eigvals) < 1) {
    cond <- NA_real_
  } else {
    cond <- max(eigvals) / min(eigvals)
  }
  
  list(cond = cond, dim = dim_unconstrained)
}




# ----------------------------------------
# Main loop
for (post_name in posterior_ids) {
  
  # Skip already processed posteriors
  if (post_name %in% names(results)) {
    cat("Skipping:", post_name, "\n")
    next
  }
  cat("==========\n")
  cat("Running:", post_name, "\n")
  
  # Initialize variables
  model_name <- NA
  data_name <- NA
  dim_unconstrained <- NA_integer_
  cond_laplace <- NA_real_
  cond_ref <- NA_real_
  has_unc_ref <- 0
  
  # Initialize failure/timeout flags
  optimize_failed <- 0
  optimize_timeout <- 0
  laplace_failed <- 0
  laplace_timeout <- 0
  
  # initialize error message
  optimize_error_msg <- NA_character_
  laplace_error_msg <- NA_character_
  ref_error_msg <- NA_character_
  outer_error_msg <- NA_character_
  
  tryCatch({
    
    # Load posterior
    po <- posterior(post_name, pdb)
    
    # Get model and data file paths
    model_path <- stan_code_file_path(po)
    data_path  <- data_file_path(po)
    model_name <- basename(model_path)
    data_name  <- basename(data_path)
    
    # Compile Stan model
    model <- cmdstan_model(model_path, force_recompile = TRUE)
    
    # Run optimization
    fit_opt <- tryCatch(
      withTimeout(
        model$optimize(
          data = data_path,
          jacobian = TRUE,
          iter = 10000
        ),
        timeout = 300   # 180 seconds
      ),
      
      # If timeout, marked optimize_timeout as 1
      TimeoutException = function(e) {
        optimize_timeout <<- 1
        optimize_error_msg <<- "Optimization timed out"
        NULL
      },
      
      # If has error, marked optimize_failed as 1
      error = function(e) {
        optimize_failed <<- 1
        optimize_error_msg <<- conditionMessage(e)
        NULL
      }
    )
    
    # Run Laplace approximation
    if (!is.null(fit_opt)) {
      fit_laplace <- tryCatch(
        withTimeout(
          model$laplace(
            data = data_path,
            mode = fit_opt,
            draws = 5000
          ),
          timeout = 300   # 300 seconds
        ),
        
        # If timeout, marked laplace_timeout as 1
        TimeoutException = function(e) {
          laplace_timeout <<- 1
          laplace_error_msg <<- "Laplace timed out"
          NULL
        },
        
        # If has error, marked optimize_failed as 1
        error = function(e) {
          laplace_failed <<- 1
          laplace_error_msg <<- conditionMessage(e)
          NULL
        }
      )
    } else {
      fit_laplace <- NULL
    }
    
    # Laplace condition number
    
    if (!is.null(fit_laplace)) {
      unc_draws_lap <- tryCatch(
        fit_laplace$unconstrain_draws(format = "matrix"),
        error = function(e) {
          laplace_error_msg <<- conditionMessage(e)
          NULL
        }
      )
      if (!is.null(unc_draws_lap)) {
        lap_info <- compute_cond_dim(unc_draws_lap)
        cond_laplace <- lap_info$cond
        dim_unconstrained <- lap_info$dim
      } else {
        cond_laplace <- NA_real_
        dim_unconstrained <- NA_integer_
      }
    }
    
    # Reference posterior condition number
    draws_ref <- tryCatch(
      reference_posterior_draws(po),
      error = function(e) {
        ref_error_msg <<- conditionMessage(e)
        NULL
      }
    )
    
    # If there is no reference posteior, marked as NA
    if (!is.null(draws_ref) && !is.null(fit_laplace)) {
      unc_draws_ref <- tryCatch(
        fit_laplace$unconstrain_draws(draws = draws_ref, format = "matrix"),
        error = function(e) {
          ref_error_msg <<- conditionMessage(e)
          NULL
        }
      )
      
      if (!is.null(unc_draws_ref)) {
        ref_info <- compute_cond_dim(unc_draws_ref)
        cond_ref <- ref_info$cond
        has_unc_ref <- 1
      } else {
        cond_ref <- NA_real_
        has_unc_ref <- 0
      }
    }
    
  }, error = function(e) {
    outer_error_msg <<- conditionMessage(e)
  })
  
  # Store results
  results[[post_name]] <- data.frame(
    posterior = post_name,
    model = model_name,
    data = data_name,
    dim = dim_unconstrained,
    cond_laplace = cond_laplace,
    cond_reference = cond_ref,
    has_unc_ref = has_unc_ref,
    optimize_failed = optimize_failed,
    optimize_timeout = optimize_timeout,
    laplace_failed = laplace_failed,
    laplace_timeout = laplace_timeout,
    optimize_error_msg = optimize_error_msg,
    laplace_error_msg = laplace_error_msg,
    ref_error_msg = ref_error_msg,
    outer_error_msg = outer_error_msg,
    stringsAsFactors = FALSE
  )
  
  # Save results for each posterior
  saveRDS(results, results_file)
  
  cat("dim:", dim_unconstrained, "\n")
  cat("Laplace:", cond_laplace, "\n")
  cat("Reference:", cond_ref, "\n")
}




# ----------------------------------------
# Summarise the result
summary_df <- do.call(rbind, results)
print(summary_df)




# ----------------------------------------
# Compare Laplace vs Reference condition numbers
# Keep only posteriors with reference samples
df_compare <- subset(summary_df, has_unc_ref != 0)

# Extract relevant columns
df_compare_cn <- df_compare[, c("cond_laplace", "cond_reference", "laplace_failed")]

# Compute ranks
df_compare_cn$rank_laplace <- rank(df_compare_cn$cond_laplace, na.last = "keep")
df_compare_cn$rank_ref <- rank(df_compare_cn$cond_reference, na.last = "keep")

# Sort by Laplace condition number
df_compare_cn <- df_compare_cn[order(df_compare_cn$cond_laplace), ]

# log version of condition number
df_compare_cn$log_laplace_cn <- log(df_compare_cn$cond_laplace)
df_compare_cn$log_ref_cn <- log(df_compare_cn$cond_reference)




# ----------------------------------------
# correlation
# raw condition number values
cor(df_compare_cn$cond_laplace, df_compare_cn$cond_reference, use = "complete.obs")
# log condition number values
cor(df_compare_cn$log_laplace_cn, df_compare_cn$log_ref_cn, use = "complete.obs")
# Rank
cor(df_compare_cn$rank_laplace, df_compare_cn$rank_ref, use = "complete.obs")
# How many NA in cond_laplace
sum(is.na(summary_df$cond_laplace))
sum(summary_df$has_unc_ref!=0)

# check if has referense posteior
summary_df$posterior <- as.character(summary_df$posterior)
summary_df$has_ref <- sapply(summary_df$posterior, function(post_name) {
  out <- tryCatch({
    po <- posterior(post_name, pdb)
    ref <- reference_posterior_draws(po)
    !is.null(ref)
  }, error = function(e) {
    FALSE
  })
  
  as.integer(out)
})
sum(summary_df$has_ref)  # 47, correct




# keep NA as the most complex case
summary_df$complexity_withNA <- NA_character_

x <- summary_df$cond_laplace
x_for_sort <- ifelse(is.na(x), Inf, x)  # see NA as Inf

ord <- order(x_for_sort)
n <- length(ord)

ranks_all <- rank(x_for_sort, ties.method = "first")

summary_df$complexity_withNA <- cut(
  ranks_all,
  breaks = c(0, n/3, 2*n/3, n),
  labels = c("simple", "moderate", "complex"),
  include.lowest = TRUE
)
table(summary_df$complexity_withNA, useNA = "ifany")



library(dplyr)
summary_df <- summary_df %>%
  select(posterior, cond_laplace, cond_reference, has_unc_ref, has_ref, complexity_withNA, everything())




# ----------------------------------------
# plot condition numbers
plot(df_compare_cn$log_laplace_cn,
     df_compare_cn$log_ref_cn,
     xlab = "log10 Laplace condition number",
     ylab = "log10 Reference condition number",
     main = "Laplace vs Reference (log10)",
     pch = 16)
abline(a = 0, b = 1, col = "red", lty = 2)

