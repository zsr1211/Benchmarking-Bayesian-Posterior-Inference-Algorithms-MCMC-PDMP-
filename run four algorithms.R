# This script is about running all posteirors in posteriordb with 4 methods 
# (1 MCMC and 3 PDMP) and 10 seeds, and saving the running reuslts as .rds files.

# PDMPSamplersR::pdmpsamplers_update()
# load in packages
library(posteriordb)
library(cmdstanr)
library(PDMPSamplersR)
library(posterior)
library(bridgestan)
# JuliaCall::julia_command("Pkg.status()")


# set posteriordb and select some posteriors
pdb <- pdb_local("/Users/zhangsr/Desktop/UvA/Thesis/posteriorDB/posteriordb")

# Posterior groups
posterior_groups <- list(
  simple = summary_df$posterior[summary_df$complexity_withNA == "simple"],
  moderate = summary_df$posterior[summary_df$complexity_withNA == "moderate"],
  complex = summary_df$posterior[summary_df$complexity_withNA == "complex"]
)
# here, summary_df is the posterior classification result from laplace condition number
# which comes from the result of script "condition number check and visulization.R"

# check
posterior_groups

# posterior_groups <- list(
#   simple = c(
#     "arK-arK",
#     "eight_schools-eight_schools_noncentered",
#     "dogs-dogs"
#   ),
#   moderate = c(
#     "diamonds-diamonds",
#     "ecdc0501-covid19imperial_v2",
#     "radon_all-radon_county_intercept"
#   ),
#   complex = c(
#     "earnings-earn_height",
#     "eight_schools-eight_schools_centered",
#     "hmm_gaussian_simulated-hmm_gaussian"
#   )
# )
# posterior_groups <- list(
#   simple = c(
#     "arK-arK")
# )


# fixed seed sets
set.seed(1234)
seeds <- sample(1:10000, 10, replace = FALSE) 
seeds


# different methods
methods <- c("MCMC", "ZigZag", "BPS", "Boomerang")


# set output folder and error file
results_dir <- "/Users/zhangsr/Desktop/UvA/Thesis/programming/Pipeline/benchmark_runs_final"
dir.create(results_dir, showWarnings = FALSE)

error_dir <- file.path(results_dir, "error_logs")
dir.create(error_dir, showWarnings = FALSE, recursive = TRUE)


# main loops
# run each posterior groups
for (group_name in names(posterior_groups)) {
  
  # run each posterior
  for (post_name in posterior_groups[[group_name]]) {
    
    cat("\n==========\n")
    cat("Group:", group_name, "Posterior:", post_name, "\n")
    
    # get posterior
    po <- posterior(post_name, pdb)
    
    # get stan model and data path
    model_path <- stan_code_file_path(po)
    data_path <- data_file_path(po)
    
    cat("Model:", model_path, "\n")
    cat("Data:", data_path, "\n")
    
    # compile model
    mod <- cmdstan_model(model_path, force_recompile = TRUE)
    
    # run each seed
    for (seed in seeds) {
      
      set.seed(seed)
      
      # run each method (algorithm)
      for (method in methods) {
        
        cat("\nMethod:", method, "Seed:", seed, "\n")
        
        # make posterior name safe for file paths
        safe_post_name <- gsub("[^A-Za-z0-9_\\-]", "_", post_name)
        
        # output file for main result
        file_out <- file.path(
          results_dir,
          paste0(
            group_name, "_", safe_post_name, "_", method, "_seed", seed, ".rds"
          )
        )
        
        # output file for MCMC unconstrained draws
        file_draws <- file.path(
          results_dir,
          paste0(
            group_name, "_", safe_post_name, "_", method, "_seed", seed,
            "_unconstrained_matrix.rds"
          )
        )
        
        # error log file
        error_out <- file.path(
          error_dir,
          paste0(
            group_name, "_", safe_post_name, "_", method, "_seed", seed,
            "_ERROR.txt"
          )
        )
        
        # skip logic
        if (method == "MCMC") {
          if (file.exists(file_out) && file.exists(file_draws)) {
            cat("MCMC result and unconstrained draws already exist → skip\n")
            next
          } else if (file.exists(file_out) && !file.exists(file_draws)) {
            cat("MCMC result exists but unconstrained draws missing → rerun MCMC\n")
          }
        } else {
          if (file.exists(file_out)) {
            cat("Result already exists → skip\n")
            next
          }
        }
        
        tryCatch({
          
          # settings for each method
          if (method == "MCMC") {
            
            fit <- mod$sample(
              data = data_path,
              chains = 4,
              parallel_chains = 4,
              iter_warmup = 2000, 
              iter_sampling = 2000,
              seed = seed,
              adapt_delta = 0.9
            )
            
            result <- fit
            
            # save MCMC result first
            saveRDS(result, file_out)
            cat("Saved MCMC result:", file_out, "\n")
            
            # extract and save unconstrained draws
            unconstrained_matrix <- fit$unconstrain_draws(format = "matrix")
            saveRDS(unconstrained_matrix, file_draws)
            cat("Saved unconstrained draws:", file_draws, "\n")
          }
          
          if (method == "ZigZag") {
            set.seed(seed)
            result <- pdmp_sample_from_stanmodel(
              model_path,
              data_path,
              flow = "PreconditionedZigZag",
              algorithm = "GridThinningStrategy",
              T = 2000,
              grid_n = 30,
              t_warmup = 1000,
              sticky = FALSE,
              adaptive_scheme = "full",
              seed = seed
            )
            
            saveRDS(result, file_out)
            cat("Saved:", file_out, "\n")
          }
          
          if (method == "BPS") {
            set.seed(seed)
            result <- pdmp_sample_from_stanmodel(
              model_path,
              data_path,
              flow = "PreconditionedBPS",
              algorithm = "GridThinningStrategy",
              T = 2000,
              grid_n = 30,
              t_warmup = 1000,
              sticky = FALSE,
              adaptive_scheme = "full",
              seed = seed
            )
            
            saveRDS(result, file_out)
            cat("Saved:", file_out, "\n")
          }
          
          if (method == "Boomerang") {
            set.seed(seed)
            result <- pdmp_sample_from_stanmodel(
              model_path,
              data_path,
              flow = "AdaptiveBoomerang",
              algorithm = "GridThinningStrategy",
              T = 2000,
              grid_n = 30,
              t_warmup = 1000,
              sticky = FALSE,
              adaptive_scheme = "full",
              seed = seed
            )
            
            saveRDS(result, file_out)
            cat("Saved:", file_out, "\n")
          }
          
          # if successful, remove old error log if it exists
          if (file.exists(error_out)) {
            file.remove(error_out)
          }
          
        }, error = function(e) {
          
          cat("ERROR:", e$message, "\n")
          
          writeLines(
            c(
              paste("time:", Sys.time()),
              paste("group:", group_name),
              paste("posterior:", post_name),
              paste("method:", method),
              paste("seed:", seed),
              paste("model_path:", model_path),
              paste("data_path:", data_path),
              paste("file_out:", file_out),
              paste("file_draws:", file_draws),
              paste("error:", e$message)
            ),
            error_out
          )
          
          cat("Saved error log:", error_out, "\n")
        })
      }
    }
  }
}


try <- readRDS("/Users/zhangsr/Desktop/UvA/Thesis/programming/Pipeline/benchmark_runs_try/simple_arK-arK_MCMC_seed7452.rds")
try$diagnostic_summary()
try$summary()




