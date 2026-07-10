# This script is about transforming reference posterior draws 
# from constrained space to unconstrained space.


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





###############
# Get unconstrained reference posterior draws -----------------------------
###############

# load in posteriordb
pdb <- pdb_local("/Users/zhangsr/Desktop/UvA/Thesis/posteriorDB/posteriordb")





# normalized the varianbels name in reference posterior
normalize_ref_name <- function(x) {
  x %>%
    str_replace("\\[.*\\]$", "")        # theta[1] -> theta
}





# function to extract the parameters from Stan model
extract_model_parameters <- function(posterior_name, pdb) {
  
  # get posterior object
  po <- posterior(posterior_name, pdb)
  
  # get Stan model code as a character
  stan_code_text <- stan_code(po)
  
  # remove comments in the Stan code
  code <- str_replace_all(stan_code_text, "//.*", "")     # remove //
  code <- str_replace_all(code, "/\\*.*?\\*/", "")        # remove /* ... */ 
  
  # exatrct the parameters block in stan code
  m <- str_match(code, "parameters\\s*\\{([\\s\\S]*?)\\n\\}")
  if (is.na(m[1, 2])) {
    return(character(0))
  }
  block <- m[1, 2]
  
  # Split into individual parameter declarations using ";"
  declarations <- str_split(block, ";", simplify = FALSE)[[1]]
  declarations <- str_trim(declarations)  # Trim leading and trailing spaces from each declaration
  declarations <- declarations[declarations != ""] # Remove empty strings (possibly caused by the ";")
  
  # remove the constrain, like <lower=0>
  declarations <- str_replace_all(declarations, "<[^>]*>", "")
  declarations <- str_squish(declarations)
  
  # extract parameter names
  param_names <- map_chr(declarations, function(x) {
    # Split each declaration by spaces into tokens and take the last token, which is the variable name
    tokens <- str_split(x, "\\s+", simplify = TRUE)
    last <- tokens[length(tokens)]
    # Remove any dimension indices, e.g., from "theta1[1]" -> "theta1"
    str_replace(last, "\\[.*\\]$", "")
  })
  
  # return the unique parameter names
  unique(param_names)
}





# function to check one posterior if their reference draws variables and base parameters are matched
check_one_posterior <- function(posterior_name, pdb, ref_override = NULL) {
  
  # get posterior object
  po <- posterior(posterior_name, pdb)
  
  # get the reference draws
  if (!is.null(ref_override)) {
    ref <- ref_override
  } else {
    ref <- tryCatch(
      reference_posterior_draws(po), 
      error = function(e) NULL)
  }

  # if no reference draws, return NULL
  if (is.null(ref)) {
    return(NULL)
  }
  
  # get model parameters
  model_params <- extract_model_parameters(posterior_name, pdb)
  
  # get reference varaible name and normaolize them
  ref_params <- if (inherits(ref, "draws") || 
                    inherits(ref, "draws_list") || 
                    inherits(ref, "draws_df")) {
    posterior::variables(ref)
  } else {
    unique(unlist(lapply(ref, names)))
  } 
  ref_params <- ref_params %>%
    normalize_ref_name() %>%
    unique()

  # summarize and return the result
  tibble(
    posterior_name = posterior_name,
    match = setequal(model_params, ref_params),       # check if ref_params match model_params
    model_params = paste(model_params, collapse = ", "),
    ref_params = paste(ref_params, collapse = ", "),
    model_only = paste(setdiff(model_params, ref_params), collapse = ", "),
    ref_only = paste(setdiff(ref_params, model_params), collapse = ", ")
  )
}





# run checks for all posteriors
result_df <- map_dfr(
  posterior_names(pdb),
  ~ check_one_posterior(.x, pdb, ref_override = NULL)
)

# check if "match" is FALSE
result_df$posterior_name[result_df$match == FALSE]
result_df[result_df$match == FALSE, ]
# there are two posteriors has "match" as FALSE: 
# "eight_schools-eight_schools_noncentered" and "gp_pois_regr-gp_pois_regr"





# first, manually transform varaibles to the base parameters
# for "eight_schools-eight_schools_noncentered", transform "theta" to "theta_trans"
po <- posterior("eight_schools-eight_schools_noncentered", pdb)  # load in the posterior
ref_ei <- reference_posterior_draws(po)                          # get reference draws (10 chains)


# lapply transformation for all the chains 
ref_base_ei <- lapply(ref_ei, function(chain) {
  
  theta_mat <- do.call(cbind, chain[grep("^theta\\[", names(chain))])
  mu_vec <- chain$mu
  tau_vec <- chain$tau
  
  # get theta_trans
  theta_trans_mat <- t( (t(theta_mat) - mu_vec) / tau_vec ) 
  
  # save it to chain
  for (j in seq_len(ncol(theta_trans_mat))) {
    chain[[paste0("theta_trans[", j, "]")]] <- theta_trans_mat[, j]
  }
  
  # remove theta[1..J]
  chain <- chain[!grepl("^theta\\[", names(chain))]
  chain
})


# change the order of parameters
ref_base_ei <- lapply(ref_base_ei, function(chain) {
  # get theta_trans column names，order it
  theta_names <- grep("^theta_trans\\[", names(chain), value = TRUE)
  theta_names <- theta_names[order(as.numeric(str_extract(theta_names, "\\d+")))]
  
  # other parameters
  other_names <- setdiff(names(chain), theta_names)
  
  # organized in the new order
  chain[c(theta_names, other_names)]
})


# check if the variables in reference draws are the same as the base parameters
check_one_posterior("eight_schools-eight_schools_noncentered", pdb, ref_override = ref_base_ei)


# transform to matrix
ref_base_ei_matrix <- posterior::as_draws_matrix(ref_base_ei)
colnames(ref_base_ei_matrix)
dim(ref_base_ei_matrix)
class(ref_base_ei_matrix)





# for "gp_pois_regr-gp_pois_regr", 
# manually transofrom the variable "f" to the base parameter "f_tilde"
po <- posterior("gp_pois_regr-gp_pois_regr", pdb)
ref_gp <- reference_posterior_draws(po)
dat <- pdb_data(po)    # get data
x <- dat$x
N <- length(x) 


# transform to f_tilde
ref_base_gp <- lapply(ref_gp, function(chain) {
  # chain is a draws_list，each parameter is vector
  f_mat <- do.call(cbind, chain[grep("^f\\[", names(chain))])  # iterations x N
  rho_vec <- chain$rho
  alpha_vec <- chain$alpha
  
  n_iter <- nrow(f_mat)
  
  # pre-allocate f_tilde matrix
  f_tilde_mat <- matrix(NA_real_, nrow = n_iter, ncol = N)
  
  for (i in seq_len(n_iter)) {
    cov_mat <- outer(x, x, function(xi, xj) alpha_vec[i]^2 * exp(-(xi - xj)^2 / (2 * rho_vec[i]^2)))
    cov_mat <- cov_mat + diag(1e-10, N)
    
    U <- chol(cov_mat)         # upper trizngular matrix
    L <- t(U)                  # lower triangular matrix
    
    f_tilde_mat[i, ] <- forwardsolve(L, f_mat[i, ])  # get f_tilde
  }
  
  # save to chain
  for (j in seq_len(N)) {
    chain[[paste0("f_tilde[", j, "]")]] <- f_tilde_mat[, j]
  }
  
  # remove the transformed f[1..N]
  chain <- chain[!grepl("^f\\[", names(chain))]
  
  chain
})


# check again if match 
check_one_posterior("gp_pois_regr-gp_pois_regr", pdb, ref_override = ref_base_gp)


# transfom to matrix
ref_base_gp_matrix <- posterior::as_draws_matrix(ref_base_gp)
colnames(ref_base_gp_matrix)
dim(ref_base_gp_matrix)
class(ref_base_gp_matrix)






# get the unconstrained draws
# some special posteriors (need to be manually transformed)
special_posteriors <- result_df$posterior_name[result_df$match == FALSE]
special_posteriors
ref_override_list <- list(
  "eight_schools-eight_schools_noncentered" = ref_base_ei_matrix,
  "gp_pois_regr-gp_pois_regr" = ref_base_gp_matrix
)






# functions to get the unconstrained draws
get_ref_uc <- function(posterior_name, pdb, 
                       ref_override_list = NULL,
                       save_dir = ".",
                       save_rds = FALSE,
                       return_draws = TRUE,
                       verbose = TRUE) {
  
  # get the file saving path
  rds_path <- file.path(save_dir, paste0(posterior_name, "_unconstrained_ref.rds"))
  
  # if file exists, then skip it and print messages
  if (file.exists(rds_path)) {
    if (verbose) message("Posterior '", posterior_name, "' RDS already exists, skipping save.")
    if (return_draws) return(readRDS(rds_path)) else return(invisible(NULL))
  }
  
  # make sure the directory exists
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  # get posterior objects
  po <- posterior(posterior_name, pdb)
  
  # special posterior need to using the manually transformed unconstrained reference
  if (!is.null(ref_override_list) && posterior_name %in% names(ref_override_list)) {
    u_ref <- ref_override_list[[posterior_name]]
    if (save_rds) saveRDS(u_ref, rds_path)
    if (verbose) message("Posterior '", posterior_name, "' uses ref_override from manual list and saved to: ", rds_path)
    if (return_draws) return(u_ref) else return(invisible(NULL))
  }
  
  # try to get reference draws
  ref <- tryCatch(
    reference_posterior_draws(po),
    error = function(e) NULL
  )
  if (is.null(ref)) return(invisible(NULL))  # if no reference draws，return NULL
  
  # compile model and generate fit_dummy (because unconstrain_draws() need fit)
  fit_dummy <- tryCatch({
    model_path <- stan_code_file_path(po)
    data_path  <- data_file_path(po)
    mod <- cmdstan_model(model_path, force_recompile = TRUE)
    mod$sample(
      data = data_path,
      chains = 1,
      iter_warmup = 1,
      iter_sampling = 1,
      refresh = 0
    )
  }, error = function(e) {
    # return error information
    list(error = paste0("fit generation failed: ", e$message))
  })
  
  # if fit_dummy is error list，return the error list
  if (is.list(fit_dummy) && "error" %in% names(fit_dummy)) {
    if (verbose) message("Posterior '", posterior_name, "' fit generation error: ", fit_dummy$error)
    return(fit_dummy)
  }
  
  # transform to unconstrained draws
  u_ref <- tryCatch({
    fit_dummy$unconstrain_draws(draws = ref, format = "matrix")
  }, error = function(e) {
    list(error = paste0("unconstrain_draws failed: ", e$message))
  })
  
  # if u_ref is error list，return the error list
  if (is.list(u_ref) && "error" %in% names(u_ref)) {
    if (verbose) message("Posterior '", posterior_name, "' unconstrain_draws error: ", u_ref$error)
    return(u_ref)
  }
  
  # save the unconstrained draws
  if (save_rds) {
    saveRDS(u_ref, rds_path)
    if (verbose) message("Posterior '", posterior_name, "' unconstrained draws saved to: ", rds_path)
  }
  
  # based on argument"return_draws", return draws or not
  if (return_draws) return(u_ref) else return(invisible(NULL))
}


# try if this function works
# check one of the special posteriors
try <- get_ref_uc("eight_schools-eight_schools_noncentered", pdb, ref_override_list = ref_override_list)

# check one of the normal posteriors
try <- get_ref_uc("eight_schools-eight_schools_centered", pdb, ref_override_list = ref_override_list)

# get the constrained reference draws
po <- posterior("eight_schools-eight_schools_centered", pdb)
ref <- reference_posterior_draws(po)
class(ref)
mod <- cmdstan_model(model_path, force_recompile = TRUE)
model_path <- stan_code_file_path(po)
data_path <- data_file_path(po)
fit <- mod$sample(
  data = data_path,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1,
  iter_sampling = 1,
  seed = 1004,
  adapt_delta = 0.9
)
# get the uncosntrained draws
u_ref <- fit$unconstrain_draws(draws = ref, format = "matrix")
# check if the function result "try" and the manually result "u_ref" is identical
identical(try, u_ref)
# # check if the saved rds file is the same with u_ref
# saveRDS(try, file = "try_uref_eight_schools_noncentered.rds")
# try_file <- readRDS("try_uref_eight_schools_noncentered.rds")
# class(try_file)
# identical(try_file, u_ref)






# get all unconstrained reference draws
posterior_names_all <- posterior_names(pdb)

# set saving directory
save_dir <- "/Users/zhangsr/Desktop/UvA/Thesis/programming/Pipeline/unconstrained reference draws"
# if directory does not exist, create one
if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)


# for loops to get all unconstrained reference draws
for (pname in posterior_names_all) {
  get_ref_uc(
    posterior_name = pname,
    pdb = pdb,
    ref_override_list = ref_override_list,
    save_dir = save_dir,
    save_rds = TRUE,
    return_draws = FALSE,  # do not return draws
    verbose = TRUE
  )
}

# check one of the saved .rds file
saved_ref <- readRDS("/Users/zhangsr/Desktop/UvA/Thesis/programming/Pipeline/unconstrained reference draws/eight_schools-eight_schools_centered_unconstrained_ref.rds")
# if they are identical to "u_ref"
identical(saved_ref, u_ref)








