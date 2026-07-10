
# load in data
results_file <- "results12345_2.rds"
results <- readRDS(results_file)
summary_df <- do.call(rbind, results)




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

