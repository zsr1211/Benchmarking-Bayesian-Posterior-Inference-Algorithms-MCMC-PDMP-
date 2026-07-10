# This script runs the exploratory random forest analysis.
#
# The input table was created by:
# compute_reference_geometry_and_performance.R
#
# The goal is to explore which posterior geometry features are 
# associated with poor performance of specific sampler.
#
# This is an exploratory analysis, not a predictive model for new posteriors.





# ============================================================
# Load in packages and data
# ============================================================

# Load packages
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(ranger)
library(ggplot2)


# Directory with saved data
output_dir <- "/Users/zhangsr/Desktop/UvA/Thesis/programming/Pipeline"


# Load posterior-level performance and geometry table
posterior_performance <- readRDS(
  file.path(output_dir, "posterior_algorithm_performance_with_geometry.rds")
)





# ============================================================
# Prepare data for the random forest analysis
# ============================================================

# The table has one row for each posterior and algorithm.
#
# Outcomes:
# bad_efficiency: -log(median ESS per gradient evaluation): 
# Larger values mean worse sampling efficiency.

# bad_rmse: log(median RMSE)
# Larger values mean larger posterior mean error.
#
# Geometry features:
#
# dim:
# Number of unconstrained parameters.
# This is used as a control variable.
#
# tail_ratio:
# Mean ratio of the 1%-99% quantile range to the IQR.
# Larger values mean heavier tails or more extreme marginal spread.
#
# mean_abs_correlation:
# Mean absolute off-diagonal posterior correlation.
# Larger values mean stronger average posterior dependence.
#
# log_condition_number:
# Log condition number of the posterior covariance matrix.
# Larger values mean stronger posterior anisotropy.

rf_data <- posterior_performance %>%
  dplyr::select(
    posterior,
    group,
    algorithm,
    bad_efficiency,
    bad_rmse,
    dim,
    mean_abs_skewness,
    tail_ratio,
    mean_abs_correlation,
    log_condition_number
  ) %>%
  dplyr::filter(
    is.finite(bad_efficiency),
    is.finite(bad_rmse),
    is.finite(dim),
    is.finite(mean_abs_skewness),
    is.finite(tail_ratio),
    is.finite(mean_abs_correlation),
    is.finite(log_condition_number)
  ) %>%
  dplyr::mutate(
    algorithm = as.factor(algorithm)
  )





# ============================================================
# Check correlations between geometry features
# ============================================================

# This check is done at posterior level.

geometry_correlation <- rf_data %>%
  dplyr::distinct(
    posterior,
    dim,
    log_condition_number,
    mean_abs_skewness,
    tail_ratio,
    mean_abs_correlation
  ) %>%
  dplyr::select(
    dim,
    log_condition_number,
    mean_abs_skewness,
    tail_ratio,
    mean_abs_correlation
  ) %>%
  stats::cor(use = "pairwise.complete.obs")


geometry_correlation

# mean_abs_skewness is not included because it is almost perfectly
# correlated with tail_ratio.





# ============================================================
# Fit random forest models
# ============================================================

# Two random forest models are fitted:
#
# 1. bad_efficiency model
# 2. bad_rmse model
#
# The models include algorithm as a predictor 
# and include dim as a control variable.

set.seed(1234)

rf_efficiency <- ranger::ranger(
  bad_efficiency ~ algorithm + dim + tail_ratio +
    mean_abs_correlation + log_condition_number,
  data = rf_data,
  num.trees = 1000,
  importance = "permutation",
  respect.unordered.factors = "order",
  seed = 1234
)


rf_rmse <- ranger::ranger(
  bad_rmse ~ algorithm + dim + tail_ratio +
    mean_abs_correlation + log_condition_number,
  data = rf_data,
  num.trees = 1000,
  importance = "permutation",
  respect.unordered.factors = "order",
  seed = 1234
)





# ============================================================
# Overall permutation importance
# ============================================================

# These importance values are overall model-level importance.

# importance_efficiency <- tibble::tibble(
#   feature = names(rf_efficiency$variable.importance),
#   importance = as.numeric(rf_efficiency$variable.importance),
#   outcome = "bad_efficiency"
# ) %>%
#   dplyr::arrange(dplyr::desc(importance))
# 
# 
# importance_rmse <- tibble::tibble(
#   feature = names(rf_rmse$variable.importance),
#   importance = as.numeric(rf_rmse$variable.importance),
#   outcome = "bad_rmse"
# ) %>%
#   dplyr::arrange(dplyr::desc(importance))
# 
# 
# importance_efficiency
# importance_rmse





# ============================================================
# Algorithm-specific permutation importance
# ============================================================

# This function computes permutation importance separately within each algorithm.
#
# The random forest is still a joint model trained on all algorithms.
# Then, for each algorithm, each geometry feature is permuted only
# within rows for that algorithm.
#
# The result shows which geometry feature is most important for
# explaining poor performance for each sampler.
#
# Importance is normalized to sum to 100% within each algorithm.

compute_algorithm_specific_importance <- function(model,
                                                  data,
                                                  outcome_name,
                                                  features,
                                                  algorithms) {
  
  result_list <- list()
  
  for (alg in algorithms) {
    
    data_alg <- data %>%
      dplyr::filter(algorithm == alg)
    
    y_true <- data_alg[[outcome_name]]
    pred_base <- predict(model, data = data_alg)$predictions
    
    baseline_mse <- mean((y_true - pred_base)^2, na.rm = TRUE)
    
    for (feat in features) {
      
      data_perm <- data_alg
      data_perm[[feat]] <- sample(data_perm[[feat]])
      
      pred_perm <- predict(model, data = data_perm)$predictions
      
      perm_mse <- mean((y_true - pred_perm)^2, na.rm = TRUE)
      
      raw_importance <- perm_mse - baseline_mse
      
      result_list[[length(result_list) + 1]] <- tibble::tibble(
        algorithm = alg,
        feature = feat,
        raw_importance = raw_importance
      )
    }
  }
  
  dplyr::bind_rows(result_list) %>%
    dplyr::group_by(algorithm) %>%
    dplyr::mutate(
      raw_importance = pmax(raw_importance, 0),
      importance_percent = 100 * raw_importance / sum(raw_importance, na.rm = TRUE)
    ) %>%
    dplyr::ungroup()
}





# ============================================================
# Compute algorithm-specific importance
# ============================================================

# These are the geometry features shown in the heatmap.
# dim is included in the random forest as a control variable,
# but it is not included in the geometry importance heatmap.

geometry_features <- c(
  "tail_ratio",
  "mean_abs_correlation",
  "log_condition_number"
)


algorithms <- c("MCMC", "ZigZag", "BPS", "Boomerang")


# Importance for poor efficiency
importance_efficiency_by_alg <- compute_algorithm_specific_importance(
  model = rf_efficiency,
  data = rf_data,
  outcome_name = "bad_efficiency",
  features = geometry_features,
  algorithms = algorithms
) %>%
  dplyr::mutate(outcome = "bad_efficiency")


# Importance for RMSE
importance_rmse_by_alg <- compute_algorithm_specific_importance(
  model = rf_rmse,
  data = rf_data,
  outcome_name = "bad_rmse",
  features = geometry_features,
  algorithms = algorithms
) %>%
  dplyr::mutate(outcome = "bad_rmse")


# Combine both outcomes
importance_by_alg <- dplyr::bind_rows(
  importance_efficiency_by_alg,
  importance_rmse_by_alg
)





# ============================================================
# Plot algorithm-specific importance heatmap
# ============================================================

# Labels used in the plot
feature_labels <- c(
  tail_ratio = "Tail heaviness",
  mean_abs_correlation = "Correlation",
  log_condition_number = "Anisotropy"
)


outcome_labels <- c(
  bad_efficiency = "Poor efficiency: -log(ESS/GE)",
  bad_rmse = "High RMSE: log(RMSE)"
)


importance_by_alg_plot <- importance_by_alg %>%
  dplyr::mutate(
    feature_label = feature_labels[feature],
    outcome_label = outcome_labels[outcome],
    algorithm = factor(
      algorithm,
      levels = c("MCMC", "ZigZag", "BPS", "Boomerang")
    ),
    feature_label = factor(
      feature_label,
      levels = c("Tail heaviness", "Correlation", "Anisotropy")
    )
  )


importance_heatmap <- ggplot(
  importance_by_alg_plot,
  aes(
    x = algorithm,
    y = feature_label,
    fill = importance_percent
  )
) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = sprintf("%.1f", importance_percent)),
    size = 3
  ) +
  facet_wrap(~ outcome_label) +
  scale_fill_gradientn(
    colours = c("#F7FBFF", "#9ECAE1", "#3182BD"),
    limits = c(0, 100),
    name = "Importance (%)"
  ) +
  labs(
    x = "Sampler",
    y = "Posterior geometry feature"
  ) +
  theme_minimal()
importance_heatmap





# ============================================================
# Direction plots for log condition number
# ============================================================

# The heatmap shows that log condition number is the dominant feature. 
# These plots show the direction of the relationship.

condition_efficiency_plot <- ggplot2::ggplot(
  rf_data,
  ggplot2::aes(
    x = log_condition_number,
    y = bad_efficiency
  )
) +
  ggplot2::geom_point(alpha = 0.6) +
  ggplot2::geom_smooth(method = "loess", se = FALSE) +
  ggplot2::facet_wrap(~ algorithm) +
  ggplot2::labs(
    x = "Log condition number",
    y = "-log(ESS/GE)"
  ) +
  ggplot2::theme_minimal()
condition_efficiency_plot


condition_rmse_plot <- ggplot2::ggplot(
  rf_data,
  ggplot2::aes(
    x = log_condition_number,
    y = bad_rmse
  )
) +
  ggplot2::geom_point(alpha = 0.6) +
  ggplot2::geom_smooth(method = "loess", se = FALSE) +
  ggplot2::facet_wrap(~ algorithm) +
  ggplot2::labs(
    x = "Log condition number",
    y = "log(RMSE)"
  ) +
  ggplot2::theme_minimal()
condition_rmse_plot

# The smooth line goes upward, 
# indicating that larger log condition number is associated with worse performance.





# ============================================================
# Check ranges of geometry features
# ============================================================

# This gives a simple summary of the feature ranges.

geometry_feature_ranges <- rf_data %>%
  dplyr::distinct(
    posterior,
    tail_ratio,
    mean_abs_correlation,
    log_condition_number
  ) %>%
  dplyr::summarise(
    tail_min = min(tail_ratio),
    tail_max = max(tail_ratio),
    corr_min = min(mean_abs_correlation),
    corr_max = max(mean_abs_correlation),
    cond_min = min(log_condition_number),
    cond_max = max(log_condition_number)
  )
geometry_feature_ranges

# condition number has a very wide range.
