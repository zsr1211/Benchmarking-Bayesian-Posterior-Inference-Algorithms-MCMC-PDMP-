# This scrip is about visualising benchmark results

# ============================================================
# Load packages
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)





# ============================================================
# Load benchmark result files
# ============================================================

# files_info2:
# RMSE is computed directly from the mean of each parameter in the PDMP output.
# W1, W2, and KL are computed from discretized PDMP traces with dt = 0.01.

files_info2 <- readRDS(
  "/Users/zhangsr/Desktop/UvA/Thesis/programming/Pipeline/files_info2.rds")

# files_info1:
# RMSE, W1, W2, and KL are computed from discretized PDMP traces with dt = 0.1.

files_info1 <- readRDS(
  "/Users/zhangsr/Desktop/UvA/Thesis/programming/Pipeline/files_info1.rds")


# Combine the two result files 
# Add a result_version column to keep track of which file each row comes from.
files_all <- bind_rows(
  files_info1 %>%
    mutate(result_version = "file1_dt0.1_discretized"),
  
  files_info2 %>%
    mutate(result_version = "file2_dt0.01_continuous_RMSE")
)






# ============================================================
# Prepare dat for plotting
# ============================================================

# Define metrics to plot:
# accuracy measures: RMSE, W1_relative, W2_relative, and KL.
# efficiency measures: ESS_per_GE and ESS_per_time.

metrics <- c(
  "RMSE",
  "W1_relative",
  "W2_relative",
  "KL",
  "ESS_per_GE",
  "ESS_per_time"
)


# Convert data to long format
# Each row in plot_long represents:
# one result_version × posterior × algorithm × seed × metric.
# Missing values are removed for each metric separately.

plot_long <- files_all %>%
  select(
    result_version,
    posterior,
    group,
    algorithm,
    seed,
    all_of(metrics)
  ) %>%
  pivot_longer(
    cols = all_of(metrics),
    names_to = "metric",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  mutate(
    group = factor(group),
    posterior_label = paste(group, posterior, sep = " / ")
  )





# Set posterior order on the x-axis
# The posterior labels are ordered by group first, and then by posterior name.

posterior_levels <- plot_long %>%
  distinct(group, posterior, posterior_label) %>%
  arrange(group, posterior) %>%
  pull(posterior_label)

plot_long <- plot_long %>%
  mutate(
    posterior_label = factor(posterior_label, levels = posterior_levels)
  )


# Summarise metric values across seeds
# For each result version, posterior, algorithm, and metric,
# I summarise the values across all seeds.
#
# The median is used as the main point estimate.
# The 25% and 75% quantiles are used as the interval.

plot_summary <- plot_long %>%
  group_by(result_version, group, posterior, posterior_label, algorithm, metric) %>%
  summarise(
    n_seed = n(),
    center_median = median(value, na.rm = TRUE),
    center_min = min(value, na.rm = TRUE),
    q25 = quantile(value, 0.25, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE),
    .groups = "drop"
  )





# ============================================================
# Function to plot one metric
# ============================================================

# This function makes one plot for one metric and one result version.
# In the plot:
# - points show the selected center value across seeds
# - error bars show the 25% to 75% interval across seeds
# - colors show different algorithms

plot_metric <- function(data,
                        metric_name,
                        result_version_name,
                        center = c("median", "min"),
                        log_y = FALSE) {
  
  center <- match.arg(center)
  
  center_col <- if (center == "median") {
    "center_median"
  } else {
    "center_min"
  }
  
  df <- data %>%
    filter(
      metric == metric_name,
      result_version == result_version_name
    )
  
  p <- ggplot(
    df,
    aes(
      x = posterior_label,
      y = .data[[center_col]],
      color = algorithm,
      group = algorithm
    )
  ) +
    geom_errorbar(
      aes(ymin = q25, ymax = q75),
      position = position_dodge(width = 0.6),
      width = 0.25,
      linewidth = 0.6
    ) +
    geom_point(
      position = position_dodge(width = 0.6),
      size = 2.5
    ) +
    labs(
      x = "Posterior (group / name)",
      y = metric_name,
      color = "Algorithm",
      title = paste0(metric_name, " by posterior and algorithm"),
      subtitle = paste0(
        result_version_name,
        " | Point = ", center,
        ", interval = 25%–75% across seeds"
      )
    ) +
    theme_bw(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
      legend.position = "top"
    )
  
  # Use a log10 y-axis when values are very spread out.
  # The axis labels still show the original metric values.
  if (log_y) {
    p <- p + scale_y_log10()
  }
  
  p
}


# ============================================================
# Plots for files_info1
# ============================================================

# files_info1 uses discretized PDMP traces with dt = 0.1.
# RMSE, W1, W2, and KL are computed from these discretized traces.

plot_RMSE_1 <- plot_metric(
  data = plot_summary,
  metric_name = "RMSE",
  result_version_name = "file1_dt0.1_discretized",
  center = "median",
  log_y = TRUE
)

plot_W1_1 <- plot_metric(
  data = plot_summary,
  metric_name = "W1_relative",
  result_version_name = "file1_dt0.1_discretized",
  center = "median",
  log_y = TRUE
)

plot_W2_1 <- plot_metric(
  data = plot_summary,
  metric_name = "W2_relative",
  result_version_name = "file1_dt0.1_discretized",
  center = "median",
  log_y = TRUE
)

plot_KL_1 <- plot_metric(
  data = plot_summary,
  metric_name = "KL",
  result_version_name = "file1_dt0.1_discretized",
  center = "median",
  log_y = TRUE
)

plot_ESS_GE_1 <- plot_metric(
  data = plot_summary,
  metric_name = "ESS_per_GE",
  result_version_name = "file1_dt0.1_discretized",
  center = "median",
  log_y = TRUE
)

plot_ESS_time_1 <- plot_metric(
  data = plot_summary,
  metric_name = "ESS_per_time",
  result_version_name = "file1_dt0.1_discretized",
  center = "median",
  log_y = TRUE
)





# ============================================================
# Plots for files_info2
# ============================================================

# files_info2 uses a different way to compute the metrics.
# RMSE is computed directly from the mean of each parameter in the PDMP output.
# W1, W2, and KL are computed from discretized PDMP traces with dt = 0.01.

plot_RMSE_2 <- plot_metric(
  data = plot_summary,
  metric_name = "RMSE",
  result_version_name = "file2_dt0.01_continuous_RMSE",
  center = "median",
  log_y = TRUE
)

plot_W1_2 <- plot_metric(
  data = plot_summary,
  metric_name = "W1_relative",
  result_version_name = "file2_dt0.01_continuous_RMSE",
  center = "median",
  log_y = TRUE
)

plot_W2_2 <- plot_metric(
  data = plot_summary,
  metric_name = "W2_relative",
  result_version_name = "file2_dt0.01_continuous_RMSE",
  center = "median",
  log_y = TRUE
)

plot_KL_2 <- plot_metric(
  data = plot_summary,
  metric_name = "KL",
  result_version_name = "file2_dt0.01_continuous_RMSE",
  center = "median",
  log_y = TRUE
)

plot_ESS_GE_2 <- plot_metric(
  data = plot_summary,
  metric_name = "ESS_per_GE",
  result_version_name = "file2_dt0.01_continuous_RMSE",
  center = "median",
  log_y = TRUE
)

plot_ESS_time_2 <- plot_metric(
  data = plot_summary,
  metric_name = "ESS_per_time",
  result_version_name = "file2_dt0.01_continuous_RMSE",
  center = "median",
  log_y = TRUE
)





# ============================================================
# Check plots
# ============================================================

# RMSE
library(patchwork)
plot_RMSE_1 + plot_RMSE_2
# RMSE does not change much between the two PDMP mean calculation methods.
# In general, complex posteriors have larger RMSE.
# The difference between moderate and simple posteriors is not very clear.
# Compared with the other algorithms, BPS tends to have larger RMSE.
# MCMC tends to have smaller RMSE.





# W1 and W2
plot_W1_1 + plot_W1_2

plot_W2_1 + plot_W2_2

# Notes:
# W1 and W2 change very little when using different dt values
# for PDMP discretization.
#
# W1 and W2 are scaled using the MCMC result as the reference.
# This makes the values more comparable across different posteriors.
# Because of this scaling, the MCMC result is around 1.
#
# In general, complex posteriors have larger W1 and W2.
# The difference between moderate and simple posteriors is not very clear.
# Compared with the other algorithms, BPS tends to have larger W1 and W2.





# Check KL plots

plot_KL_1 + plot_KL_2

# Notes:
# KL changes more than W1 and W2 when using different dt values,
# but the change is still not very large.
#
# For the posteriors "seeds_data-seeds_centered_model"
# and "surgical_data-surgical_model", the interval of BPS is relatively large.
# This suggests that BPS is less stable for these posteriors.





# Check ESS efficiency plots

plot_ESS_GE_1

plot_ESS_time_1

# Notes:
# In general, MCMC has better efficiency values,
# because it often has larger ESS/GE or larger ESS/time.
#
# For some posteriors, the intervals are large.
# This means the efficiency can vary a lot across seeds.




