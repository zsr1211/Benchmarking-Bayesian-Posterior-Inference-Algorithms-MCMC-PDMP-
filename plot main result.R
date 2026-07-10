# this script is for visulize the main result of the benchmark

library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(legendry)
library(knitr)
library(kableExtra)

# --------------------------------------------------
# 0. Read raw data
# --------------------------------------------------

files_info1_raw <- readRDS(
  "/Users/zhangsr/Desktop/UvA/Thesis/programming/Pipeline/files_info1.rds"
)

# Decide the order ONCE
group_order <- c("complex", "moderate", "simple")
# If you prefer simple first, use:
# group_order <- c("simple", "moderate", "complex")


# --------------------------------------------------
# 1. Create ONE lookup table
# --------------------------------------------------

posterior_lookup <- files_info1_raw %>%
  mutate(
    data = str_extract(posterior, "^[^-]+"),
    model = str_remove(posterior, "^[^-]+-"),
    has_reference = if_else(has_ref, "Yes", "No"),
    group = factor(group, levels = group_order)
  ) %>%
  distinct(posterior, data, model, group, has_reference, has_ref) %>%
  arrange(group, desc(has_ref), posterior) %>%
  mutate(
    model_abbr = paste0("M", row_number())
  )

# Check the lookup
posterior_lookup %>%
  select(model_abbr, posterior, data, model, group, has_reference)

files_info1 <- files_info1_raw %>%
  left_join(
    posterior_lookup %>%
      select(posterior, data, model, model_abbr),
    by = "posterior"
  ) %>%
  mutate(
    group = factor(group, levels = group_order),
    model_abbr = factor(model_abbr, levels = posterior_lookup$model_abbr)
  )

posterior_model_table <- posterior_lookup %>%
  select(
    group,
    model_abbr,
    data,
    model,
    has_reference
  )

posterior_model_table_latex <- posterior_model_table %>%
  rename(
    Group = group,
    Abbreviation = model_abbr,
    Data = data,
    Model = model,
    `Reference draws` = has_reference
  ) %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    caption = "Overview of posterior models included in the benchmark.",
    label = "tab:posterior-models",
    escape = TRUE
  ) %>%
  kable_styling(
    latex_options = c("repeat_header"),
    font_size = 8
  ) %>%
  column_spec(1, width = "1.8cm") %>%
  column_spec(2, width = "2.0cm") %>%
  column_spec(3, width = "3.0cm") %>%
  column_spec(4, width = "5.0cm") %>%
  column_spec(5, width = "2.0cm") %>%
  collapse_rows(
    columns = 1,
    valign = "top"
  )

writeLines(
  posterior_model_table_latex,
  "posterior_model_table.tex"
)

metrics <- c(
  "RMSE",
  "W1_relative",
  "W2_relative",
  "KL",
  "ESS_per_GE",
  "ESS_per_time"
)

plot_long <- files_info1 %>%
  select(
    posterior,
    model_abbr,
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
  filter(!is.na(value))

plot_summary <- plot_long %>%
  group_by(group, posterior, model_abbr, algorithm, metric) %>%
  summarise(
    n_seed = n(),
    center_median = median(value, na.rm = TRUE),
    center_min = min(value, na.rm = TRUE),
    q25 = quantile(value, 0.25, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

posterior_lookup_ref <- posterior_lookup %>%
  filter(has_ref)

x_levels_ref <- posterior_lookup_ref$model_abbr

posterior_lookup_ref <- posterior_lookup_ref %>%
  mutate(
    x_pos = row_number()
  )

plot_summary_ref <- plot_summary %>%
  filter(model_abbr %in% x_levels_ref) %>%
  left_join(
    posterior_lookup_ref %>%
      select(model_abbr, x_pos),
    by = "model_abbr"
  )

group_ranges_ref <- posterior_lookup_ref %>%
  group_by(group) %>%
  summarise(
    start = min(x_pos),
    end = max(x_pos),
    name = as.character(first(group)),
    .groups = "drop"
  )

plot_metric <- function(data,
                        metric_name,
                        lookup,
                        group_ranges,
                        center = c("median", "min"),
                        log_y = FALSE) {
  
  center <- match.arg(center)
  
  center_col <- if (center == "median") {
    "center_median"
  } else {
    "center_min"
  }
  
  df <- data %>%
    filter(metric == metric_name)
  
  p <- ggplot(
    df,
    aes(
      x = x_pos,
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
    scale_x_continuous(
      breaks = lookup$x_pos,
      labels = lookup$model_abbr,
      guide = guide_axis_nested(
        key = key_range_manual(
          start = group_ranges$start,
          end = group_ranges$end,
          name = group_ranges$name
        ),
        regular_key = key_manual(
          lookup$x_pos,
          label = lookup$model_abbr
        ),
        type = "bracket",
        angle = 0,
        pad_discrete = 0
      )
    ) +
    labs(
      x = NULL,
      y = metric_name,
      color = "Algorithm",
      title = paste0(metric_name, " by posterior and algorithm"),
      subtitle = paste0(
        "Point = ", center,
        ", interval = 25%–75% across seeds"
      )
    ) +
    theme_bw(base_size = 16) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 13),
      axis.title.y = element_text(size = 15),
      plot.title = element_text(size = 17),
      plot.subtitle = element_text(size = 13),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13)
    )
  
  if (log_y) {
    p <- p + scale_y_log10()
  }
  
  p
}

plot_RMSE <- plot_metric(
  data = plot_summary_ref,
  metric_name = "RMSE",
  lookup = posterior_lookup_ref,
  group_ranges = group_ranges_ref,
  center = "median",
  log_y = TRUE
)

plot_RMSE


plot_W1 <- plot_metric(
  data = plot_summary_ref,
  metric_name = "W1_relative",
  lookup = posterior_lookup_ref,
  group_ranges = group_ranges_ref,
  center = "median",
  log_y = TRUE
)

plot_W1


plot_W2 <- plot_metric(
  data = plot_summary_ref,
  metric_name = "W2_relative",
  lookup = posterior_lookup_ref,
  group_ranges = group_ranges_ref,
  center = "median",
  log_y = TRUE
)

plot_W2



# make_nested_axis_data <- function(lookup, plot_summary_data) {
#   
#   lookup_sub <- lookup %>%
#     arrange(group, model_abbr) %>%
#     mutate(
#       x_pos = row_number()
#     )
#   
#   plot_data_sub <- plot_summary_data %>%
#     filter(model_abbr %in% lookup_sub$model_abbr) %>%
#     mutate(
#       model_abbr = as.character(model_abbr)
#     ) %>%
#     left_join(
#       lookup_sub %>%
#         select(model_abbr, x_pos),
#       by = "model_abbr"
#     )
#   
#   group_ranges_sub <- lookup_sub %>%
#     group_by(group) %>%
#     summarise(
#       start = min(x_pos),
#       end = max(x_pos),
#       name = as.character(first(group)),
#       .groups = "drop"
#     )
#   
#   list(
#     lookup = lookup_sub,
#     plot_data = plot_data_sub,
#     group_ranges = group_ranges_sub
#   )
# }

make_nested_axis_data <- function(lookup, plot_summary_data) {
  
  lookup_sub <- lookup %>%
    mutate(
      model_abbr = as.character(model_abbr),
      model_num = as.integer(sub("^M", "", model_abbr))
    ) %>%
    arrange(group, model_num) %>%
    mutate(
      x_pos = row_number()
    )
  
  plot_data_sub <- plot_summary_data %>%
    mutate(
      model_abbr = as.character(model_abbr)
    ) %>%
    filter(model_abbr %in% lookup_sub$model_abbr) %>%
    left_join(
      lookup_sub %>%
        select(model_abbr, x_pos),
      by = "model_abbr"
    )
  
  group_ranges_sub <- lookup_sub %>%
    group_by(group) %>%
    summarise(
      start = min(x_pos),
      end = max(x_pos),
      name = as.character(first(group)),
      .groups = "drop"
    )
  
  list(
    lookup = lookup_sub,
    plot_data = plot_data_sub,
    group_ranges = group_ranges_sub
  )
}

posterior_lookup_noref <- posterior_lookup %>%
  filter(!has_ref)

nested_noref <- make_nested_axis_data(
  lookup = posterior_lookup_noref,
  plot_summary_data = plot_summary
)

plot_KL_noref <- plot_metric(
  data = nested_noref$plot_data,
  metric_name = "KL",
  lookup = nested_noref$lookup,
  group_ranges = nested_noref$group_ranges,
  center = "median",
  log_y = TRUE
)

plot_KL_noref


nested_all <- make_nested_axis_data(
  lookup = posterior_lookup,
  plot_summary_data = plot_summary
)

plot_ESS_GE_all <- plot_metric(
  data = nested_all$plot_data,
  metric_name = "ESS_per_GE",
  lookup = nested_all$lookup,
  group_ranges = nested_all$group_ranges,
  center = "median",
  log_y = TRUE
)

plot_ESS_GE_all


plot_ESS_time_all <- plot_metric(
  data = nested_all$plot_data,
  metric_name = "ESS_per_time",
  lookup = nested_all$lookup,
  group_ranges = nested_all$group_ranges,
  center = "median",
  log_y = TRUE
)

plot_ESS_time_all




plot1 <- plot_RMSE + plot_KL
plot2 <- plot_W1 + plot_W2
accuracy_combined <- (plot1 / plot2) &
  theme(
    legend.position = "top",
    text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(size = 11, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(size = 17),
    plot.subtitle = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14)
  )
accuracy_combined

ggsave(
  filename = "accuracy_combined.png",
  plot = accuracy_combined,
  width = 16,
  height = 10,
  dpi = 600
)



plot <- plot_ESS_GE_all + plot_ESS_time_all
ggsave(
  filename = "efficie_combined.png",
  plot = plot,
  width = 16,
  height = 5,
  dpi = 600
)

