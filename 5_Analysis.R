rm(list = ls())
graphics.off()

{
  library(ggplot2)
  library(tidyverse)
  library(effectsize) # rank_biserial()
  library(overlapping) # overlap(), boot.overlap()
  library(DescTools)
  library(cowplot)
  library(ggh4x)
  library(car)
}

# READ IN THE DATA, COMBINE, AND EDIT: ####
# ==================================== #

linear <- read.csv("./Data/Final_data/Linear_truth_comparison_metrics.csv")
non <- read.csv("./Data/Final_data/Non_linear_truth_comparison_metrics.csv")
semi <- read.csv(
  "./Data/Final_data/Semi_non_linear_truth_comparison_metrics.csv"
)

# Combine the three data frames into one and edit the sim column to be more descriptive:
metrics <- rbind(linear, semi, non) |>
  dplyr::mutate(
    sim_mod_id = paste0(sim, "_", mod_id),
    sim = dplyr::case_when(
      sim == "L" ~ "Linear",
      sim == "SNL" ~ "Semi-Linear",
      sim == "NL" ~ "Non-Linear"
    ),
    sim = factor(sim, levels = c("Linear", "Semi-Linear", "Non-Linear"))
  )

rm(linear, semi, non)


# Create a list of data frames, one for each metric, in wide format (LR and RF values side by side)
# ------------------------------------------------------------------------------------------------- #

# Define the metrics and simulation levels to loop over
metric_names <- c("rmse", "bias", "cor") # add "auc", "tss" here once computed
sim_levels <- c("Linear", "Semi-Linear", "Non-Linear") # ORDER MATTERS for the trend test

# Create a list to hold the wide-format data frames
wide_list <- list()

# Loop over each metric and create a wide-format data frame for it
for (m in metric_names) {
  lr_sub <- metrics[
    metrics$mod_type == "Logistic Regression",
    c("sim", "sim_mod_id", m)
  ]
  rf_sub <- metrics[metrics$mod_type == "Random Forest", c("sim_mod_id", m)]

  names(lr_sub)[3] <- "lr_val"
  names(rf_sub)[2] <- "rf_val"

  merged <- merge(lr_sub, rf_sub, by = "sim_mod_id")
  merged$diff <- merged$rf_val - merged$lr_val # positive = RF had the higher value

  wide_list[[m]] <- merged
}

rm(lr_sub, rf_sub, merged, m)


# ======================================================================================================== #
# ======================================================================================================== #

# ============= #
# OVERLAP TEST: ####
# ============= #

results_overlap <- data.frame(
  metric = character(),
  sim = character(),
  overlap_estimate = numeric(),
  overlap_ci_low = numeric(),
  overlap_ci_high = numeric(),
  stringsAsFactors = FALSE
)

for (m in metric_names) {
  cat("Calculating overlap for metric: ", m, "\n")
  dat <- wide_list[[m]]

  for (s in sim_levels) {
    cat("  Simulation: ", s, "\n")

    sub <- dat[dat$sim == s, ]
    ov_input <- list(LR = sub$lr_val, RF = sub$rf_val)

    ov_est <- overlap(ov_input)$OV
    boot_out <- boot.overlap(ov_input, B = 1000)
    ci <- quantile(as.vector(boot_out$OVboot_dist), probs = c(0.025, 0.975))

    results_overlap <- rbind(
      results_overlap,
      data.frame(
        metric = m,
        sim = s,
        overlap_estimate = ov_est,
        overlap_ci_low = unname(ci[1]),
        overlap_ci_high = unname(ci[2])
      )
    )
  }
}

print(results_overlap)

rm(m, dat, s, sub, ov_input, ov_est, boot_out, ci)

# ------------------------------------------------------------------------------------------------------- #

# Graphing the overlap results:

mod_names <- c("Logistic Regression", "Random Forest")

for (n in mod_names) {
  df <- metrics[metrics$mod_type == n, ]

  for (m in metric_names) {
    plot <- ggplot(df, aes(x = .data[[m]], color = sim)) +
      geom_density(linewidth = 1.2) +
      geom_rug(sides = "b", linewidth = 0.5, alpha = 0.7) +
      scale_color_manual(
        name = "Simulation",
        values = c(
          "Linear" = "#6a040f",
          "Semi-Linear" = "#dc2f02",
          "Non-Linear" = "#faa307"
        )
      ) +
      labs(x = m, y = "Density") +
      theme_classic() +
      theme(
        # legend.position = c(0.2, 0.2),
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10, face = "bold"),
        panel.grid = element_blank(),
        axis.line = element_line(linewidth = 1),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10, col = "black", colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black")
      )

    ggsave(
      paste0("./Data/Figures/Density_Plots/", n, "_", m, "_density_plot.png"),
      plot,
      width = 16,
      height = 12,
      units = "cm",
      dpi = 300
    )
  }
}


rf_data <- metrics[metrics$mod_type == "Logistic Regression", ]

ggplot(rf_data, aes(x = rmse, color = sim)) +
  geom_density(linewidth = 1.2) +
  geom_rug(sides = "b", linewidth = 0.5, alpha = 0.7) +
  scale_color_manual(
    name = "Simulation",
    values = c(
      "Linear" = "#6a040f",
      # "Random Forest" = "#6a040f",
      "Semi-Linear" = "#dc2f02",
      "Non-Linear" = "#faa307"
      # "Logistic Regression" = "#faa307"
    )
  ) +
  labs(x = "RMSE", y = "Density") +
  theme_classic()

rf_data <- metrics[metrics$sim == "Non-Linear", ]

ggplot(rf_data, aes(x = rmse, color = mod_type)) +
  geom_density(linewidth = 1.2, key_glyph = "path") +
  geom_rug(sides = "b", linewidth = 0.5, alpha = 0.7, show.legend = FALSE) +
  scale_color_manual(
    name = "Model",
    values = c(
      # "Linear" = "#6a040f",
      "Random Forest" = "#6a040f",
      # "Semi-Linear" = "#dc2f02",
      # "Non-Linear" = "#faa307"
      "Logistic Regression" = "#faa307"
    )
  ) +
  labs(x = "RMSE", y = "Density") +
  theme_classic()

# ======================================================================================================== #
# ======================================================================================================== #

for (i in metric_names) {
  if (i == "rmse") {
    x_label <- "RMSE"
  } else if (i == "bias") {
    x_label <- "Bias"
  } else if (i == "cor") {
    x_label <- "Correlation"
  } else {
    x_label <- i
  }

  plot <- ggplot(metrics, aes(x = .data[[i]], colour = mod_type)) +
    geom_density(linewidth = 1.2, key_glyph = "path") +
    # geom_rug(sides = "b", linewidth = 0.5, alpha = 0.7, show.legend = FALSE) +
    scale_color_manual(
      name = "Model",
      values = c(
        "Random Forest" = "#1B3A5C",
        "Logistic Regression" = "#3f8d5d"
      )
    ) +
    labs(x = x_label, y = "Density") +
    facet_wrap(~sim, nrow = 2) +
    theme_classic() +
    theme(
      legend.position = c(0.8, 0.2),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10, face = "bold"),
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 0.6),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 10, col = "black", colour = "black"),
      axis.text.x = element_text(size = 12, colour = "black")
    )

  ggsave(
    paste0("./Data/Figures/Density_Plots/", i, "_density_plot.png"),
    plot,
    width = 16,
    height = 12,
    units = "cm",
    dpi = 300
  )
}

test <- metrics |>
  dplyr::select(-mae) |>
  pivot_longer(
    cols = c(rmse, bias, cor),
    names_to = "metric",
    values_to = "value"
  ) |>
  mutate(
    metric = case_when(
      metric == "rmse" ~ "RMSE",
      metric == "bias" ~ "Bias",
      metric == "cor" ~ "Correlation"
    ),
    metric = factor(metric, levels = c("RMSE", "Bias", "Correlation")),
    sim = factor(sim, levels = c("Linear", "Semi-Linear", "Non-Linear"))
  ) |>
  ggplot(aes(x = value, colour = mod_type)) +
  geom_density(linewidth = 1.2, key_glyph = "path") +
  # geom_rug(sides = "b", linewidth = 0.5, alpha = 0.7, show.legend = FALSE) +
  scale_color_manual(
    name = "Model",
    values = c(
      "Random Forest" = "#1B3A5C",
      "Logistic Regression" = "#3f8d5d"
    )
  ) +
  labs(x = NULL, y = "Density") +
  # facet_wrap(metric ~ sim, scales = "free") +
  facet_grid2(
    metric ~ sim,
    scales = "free",
    independent = "all",
    axes = "all"
  ) +
  theme_classic() +
  theme(
    # legend.position = c(0.8, 0.2),
    legend.position = "none",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.6),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10, col = "black", colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black"),
    strip.text = element_text(size = 11, face = "bold", colour = "black")
  )
test


# ======================================================================================================== #
# ======================================================================================================== #

# ======================================================================================================== #
# ======================================================================================================== #

# ======================================================================================================== #
# ======================================================================================================== #

# ===================== #
# PAIRED WILCOXON TEST: ####
# ===================== #

# Create a data frame to hold the results of the paired Wilcoxon tests and rank-biserial correlations
results_within <- data.frame(
  metric = character(),
  sim = character(),
  n_pairs = integer(),
  wilcoxon_p = numeric(),
  hodges_lehmann = numeric(), # paired median difference (RF - LR), in original units
  rank_biserial_r = numeric(),
  rb_ci_low = numeric(),
  rb_ci_high = numeric(),
  stringsAsFactors = FALSE
)

for (m in metric_names) {
  df <- wide_list[[m]]

  for (s in sim_levels) {
    sub <- df[df$sim == s, ]

    wt <- wilcox.test(sub$rf_val, sub$lr_val, paired = TRUE, conf.int = TRUE)
    rb <- rank_biserial(sub$rf_val, sub$lr_val, paired = TRUE)

    results_within <- rbind(
      results_within,
      data.frame(
        metric = m,
        sim = s,
        n_pairs = nrow(sub),
        wilcoxon_p = wt$p.value,
        hodges_lehmann = unname(wt$estimate),
        rank_biserial_r = rb$r_rank_biserial,
        rb_ci_low = rb$CI_low,
        rb_ci_high = rb$CI_high
      )
    )
  }
}

# Benjamini-Hochberg correction across the 15 tests (per metric, or across all —
# adjust the grouping below to whichever scope you decide on)
results_within$wilcoxon_p_adj <- p.adjust(
  results_within$wilcoxon_p,
  method = "BH"
)

print(results_within)

rm(m, df, s, sub, wt, rb)

# ======================================================================================================== #
# ======================================================================================================== #

# =============================== #
# JONCKHEERE-TERPSTRA TREND TEST: ####
# =============================== #

results_trend <- data.frame(
  metric = character(),
  jt_statistic = numeric(),
  jt_p = numeric(),
  stringsAsFactors = FALSE
)

for (m in metric_names) {
  cat("Performing Jonckheere-Terpstra test for metric: ", m, "\n")

  dat <- wide_list[[m]]
  dat$sim <- factor(dat$sim, levels = sim_levels, ordered = TRUE)

  jt <- JonckheereTerpstraTest(
    dat$diff,
    dat$sim,
    alternative = "two.sided",
    nperm = 10000
  )

  results_trend <- rbind(
    results_trend,
    data.frame(
      metric = m,
      jt_statistic = unname(jt$statistic),
      jt_p = jt$p.value
    )
  )
}

print(results_trend)

results_trend_by_model <- data.frame(
  metric = character(),
  model = character(),
  jt_statistic = numeric(),
  jt_p = numeric(),
  stringsAsFactors = FALSE
)

for (m in metric_names) {
  dat <- wide_list[[m]]
  dat$sim <- factor(dat$sim, levels = sim_levels, ordered = TRUE)

  for (mod in c("rf_val", "lr_val")) {
    cat(sprintf(
      "JT trend test for %s (%s), across Linear < Semi-Linear < Non-Linear\n",
      mod,
      m
    ))

    jt <- JonckheereTerpstraTest(
      dat[[mod]],
      dat$sim,
      alternative = "two.sided",
      nperm = 10000
    )

    results_trend_by_model <- rbind(
      results_trend_by_model,
      data.frame(
        metric = m,
        model = mod,
        jt_statistic = unname(jt$statistic),
        jt_p = jt$p.value
      )
    )
  }
}

print(results_trend_by_model)

results_tau <- data.frame(
  metric = character(),
  model = character(),
  tau = numeric(),
  stringsAsFactors = FALSE
)

for (m in metric_names) {
  dat <- wide_list[[m]]
  dat$sim <- factor(dat$sim, levels = sim_levels, ordered = TRUE)

  for (mod in c("rf_val", "lr_val")) {
    kt <- cor.test(dat[[mod]], as.numeric(dat$sim), method = "kendall")

    results_tau <- rbind(
      results_tau,
      data.frame(metric = m, model = mod, tau = unname(kt$estimate))
    )
  }
}

print(results_tau)


# ======================================================================================================== #
# ======================================================================================================== #

# ============= #
# OVERLAP TEST: ####
# ============= #

results_overlap <- data.frame(
  metric = character(),
  sim = character(),
  overlap_estimate = numeric(),
  overlap_ci_low = numeric(),
  overlap_ci_high = numeric(),
  stringsAsFactors = FALSE
)

for (m in metric_names) {
  dat <- wide_list[[m]]

  for (s in sim_levels) {
    sub <- dat[dat$sim == s, ]
    ov_input <- list(LR = sub$lr_val, RF = sub$rf_val)

    ov_est <- overlap(ov_input)$OV
    boot_out <- boot.overlap(ov_input, B = 1000)
    ci <- quantile(as.vector(boot_out$OVboot_dist), probs = c(0.025, 0.975))

    results_overlap <- rbind(
      results_overlap,
      data.frame(
        metric = m,
        sim = s,
        overlap_estimate = ov_est,
        overlap_ci_low = unname(ci[1]),
        overlap_ci_high = unname(ci[2])
      )
    )
  }
}

print(results_overlap)

rm(m, dat, s, sub, ov_input, ov_est, boot_out, ci)

# ======================================================================================================== #
# ======================================================================================================== #

# BOXPLOTS: ####
# ========= #

# RMSE plot:
# ---------- #

rmse_bp <- ggplot(metrics, aes(x = sim, y = rmse, fill = mod_type)) +
  geom_boxplot() +
  labs(
    x = NULL,
    y = "RMSE",
    fill = "Model Type"
  ) +
  scale_fill_manual(
    values = c("Logistic Regression" = "#535E31", "Random Forest" = "#A87B38")
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.2, 0.82),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 1),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10, col = "black", colour = "black"),
    axis.text.x = element_text(size = 12, face = "bold", colour = "black")
  )

rmse_bp

ggsave(
  "./Data/Figures/Boxplots/RMSE_boxplot.png",
  rmse_bp,
  width = 16,
  height = 12,
  units = "cm",
  dpi = 300
)


# Bias plot:
# ---------- #

bias_bp <- ggplot(metrics, aes(x = sim, y = bias, fill = mod_type)) +
  geom_boxplot() +
  labs(
    x = NULL,
    y = "Bias",
    fill = "Model Type"
  ) +
  scale_fill_manual(
    values = c("Logistic Regression" = "#535E31", "Random Forest" = "#A87B38")
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.83, 0.82),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 1),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10, col = "black", colour = "black"),
    axis.text.x = element_text(size = 12, face = "bold", colour = "black")
  )

bias_bp

ggsave(
  "./Data/Figures/Boxplots/Bias_boxplot.png",
  bias_bp,
  width = 16,
  height = 12,
  units = "cm",
  dpi = 300
)

# Correlation plot:
# ----------------- #

cor_bp <- ggplot(metrics, aes(x = sim, y = cor, fill = mod_type)) +
  geom_boxplot() +
  labs(
    x = NULL,
    y = "Correlation",
    fill = "Model Type"
  ) +
  scale_fill_manual(
    values = c("Logistic Regression" = "#535E31", "Random Forest" = "#A87B38")
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.2, 0.2),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 1),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10, col = "black", colour = "black"),
    axis.text.x = element_text(size = 12, face = "bold", colour = "black")
  )

cor_bp

ggsave(
  "./Data/Figures/Boxplots/Correlation_boxplot.png",
  cor_bp,
  width = 16,
  height = 12,
  units = "cm",
  dpi = 300
)

# ======================================================================================================== #
# ======================================================================================================== #

# SMALL BOXPLOTS: ####
# ========= #

# RMSE plot:
# ---------- #

rmse_bp <- ggplot(metrics, aes(x = mod_type, y = rmse, fill = sim)) +
  geom_boxplot(
    median.linewidth = 0.4,
    box.linewidth = 0.4,
    outlier.size = 0.3,
    colour = "black"
  ) +
  labs(
    x = NULL,
    y = "RMSE"
  ) +
  scale_fill_manual(
    values = c(
      "Linear" = "#6a040f",
      "Semi-Linear" = "#dc2f02",
      "Non-Linear" = "#faa307"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.6),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 9, col = "black", colour = "black"),
    axis.text.x = element_text(size = 9, face = "bold", colour = "black")
  )

rmse_bp

ggsave(
  "./Data/Figures/Boxplots/RMSE_boxplot_model_s.png",
  rmse_bp,
  width = 8,
  height = 6,
  units = "cm",
  dpi = 300
)


# Bias plot:
# ---------- #

bias_bp <- ggplot(metrics, aes(x = mod_type, y = bias, fill = sim)) +
  geom_boxplot(
    median.linewidth = 0.4,
    box.linewidth = 0.4,
    outlier.size = 0.3,
    colour = "black"
  ) +
  labs(
    x = NULL,
    y = "Bias"
  ) +
  scale_fill_manual(
    values = c(
      "Linear" = "#6a040f",
      "Semi-Linear" = "#dc2f02",
      "Non-Linear" = "#faa307"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.6),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 9, col = "black", colour = "black"),
    axis.text.x = element_text(size = 9, face = "bold", colour = "black")
  )

bias_bp

ggsave(
  "./Data/Figures/Boxplots/Bias_boxplot_model_s.png",
  bias_bp,
  width = 8,
  height = 6,
  units = "cm",
  dpi = 300
)

# Correlation plot:
# ----------------- #

cor_bp <- ggplot(metrics, aes(x = mod_type, y = cor, fill = sim)) +
  geom_boxplot(
    median.linewidth = 0.4,
    box.linewidth = 0.4,
    outlier.size = 0.3,
    colour = "black"
  ) +
  labs(
    x = NULL,
    y = "Correlation"
  ) +
  scale_fill_manual(
    values = c(
      "Linear" = "#6a040f",
      "Semi-Linear" = "#dc2f02",
      "Non-Linear" = "#faa307"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.6),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 9, col = "black", colour = "black"),
    axis.text.x = element_text(size = 9, face = "bold", colour = "black")
  )

cor_bp

ggsave(
  "./Data/Figures/Boxplots/Correlation_boxplot_model_s.png",
  cor_bp,
  width = 8,
  height = 6,
  units = "cm",
  dpi = 300
)

cor_bp <- ggplot(metrics, aes(x = mod_type, y = cor, fill = sim)) +
  geom_boxplot(
    median.linewidth = 0.4,
    box.linewidth = 0.4,
    outlier.size = 0.3,
    colour = "black"
  ) +
  labs(
    x = NULL,
    y = "Correlation",
    fill = "Level of Linearity"
  ) +
  scale_fill_manual(
    values = c(
      "Linear" = "#6a040f",
      "Semi-Linear" = "#dc2f02",
      "Non-Linear" = "#faa307"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.6),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 9, col = "black", colour = "black"),
    axis.text.x = element_text(size = 9, face = "bold", colour = "black")
  )

cor_bp

legend <- get_legend(cor_bp)
plot(legend)

ggsave(
  "./Data/Figures/Boxplots/legend_sims.png",
  legend,
  width = 8,
  height = 6,
  units = "cm",
  dpi = 300
)

# ======================================================================================================== #
# ======================================================================================================== #

{
  # ---- Editable style parameters ----
  line_thickness_curves <- 2.5 # thickness of the LR/RF density outlines
  line_thickness_axis <- 1 # thickness of the x and y axis lines
  font_size_axis_title <- 20 # size of the x-axis title text
  font_size_axis_text <- 18 # size of the axis tick label text
  overlap_fill_color <- "#34a0a4" # color of the shaded overlap region

  # ---- Compute the two densities manually ----
  # (assumes `metrics` has columns `rmse` and `mod_type`, with mod_type levels
  #  "Logistic Regression" and "Random Forest")
  dens_l <- density(metrics$rmse[
    metrics$sim == "Linear" & metrics$mod_type == "Logistic Regression"
  ])
  dens_sl <- density(metrics$rmse[
    metrics$sim == "Semi-Linear" & metrics$mod_type == "Logistic Regression"
  ])
  dens_nl <- density(metrics$rmse[
    metrics$sim == "Non-Linear" & metrics$mod_type == "Logistic Regression"
  ])

  # Put both curves on a common, shared x-grid so they can be compared point-by-point
  x_common <- seq(
    min(dens_l$x, dens_sl$x, dens_nl$x),
    max(dens_l$x, dens_sl$x, dens_nl$x),
    length.out = 512
  )

  y_l <- approx(dens_l$x, dens_l$y, xout = x_common, rule = 2)$y
  y_sl <- approx(dens_sl$x, dens_sl$y, xout = x_common, rule = 2)$y
  y_nl <- approx(dens_nl$x, dens_nl$y, xout = x_common, rule = 2)$y

  # Overlap = the lower of the two densities at each x value
  overlap_df <- data.frame(
    x = x_common,
    y_overlap = pmin(y_l, y_sl, y_nl)
  )

  # Long-format data frame for the two outline curves
  curve_df <- data.frame(
    x = rep(x_common, 2),
    y = c(y_l, y_sl, y_nl),
    sim = rep(
      c("Linear", "Semi-Linear", "Non-Linear"),
      each = length(x_common)
    )
  )

  # ---- Build the plot ----
  rmse_plot <- ggplot() +
    # geom_ribbon(
    #   data = overlap_df,
    #   aes(x = x, ymin = 0, ymax = y_overlap),
    #   fill = "#34a0a4",
    #   color = NA
    # ) +
    geom_line(
      data = curve_df,
      aes(x = x, y = y, color = sim),
      linewidth = line_thickness_curves
    ) +
    scale_color_manual(
      values = c("Logistic Regression" = "#3f8d5d", "Random Forest" = "#1B3A5C")
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = line_thickness_axis),
      axis.title.x = element_text(size = font_size_axis_title, face = "bold"),
      axis.title.y = element_blank(),
      axis.text = element_text(size = font_size_axis_text, face = "bold"),
      panel.background = element_blank(),
      plot.background = element_blank()
    )

  rmse_plot
}


densityPlot(
  metrics$rmse[metrics$mod_type == "Random Forest"] ~ metrics$sim[
    metrics$mod_type == "Random Forest"
  ],
  legend = list(title = "Simulation"),
  xlab = "RMSE",
  lwd = 4,
  grid = FALSE
)

?densityPlot


rf_data <- metrics[metrics$mod_type == "Logistic Regression", ]

ggplot(rf_data, aes(x = rmse, color = sim)) +
  geom_density(linewidth = 1.2) +
  geom_rug(sides = "b", linewidth = 0.5, alpha = 0.7) +
  scale_color_manual(
    name = "Simulation",
    values = c(
      "Linear" = "#6a040f",
      # "Random Forest" = "#6a040f",
      "Semi-Linear" = "#dc2f02",
      "Non-Linear" = "#faa307"
      # "Logistic Regression" = "#faa307"
    )
  ) +
  labs(x = "RMSE", y = "Density") +
  theme_classic()

rf_data <- metrics[metrics$sim == "Non-Linear", ]

ggplot(rf_data, aes(x = rmse, color = mod_type)) +
  geom_density(linewidth = 1.2, key_glyph = "path") +
  geom_rug(sides = "b", linewidth = 0.5, alpha = 0.7, show.legend = FALSE) +
  scale_color_manual(
    name = "Model",
    values = c(
      # "Linear" = "#6a040f",
      "Random Forest" = "#6a040f",
      # "Semi-Linear" = "#dc2f02",
      # "Non-Linear" = "#faa307"
      "Logistic Regression" = "#faa307"
    )
  ) +
  labs(x = "RMSE", y = "Density") +
  theme_classic()


rmse <- wide_list[["rmse"]]

idk <- ggplot(rmse, aes(x = diff, colour = sim)) +
  geom_density(linewidth = 1.2, key_glyph = "path") +
  geom_rug(sides = "b", linewidth = 0.5, alpha = 0.7, show.legend = FALSE) +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.8) +
  scale_color_manual(
    name = "Model",
    values = c(
      "Linear" = "#6a040f",
      "Semi-Linear" = "#dc2f02",
      "Non-Linear" = "#faa307"
    )
  ) +
  labs(x = "RMSE", y = "Density") +
  theme_classic() +
  theme(
    # legend.position = c(0.2, 0.2),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 1),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10, col = "black", colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black")
  )

idk

ggsave(
  paste0("./Data/Figures/Density_Plots/rmse_diff_density_plot.png"),
  plot = idk,
  width = 16,
  height = 12,
  units = "cm",
  dpi = 300
)


for (s in sim_levels) {
  rmse <- wide_list[["rmse"]] |>
    filter(sim == s)

  plot <- ggplot(rmse, aes(x = diff)) +
    geom_density(colour = "#dc2f02", linewidth = 1.2, key_glyph = "path") +
    geom_vline(xintercept = 0, colour = "black", linewidth = 0.8) +
    geom_rug(
      sides = "b",
      linewidth = 0.5,
      alpha = 0.7,
      show.legend = FALSE,
      colour = "#dc2f02"
    ) +
    labs(x = "RMSE", y = "Density") +
    theme_classic() +
    theme(
      # legend.position = c(0.2, 0.2),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10, face = "bold"),
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 1),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 10, col = "black", colour = "black"),
      axis.text.x = element_text(size = 12, colour = "black")
    )

  ggsave(
    paste0("./Data/Figures/Density_Plots/rmse_diff_density_plot_", s, ".png"),
    plot,
    width = 16,
    height = 12,
    units = "cm",
    dpi = 300
  )
}


rmse <- wide_list[["rmse"]] |>
  filter(sim == "Non-Linear")

ggplot(rmse, aes(x = diff)) +
  geom_density(colour = "#dc2f02", linewidth = 1.2, key_glyph = "path") +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.8) +
  geom_rug(
    sides = "b",
    linewidth = 0.5,
    alpha = 0.7,
    show.legend = FALSE,
    colour = "#dc2f02"
  ) +
  labs(x = "RMSE", y = "Density") +
  theme_classic() +
  theme(
    # legend.position = c(0.2, 0.2),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 1),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10, col = "black", colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black")
  )


rmse <- wide_list[["rmse"]] |>
  filter(sim == "Semi-Linear")

rf_better <- (nrow(rmse[rmse$diff < 0, ]) / nrow(rmse)) * 100
lr_better <- (nrow(rmse[rmse$diff > 0, ]) / nrow(rmse)) * 100

rf_better
lr_better

# ======================================================================================================== #
# ======================================================================================================== #
