# analyze_coloc.R
#
# Loads per-image, per-channel-pair colocalization results from the
# Fiji batch macro, joins them to sample metadata (condition +
# biological replicate), aggregates to replicate level (separately
# for each channel pair), plots box plots per channel pair, and runs
# a stats test on the replicate-level means (NOT the per-image
# values, to avoid pseudoreplication).
#
# Expects two input files:
#
#   coloc_results.csv     (output of batch_coloc2.ijm)
#     filename, channel_1, channel_2, M1, M2, pearson, ...,
#     costes_pvalue, costes_ratio_rand_ge_actual,
#     warning_y_intercept_far_from_zero, warning_ch1_threshold_too_high,
#     warning_ch2_threshold_too_high
#
#   sample_metadata.csv    (you create this once, by hand)
#     filename, condition, replicate

library(tidyverse)
library(ggsignif)

coloc <- read_csv("/Users/jlevendis/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/m6A\ paper/Colocalisation/results/coloc_results.csv")
meta  <- read_csv("/Users/jlevendis/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/m6A\ paper/Colocalisation/sample_metadata.csv")

coloc <- coloc %>% mutate(filename = trimws(filename))
meta  <- meta  %>% mutate(filename = trimws(filename))

channel_labels <- c(
  "C1" = "mCherry",
  "C2" = "DAPI",
  "C3" = "GFP"
)

coloc <- coloc %>%
  mutate(
    channel_1 = trimws(channel_1),
    channel_2 = trimws(channel_2),
    channel_1_code = str_extract(channel_1, "^C[0-9]+"),
    channel_2_code = str_extract(channel_2, "^C[0-9]+"),
    channel_1_label = recode(channel_1_code, !!!channel_labels),
    channel_2_label = recode(channel_2_code, !!!channel_labels),
    channel_pair = paste0(channel_1_label, "_vs_", channel_2_label)
  )

df <- coloc %>%
  left_join(meta, by = "filename")

# sanity check: flag any images that didn't match metadata, or vice versa
unmatched <- df %>% filter(is.na(condition))
if (nrow(unmatched) > 0) {
  warning(nrow(unmatched), " row(s) in coloc_results.csv had no matching row in sample_metadata.csv:\n",
          paste(unique(unmatched$filename), collapse = "\n"))
  print(unique(unmatched$filename))
}

df <- df %>% filter(!is.na(condition))

# --- flag unreliable rows based on Coloc2's own warnings, rather than
#     silently dropping them - decide per-analysis whether to exclude ---
df <- df %>%
  mutate(
    any_warning = warning_y_intercept_far_from_zero == "TRUE" |
      warning_ch1_threshold_too_high == "TRUE" |
      warning_ch2_threshold_too_high == "TRUE"
  )

n_flagged <- sum(df$any_warning, na.rm = TRUE)
if (n_flagged > 0) {
  message(n_flagged, " row(s) flagged with at least one Coloc2 warning (y-intercept or threshold-too-high). ",
          "These are NOT automatically excluded - inspect df %>% filter(any_warning) before deciding.")
}

# --- aggregate to replicate level, PER CHANNEL PAIR (this is your real
#     unit of statistical replication, kept separate per pair since M1/M2/
#     pearson are not comparable across different channel combinations) ---
df_replicate <- df %>%
  group_by(condition, replicate, channel_pair) %>%
  summarise(
    mean_M1              = mean(M1, na.rm = TRUE),
    mean_M2              = mean(M2, na.rm = TRUE),
    mean_pearson          = mean(pearson, na.rm = TRUE),
    mean_costes_pvalue    = mean(costes_pvalue, na.rm = TRUE),
    mean_costes_ratio     = mean(costes_ratio_rand_ge_actual, na.rm = TRUE),
    n_images              = n(),
    n_flagged             = sum(any_warning, na.rm = TRUE),
    .groups = "drop"
  )

# --- run stats + plot separately for each channel pair ---
conditions <- unique(df_replicate$condition)
channel_pairs <- unique(df_replicate$channel_pair)

if (length(conditions) != 2) {
  message("More than two conditions detected (", paste(conditions, collapse = ", "),
          ") - use a Kruskal-Wallis test or ANOVA instead of wilcox.test. ",
          "Significance brackets below assume exactly two conditions and will be skipped.")
}

# --- condition color mapping ---
condition_colors <- c(
  "control"       = "#1bb6bb",
  "knocksideways" = "#f46b64"
)

for (cp in channel_pairs) {
  
  cat("\n========================================\n")
  cat("Channel pair:", cp, "\n")
  cat("========================================\n")
  
  df_cp <- df_replicate %>% filter(channel_pair == cp)
  df_cp_images <- df %>% filter(channel_pair == cp)
  
  if (length(conditions) == 2) {
    test_result_m1      <- wilcox.test(M1 ~ condition, data = df_cp_images, exact = FALSE)
    test_result_m2       <- wilcox.test(M2 ~ condition, data = df_cp_images, exact = FALSE)
    test_result_pearson  <- wilcox.test(pearson ~ condition, data = df_cp_images, exact = FALSE)
    
    cat("\n--- M1 test ---\n");      print(test_result_m1)
    cat("\n--- M2 test ---\n");      print(test_result_m2)
    cat("\n--- Pearson's R test ---\n"); print(test_result_pearson)
  }
  
  # --- Pearson's R violin plot: raw per-image distribution, colored by
  #     condition (control = red, knocksideways = green) ---
  
  p_pearson <- ggplot(df_cp_images, aes(x = condition, y = pearson, fill = condition)) +
    geom_violin(trim = FALSE, color = "grey30", alpha = 0.5) +
    geom_jitter(width = 0.1, size = 1.5, alpha = 0.8, aes(color = condition)) +
    scale_fill_manual(values = condition_colors) +
    scale_color_manual(values = condition_colors) +
    ylim(-0.5, 1) +
    labs(
      x = "Condition",
      y = "Pearson's R",
      title = paste0("Colocalization by condition (Pearson's R) - ", cp)
    ) +
    theme_classic(base_size = 13) +
    theme(legend.position = "none")   # colors already shown via x-axis labels; legend is redundant
  
  if (length(conditions) == 2) {
    
    # Find the highest observed data point
    max_pearson <- max(df_cp_images$pearson, na.rm = TRUE)
    
    # Place significance bar just above it
    signif_y <- max_pearson + 0.3
    if (signif_y > 0.9) {
      signif_y <- 0.9
    }
    
    p_pearson <- p_pearson +
      geom_signif(
        comparisons = list(conditions),
        annotations = paste0("p = ", signif(test_result_pearson$p.value, 3)),
        y_position = signif_y,
        textsize = 4,
        tip_length = 0.02
      )
  }
  print(p_pearson)
  print(df_cp)
  ggsave(paste0("coloc_boxplot_pearson_", cp, ".pdf"), p_pearson, width = 6, height = 5)
}

# --- save the replicate-level summary table (all channel pairs together),
#     for reporting/methods ---
write_csv(df_replicate, "coloc_replicate_summary.csv")