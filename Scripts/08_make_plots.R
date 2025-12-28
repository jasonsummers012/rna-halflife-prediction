library(tidyverse)

# ===============================================================
# Shared feature labels (define once, use everywhere)
# ===============================================================
feature_labels = c(
  "seq_length" = "Sequence length",
  "gc_content" = "GC content",
  "au_content" = "AU content",
  "are_score_3utr" = "ARE score",
  "are_score" = "ARE score",
  "polya_signal" = "PolyA signal",
  "gc_3utr" = "GC content 3'UTR",
  "mean_pairing" = "Mean pairing probability",
  "pairing_3utr" = "3'UTR mean pairing probability",
  "orf_length" = "ORF length",
  "orf_ratio" = "ORF ratio",
  "cai" = "CAI",
  "mirna_sites" = "Number of miRNA sites",
  "mirna_density" = "miRNA density"
)

# ===============================================================
# Figure 1: Model Performance
# ===============================================================
model_data = readRDS("improved_model_predictions.rds")

actual_halflife = model_data$actual
predicted_halflife = model_data$predicted

long_lived = actual_halflife > 15

fig1_data = tibble(
  actual = actual_halflife,
  predicted = predicted_halflife,
  long_lived = long_lived
)

figure1 = ggplot(fig1_data, aes(x = actual, y = predicted)) +
  geom_point(aes(color = long_lived)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_smooth(method = "lm", se = FALSE, color = "blue4", linewidth = 1) +
  annotate(
    "text",
    x = 22.5,
    y = 10,
    label = "r = 0.227",
    color = "blue4",
    size = 4
  ) +
  coord_cartesian(xlim = c(0, 25), ylim = c(0, 25)) +
  scale_color_manual(
    name = "Transcript longevity",
    values = c("FALSE" = "gray50", "TRUE" = "red4"),
    labels = c(
      "FALSE" = "Not long-lived",
      "TRUE" = "Long-lived"
    )
  ) +
  theme_bw() +
  labs(
    x = "Actual Half-life (hours)",
    y = "Predicted Half-life (hours)"
  )

ggsave("Plots/figure1.png", plot = figure1, width = 6, height = 5, dpi = 300)

# ===============================================================
# Figure 2: Test Set Feature Importance
# ===============================================================
fig2_data = readRDS("test_set_feature_correlations.rds")

# Keep only features with labels
fig2_data = fig2_data[fig2_data$feature %in% names(feature_labels), ]

# Use correlation with predicted (or you could use corr_with_actual)
fig2_data$correlation = fig2_data$corr_with_predicted

# Replace NA correlations with 0 (if any)
fig2_data$correlation[is.na(fig2_data$correlation)] = 0

# Create correlation sign factor
fig2_data$cor_sign = factor(
  ifelse(fig2_data$correlation > 0, "Positive", "Negative"),
  levels = c("Positive", "Negative")
)

# Make feature a factor ordered by correlation (descending)
fig2_data$feature = factor(
  fig2_data$feature, 
  levels = fig2_data$feature[order(-fig2_data$correlation)]
)

figure2 = ggplot(fig2_data, aes(
  x = feature,
  y = correlation,
  fill = cor_sign
)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(
    name = "Correlation",
    values = c("Positive" = "blue4", "Negative" = "red4")
  ) +
  scale_x_discrete(labels = feature_labels) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = "Feature",
    y = "Correlation with predicted half-life"
  )

ggsave("Plots/figure2.png", plot = figure2, width = 6, height = 5, dpi = 300)

# ===============================================================
# Figure 3: CASP9 Mutation Effects
# ===============================================================
casp9_results = readRDS("caspase9_predictions.rds")

fig3_data= casp9_results %>%
  mutate(mutation_group = factor(mutation_count, levels = c(0, 10, 50, 100))) %>%
  filter(mutation_group %in% c(10, 50, 100))

wildtype_prediction = casp9_results %>%
  filter(mutation_count == 0) %>%
  pull(predicted_halflife)

figure3 = ggplot(fig3_data, aes(x = mutation_group, y = predicted_halflife)) +
  geom_boxplot() +
  geom_hline(yintercept = wildtype_prediction[1], color = "red4", linetype = "dashed", linewidth = 1) +
  theme_bw() +
  labs(
    x = "Number of Mutations",
    y = "Predicted Half-life (hours)"
  )

ggsave("Plots/figure3.png", plot = figure3, width = 6, height = 5, dpi = 300)

# ===============================================================
# Figure 4: CASP9 Feature Correlations
# ===============================================================
fig4_data = readRDS("feature_correlations.rds")

# Keep only features with labels
fig4_data = fig4_data[fig4_data$feature %in% names(feature_labels), ]

# Replace NA correlations with 0
fig4_data$correlation[is.na(fig4_data$correlation)] = 0

# Create correlation sign factor (zero goes with "Negative")
fig4_data$cor_sign = factor(
  ifelse(fig4_data$correlation > 0, "Positive", "Negative"),
  levels = c("Positive", "Negative")
)

# Make feature a factor ordered by correlation (descending)
fig4_data$feature = factor(fig4_data$feature, levels = fig4_data$feature[order(-fig4_data$correlation)])

figure4 = ggplot(fig4_data, aes(
  x = feature,
  y = correlation,
  fill = cor_sign
)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(
    name = "Correlation",
    values = c("Positive" = "blue4", "Negative" = "red4")
  ) +
  scale_x_discrete(labels = feature_labels) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = "Feature",
    y = "Correlation"
  )

ggsave("Plots/figure4.png", plot = figure4, width = 6, height = 5, dpi = 300)