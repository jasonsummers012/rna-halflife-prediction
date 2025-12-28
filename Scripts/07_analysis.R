library(tidyverse)

# ===============================================================
# Load all data
# ===============================================================
casp9_data = readRDS("caspase9_prepared_for_prediction.rds")
mutant_data = readRDS("caspase9_mutants.rds")
predictions = readRDS("caspase9_predictions.rds")

cat("=== CODON MUTATION EFFECT ANALYSIS ===\n\n")

# ===============================================================
# 1. Overall effect summary
# ===============================================================
cat("1. OVERALL EFFECTS:\n")
cat("========================================\n")

summary_by_count = predictions %>%
  filter(mutation_count > 0) %>%
  group_by(mutation_count) %>%
  summarize(
    n_mutants = n(),
    mean_halflife = mean(predicted_halflife),
    sd_halflife = sd(predicted_halflife),
    mean_change_hours = mean(change_from_wildtype),
    mean_change_pct = mean(percent_change),
    min_change = min(change_from_wildtype),
    max_change = max(change_from_wildtype)
  )

print(summary_by_count)

cat("\nWildtype predicted half-life:", 
    round(predictions$predicted_halflife[1], 3), "hours\n\n")

# ===============================================================
# 2. Feature changes that correlate with stability changes
# ===============================================================
cat("2. FEATURE CORRELATIONS WITH STABILITY:\n")
cat("========================================\n")

features = casp9_data$features
feature_names = casp9_data$feature_names

# Get wildtype features
wildtype_features = features[1, ]

# Calculate feature deltas for each mutant
mutant_indices = 2:nrow(features)
feature_deltas = sweep(features[mutant_indices, ], 2, wildtype_features, "-")

# Get stability changes
stability_changes = predictions$change_from_wildtype[mutant_indices]

# Calculate correlation for each feature
correlations = data.frame(
  feature = feature_names,
  correlation = sapply(1:ncol(feature_deltas), function(i) {
    cor(feature_deltas[, i], stability_changes, use = "complete.obs")
  }),
  mean_delta = colMeans(abs(feature_deltas))
)

correlations = correlations %>%
  arrange(desc(abs(correlation))) %>%
  mutate(
    correlation = round(correlation, 3),
    mean_delta = round(mean_delta, 4)
  )

cat("\nFeatures ranked by correlation with stability change:\n")
print(correlations)

cat("\nTop features driving stability changes:\n")
top_features = head(correlations, 5)
for(i in 1:nrow(top_features)) {
  cat(sprintf("  %d. %s (r=%.3f, mean Δ=%.4f)\n",
              i, top_features$feature[i], 
              top_features$correlation[i],
              top_features$mean_delta[i]))
}

# ===============================================================
# 2B. Feature importance across test set (for comparison)
# ===============================================================
cat("\n\n2B. FEATURE IMPORTANCE IN TEST SET:\n")
cat("========================================\n")
cat("(For comparison: what features matter for predictions overall)\n\n")

# Load test set data
cnn_data = readRDS("cnn_training_data_with_features.rds")
test_results = readRDS("improved_model_predictions.rds")

features_test = cnn_data$features_test
predictions_test = test_results$predicted
actual_test = test_results$actual

# Calculate correlation for each feature
test_feature_correlations = data.frame(
  feature = feature_names,
  corr_with_predicted = sapply(1:ncol(features_test), function(i) {
    cor(features_test[, i], predictions_test, use = "complete.obs")
  }),
  corr_with_actual = sapply(1:ncol(features_test), function(i) {
    cor(features_test[, i], actual_test, use = "complete.obs")
  }),
  mean_value = colMeans(features_test),
  sd_value = apply(features_test, 2, sd)
)

# Sort by correlation with predictions
test_feature_correlations = test_feature_correlations %>%
  arrange(desc(abs(corr_with_predicted))) %>%
  mutate(
    corr_with_predicted = round(corr_with_predicted, 3),
    corr_with_actual = round(corr_with_actual, 3),
    mean_value = round(mean_value, 4),
    sd_value = round(sd_value, 4)
  )

cat("\nFeatures ranked by correlation with PREDICTED half-life (test set):\n")
print(test_feature_correlations)

cat("\n\nTop 5 features most important for overall predictions:\n")
top_test_features = head(test_feature_correlations, 5)
for(i in 1:nrow(top_test_features)) {
  cat(sprintf("  %d. %s (r=%.3f with predicted, r=%.3f with actual)\n",
              i, 
              top_test_features$feature[i], 
              top_test_features$corr_with_predicted[i],
              top_test_features$corr_with_actual[i]))
}

cat("\n\nComparison insight:\n")
cat("Features that change with synonymous mutations (Section 2) show different\n")
cat("importance than features that drive overall predictions (Section 2B).\n")
cat("This suggests the model's robustness to synonymous changes is because\n")
cat("many mutation-sensitive features (e.g., CAI, miRNA sites) have relatively\n")
cat("modest importance for overall stability predictions.\n")

# Save test set feature correlations for plotting
saveRDS(test_feature_correlations, "test_set_feature_correlations.rds")
cat("\nSaved test set feature correlations to: test_set_feature_correlations.rds\n")

# ===============================================================
# 3. Mutation count vs stability (with visualization data)
# ===============================================================
cat("\n\n3. MUTATION COUNT EFFECTS:\n")
cat("========================================\n")

# Linear model
lm_result = lm(change_from_wildtype ~ mutation_count, 
               data = predictions %>% filter(mutation_count > 0))

cat("\nLinear regression: change ~ mutation_count\n")
cat("  Slope:", round(coef(lm_result)[2], 6), "hours per mutation\n")
cat("  R-squared:", round(summary(lm_result)$r.squared, 4), "\n")
cat("  P-value:", format.pval(summary(lm_result)$coefficients[2, 4], digits = 3), "\n")

# Test for non-linear relationship
lm_quad = lm(change_from_wildtype ~ mutation_count + I(mutation_count^2),
             data = predictions %>% filter(mutation_count > 0))

cat("\nQuadratic model: change ~ mutation_count + mutation_count²\n")
cat("  Linear term:", round(coef(lm_quad)[2], 6), "\n")
cat("  Quadratic term:", round(coef(lm_quad)[3], 8), "\n")
cat("  R-squared:", round(summary(lm_quad)$r.squared, 4), "\n")
cat("  P-value:", format.pval(anova(lm_result, lm_quad)$`Pr(>F)`[2], digits = 3), "\n")

# ===============================================================
# 4. Individual mutant analysis
# ===============================================================
cat("\n\n4. EXTREME MUTANTS:\n")
cat("========================================\n")

# Most stabilizing
most_stable = predictions %>% 
  filter(mutation_count > 0) %>%
  slice_max(change_from_wildtype, n = 3)

cat("\nMost stabilizing mutants:\n")
for(i in 1:nrow(most_stable)) {
  cat(sprintf("  %s: +%.3f hours (+%.2f%%) [%d mutations]\n",
              most_stable$sequence[i],
              most_stable$change_from_wildtype[i],
              most_stable$percent_change[i],
              most_stable$mutation_count[i]))
}

# Most destabilizing
most_unstable = predictions %>%
  filter(mutation_count > 0) %>%
  slice_min(change_from_wildtype, n = 3)

cat("\nMost destabilizing mutants:\n")
for(i in 1:nrow(most_unstable)) {
  cat(sprintf("  %s: %.3f hours (%.2f%%) [%d mutations]\n",
              most_unstable$sequence[i],
              most_unstable$change_from_wildtype[i],
              most_unstable$percent_change[i],
              most_unstable$mutation_count[i]))
}

# ===============================================================
# 5. Statistical significance tests
# ===============================================================
cat("\n\n5. STATISTICAL TESTS:\n")
cat("========================================\n")

# Are mutants significantly different from wildtype?
wildtype_hl = predictions$predicted_halflife[1]
mutant_hls = predictions$predicted_halflife[predictions$mutation_count > 0]

t_test = t.test(mutant_hls, mu = wildtype_hl)
cat("\nOne-sample t-test: Are mutants different from wildtype?\n")
cat("  Mean mutant half-life:", round(mean(mutant_hls), 4), "hours\n")
cat("  Wildtype half-life:", round(wildtype_hl, 4), "hours\n")
cat("  Mean difference:", round(mean(mutant_hls) - wildtype_hl, 4), "hours\n")
cat("  P-value:", format.pval(t_test$p.value, digits = 3), "\n")

# ANOVA: Do different mutation counts have different effects?
anova_result = aov(predicted_halflife ~ as.factor(mutation_count),
                   data = predictions %>% filter(mutation_count > 0))
cat("\nANOVA: Do mutation counts differ?\n")
cat("  P-value:", format.pval(summary(anova_result)[[1]][["Pr(>F)"]][1], digits = 3), "\n")

# ===============================================================
# 6. Practical interpretation
# ===============================================================
cat("\n\n6. BIOLOGICAL INTERPRETATION:\n")
cat("========================================\n")

effect_per_mutation = mean(mutant_hls - wildtype_hl) / mean(predictions$mutation_count[predictions$mutation_count > 0])

cat(sprintf("\nAverage effect per synonymous mutation: %.4f hours (%.3f%%)\n",
            effect_per_mutation,
            (effect_per_mutation / wildtype_hl) * 100))

cat("\nFor context:\n")
cat("  - Wildtype half-life:", round(wildtype_hl, 2), "hours\n")
cat("  - Typical range with", max(predictions$mutation_count), "mutations:",
    round(min(mutant_hls), 2), "-", round(max(mutant_hls), 2), "hours\n")
cat("  - Maximum observed change:", round(max(abs(predictions$change_from_wildtype[predictions$mutation_count > 0])), 3), 
    "hours (", round(max(abs(predictions$percent_change[predictions$mutation_count > 0])), 2), "%)\n")

# ===============================================================
# 7. Create plot data
# ===============================================================
cat("\n\n7. CREATING VISUALIZATION DATA:\n")
cat("========================================\n")

plot_data = predictions %>%
  filter(mutation_count > 0) %>%
  select(mutation_count, predicted_halflife, change_from_wildtype, percent_change)

write.csv(plot_data, "mutation_effects_for_plotting.csv", row.names = FALSE)
cat("Saved plot data to: mutation_effects_for_plotting.csv\n")

# Summary statistics for plotting
summary_for_plot = predictions %>%
  filter(mutation_count > 0) %>%
  group_by(mutation_count) %>%
  summarize(
    mean_halflife = mean(predicted_halflife),
    se_halflife = sd(predicted_halflife) / sqrt(n()),
    mean_change = mean(change_from_wildtype),
    se_change = sd(change_from_wildtype) / sqrt(n())
  )

write.csv(summary_for_plot, "mutation_summary_for_plotting.csv", row.names = FALSE)
cat("Saved summary data to: mutation_summary_for_plotting.csv\n")
saveRDS(correlations, "feature_correlations.rds")

cat("\n=== ANALYSIS COMPLETE ===\n")