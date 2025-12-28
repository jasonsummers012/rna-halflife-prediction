library(keras3)
library(tidyverse)

# ===============================================================
# Load prepared caspase 9 data
# ===============================================================
casp9_data = readRDS("caspase9_prepared_for_prediction.rds")

x_sequences = casp9_data$x_sequences
features = casp9_data$features
y_mean = casp9_data$y_mean
y_sd = casp9_data$y_sd

cat("Loaded caspase 9 data:\n")
cat("Sequences:", dim(x_sequences), "\n")
cat("Features:", dim(features), "\n")
cat("Wildtype half-life:", casp9_data$wildtype_halflife, "hours\n\n")

# ===============================================================
# DIAGNOSTIC: Check if sequences are actually different
# ===============================================================
cat("=== CHECKING SEQUENCE DIFFERENCES ===\n")

# Get actual sequence data (remove padding and mask)
wildtype_seq_full = x_sequences[1, , ]
mutant1_seq_full = x_sequences[2, , ]

# Find non-padded region (where mask = 1)
wildtype_mask = wildtype_seq_full[, 6]
non_padded_start = which(wildtype_mask == 1)[1]
non_padded_end = tail(which(wildtype_mask == 1), 1)

cat("Non-padded region: positions", non_padded_start, "to", non_padded_end, "\n")
cat("Sequence length:", non_padded_end - non_padded_start + 1, "\n\n")

# Compare actual sequences in non-padded region
wildtype_seq = wildtype_seq_full[non_padded_start:non_padded_end, 1:5]
mutant1_seq = mutant1_seq_full[non_padded_start:non_padded_end, 1:5]

# Count differences
differences = sum(rowSums(abs(wildtype_seq - mutant1_seq)) > 0)
cat("Number of positions with differences:", differences, "\n")
cat("Expected mutations: ~3-10 depending on mutant\n\n")

if(differences == 0) {
  stop("ERROR: No differences detected between wildtype and mutants! Check data preparation.")
}

# ===============================================================
# DIAGNOSTIC: Check feature differences
# ===============================================================
cat("=== CHECKING FEATURE DIFFERENCES ===\n")
wildtype_features = features[1, ]
mutant1_features = features[2, ]

feature_diff = abs(wildtype_features - mutant1_features)
cat("Features with differences > 0.01:\n")
print(data.frame(
  feature = casp9_data$feature_names,
  wildtype = round(wildtype_features, 4),
  mutant1 = round(mutant1_features, 4),
  diff = round(feature_diff, 4)
) %>% filter(diff > 0.01))
cat("\n")

# ===============================================================
# Rebuild model architecture (must match training exactly!)
# ===============================================================
cat("Rebuilding model architecture...\n")

set_random_seed(12)

# Input 1: Sequence data
seq_input = layer_input(shape = c(3500, 6), name = "sequence_input")

seq_features = layer_lambda(f = function(x) x[,,1:5])(seq_input)
seq_mask = layer_lambda(f = function(x) x[,,6])(seq_input)

seq_features_masked = layer_multiply(list(
  seq_features,
  layer_lambda(f = function(x) op_expand_dims(x, axis = -1))(seq_mask)
))

# CNN layers for sequence
seq_output = seq_features_masked %>%
  layer_conv_1d(64, 3, dilation_rate=1, padding="same", activation="relu") %>%
  layer_conv_1d(64, 3, dilation_rate=2, padding="same", activation="relu") %>%
  layer_conv_1d(128, 3, dilation_rate=4, padding="same", activation="relu") %>%
  layer_global_max_pooling_1d()

# Input 2: Additional features
feature_input = layer_input(shape = ncol(features), name = "feature_input")

feature_output = feature_input %>%
  layer_dense(128, activation = "relu") %>%
  layer_batch_normalization() %>%
  layer_dropout(0.3) %>%
  
  layer_dense(64, activation = "relu") %>%
  layer_batch_normalization() %>%
  layer_dropout(0.3) %>%
  
  layer_dense(32, activation = "relu") %>%
  layer_dropout(0.2)

# Combine both pathways
combined = layer_concatenate(list(seq_output, feature_output))

# Final output layers
output = combined %>%
  layer_dense(128, activation = "relu") %>%
  layer_dropout(0.5) %>%
  layer_dense(64, activation = "relu") %>%
  layer_dropout(0.4) %>%
  layer_dense(1)

# Create model
model = keras_model(
  inputs = list(seq_input, feature_input),
  outputs = output
)

cat("Model architecture rebuilt\n\n")

# ===============================================================
# Load trained weights
# ===============================================================
cat("Loading trained weights...\n")
model$load_weights("rna_halflife_improved_weights.weights.h5")
cat("Weights loaded successfully!\n\n")

# ===============================================================
# DIAGNOSTIC: Test predictions with dropout off
# ===============================================================
cat("=== TESTING PREDICTION CONSISTENCY ===\n")

# Make 5 predictions on wildtype to check if dropout is the issue
test_x = x_sequences[1, , , drop = FALSE]
test_f = features[1, , drop = FALSE]

predictions_test = numeric(5)
for(i in 1:5) {
  pred = predict(model, list(test_x, test_f), verbose = 0)
  predictions_test[i] = pred[1]
}

cat("5 predictions on wildtype (should be identical if dropout is off):\n")
print(round(predictions_test, 6))

if(sd(predictions_test) > 0.01) {
  warning("Predictions vary! Dropout may still be active during prediction.")
}
cat("\n")

# ===============================================================
# Make predictions on all sequences
# ===============================================================
cat("Making predictions on all sequences...\n")

# Clean data (ensure numeric)
x_clean = array(as.numeric(x_sequences), dim = dim(x_sequences))
features_clean = array(as.numeric(features), dim = dim(features))

# Predict (dropout should be off automatically in predict mode)
predictions_scaled = predict(model, list(x_clean, features_clean), verbose = 0)

# Transform back to original scale (hours)
predictions_log = predictions_scaled * y_sd + y_mean
predictions_halflife = expm1(predictions_log)

cat("Predictions complete!\n\n")

# ===============================================================
# Check if predictions vary
# ===============================================================
cat("=== PREDICTION VARIATION CHECK ===\n")
cat("Unique prediction values:", length(unique(round(predictions_halflife, 6))), "\n")
cat("Range of predictions:", round(min(predictions_halflife), 4), "to", 
    round(max(predictions_halflife), 4), "hours\n")
cat("Standard deviation:", round(sd(predictions_halflife), 4), "hours\n\n")

if(length(unique(round(predictions_halflife, 4))) == 1) {
  warning("All predictions are identical! This suggests a problem.")
  cat("\nDEBUGGING INFO:\n")
  cat("This could mean:\n")
  cat("1. The mutations aren't being detected by the model\n")
  cat("2. The model is 'saturated' and predicting the same value\n")
  cat("3. There's an issue with how the data was prepared\n\n")
}

# ===============================================================
# Organize results
# ===============================================================
results = data.frame(
  sequence = casp9_data$sequence_names,
  mutation_count = casp9_data$mutation_counts,
  predicted_halflife = as.vector(predictions_halflife),
  actual_wildtype_halflife = casp9_data$wildtype_halflife
)

# Calculate change from wildtype
wildtype_prediction = results$predicted_halflife[1]
results$change_from_wildtype = results$predicted_halflife - wildtype_prediction
results$percent_change = (results$change_from_wildtype / wildtype_prediction) * 100

# ===============================================================
# Display results
# ===============================================================
cat("========================================\n")
cat("CASPASE 9 HALF-LIFE PREDICTIONS\n")
cat("========================================\n\n")

cat("WILDTYPE:\n")
cat("  Actual half-life:", casp9_data$wildtype_halflife, "hours\n")
cat("  Predicted half-life:", round(wildtype_prediction, 2), "hours\n")
cat("  Prediction error:", round(abs(casp9_data$wildtype_halflife - wildtype_prediction), 2), "hours\n\n")

cat("MUTANTS BY MUTATION COUNT:\n")
cat("========================================\n\n")

for(n_mut in sort(unique(casp9_data$mutation_counts[-1]))) {
  mutants = results %>% filter(mutation_count == n_mut)
  
  cat(n_mut, "MUTATIONS (n =", nrow(mutants), "):\n")
  cat("  Mean predicted half-life:", round(mean(mutants$predicted_halflife), 2), "hours\n")
  cat("  Range:", round(min(mutants$predicted_halflife), 2), "-", 
      round(max(mutants$predicted_halflife), 2), "hours\n")
  cat("  Mean change from wildtype:", round(mean(mutants$change_from_wildtype), 2), "hours",
      "(", round(mean(mutants$percent_change), 1), "%)\n\n")
}

# ===============================================================
# Detailed results table
# ===============================================================
cat("\nDETAILED RESULTS (first 20):\n")
cat("========================================\n")
results_display = results %>% 
  mutate(
    predicted_halflife = round(predicted_halflife, 3),
    change_from_wildtype = round(change_from_wildtype, 3),
    percent_change = round(percent_change, 2)
  ) %>%
  select(sequence, mutation_count, predicted_halflife, change_from_wildtype, percent_change)

print(results_display)

# ===============================================================
# Statistical summary (only if predictions vary)
# ===============================================================
if(length(unique(round(predictions_halflife, 4))) > 1) {
  cat("\n\nSTATISTICAL SUMMARY:\n")
  cat("========================================\n")
  
  mutants_only = results %>% filter(mutation_count > 0)
  
  cat("Correlation between mutation count and half-life change:\n")
  cor_test = cor.test(mutants_only$mutation_count, mutants_only$change_from_wildtype)
  cat("  Pearson r =", round(cor_test$estimate, 3), "\n")
  cat("  p-value =", format.pval(cor_test$p.value, digits = 3), "\n\n")
  
  summary_stats = results %>%
    filter(mutation_count > 0) %>%
    group_by(mutation_count) %>%
    summarize(
      n = n(),
      mean_change = mean(change_from_wildtype),
      sd_change = sd(change_from_wildtype),
      mean_pct = mean(percent_change)
    )
  print(summary_stats)
}

# ===============================================================
# Save results
# ===============================================================
saveRDS(results, "caspase9_predictions.rds")
write.csv(results, "caspase9_predictions.csv", row.names = FALSE)

cat("\n\nResults saved to:\n")
cat("  - caspase9_predictions.rds\n")
cat("  - caspase9_predictions.csv\n")