library(keras3)
library(tidyverse)

# ===============================================================
# Load data with features
# ===============================================================
cnn_data = readRDS("cnn_training_data_with_features.rds")

# Sequence data
x_training = cnn_data$x_training
x_val = cnn_data$x_val
x_test = cnn_data$x_test

# Additional features
features_train = cnn_data$features_train
features_val = cnn_data$features_val
features_test = cnn_data$features_test

# Target variable
y_train_scaled = cnn_data$y_training_scaled
y_val_scaled = cnn_data$y_val_scaled
y_test_scaled = cnn_data$y_test_scaled
y_mean = cnn_data$y_mean
y_sd = cnn_data$y_sd

cat("Loaded data:\n")
cat("Sequence shape:", dim(x_training), "\n")
cat("Features shape:", dim(features_train), "\n")
cat("Feature names:", paste(cnn_data$feature_names, collapse = ", "), "\n\n")

set_random_seed(12)

# ===============================================================
# Build dual-input model with improved feature network
# ===============================================================

# Input 1: Sequence data (SAME AS ORIGINAL - this works well)
seq_input = layer_input(shape = c(3500, 6), name = "sequence_input")

seq_features = layer_lambda(f = function(x) x[,,1:5])(seq_input)
seq_mask = layer_lambda(f = function(x) x[,,6])(seq_input)

seq_features_masked = layer_multiply(list(
  seq_features,
  layer_lambda(f = function(x) op_expand_dims(x, axis = -1))(seq_mask)
))

# CNN layers for sequence (UNCHANGED)
seq_output = seq_features_masked %>%
  layer_conv_1d(64, 3, dilation_rate=1, padding="same", activation="relu") %>%
  layer_conv_1d(64, 3, dilation_rate=2, padding="same", activation="relu") %>%
  layer_conv_1d(128, 3, dilation_rate=4, padding="same", activation="relu") %>%
  layer_global_max_pooling_1d()

# Input 2: Additional features (IMPROVED - deeper with batch norm)
feature_input = layer_input(shape = ncol(features_train), name = "feature_input")

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

# Create model with dual inputs
model = keras_model(
  inputs = list(seq_input, feature_input),
  outputs = output
)

# Print model summary
summary(model)

# ===============================================================
# Compile model
# ===============================================================
model %>% compile(
  optimizer = optimizer_adam(learning_rate = 1e-4, clipnorm = 1.0),
  loss = "mse",
  metrics = c("mae")
)

# ===============================================================
# Train model
# ===============================================================
history = model %>% fit(
  x = list(x_training, features_train),
  y = y_train_scaled,
  validation_data = list(list(x_val, features_val), y_val_scaled),
  epochs = 40,
  batch_size = 32,
  callbacks = list(
    callback_early_stopping(
      monitor = "val_loss", 
      patience = 10, 
      restore_best_weights = TRUE
    ),
    callback_reduce_lr_on_plateau(
      monitor = "val_loss", 
      factor = 0.5, 
      patience = 5, 
      min_lr = 1e-6
    )
  )
)

# ===============================================================
# Evaluate on test set
# ===============================================================

# Clean data for prediction
x_test_clean = array(as.numeric(x_test), dim = dim(x_test))
features_test_clean = array(as.numeric(features_test), dim = dim(features_test))

# Make predictions
predictions_scaled = predict(model, list(x_test_clean, features_test_clean))

# Calculate metrics manually
y_pred = as.vector(predictions_scaled)
y_true = as.vector(y_test_scaled)

test_mse = mean((y_true - y_pred)^2)
test_mae = mean(abs(y_true - y_pred))

cat("\n========================================\n")
cat("TEST RESULTS\n")
cat("========================================\n")
cat("Test Loss (MSE):", test_mse, "\n")
cat("Test MAE:", test_mae, "\n\n")

# ===============================================================
# Transform back to original scale
# ===============================================================
predictions_log = predictions_scaled * y_sd + y_mean
predictions_original = expm1(predictions_log)

actual_log = y_test_scaled * y_sd + y_mean
actual_original = expm1(actual_log)

# ===============================================================
# Calculate correlations
# ===============================================================
pearson_cor = cor(actual_original, predictions_original)
spearman_cor = cor(actual_original, predictions_original, method = "spearman")

cat("CORRELATION RESULTS\n")
cat("========================================\n")
cat("Pearson correlation:", pearson_cor, "\n")
cat("Spearman correlation:", spearman_cor, "\n\n")

# ===============================================================
# Detailed analysis
# ===============================================================
results_df = data.frame(
  actual = actual_original,
  predicted = predictions_original,
  error = abs(actual_original - predictions_original)
)

# Overall performance
cat("PREDICTION SUMMARY\n")
cat("========================================\n")
cat("Mean actual half-life:", mean(actual_original), "hrs\n")
cat("Mean predicted half-life:", mean(predictions_original), "hrs\n")
cat("Median absolute error:", median(results_df$error), "hrs\n")
cat("Mean absolute error:", mean(results_df$error), "hrs\n\n")

# Performance by error ranges
cat("Predictions by error magnitude:\n")
cat("========================================\n")
cat("Excellent (error < 1 hr):", sum(results_df$error < 1), "samples (", 
    round(100*mean(results_df$error < 1), 1), "%)\n", sep="")
cat("Good (error 1-2 hrs):", sum(results_df$error >= 1 & results_df$error < 2), "samples (", 
    round(100*mean(results_df$error >= 1 & results_df$error < 2), 1), "%)\n", sep="")
cat("Fair (error 2-5 hrs):", sum(results_df$error >= 2 & results_df$error < 5), "samples (", 
    round(100*mean(results_df$error >= 2 & results_df$error < 5), 1), "%)\n", sep="")
cat("Poor (error >= 5 hrs):", sum(results_df$error >= 5), "samples (", 
    round(100*mean(results_df$error >= 5), 1), "%)\n\n", sep="")

# Performance on long-lived RNAs (> 15 hrs)
cat("PERFORMANCE ON LONG-LIVED RNAs (> 15 hrs):\n")
cat("========================================\n")
long_lived = results_df[results_df$actual > 15, ]
if(nrow(long_lived) > 0) {
  cat("Count:", nrow(long_lived), "\n")
  cat("Mean error:", round(mean(long_lived$error), 2), "hrs\n")
  cat("Median error:", round(median(long_lived$error), 2), "hrs\n")
  cat("Correlation:", round(cor(long_lived$actual, long_lived$predicted), 3), "\n")
  cat("Mean predicted:", round(mean(long_lived$predicted), 2), "hrs\n")
  cat("Mean actual:", round(mean(long_lived$actual), 2), "hrs\n\n")
}

# Save model weights and results
cat("Saving model weights and results...\n")
model$save_weights("rna_halflife_improved_weights.weights.h5")
saveRDS(results_df, "improved_model_predictions.rds")
saveRDS(list(
  weights_file = "rna_halflife_improved_weights.weights.h5",
  y_mean = y_mean,
  y_sd = y_sd,
  pearson_cor = pearson_cor,
  spearman_cor = spearman_cor,
  test_mse = test_mse,
  test_mae = test_mae,
  feature_architecture = "3 layers with batch normalization"
), "improved_model_metadata.rds")

cat("\nModel weights saved to: rna_halflife_improved_weights.weights.h5\n")
cat("Predictions saved to: improved_model_predictions.rds\n")
cat("Metadata saved to: improved_model_metadata.rds\n\n")