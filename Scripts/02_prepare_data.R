# ===============================================================
# Libraries
# ===============================================================
library(tidyverse)
library(keras3)
library(abind)
library(rsample)
library(Biostrings)
library(seqinr)

# ===============================================================
# Load grouped data
# ===============================================================
gene_data = readRDS("cnn_input_grouped.rds") %>%
  filter(!is.na(gene_id))

# ===============================================================
# FILTER BY SEQUENCE LENGTH
# ===============================================================
gene_data = gene_data %>%
  mutate(seq_len = map_int(sequence, length)) %>%
  filter(seq_len >= 500, seq_len <= 5000)

cat("Filtered to", nrow(gene_data), "genes (500-5000 bp)\n")
cat("Length distribution after filtering:\n")
print(summary(map_int(gene_data$sequence, length)))
cat("\n")

# Remove helper column before continuing
gene_data = gene_data %>% select(-seq_len)

# ===============================================================
# Top 20 most abundant human miRNA seed sequences (7-mers)
# ===============================================================
top_mirna_seeds = c(
  "UGAGUGU",  # let-7 family
  "AAGGUGC",  # miR-21
  "AUCACAG",  # miR-10
  "GCACUUU",  # miR-125
  "UACCCUA",  # miR-181
  "CCAAUUA",  # miR-26
  "AGGCAGU",  # miR-143
  "UGGCGUU",  # miR-16
  "GGACAGG",  # miR-15
  "UCAGCUU",  # miR-30
  "CUGCUGU",  # miR-200
  "ACCCGUA",  # miR-155
  "AUGCAAU",  # miR-146
  "UGUACAG",  # miR-145
  "UCUCAGG",  # miR-17
  "AUAAGCU",  # miR-34
  "AACAUUC",  # miR-29
  "GCGCCUU",  # miR-103
  "GGGUGGU",  # miR-23
  "CCGAGAG"   # miR-221
)

# ===============================================================
# Helper functions
# ===============================================================

# Find longest orf
find_longest_orf = function(sequence) {
  seq_str = paste(sequence, collapse = "")
  dna_str = chartr("U", "T", seq_str)
  
  max_length = 0
  best_start = NA
  best_end = NA
  
  # Check all 3 reading frames
  for(frame in 0:2) {
    # Look for start codons (ATG) in this frame
    starts = seq(frame + 1, nchar(dna_str) - 2, by = 3)
    
    for(start_pos in starts) {
      codon = substr(dna_str, start_pos, start_pos + 2)
      if(codon != "ATG") next  # Only start from ATG
      
      # Extend until stop codon
      for(i in seq(start_pos, nchar(dna_str) - 2, by = 3)) {
        codon = substr(dna_str, i, i + 2)
        if(codon %in% c("TAA", "TAG", "TGA")) {
          orf_length = i - start_pos + 3
          if(orf_length > max_length) {
            max_length = orf_length
            best_start = start_pos
            best_end = i + 2
          }
          break
        }
      }
    }
  }
  
  return(list(start = best_start, end = best_end, length = max_length))
}

# Create human codon weights from the provided frequency table
create_human_w = function() {
  # Human codon frequencies from your table
  codon_freqs = c(
    # TTx
    "TTT" = 0.46, "TTC" = 0.54, "TTA" = 0.08, "TTG" = 0.13,
    # TCx
    "TCT" = 0.19, "TCC" = 0.22, "TCA" = 0.15, "TCG" = 0.05,
    # TATx
    "TAT" = 0.44, "TAC" = 0.56, "TAA" = 0.30, "TAG" = 0.24,
    # TGx
    "TGT" = 0.46, "TGC" = 0.54, "TGA" = 0.47, "TGG" = 1.00,
    
    # CTx
    "CTT" = 0.13, "CTC" = 0.20, "CTA" = 0.07, "CTG" = 0.40,
    # CCx
    "CCT" = 0.29, "CCC" = 0.32, "CCA" = 0.28, "CCG" = 0.11,
    # CAx
    "CAT" = 0.42, "CAC" = 0.58, "CAA" = 0.27, "CAG" = 0.73,
    # CGx
    "CGT" = 0.08, "CGC" = 0.18, "CGA" = 0.11, "CGG" = 0.20,
    
    # ATx
    "ATT" = 0.36, "ATC" = 0.47, "ATA" = 0.17, "ATG" = 1.00,
    # ACx
    "ACT" = 0.25, "ACC" = 0.36, "ACA" = 0.28, "ACG" = 0.11,
    # AAx
    "AAT" = 0.47, "AAC" = 0.53, "AAA" = 0.43, "AAG" = 0.57,
    # AGx
    "AGT" = 0.15, "AGC" = 0.24, "AGA" = 0.21, "AGG" = 0.21,
    
    # GTx
    "GTT" = 0.18, "GTC" = 0.24, "GTA" = 0.12, "GTG" = 0.46,
    # GCx
    "GCT" = 0.27, "GCC" = 0.40, "GCA" = 0.23, "GCG" = 0.11,
    # GAx
    "GAT" = 0.46, "GAC" = 0.54, "GAA" = 0.42, "GAG" = 0.58,
    # GGx
    "GGT" = 0.16, "GGC" = 0.34, "GGA" = 0.25, "GGG" = 0.25
  )
  
  # Calculate relative adaptiveness (w) for each codon
  # w = frequency of codon / frequency of most common codon for that amino acid
  
  # Genetic code grouping
  aa_groups = list(
    Phe = c("TTT", "TTC"),
    Leu = c("TTA", "TTG", "CTT", "CTC", "CTA", "CTG"),
    Ser = c("TCT", "TCC", "TCA", "TCG", "AGT", "AGC"),
    Tyr = c("TAT", "TAC"),
    Stop = c("TAA", "TAG", "TGA"),
    Cys = c("TGT", "TGC"),
    Trp = c("TGG"),
    Pro = c("CCT", "CCC", "CCA", "CCG"),
    His = c("CAT", "CAC"),
    Gln = c("CAA", "CAG"),
    Arg = c("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"),
    Ile = c("ATT", "ATC", "ATA"),
    Met = c("ATG"),
    Thr = c("ACT", "ACC", "ACA", "ACG"),
    Asn = c("AAT", "AAC"),
    Lys = c("AAA", "AAG"),
    Val = c("GTT", "GTC", "GTA", "GTG"),
    Ala = c("GCT", "GCC", "GCA", "GCG"),
    Asp = c("GAT", "GAC"),
    Glu = c("GAA", "GAG"),
    Gly = c("GGT", "GGC", "GGA", "GGG")
  )
  
  # Calculate w values
  w = numeric(length(codon_freqs))
  names(w) = names(codon_freqs)
  
  for (aa in names(aa_groups)) {
    codons = aa_groups[[aa]]
    freqs = codon_freqs[codons]
    max_freq = max(freqs)
    
    # Calculate relative adaptiveness
    w[codons] = freqs / max_freq
  }
  
  # Stop codons should be 0 for CAI calculation
  w[c("TAA", "TAG", "TGA")] = 0
  
  return(w)
}

# Create human w values once (outside the function for efficiency)
w_human = create_human_w()

# Calculate codon adaptation index using seqinr package
calculate_cai = function(sequence, orf_info, w_values = w_human) {
  if(is.na(orf_info$start) || orf_info$length < 150) return(0)
  
  seq_str = paste(sequence, collapse = "")
  dna_str = chartr("U", "T", seq_str)
  
  # Extract CDS
  cds = substr(dna_str, orf_info$start, orf_info$end)
  
  # Convert to vector of single characters for seqinr
  cds_vector = s2c(cds)
  
  # Calculate CAI using seqinr's cai function
  tryCatch({
    cai_value = cai(cds_vector, w = w_values, numcode = 1)
    return(cai_value)
  }, error = function(e) {
    # If CAI calculation fails (e.g., very short sequence), return 0
    return(0)
  })
}

# Count miRNA binding sites
count_mirna_sites = function(sequence, seed) {
  seq_str = paste(sequence, collapse = "")
  # Count overlapping matches
  matches = str_count(seq_str, paste0("(?=", seed, ")"))
  return(matches)
}

cat("Calculating advanced features...\n")

# ===============================================================
# Calculate all features
# ===============================================================
gene_data = gene_data %>%
  mutate(
    # Basic features
    seq_length = map_int(sequence, length),
    gc_content = map_dbl(sequence, ~mean(. %in% c("G", "C"))),
    au_content = map_dbl(sequence, ~mean(. %in% c("A", "U"))),
    
    # ARE score
    are_score_3utr = map_dbl(sequence, function(seq) {
      if(length(seq) < 500) {
        tail_seq = seq
      } else {
        tail_seq = tail(seq, 500)
      }
      seq_str = paste(tail_seq, collapse = "")
      auuua_count = str_count(seq_str, "AUUUA")
      auuua_count / length(tail_seq) * 1000
    }),
    
    # Poly(A) signal
    polya_signal = map_int(sequence, function(seq) {
      if(length(seq) < 100) return(0L)
      tail_seq = paste(tail(seq, 100), collapse = "")
      signals = c("AAUAAA", "AUUAAA", "AAUACA", "CAUAAA", 
                  "AAUAUA", "AAUGAA", "ACUAAA", "AGUAAA")
      
      as.integer(any(sapply(signals, function(sig) str_detect(tail_seq, sig))))  # ← Fixed
    }),
    
    # 3' UTR features
    gc_3utr = map_dbl(sequence, function(seq) {
      if(length(seq) < 500) {
        mean(seq %in% c("G", "C"))
      } else {
        mean(tail(seq, 500) %in% c("G", "C"))
      }
    }),
    
    # Pairing features
    mean_pairing = map_dbl(pairing_probs, mean),
    pairing_3utr = map_dbl(pairing_probs, function(probs) {
      if(length(probs) < 500) {
        mean(probs)
      } else {
        mean(tail(probs, 500))
      }
    }),
    
    # ORF detection
    orf_info = map(sequence, find_longest_orf),
    orf_length = map_dbl(orf_info, ~.x$length),
    orf_ratio = orf_length / seq_length,
    
    # Codon features
    cai = map2_dbl(sequence, orf_info, calculate_cai),
    
    # miRNA binding sites (total count)
    mirna_sites = map_int(sequence, function(seq) {
      sum(sapply(top_mirna_seeds, function(seed) {
        count_mirna_sites(seq, seed)
      }))
    }),
    
    # miRNA density per 1000 bases
    mirna_density = (mirna_sites / seq_length) * 1000
  ) %>%
  select(-orf_info)  # Remove list column

cat("Feature calculation complete!\n\n")
cat("Feature summary:\n")
print(summary(gene_data %>% select(seq_length, gc_content, are_score_3utr, 
                                   orf_length, cai, mirna_sites)))

# ===============================================================
# One-hot encoding for RNA bases
# ===============================================================
encode_base = function(bases) {
  bases = toupper(bases)
  valid = c("A", "C", "G", "U")
  onehot = matrix(0, nrow = length(bases), ncol = 4)
  colnames(onehot) = valid
  for (i in seq_along(bases)) {
    if (bases[i] %in% valid) {
      onehot[i, bases[i]] = 1
    }
  }
  onehot
}

# ===============================================================
# Build tensors for each gene
# ===============================================================
tensor_input = gene_data %>%
  mutate(x = map2(sequence, pairing_probs, function(seq, probs) {
    base_mat = encode_base(seq)
    cbind(base_mat, matrix(probs, ncol = 1))
  })) %>%
  pull(x)

# ===============================================================
# Padding function
# ===============================================================
pad_matrix = function(mat, max_len) {
  n = nrow(mat)
  
  if (n > max_len) {
    mat2 = mat[(n - max_len + 1):n, , drop = FALSE]
    mask = rep(1, max_len)
  } else if (n < max_len) {
    pad = matrix(0, nrow = max_len - n, ncol = ncol(mat))
    mat2 = rbind(pad, mat)
    mask = c(rep(0, max_len - n), rep(1, n))
  } else {
    mat2 = mat
    mask = rep(1, max_len)
  }
  
  cbind(mat2, mask)
}

max_len = 3500
padded_list = lapply(tensor_input, pad_matrix, max_len = max_len)

# ===============================================================
# Combine into 3D tensor
# ===============================================================
x_array = abind::abind(padded_list, along = 3)
x_array = aperm(x_array, c(3, 1, 2))

# ===============================================================
# Train/val/test split
# ===============================================================
set.seed(12)

initial = initial_split(gene_data, prop = 0.7)
training_data = training(initial)
temp_data = testing(initial)

val_split = initial_split(temp_data, prop = 0.5)
val_data = training(val_split)
test_data = testing(val_split)

y_train = training_data$half_life
y_val = val_data$half_life
y_test = test_data$half_life

training_idx = as.integer(rownames(training_data))
val_idx = as.integer(rownames(val_data))
test_idx = as.integer(rownames(test_data))

x_training = x_array[training_idx, , ]
x_val = x_array[val_idx, , ]
x_test = x_array[test_idx, , ]

# ===============================================================
# Extract features
# ===============================================================
feature_cols = c("seq_length", "gc_content", "au_content", "are_score_3utr", 
                 "polya_signal", "gc_3utr", "mean_pairing", "pairing_3utr",
                 "orf_length", "orf_ratio", "cai", "mirna_sites", "mirna_density")

features_train = training_data %>% select(all_of(feature_cols)) %>% as.matrix()
features_val = val_data %>% select(all_of(feature_cols)) %>% as.matrix()
features_test = test_data %>% select(all_of(feature_cols)) %>% as.matrix()

# Normalize features
feature_means = colMeans(features_train)
feature_sds = apply(features_train, 2, sd)

features_train_scaled = scale(features_train, center = feature_means, scale = feature_sds)
features_val_scaled = scale(features_val, center = feature_means, scale = feature_sds)
features_test_scaled = scale(features_test, center = feature_means, scale = feature_sds)

# ===============================================================
# Normalize targets
# ===============================================================
y_train_log = log1p(y_train)
y_val_log = log1p(y_val)
y_test_log = log1p(y_test)

y_mean = mean(y_train_log)
y_sd = sd(y_train_log)

y_train_scaled = (y_train_log - y_mean) / y_sd
y_val_scaled = (y_val_log - y_mean) / y_sd
y_test_scaled = (y_test_log - y_mean) / y_sd

# ===============================================================
# Save
# ===============================================================
saveRDS(list(
  x_training = x_training,
  x_val = x_val,
  x_test = x_test,
  features_train = features_train_scaled,
  features_val = features_val_scaled,
  features_test = features_test_scaled,
  feature_names = feature_cols,
  feature_means = feature_means,
  feature_sds = feature_sds,
  y_training_scaled = y_train_scaled,
  y_val_scaled = y_val_scaled,
  y_test_scaled = y_test_scaled,
  y_mean = y_mean,
  y_sd = y_sd
), "cnn_training_data_with_features.rds")

cat("\n========================================\n")
cat("ADVANCED DATA PREPARATION COMPLETE\n")
cat("========================================\n")
cat("Sequence data shape:", dim(x_training), "\n")
cat("Features shape:", dim(features_train_scaled), "\n")
cat("Number of features:", length(feature_cols), "\n")
cat("Features included:\n")
cat(paste("  -", feature_cols), sep = "\n")
cat("\nSaved to: cnn_training_data_advanced.rds\n")

# Load training data
gene_data = readRDS("cnn_input_grouped.rds") %>%
  filter(!is.na(gene_id))

# Check for NAs in pairing probs
has_na = map_lgl(gene_data$pairing_probs, ~any(is.na(.)))
cat("Genes with NA pairing probs:", sum(has_na), "out of", nrow(gene_data), "\n")