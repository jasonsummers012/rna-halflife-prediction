# ===============================================================
# Libraries
# ===============================================================
library(tidyverse)
library(keras3)
library(abind)
library(Biostrings)
library(seqinr)

# ===============================================================
# Load training data to get feature normalization parameters
# ===============================================================
training_data = readRDS("cnn_training_data_with_features.rds")
feature_means = training_data$feature_means
feature_sds = training_data$feature_sds
y_mean = training_data$y_mean
y_sd = training_data$y_sd

cat("Loaded normalization parameters from training data\n\n")

# ===============================================================
# Load mutant data
# ===============================================================
mutant_data = readRDS("caspase9_mutants.rds")

# Convert sequences to RNA
wildtype_rna = chartr("T", "U", mutant_data$wildtype)
wildtype_seq = strsplit(wildtype_rna, "")[[1]]

# ===============================================================
# Read pairing probabilities from ViennaRNA output
# ===============================================================
files = list.files("pairing_probs_mutants", pattern = "_lunp$", full.names = TRUE)

if(length(files) == 0) {
  stop("No pairing probability files found! Did you run ViennaRNA?")
}

pairing_data = map_dfr(files, function(f) {
  filename = basename(f) %>% str_remove("_lunp$")
  
  probs = read_table(f, comment = "#", col_names = c("pos", "unpaired"), 
                     show_col_types = FALSE) %>%
    mutate(
      pairing_prob = 1 - unpaired,
      seq_id = filename
    ) %>%
    # Remove rows with NA (often position 0 or header artifact)
    filter(!is.na(pairing_prob)) %>%
    select(seq_id, pos, pairing_prob)
  
  return(probs)
})

cat("Loaded pairing probabilities for", length(unique(pairing_data$seq_id)), "sequences\n\n")

# ===============================================================
# Organize pairing data
# ===============================================================
pairing_by_seq = pairing_data %>%
  group_by(seq_id) %>%
  summarize(pairing_probs = list(pairing_prob))

# Get wildtype pairing
wildtype_pairing_raw = pairing_by_seq %>% 
  filter(str_detect(seq_id, "wildtype")) %>% 
  pull(pairing_probs) %>% 
  .[[1]]

# ===============================================================
# Helper function to match and clean pairing data
# ===============================================================
match_and_clean_pairing = function(pairing_probs, sequence) {
  seq_length = length(sequence)
  prob_length = length(pairing_probs)
  
  # Remove any remaining NAs and replace with mean
  if(any(is.na(pairing_probs))) {
    mean_val = mean(pairing_probs, na.rm = TRUE)
    pairing_probs[is.na(pairing_probs)] = mean_val
  }
  
  # Match length
  if(prob_length > seq_length) {
    return(pairing_probs[1:seq_length])
  } else if(prob_length < seq_length) {
    mean_val = mean(pairing_probs, na.rm = TRUE)
    return(c(pairing_probs, rep(mean_val, seq_length - prob_length)))
  } else {
    return(pairing_probs)
  }
}

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

# Find longest ORF
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
  # Human codon frequencies
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

# Calculate all features for a given sequence
calculate_features = function(seq, pairing_probs) {
  seq_length = length(seq)
  gc_content = mean(seq %in% c("G", "C"))
  au_content = mean(seq %in% c("A", "U"))
  
  # ARE score (3' UTR focused)
  if(length(seq) < 500) {
    tail_seq = seq
  } else {
    tail_seq = tail(seq, 500)
  }
  seq_str = paste(tail_seq, collapse = "")
  auuua_count = str_count(seq_str, "AUUUA")
  are_score_3utr = auuua_count / length(tail_seq) * 1000
  
  # Poly(A) signal (with multiple variants)
  if(seq_length < 100) {
    polya_signal = 0
  } else {
    tail_seq = paste(tail(seq, 100), collapse = "")
    signals = c("AAUAAA", "AUUAAA", "AAUACA", "CAUAAA", 
                "AAUAUA", "AAUGAA", "ACUAAA", "AGUAAA")
    polya_signal = as.integer(any(sapply(signals, function(sig) str_detect(tail_seq, sig))))
  }
  
  # 3' UTR GC content
  if(seq_length < 500) {
    gc_3utr = gc_content
  } else {
    gc_3utr = mean(tail(seq, 500) %in% c("G", "C"))
  }
  
  # Pairing probabilities
  mean_pairing = mean(pairing_probs, na.rm = TRUE)
  
  if(length(pairing_probs) < 500) {
    pairing_3utr = mean_pairing
  } else {
    pairing_3utr = mean(tail(pairing_probs, 500), na.rm = TRUE)
  }
  
  # ORF detection
  orf_info = find_longest_orf(seq)
  orf_length = orf_info$length
  orf_ratio = orf_length / seq_length
  
  # CAI calculation
  cai = calculate_cai(seq, orf_info)
  
  # miRNA binding sites
  mirna_sites = sum(sapply(top_mirna_seeds, function(seed) {
    count_mirna_sites(seq, seed)
  }))
  mirna_density = (mirna_sites / seq_length) * 1000
  
  return(c(seq_length, gc_content, au_content, are_score_3utr, polya_signal,
           gc_3utr, mean_pairing, pairing_3utr, orf_length, orf_ratio,
           cai, mirna_sites, mirna_density))
}

# ===============================================================
# Calculate features for wildtype and mutants
# ===============================================================
cat("Calculating features for wildtype and mutants...\n")

# Clean wildtype pairing
wildtype_pairing = match_and_clean_pairing(wildtype_pairing_raw, wildtype_seq)
cat("Wildtype - seq length:", length(wildtype_seq), 
    ", pairing length:", length(wildtype_pairing), "\n")

# Wildtype features
wildtype_features = calculate_features(wildtype_seq, wildtype_pairing)

# Mutant features
mutant_features_list = list()
mutant_sequences_list = list()
mutant_pairing_list = list()

for(i in 1:length(mutant_data$mutants)) {
  # Convert to RNA
  mutant_rna = chartr("T", "U", mutant_data$mutants[[i]]$sequence)
  mutant_seq = strsplit(mutant_rna, "")[[1]]
  
  # Get pairing probs
  pattern = paste0("mutant_", i, "_")
  mutant_pairing_raw = pairing_by_seq %>% 
    filter(str_detect(seq_id, pattern)) %>% 
    pull(pairing_probs) %>% 
    .[[1]]
  
  # Clean and match mutant pairing
  mutant_pairing = match_and_clean_pairing(mutant_pairing_raw, mutant_seq)
  
  # Verify for first few mutants
  if(i <= 3) {
    cat("Mutant", i, "- seq length:", length(mutant_seq), 
        ", pairing length:", length(mutant_pairing), "\n")
  }
  
  # Calculate features
  features = calculate_features(mutant_seq, mutant_pairing)
  
  mutant_sequences_list[[i]] = mutant_seq
  mutant_pairing_list[[i]] = mutant_pairing
  mutant_features_list[[i]] = features
}

cat("Feature calculation complete!\n\n")

# ===============================================================
# One-hot encoding and padding functions
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

# ===============================================================
# Build tensors
# ===============================================================
max_len = 3500

# Wildtype tensor
wildtype_base_mat = encode_base(wildtype_seq)
wildtype_mat = cbind(wildtype_base_mat, matrix(wildtype_pairing, ncol = 1))
wildtype_padded = pad_matrix(wildtype_mat, max_len)

# Mutant tensors
mutant_padded_list = list()
for(i in 1:length(mutant_sequences_list)) {
  base_mat = encode_base(mutant_sequences_list[[i]])
  mat = cbind(base_mat, matrix(mutant_pairing_list[[i]], ncol = 1))
  mutant_padded_list[[i]] = pad_matrix(mat, max_len)
}

# Combine into array
all_padded = c(list(wildtype_padded), mutant_padded_list)
x_array = abind::abind(all_padded, along = 3)
x_array = aperm(x_array, c(3, 1, 2))

cat("Sequence tensors created\n")
cat("Shape:", dim(x_array), "\n\n")

# ===============================================================
# Prepare features
# ===============================================================
features_matrix = rbind(wildtype_features, do.call(rbind, mutant_features_list))

# Normalize using training data parameters
features_scaled = scale(features_matrix, 
                        center = feature_means, 
                        scale = feature_sds)

cat("Features normalized\n")
cat("Shape:", dim(features_scaled), "\n\n")

# ===============================================================
# Save prepared data
# ===============================================================
saveRDS(list(
  x_sequences = x_array,
  features = features_scaled,
  wildtype_halflife = mutant_data$wildtype_halflife,
  mutation_counts = c(0, mutant_data$mutation_counts),  # 0 for wildtype
  sequence_names = c("wildtype", paste0("mutant_", 1:length(mutant_data$mutants))),
  y_mean = y_mean,
  y_sd = y_sd,
  feature_names = training_data$feature_names
), "caspase9_prepared_for_prediction.rds")

cat("=========================================\n")
cat("CASPASE 9 DATA PREPARATION COMPLETE\n")
cat("=========================================\n")
cat("Total sequences prepared:", nrow(features_scaled), "\n")
cat("  - 1 wildtype\n")
cat("  - ", sum(mutant_data$mutation_counts == 10), " mutants (10 mutations)\n")
cat("  - ", sum(mutant_data$mutation_counts == 50), " mutants (50 mutations)\n")
cat("  - ", sum(mutant_data$mutation_counts == 100), " mutants (100 mutations)\n")
cat("\nSaved to: caspase9_prepared_for_prediction.rds\n")