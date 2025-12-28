library(tidyverse)
library(Biostrings)

# ===============================================================
# Load Caspase 9 sequence
# ===============================================================
gene_data = readRDS("cnn_input_grouped.rds") %>%
  filter(!is.na(gene_id))

casp9 = gene_data %>% filter(gene_id == "NM_001229")

casp9_seq = casp9$sequence[[1]]
casp9_str = paste(casp9_seq, collapse = "")

cat("Caspase 9 length:", length(casp9_seq), "bases\n")
cat("Half-life:", casp9$half_life, "hours\n\n")

# ===============================================================
# Find CDS region (longest ORF)
# ===============================================================
find_orf = function(seq_str) {
  dna_str = chartr("U", "T", seq_str)
  starts = str_locate_all(dna_str, "ATG")[[1]][,1]
  
  best = list(start = NA, end = NA, length = 0)
  
  for(start_pos in starts) {
    for(i in seq(start_pos, nchar(dna_str)-2, by = 3)) {
      codon = substr(dna_str, i, i+2)
      if(codon %in% c("TAA", "TAG", "TGA")) {
        orf_len = i - start_pos + 3
        if(orf_len > best$length) {
          best = list(start = start_pos, end = i+2, length = orf_len)
        }
        break
      }
    }
  }
  return(best)
}

orf = find_orf(casp9_str)
cat("ORF found:\n")
cat("  Start:", orf$start, "\n")
cat("  End:", orf$end, "\n")
cat("  Length:", orf$length, "bp\n\n")

# ===============================================================
# Define mutation region (300 bp in middle of CDS)
# ===============================================================
mutation_start = orf$start + floor(orf$length / 2) - 50
mutation_end = mutation_start + 299

# Adjust to start at codon boundary within CDS
cds_offset = mutation_start - orf$start  # How far into CDS
remainder = cds_offset %% 3

if(remainder != 0) {
  mutation_start = mutation_start + (3 - remainder)
  mutation_end = mutation_start + 299
  cat("Adjusted to codon boundary:", mutation_start, "-", mutation_end, "\n")
}

mutation_region = substr(casp9_str, mutation_start, mutation_end)
cat(mutation_region, "\n")

cat("Mutation region:", mutation_start, "-", mutation_end, "\n")
cat("Region sequence:\n")
mutation_region = substr(casp9_str, mutation_start, mutation_end)
cat(mutation_region, "\n\n")

# ===============================================================
# Genetic code (RNA codons)
# ===============================================================
genetic_code = list(
  "UUU"="F", "UUC"="F", "UUA"="L", "UUG"="L",
  "UCU"="S", "UCC"="S", "UCA"="S", "UCG"="S",
  "UAU"="Y", "UAC"="Y", "UAA"="*", "UAG"="*",
  "UGU"="C", "UGC"="C", "UGA"="*", "UGG"="W",
  "CUU"="L", "CUC"="L", "CUA"="L", "CUG"="L",
  "CCU"="P", "CCC"="P", "CCA"="P", "CCG"="P",
  "CAU"="H", "CAC"="H", "CAA"="Q", "CAG"="Q",
  "CGU"="R", "CGC"="R", "CGA"="R", "CGG"="R",
  "AUU"="I", "AUC"="I", "AUA"="I", "AUG"="M",
  "ACU"="T", "ACC"="T", "ACA"="T", "ACG"="T",
  "AAU"="N", "AAC"="N", "AAA"="K", "AAG"="K",
  "AGU"="S", "AGC"="S", "AGA"="R", "AGG"="R",
  "GUU"="V", "GUC"="V", "GUA"="V", "GUG"="V",
  "GCU"="A", "GCC"="A", "GCA"="A", "GCG"="A",
  "GAU"="D", "GAC"="D", "GAA"="E", "GAG"="E",
  "GGU"="G", "GGC"="G", "GGA"="G", "GGG"="G"
)

# Get synonymous codons for each amino acid
synonymous_codons = list()
for(aa in unique(unlist(genetic_code))) {
  synonymous_codons[[aa]] = names(genetic_code)[genetic_code == aa]
}

# ===============================================================
# Create silent mutations
# ===============================================================
create_silent_mutation = function(seq_str, orf_start, mut_start, mut_end, n_mutations = 5) {
  
  # Extract CDS
  dna_str = chartr("U", "T", seq_str)
  
  # Positions within the mutation region (in CDS frame)
  region_in_cds_start = mut_start - orf_start + 1
  region_in_cds_end = mut_end - orf_start + 1
  
  # Make sure we're in frame
  if(region_in_cds_start %% 3 != 1) {
    region_in_cds_start = region_in_cds_start - ((region_in_cds_start - 1) %% 3)
  }
  
  # Get codons in region
  codon_positions = seq(region_in_cds_start, region_in_cds_end - 2, by = 3)
  
  # Randomly select codons to mutate
  selected_codons = sample(codon_positions, min(n_mutations, length(codon_positions)))
  
  mutated_seq = seq_str
  mutations_made = list()
  
  for(codon_pos in selected_codons) {
    # Get original codon (in full sequence coordinates)
    full_pos = orf_start + codon_pos - 1
    original_codon = substr(mutated_seq, full_pos, full_pos + 2)
    original_codon_rna = chartr("T", "U", original_codon)
    
    # Get amino acid
    aa = genetic_code[[original_codon_rna]]
    
    if(is.null(aa) || aa == "*") next  # Skip if stop codon
    
    # Get synonymous codons (excluding the original)
    syn_codons = setdiff(synonymous_codons[[aa]], original_codon_rna)
    
    if(length(syn_codons) == 0) next  # Skip if no synonyms
    
    # Pick a random synonymous codon
    new_codon_rna = sample(syn_codons, 1)
    new_codon_dna = chartr("U", "T", new_codon_rna)
    
    # Make the mutation
    substr(mutated_seq, full_pos, full_pos + 2) = new_codon_dna
    
    mutations_made[[length(mutations_made) + 1]] = list(
      position = full_pos,
      original = original_codon,
      mutated = new_codon_dna,
      amino_acid = aa
    )
  }
  
  return(list(
    sequence = mutated_seq,
    mutations = mutations_made
  ))
}

# ===============================================================
# Generate multiple mutants
# ===============================================================
set.seed(12)
mutants = list()
mutant_counter = 1

# 25 mutants with 10 mutations each
for(i in 1:25) {
  mutant = create_silent_mutation(casp9_str, orf$start, mutation_start, mutation_end, n_mutations = 10)
  
  cat("\n=== Mutant", mutant_counter, "(10 mutations) ===\n")
  cat("Mutations made:", length(mutant$mutations), "\n")
  for(mut in mutant$mutations) {
    cat(sprintf("  Pos %d: %s -> %s (AA: %s)\n", 
                mut$position, mut$original, mut$mutated, mut$amino_acid))
  }
  
  mutants[[mutant_counter]] = mutant
  mutant_counter = mutant_counter + 1
}

# 25 mutants with 50 mutations each
for(i in 1:25) {
  mutant = create_silent_mutation(casp9_str, orf$start, mutation_start, mutation_end, n_mutations = 50)
  
  cat("\n=== Mutant", mutant_counter, "(50 mutations) ===\n")
  cat("Mutations made:", length(mutant$mutations), "\n")
  for(mut in mutant$mutations) {
    cat(sprintf("  Pos %d: %s -> %s (AA: %s)\n", 
                mut$position, mut$original, mut$mutated, mut$amino_acid))
  }
  
  mutants[[mutant_counter]] = mutant
  mutant_counter = mutant_counter + 1
}

# 25 mutants with 100 mutations each
for(i in 1:25) {
  mutant = create_silent_mutation(casp9_str, orf$start, mutation_start, mutation_end, n_mutations = 100)
  
  cat("\n=== Mutant", mutant_counter, "(100 mutations) ===\n")
  cat("Mutations made:", length(mutant$mutations), "\n")
  for(mut in mutant$mutations) {
    cat(sprintf("  Pos %d: %s -> %s (AA: %s)\n", 
                mut$position, mut$original, mut$mutated, mut$amino_acid))
  }
  
  mutants[[mutant_counter]] = mutant
  mutant_counter = mutant_counter + 1
}

cat("\n\nTotal mutants created:", length(mutants), "\n")

# ===============================================================
# Save mutants
# ===============================================================
saveRDS(list(
  wildtype = casp9_str,
  wildtype_halflife = casp9$half_life,
  orf = orf,
  mutation_region = c(start = mutation_start, end = mutation_end),
  mutants = mutants,
  n_mutants = length(mutants),
  mutation_counts = c(rep(10, 25), rep(50, 25), rep(100, 25))  # Track which has how many
), "caspase9_mutants.rds")

cat("\n\nSaved", length(mutants), "mutants to caspase9_mutants.rds\n")
cat("  - 25 mutants with 10 mutations each\n")
cat("  - 25 mutants with 50 mutations each\n")
cat("  - 25 mutants with 100 mutations each\n")

# ===============================================================
# Load mutants
# ===============================================================
mutant_data = readRDS("caspase9_mutants.rds")

# ===============================================================
# Create FASTA file with wildtype + all mutants
# ===============================================================
fasta_lines = c()

# Add wildtype
fasta_lines = c(fasta_lines, 
                 ">NM_001229_wildtype",
                 mutant_data$wildtype)

# Add all mutants
for(i in 1:length(mutant_data$mutants)) {
  n_muts = mutant_data$mutation_counts[i]
  fasta_lines = c(fasta_lines,
                   paste0(">NM_001229_mutant_", i, "_", n_muts, "muts"),
                   mutant_data$mutants[[i]]$sequence)
}

# Write FASTA file
write_lines(fasta_lines, "caspase9_all_sequences.fa")

cat("Created FASTA file with", 1 + length(mutant_data$mutants), "sequences\n")
cat("Saved to: caspase9_all_sequences.fa\n\n")