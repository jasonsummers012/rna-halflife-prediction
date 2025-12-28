# ===============================================================
# Libraries
# ===============================================================
library(tidyverse)
library(readr)
library(Biostrings)


# ===============================================================
# Load and prepare main dataset
# ===============================================================
data = read_csv("rna_half_lives.csv")
ids = unique(data$`Gene ID`)


# ===============================================================
# Load and clean FASTA sequences
# ===============================================================
fasta = readDNAStringSet("all_human_rna.fna")
names(fasta) = sub(" .*", "", names(fasta))
names(fasta) = sub("\\..*", "", names(fasta))

# Filter dataset to include only genes with sequences
data = data %>%
  filter(`Gene ID` %in% names(fasta))

# Add RNA sequences
data$sequence = as.character(fasta[match(data$`Gene ID`, names(fasta))])
data$sequence = chartr("T", "U", data$sequence)


# ===============================================================
# Calculate mean half-life across samples
# ===============================================================
half_life_cols = data %>%
  select(`GM07029-A1` : GM12815)

data = data %>%
  mutate(half_life_mean = rowMeans(half_life_cols, na.rm = TRUE)) %>%
  select(`Gene ID`, sequence, half_life_mean)


# ===============================================================
# Write FASTA file for ViennaRNA input
# ===============================================================
write_lines(
  paste0(">", data$`Gene ID`, "\n", data$sequence),
  "filtered_genes.fa"
)


# ===============================================================
# Combine ViennaRNA output files (_lunp)
# ===============================================================
files = list.files("pairing_probs", pattern = "_lunp$", full.names = TRUE)

pairing_probs = map_dfr(files, function(f) {
  filename = basename(f) %>% str_remove("_lunp$")
  gene_id = str_extract(filename, "N[MR]_[0-9]+")
  
  read_table(f, comment = "#", col_names = c("pos", "unpaired")) %>%
    mutate(gene_id = gene_id, pairing_prob = 1 - unpaired)
}) %>%
  filter(!is.na(pos)) %>%
  select(-unpaired)


# ===============================================================
# Expand sequence into per-base positions
# ===============================================================
sequence_expanded = data %>%
  select(`Gene ID`, sequence, half_life_mean) %>%
  mutate(base = map(sequence, ~ strsplit(.x, "")[[1]])) %>%
  unnest(base) %>%
  group_by(`Gene ID`) %>%
  mutate(pos = row_number()) %>%
  ungroup()


# ===============================================================
# Merge pairing probabilities with sequence and metadata
# ===============================================================
merged = pairing_probs %>%
  left_join(sequence_expanded,
            by = c("gene_id" = "Gene ID", "pos" = "pos")) %>%
  select(gene_id, pos, base, pairing_prob, half_life_mean)


# ===============================================================
# Save final CNN input file
# ===============================================================
write_csv(merged, "cnn_input_with_pairing.csv")


# ===============================================================
# Save final R object
# ===============================================================
cnn_input_nested = merged %>%
  group_by(gene_id) %>%
  summarize(
    sequence = list(base),
    pairing_probs = list(pairing_prob),
    half_life = unique(half_life_mean)
  )

saveRDS(cnn_input_nested, "cnn_input_grouped.rds")