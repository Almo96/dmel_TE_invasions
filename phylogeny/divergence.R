library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
blast_output_file <- args[1]
names_file <- args[2]
TE_length <- as.numeric(args[3])
output <- args[4]

blast_output <- read_tsv(blast_output_file, col_names = c("fasta_name", "familyname", "identity", "tot_length", "a_length", "mismatch", "gap_open", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
names <- read_tsv(names_file, col_names = c("name", "fasta_name"))
blast <- inner_join(names, blast_output, by="fasta_name")

identity <- blast %>% filter(bitscore>1000) %>% mutate(bases_identical = (identity/100)*a_length) %>%
  group_by(name, tot_length) %>%
  summarise(tot_bases_identical = round((sum(bases_identical)),2)) %>%
  mutate(identity = round((tot_bases_identical/TE_length),3))

write_tsv(identity, output)

# 6 qseqid sseqid pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore