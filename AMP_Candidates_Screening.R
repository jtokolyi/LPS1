################################
### AMP Candidates Screening ###
################################

library(Biostrings)

# loading trinotate file
trinotate <- read.delim(
  "C:/Users/Retno/Downloads/Trinotate_lps2_cdhit_new.xls",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# non-empty SignalP
signalp_nonempty <- trinotate[trinotate$SignalP != ".", ]

View(signalp_nonempty)

# TmHMM : empty
signalp_no_tm <- signalp_nonempty[
  signalp_nonempty$TmHMM == ".",
]

# aa_length >= 30 - <= 150
putative_amp <- signalp_no_tm[
  signalp_no_tm$aa_length >= 30 &
    signalp_no_tm$aa_length <= 150,
]

#save putative_AMP
saveRDS(
  putative_amp,
  "C:/Users/Retno/Downloads/putative_AMP.rds"
)

# merging putative AMP with DEG result 
genes_low <- genes_lowLPS_fix

genes_high <- `genes_highLPS_fix

# AMP candidates DE in low
AMP_low <- putative_amp[
  putative_amp$transcript_id %in% genes_low$transcript_id,
]


# AMP candidates DE in high
AMP_high <- putative_amp[
  putative_amp$transcript_id %in% genes_high$transcript_id,
]

# merging low LPS
AMP_low_DE <- merge(
  AMP_low,
  genes_low,
  by = "transcript_id",
  all.x = TRUE
)

# merging high LPS
AMP_high_DE <- merge(
  AMP_high,
  genes_high,
  by = "transcript_id",
  all.x = TRUE
)

# only upregulated and down-regulated
AMP_high_DE_only <- AMP_high_DE[
  AMP_high_DE$gene_type %in% c("Up-regulated", "Down-regulated"),
]

# choose only important columns
AMP_high_DE_only[, c(
  "transcript_id",
  "peptide",
  "aa_length",
  "SignalP",
  "TmHMM",
  "log2FoldChange",
  "pvalue",
  "padj",
  "gene_type",
  "GENE",
  "FullName"
)]

AMP_low_DE_only <- AMP_low_DE[
  AMP_low_DE$gene_type %in% c("Up-regulated", "Down-regulated"),
] ## zero

AMP_high_candidates <- AMP_high_DE_only[, c(
  "transcript_id",
  "X.gene_id",
  "peptide",
  "aa_length",
  "SignalP",
  "TmHMM",
  "log2FoldChange",
  "pvalue",
  "padj",
  "gene_type",
  "UniprotID",
  "GENE",
  "FullName"
)]

# load peptide fasta sequence file

pep_file <- "C:/Users/Retno/Downloads/lps2_cdhit.Trinity.fasta.transdecoder.pep"

readLines(pep_file, n = 10)

pep_sequences <- readAAStringSet(pep_file)

## extracting first field of fasta header
protein_ids_all <- sub(" .*", "", names(pep_sequences))

#find all matching protein prediction
candidate_proteins <- protein_ids_all[
  sub("\\.p[0-9]+$", "", protein_ids_all) %in%
    AMP_high_candidates$transcript_id
]

# count how many protein each transcript are

table(
  sub("\\.p[0-9]+$", "", candidate_proteins)
)


##Find missing transcript
id_candidates <- AMP_high_candidates$transcript_id

protein_transcript_ids <- sub(
  "\\.p[0-9]+$",
  "",
  candidate_proteins
)

missing_transcripts <- setdiff(
  id_candidates,
  protein_transcript_ids
)

missing_transcripts

# retrieve prot id for DE 26 candidates (4 of them duplicated)
AMP_prot_ids <- putative_amp[
  putative_amp$transcript_id %in% AMP_high_candidates$transcript_id,
  c("transcript_id", "prot_id", "aa_length")
]

# create 22 candidate table
AMP_high_unique <- AMP_high_candidates[
  !duplicated(AMP_high_candidates$transcript_id),
]

# retrieve prot_id from putative_amp
AMP_prot_ids <- putative_amp[
  putative_amp$transcript_id %in% AMP_high_unique$transcript_id,
  c("transcript_id", "prot_id", "aa_length") 
] # 23 protein prediction, 22 unique transcrip

# extract actual peptide sequences
AMP_prot_ids$peptide_sequence <- as.character(
  pep_sequences[match_id]
)

# remove the * in the end of the sequences
AMP_prot_ids$peptide_clean <- gsub(
  "\\*",
  "",
  AMP_prot_ids$peptide_sequence
)

AMP_prot_ids$sequence_length <- nchar(
  AMP_prot_ids$peptide_clean
)

#creating sequence_length column\
AMP_prot_ids$sequence_length <- nchar(
  AMP_prot_ids$peptide_clean
)

## clean the protein id
pep_check <- readAAStringSet(pep_file)
names(pep_check) <- sub(" .*", "", names(pep_check))

# create AAString set
peptides_to_blast <- AAStringSet(
  AMP_prot_ids$peptide_clean
)

#assign the protein ID
names(peptides_to_blast) <- AMP_prot_ids$prot_id

#save fasta
blastp_file <- "C:/Users/Retno/Downloads/AMP_high_LPS_candidates_BLASTp.fasta"

writeXStringSet(
  peptides_to_blast,
  blastp_file
)

## make a candidate table
AMP_final <- AMP_prot_ids[, c(
  "transcript_id",
  "prot_id",
  "aa_length",
  "peptide_clean"
)]

#Add DE information from AMP_high_candidates
AMP_final <- merge(
  AMP_final,
  unique(
    AMP_high_candidates[, c(
      "transcript_id",
      "gene_type",
      "log2FoldChange",
      "padj"
    )]
  ),
  by = "transcript_id",
  all.x = TRUE
)

AMP_final[, c(
  "transcript_id",
  "prot_id",
  "aa_length",
  "gene_type",
  "log2FoldChange",
  "padj"
)]

#save it as csv file
write.csv(
  AMP_final,
  "C:/Users/Retno/Downloads/AMP_final.csv",
  row.names = FALSE
)

#save as rds
saveRDS(
  AMP_final,
  "C:/Users/Retno/Downloads/AMP_final.rds"
)
