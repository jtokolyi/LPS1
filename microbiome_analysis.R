#####################################
library(phyloseq)
library(dplyr)
library(tibble)
library(readr)
library(DESeq2)
library(ggrepel)
library(ggplot2)
library(forcats)


# define input/output directories
data_dir <- "C:/Users/Retno/Downloads/lps_abundance"

output_dir <- "C:/Users/Retno/Downloads/lps_microbiome"

# get the TSV files
files <- list.files(
  data_dir,
  pattern = "\\.tsv$",
  full.names = TRUE
)

# read each TSV file
sample_data_list <- lapply(files, function(f) {
  
  dat <- read.delim(
    f,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
# keep the information needed
  dat <- dat %>%
    dplyr::select(
      tax_id,
      abundance,
      species,
      genus,
      family,
      order,
      class,
      phylum,
      clade,
      superkingdom,
      subspecies,
      `species subgroup`,
      `species group`,
      `estimated counts`
    )
  
# add sample name
  dat$SampleID <- tools::file_path_sans_ext(basename(f))
  
  dat
})

# combine all samples
all_data <- bind_rows(sample_data_list)


# create OTU/count table
# Rows = taxa
# Columns = samples
# Values = estimated counts

otu_df <- all_data %>%
  dplyr::select(tax_id, SampleID, `estimated counts`) %>%
  group_by(tax_id, SampleID) %>%
  summarise(
    Count = sum(`estimated counts`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = SampleID,
    values_from = Count,
    values_fill = 0
  )

# convert to data frame
otu_df <- as.data.frame(otu_df)

# use tax_id as row names
rownames(otu_df) <- as.character(otu_df$tax_id)

# remove tax_id column
otu_df$tax_id <- NULL

# convert to matrix
otu_mat <- as.matrix(otu_df)

# make sure columns are in the desired order
sample_order <- c(
  paste0("C", 1:5),
  paste0("H", 1:5),
  paste0("L", 1:5)
)

otu_mat <- otu_mat[, sample_order, drop = FALSE]

# create OTU table
OTU <- otu_table(
  otu_mat,
  taxa_are_rows = TRUE
)

# create taxonomy table
tax_df <- all_data %>%
  dplyr::select(
    tax_id,
    species,
    genus,
    family,
    order,
    class,
    phylum,
    clade,
    superkingdom
  ) %>%
  distinct(tax_id, .keep_all = TRUE)

# make tax_id the row names
rownames(tax_df) <- as.character(tax_df$tax_id)

# remove tax_id column
tax_df$tax_id <- NULL

# taxonomy follows OTU order
tax_df <- tax_df[rownames(otu_mat), , drop = FALSE]

# convert to matrix
tax_mat <- as.matrix(tax_df)

# phyloseq taxonomy table
TAX <- tax_table(tax_mat)

# create sample metadata
sample_metadata <- data.frame(
  SampleID = sample_order,
  
  Treatment = c(
    rep("Control", 5),
    rep("High", 5),
    rep("Low", 5)
  ),
  
  Replicate = rep(1:5, 3),
  
  row.names = sample_order,
  
  stringsAsFactors = FALSE
)

# make Treatment a factor
sample_metadata$Treatment <- factor(
  sample_metadata$Treatment,
  levels = c("Control", "High", "Low")
)

# Create phyloseq sample_data object
SAM <- sample_data(sample_metadata)

# create the unfiltered phyloseq object
ps <- phyloseq(
  OTU,
  TAX,
  SAM
)

# filter taxa based on prevalence
# Retain taxa present (> 0) in at least 2 samples
presence <- apply(
  as(otu_table(ps), "matrix"),
  1,
  function(x) sum(x > 0)
)

ps_counts_filtered <- prune_taxa(
  presence >= 2,
  ps
)

# save the two phyloseq objects
saveRDS(
  ps,
  file = file.path(
    output_dir,
    "phyloseq_object.rds"
  )
)

saveRDS(
  ps_counts_filtered,
  file = file.path(
    output_dir,
    "phyloseq_filtered.rds"
  )
)


######################################
### DESEq2 (Differential Abundance ###
######################################


# Phyloseq data : ps_counts_filtered
ps_counts_filtered <- readRDS("C:/Users/Retno/Downloads/lps_microbiome/phyloseq_filtered.rds")
ps <- readRDS("C:/Users/Retno/Downloads/lps_microbiome/phyloseq_object.rds")

# Run DESeq2
dds <- phyloseq_to_deseq2(ps_counts_filtered, ~ Treatment)
dds <- DESeq(dds)

# Low vs Control 
L_vs_C <- results(dds, contrast = c("Treatment", "Low", "Control"))

# Convert DESeq2 results into a tibble, keeping ASV/OTU IDs as a column
L_vs_C_df <- as(L_vs_C, "data.frame") %>%
  rownames_to_column("ASV")

# Convert taxonomy table into a tibble, also keeping ASV/OTU IDs
tax_df <- as.data.frame(tax_table(ps)) %>%
  rownames_to_column("ASV")

# Join DESeq2 results with taxonomy by ASV IDs
L_vs_C_tax <- left_join(L_vs_C_df, tax_df, by = "ASV")


# High vs Control 
H_vs_C <- results(dds, contrast = c("Treatment", "High", "Control"))

H_vs_C_df <- as(H_vs_C, "data.frame") %>%
  rownames_to_column("ASV")

H_vs_C_tax <- left_join(H_vs_C_df, tax_df, by = "ASV")


#add Comparison column
L_vs_C_tax$Comparison <- "Low_vs_Control"
H_vs_C_tax$Comparison <- "High_vs_Control"


#combine
res_combine <- rbind(L_vs_C_tax, H_vs_C_tax)
res_combine$Comparison <- factor(
  res_combine$Comparison,
  levels = c("Low_vs_Control", "High_vs_Control"),
  labels = c("Low vs Control", "High vs Control")
)

#save 
saveRDS(res_combine, file = "C:/Users/Retno/Downloads/lps_microbiome/res_combine.rds")

#add significance column to res_combine
res_combine$Significance <- ifelse(res_combine$padj < 0.05, "Significant", "Not Significant")
res_combine$Significance <- factor(res_combine$Significance, levels = c("Not Significant", "Significant"))

# Volcano Plot
# collapse categories and handle NA
# clean significance into exactly two classes
res_combine <- res_combine %>%
  mutate(
    Significance = ifelse(!is.na(padj) & padj < 0.05, "Significant", "Not Significant"),
    Significance = factor(Significance, levels = c("Not Significant", "Significant")),
    logP = ifelse(is.na(pvalue) | pvalue <= 0, NA_real_, -log10(pvalue))
  )

# re-ensure it's a factor with correct order 
res_combine$Significance <- factor(
  res_combine$Significance,
  levels = c("Not Significant", "Significant")
)

# volcano plot
ggplot(res_combine, aes(x = log2FoldChange, y = -log10(pvalue), color = Significance)) +
  geom_point(alpha = 0.7, size = 2) +
  
  # Add threshold lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  
  # Add labels only for significant points
  geom_text_repel(
    data = subset(res_combine, Significance == "Significant"),
    aes(label = species),
    size = 4.5,
    fontface = "italic",
    max.overlaps = 20
  ) +
  
  facet_wrap(~Comparison) +
  
  # Explicitly match colors to factor levels
  scale_color_manual(
    values = c("Not Significant" = "grey", "Significant" = "red"),
    drop = TRUE
  ) +
  
  # Set x-axis left limit
  coord_cartesian(xlim = c(-15, NA)) +  # NA keeps the upper limit automatic
  
  theme_bw() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.position = "bottom",
    
    # Axis text
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    
    # Facet title (Low vs Control, High vs Control)
    strip.text = element_text(size = 14, face = "bold"),
    
    # Legend text + title
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 15),
    
    axis.line = element_blank(),
    axis.ticks = element_line(),
    plot.caption = element_text(hjust = 0, size = 10)
  )+
  
  labs(
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Significance:"
  )

ggsave("C:/Users/Retno/Downloads/Volcano_DESeq2.png", width = 10, height = 5, dpi = 300)



##################################
### Pseudomonas Lollipop Chart ###
##################################

readRDS(res_combine, file = "C:/Users/Retno/Downloads/lps_microbiome/res_combine.rds")


# define ordered species (bottom -> top)
species_order <- c(
  "Pseudomonas sp. SWI44",
  "Pseudomonas sp. SWI6",
  "Pseudomonas plecoglossicida",
  "Pseudomonas sp. MYb193",
  "Pseudomonas fluorescens",
  "Pseudomonas putida",
  "Pseudomonas sp. MPC6",
  "Pseudomonas brassicacearum"
)

# filter + clean data
lollipop_data <- res_combine %>%
  filter(species %in% species_order,
         !is.na(log2FoldChange)) %>%
  
  mutate(
    species = factor(species, levels = species_order),
    
    # Standardize direction:
    # High vs Control → negative side (left)
    # Low vs Control → positive side (right)
    log2FC_plot = case_when(
      Comparison == "High vs Control" ~ -log2FoldChange,
      Comparison == "Low vs Control"  ~  log2FoldChange
    ),
    
    Significance = ifelse(padj < 0.05, "Significant", "Not Significant")
  )

# lollipop chart
ggplot(lollipop_data,
       aes(x = log2FoldChange,
           y = species)) +
  
  # Lollipop stems
  
  geom_segment(aes(x = 0,
                   xend = log2FoldChange,
                   y = species,
                   yend = species,
                   color = Significance),
               linewidth = 0.9) +
  
  # Dots
  geom_point(aes(color = Significance,
                 size = log10(baseMean))) +
  scale_size_continuous(
    range = c(3, 10),
    name = expression(log[10](baseMean))
  ) +
  
  # Make size legend dots grey
  guides(
    size = guide_legend(
      override.aes = list(color = "grey60")
    )
  ) +
  
  # Central reference line
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "black") +
  
  # Separate Low and High side-by-side
  facet_wrap(~Comparison, nrow = 1) +
  
  # Colors
  scale_color_manual(
    values = c(
      "Significant" = "red",
      "Not Significant" = "grey60"
    )
  ) +
  
  labs(
    x = expression(log[2]~(Fold~Change)),
    y = NULL,
    color = "FDR < 0.05"
  ) +
  
  theme_bw(base_size = 13) +
  
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(
      face = "italic",
      size = 14,
      color = "black")
  )

ggsave("C:/Users/Retno/Downloads/Lolli_Pseudomonas_3.png", width = 13, height = 4, dpi = 300)

ggsave("C:/Users/Retno/Downloads/Lolli_Pseudomonas_4.png", width = 9, height = 5, dpi = 300)









