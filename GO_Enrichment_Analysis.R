#########################################
### Gene Ontology enrichment analysis ###
#########################################


library(topGO); 
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)

# load trinotate file 
trin <- read.delim(
  "C:/Users/Retno/Downloads/Trinotate_lps2_cdhit.xls",
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Clean quotation marks from column names
colnames(trin) <- gsub('"', '', colnames(trin))

# Remove empty columns
trin <- trin[, colnames(trin) != "", drop = FALSE]

# gene2go
gene2go <- trin %>%
  dplyr::select(gene_id = transcript_id, GO_raw = gene_ontology_BLASTP) %>%
  filter(!is.na(GO_raw) & GO_raw != "") %>%
  mutate(GO_list = strsplit(as.character(GO_raw), "`")) %>%
  unnest(GO_list) %>%
  mutate(GO = sub("\\^.*", "", GO_list)) %>%
  dplyr::select(gene_id, GO) %>%
  filter(GO!=".")

# Load The Deseq2 analysis result file 
# High LPS : genes_highLPS_fix
# Low LPS : genes_lowLPS_fix

# Change the file name 
resL <- genes_lowLPS_fix
resH <- genes_highLPS_fix

################################
### GO Enrichment in Low LPS ###
################################

# extract adjusted p-value
gene_pvals_L <- resL$padj
names(gene_pvals_L) <- resL$transcript_id
gene_pvals_L <- gene_pvals_L[!is.na(gene_pvals_L)]
gene_universe_L <- names(gene_pvals_L)

# Upregulated genes #
up_genes_L <- resL$transcript_id[which(resL$log2FoldChange > 0 & resL$padj < 0.05)]
up_gene_list_L <- factor(as.integer(gene_universe_L %in% up_genes_L))
names(up_gene_list_L) <- gene_universe_L

GOdata_up_L <- new("topGOdata",
                   ontology = "BP",
                   allGenes = up_gene_list_L,
                   annot = annFUN.gene2GO,
                   gene2GO = split(gene2go$GO, gene2go$gene_id))

resultFisher_up_L <- runTest(GOdata_up_L, algorithm = "weight01", statistic = "fisher")

# Downregulated genes #
down_genes_L <- resL$transcript_id[which(resL$log2FoldChange < 0 & resL$padj < 0.05)]
down_gene_list_L <- factor(as.integer(gene_universe_L %in% down_genes_L))
names(down_gene_list_L) <- gene_universe_L

GOdata_down_L <- new("topGOdata",
                     ontology = "BP",
                     allGenes = down_gene_list_L,
                     annot = annFUN.gene2GO,
                     gene2GO = split(gene2go$GO, gene2go$gene_id))

resultFisher_down_L <- runTest(GOdata_down_L, algorithm = "weight01", statistic = "fisher")

# Prepare two tables (down- and up-regulated)
go_up_tableL <- GenTable(GOdata_up_L, weightFisher = resultFisher_up_L, orderBy = "weightFisher", topNodes = 30)
go_down_tableL <- GenTable(GOdata_down_L, weightFisher = resultFisher_down_L, orderBy = "weightFisher", topNodes = 30)

# adding other columns
go_up_tableL <- go_up_tableL %>%
  mutate(
    weightFisher = as.numeric(weightFisher),
    Significant = as.numeric(Significant),
    Expected = as.numeric(Expected),
    FoldEnrichment = Significant / Expected,
    log10p = -log10(weightFisher),
    Direction = "Upregulated"
  )

go_down_tableL <- go_down_tableL %>%
  mutate(
    weightFisher = as.numeric(weightFisher),
    Significant = as.numeric(Significant),
    Expected = as.numeric(Expected),
    FoldEnrichment = Significant / Expected,
    log10p = -log10(weightFisher),
    Direction = "Downregulated"
  )

# combine the table
go_combinedL <- bind_rows(go_up_tableL, go_down_tableL)

# only choose top10
go_up_tableL_10 <- go_up_tableL[1:10, ]
go_down_tableL_10 <- go_down_tableL[1:10, ]
go_combinedL_10 <- bind_rows(go_up_tableL_10, go_down_tableL_10)

# filter by threshold
go_combinedL_10_filtered <- go_combinedL_10 %>%
  filter((Direction == "Upregulated" & log10p >= 0.6) |
           (Direction == "Downregulated" & log10p <= -0.6))

# create stacked order (y-axis)
up_termsL <- go_combinedL_10_filtered %>%
  filter(Direction == "Upregulated") %>%
  arrange(desc(log10p)) %>%
  pull(Term)

down_termsL <- go_combinedL_10_filtered %>%
  filter(Direction == "Downregulated") %>%
  arrange(log10p) %>%  # most negative first → largest absolute first
  pull(Term)

stacked_order_L <- c(up_termsL, down_termsL)

# apply factor levels
go_combinedL_10_filtered <- go_combinedL_10_filtered %>%
  mutate(Term_ordered = factor(Term, levels = stacked_order_L))

# save it
write.csv(go_combinedL_10_filtered, 
          "C:/Users/Retno/Downloads/lps_microbiome/GO_terms_L_new.csv",   # file name
          row.names = FALSE) 

# Lollipop Plot

ggplot(go_L, aes(y = Term_ordered, x = log10p, color = Direction)) +
  # Lollipop lines
  geom_segment(aes(x = 0, xend = log10p, y = Term_ordered, yend = Term_ordered), size = 1) +
  # Lollipop points
  geom_point(aes(fill = Direction), size = 3, shape = 21, stroke = 0.5) +
  # Dashed vertical lines for significance threshold
  geom_vline(xintercept = c(-1.3, 1.3), linetype = "dashed", color = "grey60") +
  # X-axis
  scale_x_continuous(
    name = "-log10(p-value) (diverging)",
    labels = function(x) abs(x),       # keep the labels positive if you want
    breaks = seq(-5, 5, by = 1),       # set ticks every 1 unit from -5 to 5
    limits = c(-5, 5)                 # optional: match the limits to the breaks
  ) +
  # Y-axis
  scale_y_discrete(expand = expansion(mult = c(0.005, 0.02))) +
  # Colors
  scale_fill_manual(values = c("Upregulated" = "#CA0020", "Downregulated" = "#0571B0")) +
  scale_color_manual(values = c("Upregulated" = "#CA0020", "Downregulated" = "#0571B0")) +
  # Labels and theme
  labs(y = "GO Term", title = "Enriched GO Terms in Low LPS") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 11),
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor.y = element_blank()
  )

ggsave("C:/Users/Retno/Downloads/lps_microbiome/GO_LowLPS2.png", width = 7, height = 2, dpi = 300)


#################################
### GO Enrichment in High LPS ###
#################################

gene_pvals_H <- resH$padj
names(gene_pvals_H) <- resH$transcript_id
gene_pvals_H <- gene_pvals_H[!is.na(gene_pvals_H)]
gene_universe_H <- names(gene_pvals_H)

# Upregulated genes
up_genes_H <- resH$transcript_id[which(resH$log2FoldChange > 0 & resH$padj < 0.05)]
up_gene_list_H <- factor(as.integer(gene_universe_H %in% up_genes_H))
names(up_gene_list_H) <- gene_universe_H

GOdata_up_H <- new("topGOdata",
                   ontology = "BP",
                   allGenes = up_gene_list_H,
                   annot = annFUN.gene2GO,
                   gene2GO = split(gene2go$GO, gene2go$gene_id))

resultFisher_up_H <- runTest(GOdata_up_H, algorithm = "weight01", statistic = "fisher")

# Downregulated genes 
down_genes_H <- resH$transcript_id[which(resH$log2FoldChange < 0 & resH$padj < 0.05)]
down_gene_list_H <- factor(as.integer(gene_universe_H %in% down_genes_H))
names(down_gene_list_H) <- gene_universe_H

GOdata_down_H <- new("topGOdata",
                     ontology = "BP",
                     allGenes = down_gene_list_H,
                     annot = annFUN.gene2GO,
                     gene2GO = split(gene2go$GO, gene2go$gene_id))

resultFisher_down_H <- runTest(GOdata_down_H, algorithm = "weight01", statistic = "fisher")

## Prepare two tables
go_up_tableH <- GenTable(GOdata_up_H, weightFisher = resultFisher_up_H, orderBy = "weightFisher", topNodes = 30)
go_down_tableH <- GenTable(GOdata_down_H, weightFisher = resultFisher_down_H, orderBy = "weightFisher", topNodes = 30)

# adding other columns

go_up_tableH <- go_up_tableH %>%
  mutate(
    weightFisher = as.numeric(weightFisher),
    Significant = as.numeric(Significant),
    Expected = as.numeric(Expected),
    FoldEnrichment = Significant / Expected,
    log10p = -log10(weightFisher),
    Direction = "Upregulated"
  )

go_down_tableH <- go_down_tableH %>%
  mutate(
    weightFisher = as.numeric(weightFisher),
    Significant = as.numeric(Significant),
    Expected = as.numeric(Expected),
    FoldEnrichment = Significant / Expected,
    log10p = -log10(weightFisher),
    Direction = "Downregulated"
  )

# combine the table
go_combinedH <- bind_rows(go_up_tableH, go_down_tableH)

## choose top15
go_up_tableH_15 <- go_up_tableH[1:15, ]
go_down_tableH_15 <- go_down_tableH[1:15, ]
go_combinedH_15 <- bind_rows(go_up_tableH_15, go_down_tableH_15)



# create diverging x-axis
go_combinedH_15 <- GO_combinedH_15 %>%
  mutate(
    # Make negative for downregulated
    log10p = ifelse(Direction == "Downregulated", -log10p, log10p)
  )

# create stacked order (y-axis)
# Upregulated first (top), order by largest log10p_div
up_terms_H <- go_combinedH_15 %>%
  filter(Direction == "Upregulated") %>%
  arrange(desc(log10p)) %>%
  pull(Term)

# downregulated next (bottom), order by largest absolute log10p_div
down_terms_H <- go_combinedH_15 %>%
  filter(Direction == "Downregulated") %>%
  arrange(log10p) %>%  # most negative first → largest absolute first
  pull(Term)

# combine order
stacked_order_H <- c(up_terms_H, down_terms_H)

# apply as factor
go_combinedH_15 <- go_combinedH_15 %>%
  mutate(Term_ordered = factor(Term, levels = stacked_order_H))

# plot
ggplot(go_combinedH_15, aes(y = Term_ordered, x = log10p, color = Direction)) +
  # Lollipop lines
  geom_segment(aes(x = 0, xend = log10p, y = Term_ordered, yend = Term_ordered), size = 1) +
  # Lollipop points
  geom_point(aes(fill = Direction), size = 3, shape = 21, stroke = 0.5) +
  # Dashed vertical lines for significance threshold
  geom_vline(xintercept = c(-1.3, 1.3), linetype = "dashed", color = "grey60") +
  # X-axis
  scale_x_continuous(
    name = "-log10(p-value) (diverging)",
    labels = function(x) abs(x),       # keep the labels positive if you want
    breaks = seq(-5, 5, by = 1),       # set ticks every 1 unit from -5 to 5
    limits = c(-5, 5)                 # optional: match the limits to the breaks
  ) +
  # Y-axis
  scale_y_discrete() +
  # Colors
  scale_fill_manual(values = c("Upregulated" = "#CA0020", "Downregulated" = "#0571B0")) +
  scale_color_manual(values = c("Upregulated" = "#CA0020", "Downregulated" = "#0571B0")) +
  # Labels and theme
  labs(y = "GO Term", title = "Enriched GO Terms in High LPS") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 11),
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor.y = element_blank()
  )

ggsave("C:/Users/Retno/Downloads/lps_microbiome/GO_HighLPS3.png", width = 7.2, height = 8, dpi = 300)

######################################################################