##############################
### Reading Phylogeny Tree ###
##############################

# loading required libraries
library(ape)
library(ggtree)
library(dplyr)


# reading  IQ-TREE output
tree_file <- "C:/Users/Retno/Downloads/phylogeny/pseudomonas_aligned.fasta.treefile"  # replace with your actual file path
tree <- read.tree(tree_file)


# root the tree using outgroup
tree_rooted <- root(tree, outgroup = "Cellvibrio_japonicus_Ueda_107", resolve.root = TRUE)


# convert bootstrap labels to numeric

tree_rooted$node.label <- as.numeric(tree_rooted$node.label)
labels <- tree_rooted$tip.label


# define groups
red_vec  <- c("Pseudomonas_plecoglossicida", "Pseudomonas_sp._SWI6", "Pseudomonas_sp._SWI44")
blue_vec <- c("Pseudomonas_sp._MPC6", "Pseudomonas_brassicaserum", "Pseudomonas_fluorescens",
              "Pseudomonas_putida", "Pseudomonas_migulae", "Pseudomonas_alcaligenes",
              "Pseudomonas_alcaliphila_JCM_10630", "Pseudomonas_sp._MYb193")
green_vec <- c("Pseudomonas_laurentiana_X11/14_1", "Pseudomonas_kielensis_M26/9/10_0")
species_groups <- data.frame(label = labels, stringsAsFactors = FALSE) %>%
  mutate(group = case_when(
    label %in% red_vec  ~ "Downregulated",
    label %in% blue_vec ~ "Not_significant",
    label %in% green_vec ~ "Antibacterial_activity_strain",
    label == "Cellvibrio_japonicus_Ueda_107"    ~ "Outgroup",
    TRUE ~ "Type_strain"
  ))


# clean bootstrap values: keep only >= 70
tree_rooted$node.label <- ifelse(
  as.numeric(tree_rooted$node.label) >= 70,
  tree_rooted$node.label,
  NA
)

# find maximum branch length
max_depth <- max(node.depth.edgelength(tree_rooted))


# plot the tree
p <- ggtree(tree_rooted, layout = "rectangular") %<+% species_groups +
  geom_tiplab(aes(color = group), size = 5, fontface = "italic") +
  geom_text2(aes(subset = !isTip & !is.na(label), label = label),
             vjust = -0.5, size = 3, color = "gray40") +
  scale_color_manual(values = c(
    "Downregulated" = "#FF0800",
    "Not_significant" = "#0000FF",
    "Antibacterial_activity_strain" = "#568203",
    "Outgroup" = "orange",
    "Type_strain" = "black"
  )) +
  xlim(0, max_depth * 1.25) +  # <-- give extra space
  theme(legend.position = "none") 

# save figure

ggsave("C:/Users/Retno/Downloads/phylogeny/Pseudomonas_tree.png", plot = p,
       width = 14, height = 10, dpi = 500)










