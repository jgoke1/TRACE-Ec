library(tidyverse)
library(tidyverse)
library(igraph)
library(ggraph)


# First column is sample ID
snp_mat <- read.delim("snp_dist_matrix.tsv")
rownames(snp_mat) <- snp_mat[,1]
snp_mat <- snp_mat[,-1]


# Read SNP distance matrix, skipping the first "snp-dists 0.8.2" line
snp_mat <- read.table("snp_dist_matrix.tsv", 
                      header = TRUE, 
                      check.names = FALSE, 
                      sep = "\t", 
                      skip = 1, 
                      row.names = 1)

# Check first few rows
head(snp_mat[, 1:5])


# Convert to long format
snp_long <- snp_mat %>%
  as.data.frame() %>%
  rownames_to_column("from") %>%
  pivot_longer(-from, names_to = "to", values_to = "distance") %>%
  filter(from != to)

# Keep only one direction
snp_long <- snp_long %>%
  rowwise() %>%
  mutate(pair = paste(sort(c(from, to)), collapse = "_")) %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(-pair)

# Apply SNP threshold <= 40
edges <- snp_long %>% filter(distance <= 40)


# Load metadata
nodes <- read.csv("mlst.csv")

# Ensure 'name' column in nodes matches sample IDs
nodes <- nodes %>% rename(name = ERR) %>% mutate(name = as.character(name))

# Filter edges to keep only ones where both nodes exist in metadata
edges <- edges %>%
  filter(from %in% nodes$name & to %in% nodes$name)


library(igraph)
library(ggraph)

# Create graph
g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)

# Define colors
source_colors <- c(
  "Human" = "#e41a1c",
  "Animal" = "#ffcc00",
  "Environment" = "#377eb8"
)

# Plot
ggraph(g, layout = "fr") +
  geom_edge_link(alpha = 0.3, color = "grey50") +
  geom_node_point(aes(color = Source), size = 4) +
  geom_node_text(aes(label = ST), repel = TRUE, size = 3) +
  scale_color_manual(values = source_colors) +
  theme_void() +
  theme(legend.position = "bottom")
