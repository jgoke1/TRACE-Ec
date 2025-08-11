# Clear the R environment
rm(list = ls())

# Load libraries
library(tidyverse)
library(plotly)
library(colorspace)
library(readxl)
library(forcats)
library(ggrepel)
library(networkD3)
library(htmlwidgets)
library(viridis)
library(circlize)

# Disable scientific notation
options(scipen = 999)

# Read data
main_data <- read_csv("/Users/seqafrica_ghru/atlas_antibiotics-subset.csv")

# Add continent column
main_data <- main_data %>%
  mutate(continent = case_when(
    Country %in% c("France", "Spain", "Belgium", "Italy", "Germany", "Ireland", "Portugal", 
                   "Greece", "United Kingdom", "Poland", "Switzerland", "Hungary", "Austria", 
                   "Finland", "Denmark", "Sweden", "Croatia", "Czech Republic", "Netherlands", 
                   "Russia", "Romania", "Latvia", "Lithuania", "Serbia", "Ukraine", "Slovenia", 
                   "Bulgaria", "Norway", "Slovak Republic", "Estonia") ~ "Europe",
    Country %in% c("Canada", "United States", "Mexico", "Guatemala", "Dominican Republic", 
                   "Costa Rica", "Jamaica", "Honduras", "Puerto Rico", "Nicaragua", "Panama", 
                   "El Salvador") ~ "North America",
    Country %in% c("Colombia", "Chile", "Venezuela", "Argentina", "Brazil") ~ "South America",
    Country %in% c("Australia", "New Zealand") ~ "Australia",
    Country %in% c("China", "Hong Kong", "Japan", "Malaysia", "Thailand", "Philippines", 
                   "Korea, South", "Taiwan", "India", "Singapore", "Vietnam", "Indonesia") ~ "Asia",
    Country %in% c("Nigeria", "Kenya", "South Africa", "Ivory Coast", "Morocco", "Cameroon", 
                   "Malawi", "Uganda", "Ghana", "Namibia", "Mauritius", "Tunisia", "Egypt") ~ "Africa",
    Country %in% c("Israel", "Kuwait", "Turkey", "Jordan", "Saudi Arabia", "Pakistan", 
                   "Lebanon", "Qatar", "Oman") ~ "Middle East",
    TRUE ~ NA_character_
  ))

# Select columns
main_data <- main_data %>%
  select(
    # Metadata
    "Isolate Id", "Study", "Species", "Family", "Country", "continent", "State", "Gender", 
    "Age Group", "Speciality", "Source", "In / Out Patient", "Year", "Phenotype",
    # Beta-lactam antibiotics
    "Amoxycillin clavulanate", "Amoxycillin clavulanate_I",
    "Ampicillin", "Ampicillin_I",
    "Penicillin", "Penicillin_I",
    "Piperacillin tazobactam", "Piperacillin tazobactam_I",
    "Ampicillin sulbactam", "Ampicillin sulbactam_I",
    "Aztreonam", "Aztreonam_I",
    "Aztreonam avibactam", "Aztreonam avibactam_I",
    "Cefepime", "Cefepime_I",
    "Cefoxitin", "Cefoxitin_I",
    "Ceftazidime", "Ceftazidime_I",
    "Ceftriaxone", "Ceftriaxone_I",
    "Cefixime", "Cefixime_I",
    "Ceftaroline", "Ceftaroline_I",
    "Ceftaroline avibactam", "Ceftaroline avibactam_I",
    "Ceftazidime avibactam", "Ceftazidime avibactam_I",
    "Ceftolozane tazobactam", "Ceftolozane tazobactam_I",
    "Cefoperazone sulbactam", "Cefoperazone sulbactam_I",
    "Meropenem", "Meropenem_I",
    "Meropenem vaborbactam", "Meropenem vaborbactam_I",
    "Doripenem", "Doripenem_I",
    "Imipenem", "Imipenem_I",
    "Ertapenem", "Ertapenem_I",
    "Oxacillin", "Oxacillin_I",
    "Cefpodoxime", "Cefpodoxime_I",
    "Ceftibuten", "Ceftibuten_I",
    "Ceftibuten avibactam", "Ceftibuten avibactam_I",
    "Tebipenem", "Tebipenem_I",
    "FOX",
    # Resistance genes
    "AMPC", "SHV", "TEM", "CTXM1", "CTXM2", "CTXM825", "CTXM9",
    "VEB", "PER", "GES", "ACC", "CMY1MOX", "CMY11", "DHA",
    "ACTMIR", "KPC", "OXA", "NDM", "IMP", "VIM", "SPM", "GIM"
  ) %>%
  select(where(~ !all(is.na(.))))

# Define metadata and gene columns
metadata_cols <- c("Isolate Id", "Study", "Species", "Family", "Country", 
                   "Gender", "Age Group", "Speciality", "Source", 
                   "In / Out Patient", "Year", "Phenotype", "continent")
all_genes <- c("AMPC", "SHV", "TEM", "CTXM1", "CTXM2", "CTXM825", "CTXM9", 
               "VEB", "PER", "GES", "ACC", "CMY1MOX", "CMY11", "DHA", "FOX", 
               "ACTMIR", "KPC", "OXA", "NDM", "IMP", "VIM", "SPM", "GIM")
gene_cols_present <- all_genes[all_genes %in% names(main_data)]
all_other_cols <- setdiff(colnames(main_data), c(metadata_cols, all_genes))
mic_cols <- all_other_cols[!str_ends(all_other_cols, "_I")]
interp_cols <- all_other_cols[str_ends(all_other_cols, "_I")]

# Debug: Check available genes
cat("Genes in main_data:", paste(colnames(main_data)[colnames(main_data) %in% all_genes], collapse = ", "), "\n")

# Create long-format data
mic_long <- main_data %>%
  mutate(across(all_of(mic_cols), as.character)) %>%
  select(all_of(metadata_cols), all_of(mic_cols)) %>%
  pivot_longer(
    cols = all_of(mic_cols),
    names_to = "antibiotic",
    values_to = "MIC"
  ) %>%
  filter(!is.na(MIC) & MIC != "")

interp_long <- main_data %>%
  mutate(across(all_of(interp_cols), as.character)) %>%
  select(all_of(metadata_cols), all_of(interp_cols)) %>%
  pivot_longer(
    cols = all_of(interp_cols),
    names_to = "antibiotic_I",
    values_to = "Interpretation"
  ) %>%
  mutate(antibiotic = str_remove(antibiotic_I, "_I")) %>%
  select(-antibiotic_I) %>%
  filter(!is.na(Interpretation) & Interpretation != "")

gene_long <- main_data %>%
  mutate(across(all_of(gene_cols_present), as.character)) %>%
  select(all_of(metadata_cols), all_of(gene_cols_present)) %>%
  pivot_longer(
    cols = all_of(gene_cols_present),
    names_to = "Gene",
    values_to = "Presence"
  ) %>%
  filter(!is.na(Presence) & Presence != "")

# Debug: Check genes in gene_long
cat("Genes in gene_long:", paste(unique(gene_long$Gene), collapse = ", "), "\n")

mic_interp_joined <- mic_long %>%
  left_join(interp_long, by = c(metadata_cols, "antibiotic"))

# Debug: Check columns in mic_interp_joined
cat("Columns in mic_interp_joined:", paste(colnames(mic_interp_joined), collapse = ", "), "\n")

gene_wide <- gene_long %>%
  pivot_wider(
    names_from = Gene,
    values_from = Presence,
    values_fill = ""
  )

full_data <- mic_interp_joined %>%
  left_join(gene_wide, by = metadata_cols)

# Merge gene_long with mic_interp_joined
mic_interp_gene_long <- mic_interp_joined %>%
  left_join(gene_long, by = metadata_cols) %>%
  filter(!is.na(Presence)) # Optional: Remove rows where gene presence is NA if desired

# Debug: Check columns in merged dataset
cat("Columns in mic_interp_gene_long:", paste(colnames(mic_interp_gene_long), collapse = ", "), "\n")

# Export data
dir.create("data_exports", showWarnings = FALSE)
write_csv(mic_interp_joined, "data_exports/mic_interp_joined.csv")
write_csv(interp_long, "data_exports/interp_long.csv")
write_csv(mic_long, "data_exports/mic_long.csv")
write_csv(gene_long, "data_exports/gene_long.csv")
write_csv(mic_interp_gene_long, "data_exports/mic_interp_gene_long.csv") # Export merged dataset

# 1. Gene Prevalence Bar Plot
gene_freq <- gene_long %>%
  filter(Gene %in% gene_cols_present) %>%
  mutate(is_positive = !str_ends(Presence, "-Neg")) %>%
  group_by(Gene) %>%
  summarise(
    total_isolates = n_distinct(`Isolate Id`),
    positive_isolates = n_distinct(`Isolate Id`[is_positive]),
    prevalence = positive_isolates / total_isolates * 100
  ) %>%
  arrange(desc(prevalence))

cat("Gene prevalence summary:\n")
print(gene_freq)

p <- ggplot(gene_freq, aes(x = reorder(Gene, prevalence), y = prevalence, fill = prevalence)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_c(option = "magma", name = "Prevalence (%)") +
  coord_flip() +
  labs(title = "Prevalence of Resistance Genes", x = "Gene", y = "Prevalence (%)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "#333333"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

print(p)
ggsave("gene_prevalence.png", p, width = 8, height = 6, dpi = 300)
cat("Generated and saved gene prevalence bar plot\n")

# 2. Sankey Diagram (Beta-Lactam Resistance with Metadata)
beta_lactam_antibiotics <- c(
  "Amoxycillin clavulanate", "Ampicillin", "Penicillin", "Piperacillin tazobactam",
  "Ampicillin sulbactam", "Aztreonam", "Aztreonam avibactam", "Cefepime",
  "Cefoxitin", "Ceftazidime", "Ceftriaxone", "Cefixime", "Ceftaroline",
  "Ceftaroline avibactam", "Ceftazidime avibactam", "Ceftolozane tazobactam",
  "Cefoperazone sulbactam", "Meropenem", "Meropenem vaborbactam", "Doripenem",
  "Imipenem", "Ertapenem", "Oxacillin", "Cefpodoxime", "Ceftibuten",
  "Ceftibuten avibactam", "Tebipenem"
)

base_data_sankey <- gene_long %>%
  filter(Gene %in% gene_cols_present, !str_ends(Presence, "-Neg")) %>%
  select(`Isolate Id`, Gene, continent, Source, Gender, `In / Out Patient`) %>%
  left_join(
    mic_interp_joined %>%
      filter(antibiotic %in% beta_lactam_antibiotics) %>%
      select(`Isolate Id`, continent, Interpretation),
    by = c("Isolate Id", "continent")
  ) %>%
  filter(!is.na(Interpretation)) %>%
  distinct()

continents <- unique(base_data_sankey$continent[!is.na(base_data_sankey$continent)])

create_sankey <- function(continent_name, metadata_col, metadata_name) {
  data_subset <- base_data_sankey %>%
    filter(continent == continent_name) %>%
    select(`Isolate Id`, Gene, !!sym(metadata_col), Interpretation) %>%
    group_by(Gene, !!sym(metadata_col), Interpretation) %>%
    summarise(count = n_distinct(`Isolate Id`)) %>%
    ungroup() %>%
    rename(Metadata = !!sym(metadata_col))
  
  if (nrow(data_subset) == 0) {
    cat("No data available for", continent_name, "with", metadata_name, "\n")
    return(NULL)
  }
  
  cat("Data for", continent_name, "with", metadata_name, ": ", nrow(data_subset), "rows\n")
  print(head(data_subset))
  
  nodes <- data.frame(name = unique(c(data_subset$Gene, data_subset$Metadata, data_subset$Interpretation)))
  links <- data_subset %>%
    group_by(Gene, Metadata) %>%
    summarise(value = sum(count)) %>%
    ungroup() %>%
    mutate(
      source = match(Gene, nodes$name) - 1,
      target = match(Metadata, nodes$name) - 1
    ) %>%
    select(source, target, value) %>%
    bind_rows(
      data_subset %>%
        group_by(Metadata, Interpretation) %>%
        summarise(value = sum(count)) %>%
        ungroup() %>%
        mutate(
          source = match(Metadata, nodes$name) - 1,
          target = match(Interpretation, nodes$name) - 1
        ) %>%
        select(source, target, value)
    )
  
  if (nrow(links) == 0) {
    cat("No valid links for", continent_name, "with", metadata_name, "\n")
    return(NULL)
  }
  
  node_colors <- c(
    viridis(length(unique(data_subset$Gene)), option = "magma", begin = 0, end = 0.4),
    viridis(length(unique(data_subset$Metadata)), option = "magma", begin = 0.5, end = 0.7),
    rep("#666666", length(unique(data_subset$Interpretation)))
  )
  
  sankey_plot <- sankeyNetwork(
    Links = links,
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "name",
    units = "isolates",
    fontSize = 12,
    nodeWidth = 30,
    fontFamily = "Arial",
    iterations = 0,
    colourScale = htmlwidgets::JS(sprintf("d3.scaleOrdinal().range(['%s'])", paste(node_colors, collapse = "','")))
  )
  
  title_html <- htmltools::tags$h3(
    paste("Resistance Genes, ", metadata_name, ", and Beta-Lactam Resistance in ", continent_name, sep = ""),
    style = "text-align: center; font-family: Arial; color: #333;"
  )
  sankey_plot <- htmlwidgets::prependContent(sankey_plot, title_html)
  
  print(sankey_plot)
  saveWidget(
    sankey_plot,
    file = paste0("sankey_", metadata_name, "_", gsub(" ", "_", continent_name), ".html"),
    title = paste("Resistance Genes, ", metadata_name, ", and Beta-Lactam Resistance in ", continent_name, sep = ""),
    selfcontained = TRUE
  )
  
  cat("Generated and saved Sankey diagram for", continent_name, "with", metadata_name, "\n")
  return(sankey_plot)
}

metadata_cols <- c("Source", "Gender", "In / Out Patient")
metadata_names <- c("Source", "Gender", "In_Out_Patient")
for (i in seq_along(metadata_cols)) {
  cat("\nGenerating Sankey diagrams for", metadata_names[i], "\n")
  for (cont in continents) {
    cat("Processing", cont, "\n")
    create_sankey(cont, metadata_cols[i], metadata_names[i])
  }
}

# 3. Chord Diagram (Gene Co-occurrence)
gene_wide <- gene_long %>%
  filter(Gene %in% gene_cols_present) %>%
  mutate(is_positive = !str_ends(Presence, "-Neg")) %>%
  select(`Isolate Id`, Gene, is_positive) %>%
  pivot_wider(
    names_from = Gene,
    values_from = is_positive,
    values_fill = FALSE
  )

# Debug: Check available genes in gene_wide
cat("Genes in gene_wide:", paste(colnames(gene_wide)[-1], collapse = ", "), "\n")

co_occurrence <- gene_wide %>%
  select(all_of(gene_cols_present)) %>%
  as.matrix() %>%
  crossprod()
diag(co_occurrence) <- 0

colors <- viridis(length(gene_cols_present), option = "magma")
grid_col <- setNames(colors, gene_cols_present)

circos.clear()
circos.par(
  gap.after = 10,
  start.degree = 90,
  clock.wise = TRUE,
  canvas.xlim = c(-1.1, 1.1),
  canvas.ylim = c(-1.1, 1.1)
)

chordDiagram(
  co_occurrence,
  grid.col = grid_col,
  transparency = 0.3,
  annotationTrack = c("name", "grid"),
  preAllocateTracks = list(
    list(track.height = 0.1),
    list(track.height = 0.05)
  ),
  annotationTrackHeight = c(0.05, 0.05),
  link.lwd = 2,
  link.lty = 1,
  link.border = "black",
  big.gap = 20
)

title(
  "Co-occurrence of Resistance Genes Across Isolates",
  cex.main = 2,
  col.main = "#333333",
  font.main = 2,
  adj = 0.5,
  line = -1
)

circos.clear()
cat("Generated chord diagram\n")

# Save chord diagram as PNG
png("chord_diagram.png", width = 800, height = 800)
circos.par(
  gap.after = 10,
  start.degree = 90,
  clock.wise = TRUE,
  canvas.xlim = c(-1.1, 1.1),
  canvas.ylim = c(-1.1, 1.1)
)
chordDiagram(
  co_occurrence,
  grid.col = grid_col,
  transparency = 0.3,
  annotationTrack = c("name", "grid"),
  preAllocateTracks = list(
    list(track.height = 0.1),
    list(track.height = 0.05)
  ),
  annotationTrackHeight = c(0.05, 0.05),
  link.lwd = 2,
  link.lty = 1,
  link.border = "black",
  big.gap = 20
)
title(
  "Co-occurrence of Resistance Genes Across Isolates",
  cex.main = 2,
  col.main = "#333333",
  font.main = 2,
  adj = 0.5,
  line = -1
)
circos.clear()
dev.off()
cat("Saved chord diagram as chord_diagram.png\n")

# 4. Heatmap (Beta-Lactam Resistance by Age Group)
continents <- unique(mic_interp_joined$continent[!is.na(mic_interp_joined$continent)])

create_heatmap <- function(continent_name) {
  heatmap_data <- mic_interp_joined %>%
    filter(continent == continent_name, antibiotic %in% beta_lactam_antibiotics) %>%
    filter(!is.na(Interpretation), Interpretation != "", !is.na(`Age Group`), `Age Group` != "") %>%
    group_by(antibiotic, `Age Group`) %>%
    summarise(
      total_isolates = n_distinct(`Isolate Id`),
      resistant_isolates = n_distinct(`Isolate Id`[Interpretation == "Resistant"]),
      proportion_resistant = resistant_isolates / total_isolates
    ) %>%
    ungroup() %>%
    complete(antibiotic = beta_lactam_antibiotics, `Age Group`, fill = list(total_isolates = 0, resistant_isolates = 0, proportion_resistant = 0))
  
  cat("Heatmap data summary for", continent_name, ": ", nrow(heatmap_data), "rows\n")
  print(head(heatmap_data))
  
  if (nrow(heatmap_data) == 0 || all(heatmap_data$total_isolates == 0)) {
    cat("No valid data for", continent_name, "\n")
    return(NULL)
  }
  
  p <- ggplot(heatmap_data, aes(x = `Age Group`, y = antibiotic, fill = proportion_resistant)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_viridis(
      option = "magma",
      name = "Proportion Resistant",
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      na.value = "grey50"
    ) +
    labs(
      title = paste("Resistance to Beta-Lactam Antibiotics by Age Group in", continent_name),
      x = "Age Group",
      y = "Antibiotic"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "#333333"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 10),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
  
  print(p)
  ggsave(
    filename = paste0("heatmap_resistance_age_group_", gsub(" ", "_", continent_name), ".png"),
    plot = p,
    width = 10,
    height = 8,
    dpi = 300
  )
  
  cat("Generated and saved heatmap for", continent_name, "\n")
  return(p)
}

for (cont in continents) {
  cat("Processing", cont, "\n")
  create_heatmap(cont)
}

# 5. Dodged Bar Chart (Beta-Lactam Resistance vs. Gene Presence)
base_data <- gene_long %>%
  filter(Gene %in% gene_cols_present) %>%
  mutate(gene_positive = !str_ends(Presence, "-Neg")) %>%
  select(`Isolate Id`, continent, Gene, gene_positive) %>%
  left_join(
    mic_interp_joined %>%
      filter(antibiotic %in% beta_lactam_antibiotics) %>%
      select(`Isolate Id`, continent, antibiotic, Interpretation),
    by = c("Isolate Id", "continent")
  ) %>%
  filter(!is.na(Interpretation), Interpretation != "") %>%
  distinct()

# Debug: Check columns in base_data
cat("Columns in base_data:", paste(colnames(base_data), collapse = ", "), "\n")
cat("Rows in base_data:", nrow(base_data), "\n")

gene_summary <- base_data %>%
  group_by(`Isolate Id`, continent) %>%
  summarise(has_gene = any(gene_positive)) %>%
  ungroup()

resistance_summary <- base_data %>%
  filter(antibiotic %in% beta_lactam_antibiotics) %>%
  group_by(`Isolate Id`, continent) %>%
  summarise(has_resistance = any(Interpretation == "Resistant")) %>%
  ungroup()

# Debug: Check column names
cat("Columns in gene_summary:", paste(colnames(gene_summary), collapse = ", "), "\n")
cat("Columns in resistance_summary:", paste(colnames(resistance_summary), collapse = ", "), "\n")

group_data <- gene_summary %>%
  left_join(resistance_summary, by = c("Isolate Id", "continent")) %>%
  mutate(
    group = case_when(
      has_resistance & !has_gene ~ "Resistant_No_Gene",
      has_gene & !has_resistance ~ "Gene_No_Resistance",
      TRUE ~ "Other"
    ),
    `Isolate Id` = `Isolate Id` # Explicitly retain Isolate Id
  ) %>%
  filter(group != "Other")

# Debug: Check columns in group_data before summarization
cat("Columns in group_data before summarization:", paste(colnames(group_data), collapse = ", "), "\n")
cat("Rows in group_data before summarization:", nrow(group_data), "\n")

group_data_summary <- group_data %>%
  group_by(continent, group) %>%
  summarise(count = n_distinct(`Isolate Id`)) %>%
  ungroup()

# Debug: Check columns in group_data_summary
cat("Columns in group_data_summary:", paste(colnames(group_data_summary), collapse = ", "), "\n")
cat("Rows in group_data_summary:", nrow(group_data_summary), "\n")
print(head(group_data_summary))

p <- ggplot(group_data_summary, aes(x = continent, y = count, fill = group)) +
  geom_bar(stat = "identity", position = "dodge", color = "white", width = 0.45) +
  scale_fill_viridis_d(
    option = "magma",
    name = "Isolate Group",
    labels = c("Gene_No_Resistance" = "Gene+, No Resistance", 
               "Resistant_No_Gene" = "Resistant, No Gene"),
    begin = 0.2, end = 0.8
  ) +
  labs(
    title = "Isolates with Beta-Lactam Resistance vs. Resistance Gene Presence by Continent",
    x = "Continent",
    y = "Number of Isolates"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "#333333"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "top",
    legend.title = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

print(p)
ggsave(
  filename = "isolate_groups_by_continent.png",
  plot = p,
  width = 10,
  height = 6,
  dpi = 300
)
cat("Generated and saved bar chart\n")

# 6. Dodged Bar Chart with Stacked Antibiotics
# Debug: Verify columns before join
cat("Columns in group_data before resistant_no_gene:", paste(colnames(group_data), collapse = ", "), "\n")
cat("Rows in group_data:", nrow(group_data), "\n")
cat("Columns in base_data for join :", paste(colnames(base_data), collapse = ", "), "\n")

resistant_no_gene <- group_data %>%
  filter(group == "Resistant_No_Gene") %>%
  left_join(
    base_data %>%
      filter(Interpretation == "Resistant") %>%
      select(`Isolate Id`, continent, antibiotic),
    by = c("Isolate Id" = "Isolate Id", "continent" = "continent")
  ) %>%
  group_by(continent, group, antibiotic) %>%
  summarise(count = n_distinct(`Isolate Id`)) %>%
  ungroup() %>%
  mutate(antibiotic = replace_na(antibiotic, "Unknown"))

# Debug: Check resistant_no_gene
cat("Columns in resistant_no_gene:", paste(colnames(resistant_no_gene), collapse = ", "), "\n")
cat("Rows in resistant_no_gene:", nrow(resistant_no_gene), "\n")
print(head(resistant_no_gene))

gene_no_resistance <- group_data %>%
  filter(group == "Gene_No_Resistance") %>%
  group_by(continent, group) %>%
  summarise(count = n_distinct(`Isolate Id`)) %>%
  ungroup() %>%
  mutate(antibiotic = "Gene_No_Resistance")

plot_data <- bind_rows(resistant_no_gene, gene_no_resistance)

cat("Plot data summary:\n")
print(head(plot_data))

p <- ggplot(plot_data, aes(x = group, y = count, fill = antibiotic)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), color = "white", linewidth = 0.5) +
  facet_wrap(~ continent, scales = "free_x", ncol = length(unique(plot_data$continent))) +
  scale_fill_viridis_d(
    option = "magma",
    name = "Antibiotic (or Group)",
    labels = c("Gene_No_Resistance" = "Gene+, No Resistance", 
               setNames(beta_lactam_antibiotics, beta_lactam_antibiotics),
               "Unknown" = "Unknown Antibiotic"),
    begin = 0, end = 0.9
  ) +
  labs(
    title = "Beta-Lactam Resistance vs. Resistance Gene Presence by Continent",
    x = "Isolate Group",
    y = "Number of Isolates"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "#333333"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.text = element_text(size = 10, face = "bold")
  )

print(p)
ggsave(
  filename = "isolate_groups_by_continent_dodged_stacked.png",
  plot = p,
  width = 12,
  height = 6,
  dpi = 300
)
cat("Generated and saved dodged bar chart\n")

# 7. Sankey Diagram (Non-Beta-Lactam Antibiotic Classes)
non_beta_lactams <- unique(mic_interp_joined$antibiotic[!mic_interp_joined$antibiotic %in% beta_lactam_antibiotics])
cat("Non-beta-lactam antibiotics found:", paste(non_beta_lactams, collapse = ", "), "\n")
if (length(non_beta_lactams) == 0) {
  cat("WARNING: No non-beta-lactam antibiotics found. Using all antibiotics as fallback.\n")
  non_beta_lactams <- unique(mic_interp_joined$antibiotic)
}

if ("antibiotic_class" %in% colnames(mic_interp_joined)) {
  cat("Using existing antibiotic_class column.\n")
  mic_data <- mic_interp_joined %>%
    filter(antibiotic %in% non_beta_lactams) %>%
    select(`Isolate Id`, continent, antibiotic, antibiotic_class, Interpretation)
} else {
  cat("No antibiotic_class column found. Assigning classes to non-beta-lactam antibiotics.\n")
  class_mapping <- tibble(
    antibiotic = non_beta_lactams,
    antibiotic_class = case_when(
      antibiotic %in% c("Gentamicin", "Amikacin", "Tobramycin") ~ "Aminoglycosides",
      antibiotic %in% c("Ciprofloxacin", "Levofloxacin") ~ "Fluoroquinolones",
      antibiotic %in% c("Tetracycline", "Doxycycline") ~ "Tetracyclines",
      antibiotic %in% c("Erythromycin", "Azithromycin", "Clarithromycin") ~ "Macrolides",
      antibiotic %in% c("Trimethoprim", "Sulfamethoxazole") ~ "Sulfonamides",
      TRUE ~ "Other"
    )
  )
  mic_data <- mic_interp_joined %>%
    filter(antibiotic %in% non_beta_lactams) %>%
    left_join(class_mapping, by = "antibiotic") %>%
    select(`Isolate Id`, continent, antibiotic, antibiotic_class, Interpretation)
}

cat("Antibiotic classes:", paste(unique(mic_data$antibiotic_class), collapse = ", "), "\n")

sankey_data <- gene_long %>%
  filter(Gene %in% gene_cols_present, !str_ends(Presence, "-Neg")) %>%
  select(`Isolate Id`, continent, Gene) %>%
  left_join(
    mic_data %>%
      filter(!is.na(Interpretation), Interpretation != ""),
    by = c("Isolate Id" = "Isolate Id", "continent" = "continent")
  ) %>%
  filter(!is.na(antibiotic_class)) %>%
  distinct() %>%
  group_by(continent, Gene, antibiotic_class, Interpretation) %>%
  summarise(count = n_distinct(`Isolate Id`)) %>%
  ungroup()

cat("Sankey data summary:", nrow(sankey_data), "rows\n")
print(head(sankey_data))

continents <- unique(sankey_data$continent[!is.na(sankey_data$continent)])
cat("Continents found:", paste(continents, collapse = ", "), "\n")

create_sankey_non_beta <- function(continent_name) {
  data_subset <- sankey_data %>%
    filter(continent == continent_name)
  
  if (nrow(data_subset) == 0) {
    cat("No data available for", continent_name, "\n")
    return(NULL)
  }
  
  cat("Data for", continent_name, ": ", nrow(data_subset), "rows\n")
  print(head(data_subset))
  
  nodes <- data.frame(name = unique(c(data_subset$Gene, data_subset$antibiotic_class, data_subset$Interpretation)))
  
  links <- data_subset %>%
    group_by(Gene, antibiotic_class) %>%
    summarise(value = sum(count)) %>%
    ungroup() %>%
    mutate(
      source = match(Gene, nodes$name) - 1,
      target = match(antibiotic_class, nodes$name) - 1
    ) %>%
    select(source, target, value) %>%
    bind_rows(
      data_subset %>%
        group_by(antibiotic_class, Interpretation) %>%
        summarise(value = sum(count)) %>%
        ungroup() %>%
        mutate(
          source = match(antibiotic_class, nodes$name) - 1,
          target = match(Interpretation, nodes$name) - 1
        ) %>%
        select(source, target, value)
    )
  
  if (nrow(links) == 0) {
    cat("No valid links for", continent_name, "\n")
    return(NULL)
  }
  
  node_colors <- c(
    viridis(length(unique(data_subset$Gene)), option = "magma", begin = 0, end = 0.4),
    viridis(length(unique(data_subset$antibiotic_class)), option = "magma", begin = 0.5, end = 0.7),
    rep("#666666", length(unique(data_subset$Interpretation)))
  )
  
  sankey_plot <- sankeyNetwork(
    Links = links,
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "name",
    units = "isolates",
    fontSize = 12,
    nodeWidth = 30,
    fontFamily = "Arial",
    iterations = 0,
    colourScale = htmlwidgets::JS(sprintf("d3.scaleOrdinal().range(['%s'])", paste(node_colors, collapse = "','")))
  )
  
  title_html <- htmltools::tags$h3(
    paste("Resistance Genes, Non-Beta-Lactam Antibiotic Classes, and Resistance in", continent_name),
    style = "text-align: center; font-family: Arial; color: #333;"
  )
  sankey_plot <- htmlwidgets::prependContent(sankey_plot, title_html)
  
  print(sankey_plot)
  saveWidget(
    sankey_plot,
    file = paste0("sankey_non_beta_lactam_classes_", gsub(" ", "_", continent_name), ".html"),
    title = paste("Resistance Genes, Non-Beta-Lactam Antibiotic Classes, and Resistance in", continent_name),
    selfcontained = TRUE
  )
  
  cat("Generated and saved Sankey diagram for", continent_name, "\n")
  return(sankey_plot)
}

for (cont in continents) {
  cat("Processing", cont, "\n")
  create_sankey_non_beta(cont)
}



##############################################################################

# 12. Box Plot (Number of Resistant Antibiotics by Source)
resistance_count <- mic_interp_joined %>%
  filter(antibiotic %in% beta_lactam_antibiotics, !is.na(Interpretation), Interpretation != "") %>%
  group_by(`Isolate Id`, Source) %>%
  summarise(resistant_count = sum(Interpretation == "Resistant")) %>%
  ungroup() %>%
  filter(!is.na(Source))

cat("Resistance count by source summary:", nrow(resistance_count), "rows\n")
print(head(resistance_count))

p <- ggplot(resistance_count, aes(x = Source, y = resistant_count, fill = Source)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_fill_viridis_d(option = "magma", name = "Source") +
  labs(
    title = "Number of Resistant Beta-Lactam Antibiotics by Isolate Source",
    x = "Source",
    y = "Number of Resistant Antibiotics"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "#333333"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

print(p)
ggsave("resistance_by_source_boxplot.png", p, width = 10, height = 6, dpi = 300)
cat("Generated and saved resistance by source box plot\n")


#################################################################################




# 3. Heatmap (Gene Co-occurrence by Continent)
# Prepare co-occurrence data stratified by continent
gene_cooccurrence_data <- gene_long %>%
  filter(Gene %in% gene_cols_present, !str_ends(Presence, "-Neg")) %>%
  select(`Isolate Id`, continent, Gene) %>%
  group_by(`Isolate Id`, continent) %>%
  summarise(genes_present = list(Gene)) %>%
  ungroup() %>%
  filter(!is.na(continent))

# Debug: Check continents and genes
cat("Continents in gene_cooccurrence_data:", paste(unique(gene_cooccurrence_data$continent), collapse = ", "), "\n")
cat("Genes in gene_cooccurrence_data:", paste(unique(unlist(gene_cooccurrence_data$genes_present)), collapse = ", "), "\n")

# Debug: Check number of genes per isolate by continent
gene_counts <- gene_cooccurrence_data %>%
  group_by(continent, `Isolate Id`) %>%
  summarise(num_genes = length(genes_present[[1]])) %>%
  ungroup()
cat("Gene counts per isolate by continent:\n")
print(gene_counts %>% group_by(continent) %>% summarise(min_genes = min(num_genes), max_genes = max(num_genes), total_isolates = n()))
cat("Isolates with fewer than 2 genes:\n")
print(gene_counts %>% filter(num_genes < 2))

# Function to compute co-occurrence matrix for a given continent
compute_cooccurrence <- function(cont_data, continent_name) {
  # Filter isolates with at least 2 positive genes
  cont_data <- cont_data %>%
    mutate(num_genes = map_int(genes_present, length)) %>%
    filter(num_genes >= 2)
  
  if (nrow(cont_data) == 0) {
    cat("No isolates with 2 or more genes for", continent_name, "\n")
    return(NULL)
  }
  
  # Get all pairs of genes for each isolate
  cooccur_pairs <- cont_data %>%
    mutate(pairs = map(genes_present, ~ combn(.x, 2, simplify = FALSE))) %>%
    unnest(pairs) %>%
    mutate(
      Gene1 = map_chr(pairs, ~ .x[1]),
      Gene2 = map_chr(pairs, ~ .x[2])
    ) %>%
    select(Gene1, Gene2)
  
  # Total isolates in this continent (with at least 2 genes)
  total_isolates <- n_distinct(cont_data$`Isolate Id`)
  
  # Compute co-occurrence counts and percentages
  cooccur_matrix <- cooccur_pairs %>%
    group_by(Gene1, Gene2) %>%
    summarise(count = n()) %>%
    ungroup() %>%
    mutate(percentage = (count / total_isolates) * 100) %>%
    select(Gene1, Gene2, percentage)
  
  # Create a complete matrix with zeros for non-co-occurring pairs
  all_genes <- unique(c(cooccur_matrix$Gene1, cooccur_matrix$Gene2))
  if (length(all_genes) == 0) {
    cat("No valid gene pairs for", continent_name, "\n")
    return(NULL)
  }
  
  cooccur_matrix <- cooccur_matrix %>%
    complete(Gene1 = all_genes, Gene2 = all_genes, fill = list(percentage = 0)) %>%
    filter(Gene1 != Gene2) %>%  # Remove diagonal
    mutate(continent = continent_name)
  
  return(cooccur_matrix)
}

# Compute co-occurrence for each continent
continents <- unique(gene_cooccurrence_data$continent[!is.na(gene_cooccurrence_data$continent)])
cooccurrence_list <- map(continents, ~ compute_cooccurrence(
  gene_cooccurrence_data %>% filter(continent == .x), .x
))
# Filter out NULL results
cooccurrence_list <- cooccurrence_list[!sapply(cooccurrence_list, is.null)]
cooccurrence_data <- bind_rows(cooccurrence_list)

# Debug: Check co-occurrence data
cat("Co-occurrence data summary:", nrow(cooccurrence_data), "rows\n")
print(head(cooccurrence_data))

# Create heatmap for each continent
create_cooccur_heatmap <- function(cont_name) {
  heatmap_data <- cooccurrence_data %>%
    filter(continent == cont_name) %>%
    mutate(
      Gene1 = factor(Gene1, levels = gene_cols_present),
      Gene2 = factor(Gene2, levels = gene_cols_present)
    ) %>%
    filter(!is.na(Gene1), !is.na(Gene2))
  
  if (nrow(heatmap_data) == 0) {
    cat("No co-occurrence data for", cont_name, "\n")
    return(NULL)
  }
  
  cat("Heatmap data for", cont_name, ":", nrow(heatmap_data), "rows\n")
  print(head(heatmap_data))
  
  p <- ggplot(heatmap_data, aes(x = Gene1, y = Gene2, fill = percentage)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f%%", percentage)), color = "white", size = 3) +
    scale_fill_viridis_c(
      option = "magma",
      name = "Co-occurrence (%)",
      limits = c(0, max(cooccurrence_data$percentage, na.rm = TRUE)),
      breaks = seq(0, max(cooccurrence_data$percentage, na.rm = TRUE), by = 10)
    ) +
    labs(
      title = paste("Resistance Gene Co-occurrence in", cont_name),
      x = "Gene 1",
      y = "Gene 2"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "#333333"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "right",
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
  
  print(p)
  ggsave(
    filename = paste0("gene_cooccurrence_heatmap_", gsub(" ", "_", cont_name), ".png"),
    plot = p,
    width = 8,
    height = 8,
    dpi = 300
  )
  cat("Generated and saved co-occurrence heatmap for", cont_name, "\n")
  return(p)
}

# Generate heatmaps for each continent
for (cont in continents) {
  cat("Processing", cont, "\n")
  create_cooccur_heatmap(cont)
}

##############################################################################

###############################################################################













