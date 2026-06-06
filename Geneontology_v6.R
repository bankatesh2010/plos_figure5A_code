# Load required libraries
library(dplyr)
library(ggplot2)
library(readxl)
library(tidyr)
library(ggpubr)

# Define file path and sheet name
# Input file
file_path <- "Gene_ontology_KO_exclusive.xlsx"
sheet_name <- "GoEnrichmentResult"

# Read the data 
go_data <- read_excel(file_path, sheet = sheet_name)

# Check column names
print(colnames(go_data))

# Rename columns for consistency
go_data <- go_data %>% rename(
  GO_Term           = `ID`,
  result_count      = `Result count`,
  Description       = `Name`,
  p_value           = `P-value`,
  Bgd_count         = `Bgd count`,
  Result_gene_list  = `Result gene list`,
  Pct_of_bgd        = `Pct of bgd`,
  Fold_enrichment   = `Fold enrichment`,
  Odds_ratio        = `Odds ratio`,  
  Benjamini         = `Benjamini`,
  Bonferroni        = `Bonferroni`
)

# Check column types
str(go_data)

# Convert necessary columns to numeric
go_data <- go_data %>%
  mutate(
    Odds_ratio      = suppressWarnings(as.numeric(Odds_ratio)),
    Fold_enrichment = suppressWarnings(as.numeric(Fold_enrichment))
  )

# Convert p-value to -log10 scale and apply correction
go_data <- go_data %>%
  mutate(
    adj_p_value         = p.adjust(p_value, method = "BH"),
    log_p_value         = -log10(adj_p_value),
    log_odds_ratio      = log2(Odds_ratio),
    log_fold_enrichment = log2(Fold_enrichment)
  )

# Filter significant GO terms
significant_go_data <- go_data %>% 
  filter(adj_p_value < 0.05, result_count >= 2)

# Dot Plot
dot_plot <- ggplot(
  significant_go_data,
  aes(
    x     = log_fold_enrichment, 
    y     = reorder(Description, result_count), 
    size  = result_count, 
    color = log_p_value
  )
) +
  geom_point() +
  scale_color_gradient(
    low    = "blue", 
    high   = "red", 
    limits = c(0, max(significant_go_data$log_p_value))
  ) +
  scale_size_continuous(
    breaks = c(5, 10, 15, 20, 30),
    range  = c(2, 10)
  ) +
  scale_y_discrete(expand = expansion(mult = c(0.03, 0.03))) +
  labs(
    title = "GO Term Enrichment",
    x     = "Log2 Fold Enrichment",
    y     = "GO Term",
    size  = "Count",
    color = expression("-log"[10]~italic(P))
  ) +
  theme_minimal() +
  theme(
    plot.title        = element_text(hjust = 0.5),
    panel.grid.major  = element_blank(),   # remove grid
    panel.grid.minor  = element_blank(),   # remove grid
    axis.line         = element_line(color = "black", linewidth = 0.6), # keep axes
    axis.ticks        = element_line(color = "black")
  )

# Display plot
print(dot_plot)

# Save plot
# Output directory
save_dir <- "output"

ggsave(
  filename = file.path(save_dir, "GO_DotPlot_new_3.png"),
  plot     = dot_plot,
  width    = 10,
  height   = 6,
  dpi      = 300
)
