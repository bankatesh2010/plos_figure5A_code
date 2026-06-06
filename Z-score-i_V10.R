# Load required libraries
library(dplyr)
library(ggplot2)
library(readxl)
library(tidyr)
library(ggpubr)

# Read the specific sheet
file_path <- "Proteomics_data_codon_count_D8.xlsx"
proteomics_data <- read_excel(file_path, sheet = "KO-D8-WTAB")

# Define save directory
save_dir  <- "output"

# Create the directory if it doesn't exist
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}

# Filter data for Up and Down groups (excluding "Not significant")
data_all <- proteomics_data %>%
  filter(Regulation %in% c("up", "down"))

# Calculate Z-score for each codon
data_all <- data_all %>%
  mutate(
    Z_score_Asp_GAU = (O_f_Asp_GAU - E_f_Asp_GAU) / sd(O_f_Asp_GAU, na.rm = TRUE),
    Z_score_Asp_GAC = (O_f_Asp_GAC - E_f_Asp_GAC) / sd(O_f_Asp_GAC, na.rm = TRUE),
    Z_score_His_CAC = (O_f_His_CAC - E_f_His_CAC) / sd(O_f_His_CAC, na.rm = TRUE),
    Z_score_His_CAU = (O_f_His_CAU - E_f_His_CAU) / sd(O_f_His_CAU, na.rm = TRUE),
    Z_score_Tyr_UAC = (O_f_Tyr_UAC - E_f_Tyr_UAC) / sd(O_f_Tyr_UAC, na.rm = TRUE),
    Z_score_Tyr_UAU = (O_f_Tyr_UAU - E_f_Tyr_UAU) / sd(O_f_Tyr_UAU, na.rm = TRUE),
    Z_score_Asn_AAC = (O_f_Asn_AAC - E_f_Asn_AAC) / sd(O_f_Asn_AAC, na.rm = TRUE),
    Z_score_Asn_AAU = (O_f_Asn_AAU - E_f_Asn_AAU) / sd(O_f_Asn_AAU, na.rm = TRUE),
    Z_score_NAU = (O_f_NAU - E_f_NAU) / sd(O_f_NAU, na.rm = TRUE),
    Z_score_NAC = (O_f_NAC - E_f_NAC) / sd(O_f_NAC, na.rm = TRUE)
  )

# Reshape data function
reshape_data <- function(data, cols, codon_names) {
  data %>%
    select(Regulation, all_of(cols)) %>%
    pivot_longer(cols = all_of(cols), names_to = "Codon", values_to = "Z_score") %>%
    mutate(Codon = recode(Codon, !!!codon_names))
}

data_Asp <- reshape_data(data_all, c("Z_score_Asp_GAU", "Z_score_Asp_GAC"), c("Z_score_Asp_GAU" = "GAU", "Z_score_Asp_GAC" = "GAC"))
data_His <- reshape_data(data_all, c("Z_score_His_CAC", "Z_score_His_CAU"), c("Z_score_His_CAC" = "CAC", "Z_score_His_CAU" = "CAU"))
data_Tyr <- reshape_data(data_all, c("Z_score_Tyr_UAC", "Z_score_Tyr_UAU"), c("Z_score_Tyr_UAC" = "UAC", "Z_score_Tyr_UAU" = "UAU"))
data_Asn <- reshape_data(data_all, c("Z_score_Asn_AAC", "Z_score_Asn_AAU"), c("Z_score_Asn_AAC" = "AAC", "Z_score_Asn_AAU" = "AAU"))
data_NAC_NAU <- reshape_data(data_all, c("Z_score_NAC", "Z_score_NAU"), c("Z_score_NAC" = "NAC", "Z_score_NAU" = "NAU"))

# Custom theme for publication-quality (removing grid lines)
custom_theme <- theme_minimal(base_size = 14, base_family = "Arial") +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.line = element_line(color = "black", linewidth = 0.5),  # **Fixed this line**
    axis.title.x = element_blank()  # Remove x-axis title
  )

# Function to get max Z-score for stat text placement
max_z_score_safe <- function(data) {
  if (any(!is.na(data$Z_score))) {
    return(max(data$Z_score, na.rm = TRUE) * 1.1)
  } else {
    return(2)  # Default value if no valid Z-score
  }
}

# Function to generate plots (showing outliers)
plot_codon <- function(data, y_label, plot_title) {
  ggplot(data, aes(x = Regulation, y = Z_score, fill = Codon)) +
    geom_boxplot(color = "black") +  # **Show outliers (default behavior)**
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkred", size = 1) +  
    labs(title = plot_title, y = y_label, fill = "Codon") +  
    custom_theme +
    scale_fill_manual(values = c("#4E79A7", "#F28E2B")) +
    scale_y_continuous(limits = c(-3, 8), breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8)) +
    scale_x_discrete(labels = c("down" = "Down-regulated", "up" = "Up-regulated")) +  
    stat_compare_means(aes(label = ..p.signif..), method = "t.test", label.y = max_z_score_safe(data), paired = TRUE, size = 5, color = "black") +  
    theme(plot.title = element_text(hjust = 0.5))  
}


# Generate plots
plot_Asp <- plot_codon(data_Asp, "Codon Z-score", "Asp")
plot_His <- plot_codon(data_His, "Codon Z-score", "His")
plot_Tyr <- plot_codon(data_Tyr, "Codon Z-score", "Tyr")
plot_Asn <- plot_codon(data_Asn, "Codon Z-score", "Asp")
plot_NAC_NAU <- plot_codon(data_NAC_NAU, "Codon Z-score", "Q-codons")

# Save images in publication quality (300 dpi)
ggsave(filename = file.path(save_dir, "Zscore_Asp.png"), plot = plot_Asp, width = 6, height = 5, dpi = 300)
ggsave(filename = file.path(save_dir, "Zscore_His.png"), plot = plot_His, width = 6, height = 5, dpi = 300)
ggsave(filename = file.path(save_dir, "Zscore_Tyr.png"), plot = plot_Tyr, width = 6, height = 5, dpi = 300)
ggsave(filename = file.path(save_dir, "Zscore_Asn.png"), plot = plot_Asn, width = 6, height = 5, dpi = 300)
ggsave(filename = file.path(save_dir, "Zscore_NAC_NAU.png"), plot = plot_NAC_NAU, width = 6, height = 5, dpi = 300)

# Display plots
print(plot_Asp)
print(plot_His)
print(plot_Tyr)
print(plot_Asn)
print(plot_NAC_NAU)
