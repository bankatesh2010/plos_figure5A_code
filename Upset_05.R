# ----------------------------------------------------------
# Load libraries
# ----------------------------------------------------------
library(readxl)
library(dplyr)
library(ggplot2)
library(ComplexUpset)

# ----------------------------------------------------------
# File path and sheet names
# ----------------------------------------------------------
file_path <- "E:/USB/2024/Project work/Proteomics/Differeniation/Result file/Proteomics_down and significant.xlsx"
sheet_down <- "KO-WTAB-down-all"
sheet_up   <- "KO-WTAB-up-all"

# ----------------------------------------------------------
# Read Excel sheets
# ----------------------------------------------------------
down_df <- read_excel(file_path, sheet = sheet_down)
up_df   <- read_excel(file_path, sheet = sheet_up)

# ----------------------------------------------------------
# Define days (in the correct order)
# ----------------------------------------------------------
days <- c("Day0", "Day8", "Day10", "Day13", "Day17")

# ----------------------------------------------------------
# Function to make binary presence/absence matrix
# ----------------------------------------------------------
make_binary_matrix <- function(df, days_vec) {
  sets <- lapply(df[days_vec], function(x) na.omit(x))
  all_proteins <- unique(unlist(sets))
  mat <- data.frame(Protein = all_proteins, stringsAsFactors = FALSE)
  for (d in days_vec) {
    mat[[d]] <- all_proteins %in% sets[[d]]
  }
  mat[, c("Protein", days_vec)]
}

# ----------------------------------------------------------
# Prepare matrices
# ----------------------------------------------------------
binary_down <- make_binary_matrix(down_df, days)
binary_up   <- make_binary_matrix(up_df,   days)
binary_all  <- unique(rbind(binary_down, binary_up))

# ----------------------------------------------------------
# Helper to make consistent UpSet plots
# ----------------------------------------------------------
make_upset_plot <- function(data, days, title, color) {
  upset(
    data,
    intersect = days,
    name = title,
    base_annotations = list(
      'Intersection size' = intersection_size(
        counts = TRUE,
        text = list(size = 5, vjust = -0.8, colour = "black"),
        mapping = aes(fill = I(color))
      ) +
        ylab("Intersection Size")
    ),
    stripes = upset_stripes(
      geom = geom_segment(linewidth = 0.6, colour = "black")
    ),
    width_ratio = 0.15,
    sort_sets = FALSE,
    sort_intersections = "descending"
  ) +
    scale_fill_identity() +
    theme(
      text = element_text(family = "Arial", size = 12, colour = "black"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5, colour = "black"),
      axis.title = element_text(size = 12, face = "bold", colour = "black"),
      axis.text = element_text(size = 12, colour = "black"),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      legend.position = "none"
    )
}

# ----------------------------------------------------------
# Generate and save plots separately
# ----------------------------------------------------------

# 1️⃣ All DE proteins
p_all <- make_upset_plot(binary_all, days, "All DE proteins", "black")
ggsave(
  "E:/USB/2024/Project work/Proteomics/Differeniation/Result file/UpSet_All_DE_Proteins.png",
  p_all, width = 10, height = 7, dpi = 300
)
print(p_all)

# 2️⃣ Up-regulated proteins
p_up <- make_upset_plot(binary_up, days, "Up-regulated proteins", "#D55E00")
ggsave(
  "E:/USB/2024/Project work/Proteomics/Differeniation/Result file/UpSet_Upregulated_Proteins.png",
  p_up, width = 10, height = 7, dpi = 300
)
print(p_up)

# 3️⃣ Down-regulated proteins
p_down <- make_upset_plot(binary_down, days, "Down-regulated proteins", "#0072B2")
ggsave(
  "E:/USB/2024/Project work/Proteomics/Differeniation/Result file/UpSet_Downregulated_Proteins.png",
  p_down, width = 10, height = 7, dpi = 300
)
print(p_down)

