# ------------------ Libraries ------------------
library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(tidyr)
library(ggrepel)

# ------------------ User settings ------------------
file_path  <- "E:/USB/2024/Project work/Proteomics/Differeniation/Result file/data/Protein_Groups__2025-10-02T19-45-01.xlsx"
sheet_name <- "Data"
save_dir   <- "E:/USB/2025/Manuscripts/data"

plot_width  <- 8
plot_height <- 6
plot_dpi    <- 300

if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

# ------------------ Step 1: Load data ------------------
prot <- read_excel(file_path, sheet = sheet_name)

# ------------------ Step 2: Identify LFQ columns ------------------
expr_cols <- grep("^PG\\.MaxLFQ \\(imp\\)", colnames(prot), value = TRUE)
if (length(expr_cols) == 0) stop("❌ No 'PG.MaxLFQ (imp...)' columns found in the dataset!")

# ------------------ Step 3: Parse sample metadata directly from column names ------------------
sample_info <- data.frame(
  Sample = expr_cols,
  Group  = str_extract(expr_cols, regex("\\b(WT|KO|AB)\\b", ignore_case = TRUE)),
  Day_raw = str_extract(expr_cols, regex("\\b(Day\\s*\\d+|D\\s*\\d+|\\b\\d+\\b)\\b", ignore = TRUE)),
  stringsAsFactors = FALSE
)

# Normalize Group capitalization and ensure factor order: WT, KO, AB
sample_info$Group <- toupper(sample_info$Group)
sample_info$Group <- factor(sample_info$Group, levels = c("WT", "KO", "AB"))

# Normalize day into "Day X"
sample_info <- sample_info %>%
  mutate(Day_num = as.numeric(str_extract(Day_raw, "\\d+")),
         Day = paste0("Day ", Day_num))

message("Detected groups and days:")
print(head(sample_info))

# ------------------ Step 4: Extract expression data ------------------
expr_all <- prot %>% select(all_of(sample_info$Sample))

# Attach rownames
if ("Protein Group" %in% colnames(prot)) {
  row_ids <- make.unique(as.character(prot$`Protein Group`))
} else {
  row_ids <- seq_len(nrow(prot))
}
rownames(expr_all) <- row_ids

# ------------------ Step 5: Transform + impute ------------------
expr_all <- as.data.frame(expr_all)
expr_all <- log2(expr_all + 1)

# Mean imputation
expr_all <- expr_all %>%
  mutate(across(everything(), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))

# ------------------ Step 6: Prepare PCA ------------------
expr_t_all <- as.data.frame(t(expr_all))
rownames(expr_t_all) <- sample_info$Sample

# Remove low variance proteins
variances <- apply(expr_t_all, 2, var, na.rm = TRUE)
expr_t_all <- expr_t_all[, variances > 1e-8, drop = FALSE]

message("PCA-ready: ", nrow(expr_t_all), " samples × ", ncol(expr_t_all), " proteins")

# ------------------ Step 7: Run PCA ------------------
pca <- prcomp(expr_t_all, center = TRUE, scale. = TRUE)
pcvar <- summary(pca)$importance[2, ]

# Prepare PCA dataframe and merge metadata
pca_df <- data.frame(pca$x[, 1:2], Sample = rownames(pca$x)) %>%
  left_join(sample_info, by = "Sample")

pca_df$Group <- factor(pca_df$Group, levels = c("WT", "KO", "AB"))
day_levels <- paste0("Day ", sort(unique(pca_df$Day_num)))
pca_df$Day <- factor(pca_df$Day, levels = day_levels)

# ------------------ Step 8: Plot PCA (framed panel, no caption) ------------------
group_colors <- c("WT" = "#3366CC", "KO" = "#CC3333", "AB" = "#99CCFF")
shapes <- c(16, 17, 15, 3, 8)
names(shapes) <- day_levels

p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, shape = Day)) +
  geom_point(size = 4, stroke = 0.5) +
  scale_color_manual(values = group_colors) +
  scale_shape_manual(values = shapes) +
  xlab(paste0("PC1 (", round(pcvar[1] * 100, 1), "% variance)")) +
  ylab(paste0("PC2 (", round(pcvar[2] * 100, 1), "% variance)")) +
  theme_pubr(base_size = 14) +
  theme(
    plot.title = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11, color = "black"),
    panel.grid.major = element_line(color = "gray85", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    # ---- FRAME ADDED HERE (top, right, bottom, left) ----
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  ) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  )

print(p)

# ------------------ Step 9: Save output ------------------
out_file <- file.path(save_dir, "PCA_All_Days_KO_WT_AB_2.png")
ggsave(out_file, p, width = plot_width, height = plot_height, dpi = plot_dpi)

message("✅ PCA plot saved: ", normalizePath(out_file, winslash = "/"))


