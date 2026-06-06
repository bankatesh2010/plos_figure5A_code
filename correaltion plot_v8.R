# Load required libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)

# ----------------------------------------------------------
# Input file and settings
# ----------------------------------------------------------

# Input file
file_path <- "Protein_Groups_D8_D13.xlsx"

sheet_name   <- "Data"
fold_col     <- "logFC (KOD8-WTD8)"
codon_column <- "O_f_NAC"

# Output settings
save_dir <- "output"
out_base <- "codon_vs_foldchange_KOD8_WTD8"
# Plot template
plot_width  <- 6    # inches
plot_height <- 5    # inches
plot_dpi    <- 300
# ————————————————————

# Create save directory if missing
if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

# Ensure the directory is writable
if (file.access(save_dir, 2) != 0) {
  stop("Save directory is not writable: ", save_dir)
}

# 1) Load data
prot <- read_excel(file_path, sheet = sheet_name)

# 2) Compute regulation status
prot <- prot %>%
  mutate(
    Regulation = ifelse(
      abs(`logFC (KOD8-ABD8)`) >= 1 & `P.Val (KOD8-ABD8)` < 0.05 &
        abs(`logFC (KOD8-WTD8)`) >= 1 & `P.Val (KOD8-WTD8)` < 0.05,
      ifelse(`logFC (KOD8-ABD8)` >= 1 & `logFC (KOD8-WTD8)` >= 1, "Up", "Down"),
      "Not significant"
    )
  )

# 3) Select & clean
sel <- prot %>%
  select(
    fold_change     = all_of(fold_col),
    codon_frequency = all_of(codon_column),
    Regulation
  ) %>%
  mutate(across(c(fold_change, codon_frequency), as.numeric)) %>%
  filter(is.finite(fold_change), is.finite(codon_frequency))

# 4) Pearson test and build plotmath label
ct    <- cor.test(sel$fold_change, sel$codon_frequency, method = "pearson")
R     <- round(ct$estimate, 3)
P_str <- format.pval(ct$p.value, digits = 3)

lab <- if (grepl("^<", P_str)) {
  paste0("r==", R, "~italic(P)", P_str)
} else {
  paste0("r==", R, "~italic(P)==", P_str)
}

# 4b) Compute horizontal center
x_center <- mean(range(sel$fold_change, na.rm = TRUE))

# 5) Build plot
p <- ggplot(sel, aes(x = fold_change, y = codon_frequency)) +
  # --- draw greys first (background) ---
  geom_point(
    data = dplyr::filter(sel, Regulation == "Not significant"),
    color = "grey80", size = 2
  ) +
  # --- then colored points on top ---
  geom_point(
    data = dplyr::filter(sel, Regulation != "Not significant"),
    aes(color = Regulation),
    size = 2
  ) +
  # regression line on top (unchanged)
  geom_smooth(method = "lm", linetype = "dashed", color = "black", se = FALSE) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not significant" = "grey80")) +
  labs(
    title = ifelse(codon_column == "O_f_NAU", "NAU", "NAC"),
    x     = "Fold change (logFC KO-D8 vs WT-D8)",
    y     = "Codon frequency"
  ) +
  annotate(
    "text",
    x     = x_center,
    y     = Inf,
    label = lab,
    parse = TRUE,
    hjust = 0.5,
    vjust = 1.05,   # slightly higher to avoid overlap
    size  = 5,
    color = "black",
    family = "sans"
  ) +
  theme_pubr(base_size = 14, base_family = "sans") +
  theme(
    axis.text       = element_text(size = 12),
    axis.title      = element_text(size = 14),
    plot.title      = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

# 6) Display and save
print(p)

# Build output filename in the requested folder
out_file <- file.path(save_dir, paste0(out_base, ".png"))

ggsave(
  filename = out_file,
  plot     = p,
  width    = plot_width,
  height   = plot_height,
  dpi      = plot_dpi,
  units    = "in"
)

# Confirm save
if (file.exists(out_file)) {
  message("Saved figure to: ", normalizePath(out_file, winslash = "/"))
} else {
  stop("Save failed; file not found at: ", out_file)
}

