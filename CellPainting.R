library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(stringr)
library(tidyverse)
library(tidymodels)
library(DESeq2)
library(GSVA)
library(msigdbr)
library(org.Mm.eg.db)
library(forcats)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(grid)

setwd("~/Research/Data/CellPaintingData/")

# Setup Output Directories
output <- "~/Research/Figures_and_Results/CellPainting"
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

#=====================#
# Import RNA-seq Data #
#=====================#
countsData <- read.delim("Salmonmerge_count_genes.txt", header=TRUE, row.names=1, sep="\t", check.names = FALSE)
countsMetadata <- read.delim("Metadata.txt", header = TRUE, sep = "\t", row.names=1)

# Remove 2D Monolayer samples
monolayer_ids <- c("ML1_S45", "ML2_S93", "ML3_S141")
countsData <- countsData[, !colnames(countsData) %in% monolayer_ids]
countsMetadata <- countsMetadata[!rownames(countsMetadata) %in% monolayer_ids, ]

#-----------------------------------------------#
# 🔧 Standardize metadata labels for consistency
#-----------------------------------------------#
countsMetadata$dose <- as.character(countsMetadata$dose)
countsMetadata$dose <- trimws(countsMetadata$dose)

# Replace "untr" with "water" for consistency
countsMetadata$dose[countsMetadata$dose == "untr"] <- "water"

# Define dose as an ordered factor for plotting and integration
countsMetadata$dose <- factor(
  countsMetadata$dose,
  levels = c("control", "water", "dmso1", "dmso2", "Low", "Mid", "High")
)

# Reorder with so sample IDs align
countsData    <- countsData[, order(colnames(countsData))]
countsMetadata <- countsMetadata[order(rownames(countsMetadata)), ]

# Sanity check
stopifnot(all(colnames(countsData) == rownames(countsMetadata)))

# Filter low coverage Genes
# Compute CPM 
lib_sizes <- colSums(countsData)
cpm_vals <- t(t(countsData) / lib_sizes * 1e6)

# Filter: keep genes with CPM > 1 in at least 75% of samples
keep <- rowMeans(cpm_vals > 1) > 0.75
countsData <- countsData[keep, ]

 # Barplot summarising number of genes removed with CPM filter
statsbarplot <- tibble(
  Category = factor(c("Before filter", "Removed", "After filter"),
                    levels = c("Before filter", "Removed", "After filter")),
  Genes = c(nrow(cpm_vals), nrow(cpm_vals) - sum(keep), sum(keep))
)

ggplot(statsbarplot, aes(x = Category, y = Genes, fill = Category)) +
  geom_col() +
  geom_text(aes(label = Genes), vjust = -0.5) +
  scale_fill_manual(values = c("Before filter" = "blue",
                               "Removed" = "red",
                               "After filter" = "darkgreen")) +
  theme_minimal() +
  labs(title = "Gene Filtering Summary (CPM > 1 in 75% of samples)",
       x = "", y = "Number of genes") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

ggsave(file.path(output, "Gene_Filtering_Summary.pdf"), width = 8, height = 6, dpi = 1000, device = "pdf")

# PCA of Raw counts before normalization & transformation
# Transpose so samples are rows
pca_raw <- prcomp(t(countsData), scale. = TRUE)

# Put PCA coordinates into a data frame
pca_raw_df <- as.data.frame(pca_raw$x)
pca_raw_df$compound <- countsMetadata$Label
pca_raw_df$dose <- countsMetadata$dose

# Variance explained
var_explained <- (pca_raw$sdev^2) / sum(pca_raw$sdev^2) * 100
pc1_label <- paste0("PC1 (", round(var_explained[1], 1), "%)")
pc2_label <- paste0("PC2 (", round(var_explained[2], 1), "%)")

# PCA plot
pca_raw_plot <- ggplot(pca_raw_df, aes(PC1, PC2, color = dose)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA of Raw RNA-seq Counts (Without 2D samples)",
    x = pc1_label,
    y = pc2_label,
    color = "Dose"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

pca_raw_plot

ggsave(
  file.path(output, "PCA_RawCounts_Unnormalized.pdf"),
  pca_raw_plot,
  width = 8, height = 6, dpi = 1000, device = "pdf"
)

# Visualize replicate clustering
# Add replicate IDs (rownames of countsMetadata)
pca_raw_df$sample_id <- rownames(countsMetadata)

# Plot by compound and dose, label each replicate
pca_replicate_plot <- ggplot(pca_raw_df, aes(PC1, PC2, color = compound, shape = dose)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text(aes(label = sample_id), size = 2, vjust = -1, check_overlap = TRUE) +
  scale_shape_manual(
    values = c(
      "control" = 16,  # solid circle
      "water"   = 17,  # triangle
      "dmso1"   = 15,  # square
      "dmso2"   = 18,  # diamond
      "Low"     = 3,   # plus
      "Mid"     = 4,   # cross
      "High"    = 8    # star
    )
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA of Raw Counts — Replicate Clustering Check",
    x = pc1_label,
    y = pc2_label,
    shape = "Dose"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

pca_replicate_plot

ggsave(
  filename = file.path(output, "PCA_RawCounts_Replicate_Clustering.pdf"),
  plot = pca_replicate_plot,
  width = 15,
  height = 12,
  dpi = 1000,
  device = "pdf"
)

# Quantify replicate similarity
# Replicates for the same compound–dose should have correlation > 0.9.
# If one replicate shows much lower correlation (< 0.7–0.8), it’s suspicious.

# Compute pairwise Pearson correlations across all samples
sample_cor <- cor(countsData, method = "pearson")
dim(sample_cor)

# Quantify replicate correlation per condition
# Add a simplified condition label for grouping
countsMetadata$condition <- countsMetadata$Label

# Compute within-condition correlations
rep_cor_summary <- lapply(unique(countsMetadata$condition), function(cond) {
  samples <- rownames(countsMetadata[countsMetadata$condition == cond, ])
  sub_cor <- sample_cor[samples, samples, drop = FALSE]
  
  tibble(
    condition = cond,
    n_replicates = nrow(sub_cor),
    mean_cor = mean(sub_cor[lower.tri(sub_cor)]),
    min_cor = min(sub_cor[lower.tri(sub_cor)])
  )
}) %>%
  bind_rows()

rep_cor_summary <- rep_cor_summary %>%
  arrange(mean_cor)

rep_cor_summary

# Identify which replicate(s) to remove
# Identify outliers: replicate correlation < 0.8 within its group
outliers <- list()

for (cond in unique(countsMetadata$condition)) {
  samples <- rownames(countsMetadata[countsMetadata$condition == cond, ])
  sub_cor <- sample_cor[samples, samples, drop = FALSE]
  
  if (nrow(sub_cor) > 2) {
    mean_by_sample <- colMeans(sub_cor)
    bad <- samples[mean_by_sample < 0.8]
    if (length(bad) > 0) outliers[[cond]] <- bad
  }
}

outliers

# Remove outliers
countsData <- countsData[, !colnames(countsData) %in% unlist(outliers)]
countsMetadata <- countsMetadata[!rownames(countsMetadata) %in% unlist(outliers), ]

#==============================# 
# Deseq2 default normalization #
#==============================#
# Ensure sample names match
stopifnot(all(colnames(countsData) == rownames(countsMetadata)))

# Check if all counts are integers
summary(countsData %% 1 == 0)

# If they are not, round to nearest integer 
countsData <- round(countsData)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = countsData,
  colData   = countsMetadata,
  design    = ~ 1   # no design yet, just normalization
)

# Estimate size factors (normalization)
dds <- estimateSizeFactors(dds)

# Extract normalized counts
norm_counts <- counts(dds, normalized=TRUE)

# Variance-stabilizing transform 
vsd <- vst(dds, blind = TRUE) # Blind so that it doesn’t look at design formula (compound, dose)

# Matrix: genes x samples (same as counts, but stabilized)
vsd_mat <- assay(vsd)

dim(vsd_mat)
vsd_mat[1:5, 1:5]

countsData <- vsd_mat

#  PCA
# Make sure sample order matches columns of vsd_mat
stopifnot(all(colnames(vsd_mat) == rownames(countsMetadata)))

# PCA
rna_pca <- prcomp(t(countsData), scale. = TRUE)
rna_pca_df <- as.data.frame(rna_pca$x)

# Add metadata
rna_pca_df$compound <- countsMetadata$Label
rna_pca_df$dose <- countsMetadata$dose  

# Variance explained
var_explained <- (rna_pca$sdev^2) / sum(rna_pca$sdev^2) * 100
pc1_label <- paste0("PC1 (", round(var_explained[1], 1), "%)")
pc2_label <- paste0("PC2 (", round(var_explained[2], 1), "%)")

# PCA plot — color by compound (label)
rna_pca_plot <- ggplot(rna_pca_df, aes(x = PC1, y = PC2, color = compound)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA of Normalized VST RNA-seq Data",
    x = pc1_label,
    y = pc2_label,
    color = "Compound"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

rna_pca_plot

#================================#
# Calculating cytotoxicity Index #
#================================#
#Pathway-based Cytotoxicity Index (via enrichment scores)
# Retrieve mouse-native pathway collections
# Hallmark (MH), Canonical Pathways (MC2), and GO Biological Processes (MC5)
msig_h  <- msigdbr(db_species = "MM", species = "Mus musculus", collection = "MH")
msig_c2 <- msigdbr(db_species = "MM", species = "Mus musculus", collection = "M2")
msig_c5 <- msigdbr(db_species = "MM", species = "Mus musculus", collection = "M5")

# Combine them
gene_sets_full <- bind_rows(msig_h, msig_c2, msig_c5)

# Select only the four key cell death pathways (Apoptosis, Ferroptosis, Macroautophagy, & Anoikis)
# Filter to core apoptosis, ferroptosis, macroautophagy, and anoikis
selected_sets <- gene_sets_full %>%
  filter(
    gs_name %in% c(
      "HALLMARK_APOPTOSIS",
      "REACTOME_APOPTOSIS",
      "WP_APOPTOSIS",
      "GOBP_FERROPTOSIS",
      "GOBP_MACROAUTOPHAGY",
      "REACTOME_MACROAUTOPHAGY",
      "GOBP_ANOIKIS",
      "GOBP_REGULATION_OF_ANOIKIS",
      "GOBP_POSITIVE_REGULATION_OF_ANOIKIS"
    )
  ) %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%  # remove duplicate genes
  split(x = .$gene_symbol, f = .$gs_name)

cat("✅ Final selected pathways:", length(selected_sets), "\n")
print(names(selected_sets))

# Merge all Apoptosis-related pathways into one set
apoptosis_genes <- gene_sets_full %>%
  filter(gs_name %in% c("HALLMARK_APOPTOSIS", "REACTOME_APOPTOSIS", "WP_APOPTOSIS")) %>%
  pull(gene_symbol) %>%
  unique()

# Merge all Macroautophagy-related pathways into one set
macroautophagy_genes <- gene_sets_full %>%
  filter(gs_name %in% c("GOBP_MACROAUTOPHAGY", "REACTOME_MACROAUTOPHAGY")) %>%
  pull(gene_symbol) %>%
  unique()

# Merge all Anoikis-related pathways into one set
anoikis_genes <- gene_sets_full %>%
  filter(gs_name %in% c("GOBP_ANOIKIS", 
                        "GOBP_REGULATION_OF_ANOIKIS",
                        "GOBP_POSITIVE_REGULATION_OF_ANOIKIS")) %>%
  pull(gene_symbol) %>%
  unique()

# Build final unified pathway list
selected_sets <- list(
  APOPTOSIS = apoptosis_genes,
  FERROPTOSIS = unique(gene_sets_full$gene_symbol[gene_sets_full$gs_name == "GOBP_FERROPTOSIS"]),
  MACROAUTOPHAGY = macroautophagy_genes,
  ANOIKIS = anoikis_genes
)

cat("✅ Final unified pathways:", length(selected_sets), "\n")
print(names(selected_sets))

#==============================#
# 📊 Barplot: Genes per pathway
#==============================#

# Count number of genes in each pathway
pathway_sizes <- sapply(selected_sets, length)
pathway_df <- data.frame(
  Pathway = names(pathway_sizes),
  GeneCount = as.numeric(pathway_sizes)
)

# Define consistent custom colors (from your established palette)
pathway_colors <- c(
  "APOPTOSIS" = "#D55E00",      # orange/red tone
  "FERROPTOSIS" = "#0072B2",    # blue tone
  "MACROAUTOPHAGY" = "#009E73", # green tone
  "ANOIKIS" = "#E69F00"         # amber tone
)

# Plot
p_pathway_bar <- ggplot(pathway_df, aes(x = reorder(Pathway, -GeneCount), y = GeneCount, fill = Pathway)) +
  geom_col(width = 0.7, color = "black") +
  geom_text(aes(label = GeneCount), vjust = -0.5, size = 5) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = pathway_colors) +
  labs(
    #title = "Number of Genes per Pathway",
    x = "Pathway",
    y = "Number of Genes"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 13)
  )

# Show plot
p_pathway_bar

# Save
ggsave(
  filename = file.path(output, "Gene_Counts_Per_Pathway.pdf"),
  plot = p_pathway_bar,
  width = 8,
  height = 6,
  dpi = 1000,
  device = "pdf"
)
# ========== End of Bar Plot ===================#

# Convert countsData Ensembl IDs to gene symbols
ensembl_ids <- rownames(countsData)
gene_map <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = ensembl_ids,
  keytype = "ENSEMBL",
  columns = "SYMBOL"
)

# Remove duplicates and NAs
gene_map <- gene_map %>%
  distinct(ENSEMBL, .keep_all = TRUE) %>%
  filter(!is.na(SYMBOL))

# Replace Ensembl gene IDs with gene symbols
countsData_symbol <- countsData[gene_map$ENSEMBL, ]
rownames(countsData_symbol) <- gene_map$SYMBOL

# Prepare numeric Matrix for GSVA & Heatmap plotting
expr_mat <- as.matrix(countsData_symbol)

# Remove duplicated gene symbols by keeping the first occurrence
expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), ]

# Ensure sample IDs match
stopifnot(all(colnames(expr_mat) %in% rownames(countsMetadata)))

# Rename columns using 'Label' to show compound_dose
colnames(expr_mat) <- countsMetadata[colnames(expr_mat), "Label"]

# Verify
head(colnames(expr_mat))

# Confirm all rownames are unique
cat("✅ Unique row names:", length(unique(rownames(expr_mat))), "/", nrow(expr_mat), "\n")

# Run ssGSEA to calculate enrichment scores
# Ensure GSVA uses unique IDs 
colnames(expr_mat) <- rownames(countsMetadata)
params <- ssgseaParam(exprData = expr_mat, geneSets = selected_sets)
gsva_result <- gsva(params) #calculates a pathway activity score for every pathway in each sample

# Extract numeric matrix of scores
if (is.matrix(gsva_result)) {
  gsva_scores <- gsva_result
} else if (inherits(gsva_result, "SummarizedExperiment")) {
  gsva_scores <- SummarizedExperiment::assay(gsva_result)
} else {
  stop("Unexpected gsva_result type: ", class(gsva_result))
}

cat("✅ GSVA score matrix dimensions:", dim(gsva_scores), "\n")

head(gsva_scores)
summary(as.vector(gsva_scores))

all(colnames(gsva_scores) == rownames(countsMetadata))

# Summarize sample names by compound–dose
sample_map <- countsMetadata[, c("Label", "dose")]

# Average replicates per compound–dose
heatmap_mat <- sapply(unique(sample_map$Label), function(lbl) {
  samples <- rownames(sample_map)[sample_map$Label == lbl]
  rowMeans(gsva_scores[, samples, drop = FALSE], na.rm = TRUE)
})

# Z-score scale and visualize
gsva_scaled <- t(scale(t(heatmap_mat)))
ssGSEA_Heatmap <- Heatmap(
  gsva_scaled,
  name = "Pathway Activity",
  col = colorRamp2(c(-2, 0, 2), c("#4575B4", "white", "#D73027")),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_rot = 45,
  heatmap_width = unit(14, "in"),
  heatmap_height = unit(8, "in"),
  #column_title = "Compound–Dose",
  row_title = "Pathway"
)

ssGSEA_Heatmap

# Save heatmap
pdf(file.path(output, "ssGSEA_Pathway_Activity_Heatmap.pdf"),
    width = 18, height = 9)
draw(ssGSEA_Heatmap, heatmap_legend_side = "right")
dev.off()

#====================#
# PCA of GSVA Scores #
#====================#
pca_gsva <- prcomp(t(gsva_scaled), scale. = TRUE)

# Create a tidy dataframe for plotting
pca_gsva_df <- as.data.frame(pca_gsva$x)

# Extract variance explained
var_explained <- (pca_gsva$sdev^2) / sum(pca_gsva$sdev^2) * 100
pc1_label <- paste0("PC1 (", round(var_explained[1], 1), "%)")
pc2_label <- paste0("PC2 (", round(var_explained[2], 1), "%)")

# Add metadata for each compound-dose
pca_gsva_df$Label <- colnames(gsva_scaled)
pca_gsva_df <- pca_gsva_df %>%
  left_join(sample_map, by = "Label")

# Plot PCA
p_gsva_pca <- ggplot(pca_gsva_df, aes(x = PC1, y = PC2, color = dose)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text(aes(label = Label), size = 3, vjust = -0.7, check_overlap = TRUE) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA of GSVA Pathway Activity (Z-scored)",
    x = pc1_label,
    y = pc2_label,
    color = "Dose"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "right"
  )

p_gsva_pca

# Save high-resolution PCA figure
ggsave(
  filename = file.path(output, "PCA_GSVA_PathwayActivity.pdf"),
  plot = p_gsva_pca,
  width = 10,
  height = 8,
  dpi = 1000,
  device = "pdf"
)

# ======= End of PCA ===================#

# Compute the Cytotoxicity Index
cyto_index_raw <- colMeans(gsva_scores, na.rm = TRUE)      # absolute stress activity
cyto_index_scaled <- colMeans(t(scale(t(gsva_scores))), na.rm = TRUE)  # relative stress activity

# Add to metadata
countsMetadata$Cytotoxicity_Index_raw <- cyto_index_raw[rownames(countsMetadata)]
countsMetadata$Cytotoxicity_Index_scaled <- cyto_index_scaled[rownames(countsMetadata)]

#=======================================#
# Average replicate cytotoxicity scores #
#=======================================#
#--- Ensure dose is ordered properly ---
dose_order <- c("control", "water", "dmso1", "dmso2", "Low", "Mid", "High")
countsMetadata$dose <- factor(countsMetadata$dose, levels = dose_order)

#--- Compute averages per condition ---
Cytotox_summary <- countsMetadata %>%
  dplyr::group_by(Label, dose) %>%
  dplyr::summarise(
    mean_CI_raw = mean(Cytotoxicity_Index_raw, na.rm = TRUE),
    sd_CI_raw   = sd(Cytotoxicity_Index_raw, na.rm = TRUE),
    mean_CI_scaled = mean(Cytotoxicity_Index_scaled, na.rm = TRUE),
    sd_CI_scaled   = sd(Cytotoxicity_Index_scaled, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"   # removes the "grouped output" message
  )

# Use scaled CI for visualization == #
countsMetadata_avg <- Cytotox_summary %>%
  dplyr::select(Label, dose, Cytotoxicity_Index = mean_CI_scaled)

#---------------------------------------#
# 🔹 Cytotoxicity Index Visualization  #
#---------------------------------------#

# Define order and colors
dose_colors <- c(
  "control" = "#999999",
  "water"   = "#0072B2",
  "dmso1"   = "#009E73",
  "dmso2"   = "#56B4E9",
  "Low"     = "#F0E442",
  "Mid"     = "#E69F00",
  "High"    = "#D55E00"
)

# Clean and structure compound labels
countsMetadata_avg <- countsMetadata_avg %>%
  mutate(
    compound = str_remove(Label, "_(Low|Mid|High)$"),  # remove dose suffix
    dose_level = factor(dose, levels = dose_order),
    compound_dose = paste0(compound, "_", dose_level),
    # Remove redundant names for controls
    compound_dose = case_when(
      compound_dose %in% c("control_control", "dmso1_dmso1",
                           "dmso2_dmso2", "water_water") ~
        str_remove(compound_dose, "_.*$"),
      TRUE ~ compound_dose
    )
  )

# Define x-axis order: controls first, then compounds (Low→Mid→High)
control_first <- c("control", "water", "dmso1", "dmso2")
compound_levels <- unique(countsMetadata_avg$compound[!countsMetadata_avg$compound %in% control_first])

x_axis_order <- c(
  control_first,
  unlist(lapply(compound_levels, function(comp) paste0(comp, "_", c("Low", "Mid", "High"))))
)
x_axis_order <- x_axis_order[x_axis_order %in% countsMetadata_avg$compound_dose]
countsMetadata_avg$compound_dose <- factor(countsMetadata_avg$compound_dose, levels = x_axis_order)

#--- Plot Cytotoxicity Index ---
p_CI <- ggplot(countsMetadata_avg, aes(x = compound_dose, y = Cytotoxicity_Index, color = dose_level)) +
  geom_line(aes(group = compound), linewidth = 0.8, alpha = 0.8) +   # connect Low–Mid–High per compound
  geom_point(size = 3) +
  scale_color_manual(values = dose_colors, breaks = dose_order) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 9),
    panel.grid.major.x = element_line(color = "grey90"),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  labs(
    #title = "Dose–Response Cytotoxicity Trends Across Compounds (RNA-seq)",
    x = "Compound and Dose Level",
    y = "Cytotoxicity Index",
    color = "Dose Level"
  )

p_CI

#--- Save Cytotoxicity Index plot to output folder ---
ggsave(
  filename = file.path(output, "Cytotoxicity_Index_Dose_Response.pdf"),
  plot = p_CI,
  width = 12,
  height = 8,
  dpi = 1000,
  device = "pdf"
)

# Summarise CI across all doses
CI_summary_rank <- countsMetadata_avg %>%
  arrange(desc(Cytotoxicity_Index))

write.csv(CI_summary_rank, file.path(output, "Cytotoxicity_Index_Ranking.csv"), row.names = FALSE)

#==========================================#
# Save environment for efficient reloading #
#==========================================#
rdata_dir <- file.path(output, "RDS_Files")
if (!dir.exists(rdata_dir)) dir.create(rdata_dir, recursive = TRUE)

save.image(file = file.path(rdata_dir, "Cytotoxicity_Index_workspace.RData"))

#============================#
# Quick load of RNA-seq Data #
#============================#
rdata_path <- file.path(output, "RDS_Files", "Cytotoxicity_Index_workspace.RData")
load(rdata_path)
#=============================#

# Visualize global cytotoxicty Index Distribution
ggplot(countsMetadata_avg, aes(x = dose, y = Cytotoxicity_Index, color = dose)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2) +
  theme_minimal() +
  labs(y = "Cytotoxicity Index", x = "Dose", title = "Global Stress Activation per Dose Level")

#====================# 
#  Import image Data #
#====================#
ImageData <- read.delim("merged_EB_clean.csv", header=TRUE, sep=",", check.names = FALSE)

# Exclude all the rows with NA from the Image Data
ImageData <- ImageData[!(is.na(ImageData$ImageNumber) | grepl("^NA", trimws(ImageData$ImageNumber))), ]

# Remove non-feature and NA columns
ImageData <- ImageData %>%
  dplyr::select(
    -c(ImageNumber, ObjectNumber,          # Drop ID/index columns
       Location_Center_Z, Number_Object_Number, AreaShape_EulerNumber)
  ) %>%
  dplyr::select(where(~ !any(is.na(.))))  # Remove columns that contain any NA values

# Remove 'AreaShape_' from column names
colnames(ImageData) <- gsub("^AreaShape_", "", colnames(ImageData))

# Summarize all numeric features by median per treatment
ImageData <- ImageData %>%
  group_by(FileName_OrigGreen) %>%
  summarise(
    across(where(is.numeric), \(x) median(x, na.rm = TRUE))
  ) %>%
  ungroup()

# Scale all features (z-score standardization) 
ImageData <- ImageData %>%
  mutate(across(where(is.numeric), ~ as.numeric(scale(.x))))

# Extract condition and numeric concentration to match dose column in RNA-seq counts Metadata
ImageData <- ImageData %>%
  mutate(
    compound = str_trim(str_extract(FileName_OrigGreen, "^[^(]+")),
    conc_value = as.numeric(str_extract(FileName_OrigGreen, "(?<=\\()[0-9\\.]+")),
    conc_value = case_when(
      is.na(conc_value) & compound == "Carrier 1" ~ "dmso1",
      is.na(conc_value) & compound == "Carrier 2" ~ "dmso2",
      is.na(conc_value) & compound == "Carrier 3" ~ "water",
      is.na(conc_value) & compound == "Control"   ~ "control",
      TRUE ~ as.character(conc_value)
    )
  )

# Assign Low / Mid / High per compound
ImageData <- ImageData %>%
  group_by(compound) %>%
  mutate(
    dose_level = case_when(
      all(is.na(conc_value)) & compound == "Carrier 1" ~ "dmso1",
      all(is.na(conc_value)) & compound == "Carrier 2" ~ "dmso2",
      all(is.na(conc_value)) & compound == "Carrier 3" ~ "water",
      all(is.na(conc_value)) & compound == "Control"   ~ "control",
      conc_value == min(conc_value, na.rm = TRUE) ~ "Low",
      conc_value == max(conc_value, na.rm = TRUE) ~ "High",
      TRUE ~ "Mid"
    )
  ) %>%
  ungroup()

# Standardize carrier compound names to match RNA-seq metadata
ImageData <- ImageData %>%
  mutate(
    compound = case_when(
      compound == "Carrier 1" ~ "dmso1",
      compound == "Carrier 2" ~ "dmso2",
      compound == "Carrier 3" ~ "water",
      compound == "Control"   ~ "control",
      TRUE ~ compound
    ),
    dose_level = case_when(
      compound %in% c("control", "water", "dmso1", "dmso2") ~ NA_character_,  # remove dose
      TRUE ~ dose_level
    )
  )

# Image Data PCA
img_pca <- prcomp(dplyr::select(ImageData, where(is.numeric)), scale. = FALSE)

pca_df <- data.frame(
  img_pca$x[, 1:2],
  compound = ImageData$compound,
  dose_level = ImageData$dose_level
)

# Clean descriptive labels (omit dose for controls)
control_labels <- c("control", "water", "dmso1", "dmso2")
pca_df$Label <- ifelse(
  pca_df$compound %in% control_labels,
  pca_df$compound,
  paste0(pca_df$compound, " (", pca_df$dose_level, ")")
)

Image_pca_plot <- ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(color = "#1b9e77", size = 3, alpha = 0.9) +
  ggrepel::geom_text_repel(
    aes(label = Label),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = NA
  ) +
  theme_minimal(base_size = 14) +
  labs(
    #title = "PCA of Morphological Features",
    x = paste0("PC1 (", round(100*summary(img_pca)$importance[2,1], 1), "% variance)"),
    y = paste0("PC2 (", round(100*summary(img_pca)$importance[2,2], 1), "% variance)")
  ) +
  theme(plot.title = element_text(hjust = 0.5))

Image_pca_plot

#====================================================#
# 🧩 Harmonize and Integrate Image with RNA-seq Data #
#====================================================#
# Clean and standardize compound and dose naming in ImageData
ImageData <- ImageData %>%
  mutate(
    compound = str_to_lower(compound),
    compound = str_replace_all(compound, " ", "_"),
    dose_level = case_when(
      compound %in% c("control", "dmso1", "dmso2", "water") ~ compound,
      TRUE ~ str_to_lower(dose_level)
    )
  )

# Standardize in RNA-seq data
countsMetadata_avg <- countsMetadata_avg %>%
  mutate(
    compound = str_to_lower(compound),
    compound = str_replace_all(compound, " ", "_"),
    dose = str_to_lower(dose)
  )

# Merge Image features with RNA-seq Cytotoxicity Index
IntegratedData <- inner_join(
  ImageData,
  countsMetadata_avg,
  by = c("compound" = "compound", "dose_level" = "dose")
)

cat("✅ Integrated dataset:", nrow(IntegratedData), "samples ×", ncol(IntegratedData), "features\n")

# Confirm One-to-One Matching of integrated Data
# This ensures that each condition in RNA-seq has a corresponding image-derived phenotype
table(IntegratedData$compound, IntegratedData$dose_level)

# RQ: How do morphological features change with increasing cytotoxicity
# Compute feature-wise correlations between morphological features and the CI
# Select only numeric image features and CI
numeric_cols <- sapply(IntegratedData, is.numeric)
numeric_features <- IntegratedData[, numeric_cols, drop = FALSE]
numeric_features <- numeric_features[, setdiff(names(numeric_features), "Cytotoxicity_Index"), drop = FALSE]

# Compute correlations for each numeric feature
cor_results <- apply(numeric_features, 2, function(x)
  cor(x, IntegratedData$Cytotoxicity_Index, use = "complete.obs", method = "pearson")
)

# Convert to a tidy data frame
cor_summary <- data.frame(
  Feature = names(cor_results),
  Correlation = as.numeric(cor_results)
) %>%
  arrange(desc(abs(Correlation)))  # sort by strength of association

# Inspect the top 10 features
head(cor_summary, 10)

# Heatmap
pheatmap(
  matrix(cor_summary$Correlation, nrow = 1,
         dimnames = list("Correlation", cor_summary$Feature)),
  color = colorRampPalette(c("blue", "white", "red"))(100),
  cluster_cols = FALSE,     # 🔹 disable clustering (fixes the error)
  cluster_rows = FALSE,
  main = "Correlation of Morphological Features with Cytotoxicity Index",
  fontsize_col = 6,
  border_color = NA
)

# Barplot
ggplot(cor_summary, aes(x = reorder(Feature, Correlation), y = Correlation,
                        fill = Correlation > 0)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  theme_minimal(base_size = 13) +
  scale_fill_manual(values = c("#4575B4", "#D73027")) +
  labs(
    title = "Morphological Features Correlated with Cytotoxicity Index",
    x = "Feature",
    y = "Pearson Correlation"
  )

# Correlation plots
# Select top 20 features most correlated with CI
top_features <- cor_summary$Feature[1:20]

# Generate scatter plots for each feature
plot_list <- list()

for (i in seq_along(top_features)) {
  f <- top_features[i]

  # Compute correlation and p-value
  r_val <- cor(IntegratedData[[f]], IntegratedData$Cytotoxicity_Index, use = "complete.obs")
  p_val <- cor.test(IntegratedData[[f]], IntegratedData$Cytotoxicity_Index)$p.value

  # Label for annotation
  label_text <- paste0("r = ", round(r_val, 2), 
                       ", p = ", format.pval(p_val, digits = 2, eps = .001))

  plot_list[[i]] <- ggplot(IntegratedData, aes_string(x = "Cytotoxicity_Index", y = f)) +
    geom_point(alpha = 0.7, color = "#1b9e77", size = 1.8) +
    geom_smooth(method = "lm", se = FALSE, color = "#d95f02", linewidth = 0.6) +
    theme_minimal(base_size = 11) +
    labs(
      title = paste0(LETTERS[i], ". ", f),
      x = "Cytotoxicity Index (CI)",
      y = f
    ) +
    annotate("text", x = Inf, y = Inf, label = label_text, 
             hjust = 1.1, vjust = 1.5, size = 2.8, color = "black") +
    theme(
      plot.title = element_text(face = "bold", size = 9, hjust = 0),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7),
      plot.margin = margin(2, 2, 2, 2)
    )
}

# Combine into 4 × 5 panel layout
Feature_CI_Top20 <- wrap_plots(plotlist = plot_list, ncol = 4) +
  plot_annotation(
    title = "Top 20 Morphological Features Correlated with Cytotoxicity Index",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

# Display
Feature_CI_Top20

# Save high-quality figure
ggsave(
  filename = file.path(output, "Figure_Morphology_vs_Cytotoxicity_Top20.pdf"),
  plot = Feature_CI_Top20,
  width = 14,
  height = 12,
  dpi = 1000,
  device = "pdf"
)

# ===== Figure 2 ==========
# Correlation of selected morphological features
selected_features <- c(
  "MeanRadius",
  "Area",
  "Zernike_2_2",
  "Eccentricity",
  "FormFactor",
  "Zernike_6_0",
  "Zernike_8_0"
)

plot_list <- list()

for (i in seq_along(selected_features)) {
  
  f <- selected_features[i]

  r_val <- cor(
    IntegratedData[[f]],
    IntegratedData$Cytotoxicity_Index,
    use = "complete.obs"
  )
  
  p_val <- cor.test(
    IntegratedData[[f]],
    IntegratedData$Cytotoxicity_Index
  )$p.value

  label_text <- paste0(
    "r = ", round(r_val, 2),
    ", p = ", format.pval(p_val, digits = 2, eps = .001)
  )

  plot_list[[i]] <- ggplot(
    IntegratedData,
    aes_string(x = "Cytotoxicity_Index", y = f)
  ) +
    geom_point(alpha = 0.7,
               color = "#1b9e77",
               size = 1.8) +
    geom_smooth(method = "lm",
                se = FALSE,
                color = "#d95f02",
                linewidth = 0.6) +
    theme_minimal(base_size = 11) +
    labs(
      x = "Cytotoxicity Index (CI)",
      y = f
    ) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = label_text,
      hjust = 1.1,
      vjust = 1.5,
      size = 2.8
    ) +
    theme(
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7),
      plot.margin = margin(4, 4, 4, 4)
    )
}

# ---------------------------
# Two-row layout (4 + 3)
# ---------------------------

Feature_CI_Selected <-
  (wrap_plots(plot_list[1:4], ncol = 4)) /
  (wrap_plots(plot_list[5:7], ncol = 3)) +
  plot_annotation(
    tag_levels = "A"   # Adds A–G only
  ) &
  theme(
    plot.tag = element_text(
      size = 14,
      face = "bold"
    ),
    plot.tag.position = c(0, 1)  # top-left
  )

Feature_CI_Selected


# High-resolution export
ggsave(
  filename = file.path(output, "Figure_2.pdf"),
  plot = Feature_CI_Selected,
  width = 14,
  height = 8,
  dpi = 1000,
  device = "pdf"
)

#==== End of Figure 2 ===========#

#======================================================#
# Compute cytotoxic threshold (Estimate from controls) #
#======================================================#

# Extract replicate-level control CI values
control_vals <- countsMetadata %>%
  dplyr::filter(dose == "control") %>%
  dplyr::pull(Cytotoxicity_Index_scaled)

# Compute threshold: mean + 2 SD
mean_ctrl <- mean(control_vals, na.rm = TRUE)
sd_ctrl   <- sd(control_vals, na.rm = TRUE)
cytotox_threshold <- mean_ctrl + 2 * sd_ctrl

cat("Control mean CI:", round(mean_ctrl, 3), "\n")
cat("Control SD:", round(sd_ctrl, 3), "\n")
cat("Cytotoxic threshold (mean + 2SD):",
    round(cytotox_threshold, 3), "\n")

# classify compound–dose conditions
countsMetadata_avg <- countsMetadata_avg %>%
  dplyr::mutate(
    Cytotox_Class =
      ifelse(Cytotoxicity_Index > cytotox_threshold,
             "Cytotoxic",
             "Non-cytotoxic")
  )

table(countsMetadata_avg$Cytotox_Class)

p_CI_threshold <- p_CI +
  geom_hline(yintercept = cytotox_threshold,
             linetype = "dashed",
             color = "red",
             linewidth = 1)

p_CI_threshold

#===================================================================================# 
# Can morphological features alone identify transcriptionally cytotoxic conditions? #
#===================================================================================#
# PCA of morphology features colored by Cytotox_Class  

# Attach Cytotox_Class to IntegratedData (join by compound and dose)
IntegratedData_class <- IntegratedData %>%
  dplyr::left_join(
    countsMetadata_avg %>% dplyr::select(compound, dose, Cytotox_Class),
    by = c("compound" = "compound", "dose_level" = "dose")
  )

# Subset to selected morphological features only
morph_mat <- IntegratedData_class %>%
  dplyr::select(dplyr::all_of(selected_features)) %>%
  as.matrix()

morph_pca <- prcomp(morph_mat, scale. = TRUE)

# Variance explained
var_explained <- (morph_pca$sdev^2) / sum(morph_pca$sdev^2) * 100
pc1_label <- paste0("PC1 (", round(var_explained[1], 1), "%)")
pc2_label <- paste0("PC2 (", round(var_explained[2], 1), "%)")

# Build plotting data frame
morph_pca_df <- as.data.frame(morph_pca$x[, 1:2]) %>%
  dplyr::mutate(
    compound      = IntegratedData_class$compound,
    dose_level    = IntegratedData_class$dose_level,
    Cytotox_Class = IntegratedData_class$Cytotox_Class,
    Label = ifelse(
      compound %in% c("control", "water", "dmso1", "dmso2"),
      compound,
      paste0(compound, " (", dose_level, ")")
    )
  )

# Plot PCA colored by Cytotox_Class
p_morph_pca_class <- ggplot(
  morph_pca_df,
  aes(x = PC1, y = PC2, color = Cytotox_Class)
) +
  geom_point(size = 3.5, alpha = 0.9) +
  ggrepel::geom_text_repel(
    aes(label = Label),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "grey70",
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c("Cytotoxic" = "#D73027", "Non-cytotoxic" = "#4575B4"),
    na.value = "grey60"
  ) +
  theme_minimal(base_size = 14) +
  labs(
    #title = "PCA of Selected Morphological Features",
    x = pc1_label,
    y = pc2_label,
    color = "Cytotoxicity class"
  ) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

p_morph_pca_class

#==============================
# Supplementary Figure3c
#========================
# ==========================================================#
# QUESTION: Is morphological magnitude a reliable indicator #
# of transcriptional cytotoxicity?
# ==========================================================#

# --- 1. Define control set and compute centroid in z-scored morphology space ---
control_ids <- c("water", "dmso1", "dmso2", "control")

ctrl_rows <- IntegratedData %>%
  filter(tolower(compound) %in% control_ids)

ctrl_centroid <- ctrl_rows %>%
  dplyr::select(all_of(selected_features)) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
  as.numeric()

# --- 2. Euclidean distance of every condition from the control centroid ---
feat_mat   <- as.matrix(IntegratedData[, selected_features])
morph_dist <- sqrt(rowSums(sweep(feat_mat, 2, ctrl_centroid, "-")^2))

QuestionA <- IntegratedData %>%
  mutate(
    morph_dist = morph_dist,
    label      = paste0(tolower(compound), "_", tolower(dose_level))
  )

# --- 3. Treatment-only set and correlation against CI ---
treat <- QuestionA %>% filter(!tolower(compound) %in% control_ids)

pear  <- cor.test(treat$morph_dist, treat$Cytotoxicity_Index, method = "pearson")
spear <- cor.test(treat$morph_dist, treat$Cytotoxicity_Index, method = "spearman")

cat("Pearson  r   =", round(pear$estimate, 3),
    " p =", signif(pear$p.value, 3), "\n")
cat("Spearman rho =", round(spear$estimate, 3),
    " p =", signif(spear$p.value, 3), "\n")

# --- 4. Flag the PC1-driver outliers ---
flag_ids <- c("retinoic_acid_mid", "retinoic_acid_high",
              "hydroxyurea_high",   "glycolic_acid_high")

treat <- treat %>%
  mutate(is_flag = label %in% flag_ids)

# --- 5. CI threshold from controls (same rule used for Cytotox_Class) ---
threshold_ci <- mean(ctrl_rows$Cytotoxicity_Index, na.rm = TRUE) +
                2 * sd(ctrl_rows$Cytotoxicity_Index,  na.rm = TRUE)

# Recompute Cytotox_Class directly on treat so the plot has the column
treat <- treat %>%
  mutate(
    Cytotox_Class = ifelse(Cytotoxicity_Index >= threshold_ci,
                           "Cytotoxic", "Non-cytotoxic"),
    Cytotox_Class = factor(Cytotox_Class,
                           levels = c("Non-cytotoxic", "Cytotoxic"))
  )
# --- 6. Plot: morphological distance vs CI ---
p_qA <- ggplot(treat, aes(x = morph_dist, y = Cytotoxicity_Index)) +
  geom_hline(yintercept = threshold_ci, linetype = "dashed", color = "grey50") +
  geom_smooth(method = "lm", se = TRUE, color = "black", alpha = 0.12) +
  geom_point(aes(color = Cytotox_Class), size = 3, alpha = 0.85) +
  geom_text_repel(
    data        = filter(treat, is_flag),
    aes(label   = label),
    box.padding = 0.6, max.overlaps = Inf, size = 3.5
  ) +
  scale_color_manual(values = c("Cytotoxic"     = "#D62728",
                                "Non-cytotoxic" = "#1F77B4")) +
  annotate("text",
           x = Inf, y = -Inf, hjust = 1.05, vjust = -0.6, size = 3.3,
           label = sprintf("Pearson r = %.2f (p = %.2g)\nSpearman ρ = %.2f (p = %.2g)",
                           pear$estimate,  pear$p.value,
                           spear$estimate, spear$p.value)) +
  labs(x = "Morphological distance from control centroid",
       y = "Cytotoxicity Index (CI)",
       color = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position  = "top",
        panel.grid.minor = element_blank())

print(p_qA)
ggsave("QuestionA_morph_distance_vs_CI.pdf", p_qA, width = 6, height = 5)

#===============================================#
