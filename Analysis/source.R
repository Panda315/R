# Specify the file path of the CSV file
file_path <- "/Users/panda/Desktop/Hypothalamus_duplicates-removed-1.csv"

# Read the CSV file into a data frame
data <- read.csv(file_path)

# Check the structure of the imported data
str(data)

# Creating metadata
metadata <- data.frame(
  Samples = c(
    "HYP_CD_LD024_1", "HYP_CD_LD024_2", "HYP_CD_LD024_3",
    "HYP_CD_LD816_1", "HYP_CD_LD816_2", "HYP_CD_LD816_3",
    "HYP_CD_LD1212_1", "HYP_CD_LD1212_2", "HYP_CD_LD1212_3",
    "HYP_CD_LD168_1", "HYP_CD_LD168_2", "HYP_CD_LD168_3",
    "HYP_CD_LD240_1", "HYP_CD_LD240_2", "HYP_CD_LD240_3",
    "HYP_ABX_LD024_1", "HYP_ABX_LD024_2", "HYP_ABX_LD024_3",
    "HYP_ABX_LD816_1", "HYP_ABX_LD816_2", "HYP_ABX_LD816_3",
    "HYP_ABX_LD1212_1", "HYP_ABX_LD1212_2", "HYP_ABX_LD1212_3",
    "HYP_ABX_LD168_1", "HYP_ABX_LD168_2", "HYP_ABX_LD168_3",
    "HYP_ABX_LD240_1", "HYP_ABX_LD240_2", "HYP_ABX_LD240_3"
  ),
  LD_Cycle = c(
    "LD0/24", "LD0/24", "LD0/24",
    "LD8/16", "LD8/16", "LD8/16",
    "LD12/12", "LD12/12", "LD12/12",
    "LD16/8", "LD16/8", "LD16/8",
    "LD24/0", "LD24/0", "LD24/0",
    "LD0/24", "LD0/24", "LD0/24",
    "LD8/16", "LD8/16", "LD8/16",
    "LD12/12", "LD12/12", "LD12/12",
    "LD16/8", "LD16/8", "LD16/8",
    "LD24/0", "LD24/0", "LD24/0"
  ),
  Treatment = c(
    rep("CD", 15),
    rep("ABX", 15)
  ),
  group = c(
    "CD_LD0/24", "CD_LD0/24", "CD_LD0/24", "CD_LD8/16", "CD_LD8/16",
    "CD_LD8/16", "CD_LD12/12", "CD_LD12/12", "CD_LD12/12", "CD_LD16/8",
    "CD_LD16/8", "CD_LD16/8", "CD_LD24/0", "CD_LD24/0", "CD_LD24/0",
    "ABX_LD0/24", "ABX_LD0/24", "ABX_LD0/24", "ABX_LD8/16", "ABX_LD8/16",
    "ABX_LD8/16", "ABX_LD12/12", "ABX_LD12/12", "ABX_LD12/12", "ABX_LD16/8",
    "ABX_LD16/8", "ABX_LD16/8", "ABX_LD24/0", "ABX_LD24/0", "ABX_LD24/0"
  )
)

# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)

# Exclude the first column from data
countData <- data[, -1]

# Convert variables to factors
metadata$LD_Cycle <- factor(metadata$LD_Cycle)
metadata$Treatment <- factor(metadata$Treatment)

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = metadata,
                              design = ~ LD_Cycle + Treatment)

# Pre-filtering to keep raw counts above 10
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Perform normalization
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

# Extract normalized counts
normalizedmain_data <- counts(dds, normalized = TRUE)
dim(normalizedmain_data)
head(normalizedmain_data)

# Performing log transformation
log <- vst(dds, blind = TRUE)

# Perform PCA
pca <- prcomp(t(assay(log)))

# Create a data frame for the PCA plot
pca_data <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], LD_Cycle = metadata$LD_Cycle, Treatment = metadata$Treatment)

# Create the PCA plot
p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = LD_Cycle, shape = Treatment)) +
  geom_point() +
  geom_text(aes(label = rownames(pca_data)), size = 3, nudge_x = 0.1, nudge_y = 0.1) +
  labs(x = "PC1", y = "PC2", color = "LD Cycle", shape = "Treatment") +
  theme_minimal()

# Print the PCA plot
print(p)

# Running differential expression analysis
dds <- DESeq(dds)
resultsNames(dds)

# Exploring Log2 fold change
l2f <- results(dds, contrast = c("Treatment", "ABX", "CD"), alpha = 0.05)
summary(l2f)

# MA plot
plotMA(l2f, ylim = c(-8, 8))

# Assign row names to the gene names column
rownames(l2f) <- l2f$gene_name
# Remove the redundant gene name column
l2f$gene_name <- NULL

# Create a new column for differential expression status
l2f$diffexpressed <- "NO CHANGE"

# Update differential expression status based on log2 fold change and adjusted p-value cutoffs
l2f$diffexpressed[l2f$log2FoldChange > 0.32 & l2f$padj < 0.05] <- "UP"
l2f$diffexpressed[l2f$log2FoldChange < -0.32 & l2f$padj < 0.05] <- "DOWN"

# Subset significantly differentially expressed genes only, excluding NAs
l2f_sig <- l2f[complete.cases(l2f) & l2f$padj < 0.05 & abs(l2f$log2FoldChange) > 0.32, ]

# Subsetting normalized counts to only significant DE genes
sig_norm_counts <- normalizedmain_data[rownames(l2f_sig), ]
sig_norm_counts <- t(sig_norm_counts)

# Convert normalized counts to a matrix
sig_norm_counts <- as.matrix(sig_norm_counts)

# Heatmap
heat_colors <- brewer.pal(6, "YlOrRd")
print(dim(normalizedmain_data))
head(normalizedmain_data)
pheatmap(sig_norm_counts, color = heat_colors, cluster_rows = TRUE, show_rownames = FALSE, annotation = dplyr::select(metadata, LD_Cycle), scale = "row")
