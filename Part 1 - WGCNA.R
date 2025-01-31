# Script: WGCNA pancreatic cancer - body composition
# Date: 2024-01-17
# Author: mkutmon, arina

#===============================================
# Set up environment
#===============================================
# Clear workspace and set string as factors to false
options(stringsAsFactors = F)

#===============================================
# Install and load required packages
#===============================================
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!("WGCNA" %in% installed.packages())) { BiocManager::install("WGCNA") }
if (!("dplyr" %in% installed.packages())) { BiocManager::install("dplyr") }
if (!("rstudioapi" %in% installed.packages())) { install.packages("rstudioapi") }  
if (!("biomaRt" %in% installed.packages())) { BiocManager::install("biomaRt") }
if (!("DESeq2" %in% installed.packages())) { BiocManager::install("DESeq2") }
if (!("ggplot2" %in% installed.packages())) { BiocManager::install("ggplot2") }
if (!("gridExtra" %in% installed.packages())) { BiocManager::install("gridExtra") }
if (!("CorLevelPlot" %in% installed.packages())) { devtools::install_github("kevinblighe/CorLevelPlot") }
if (!("purrr" %in% installed.packages())) { devtools::install_github("purrr") }
if (!("RCy3" %in% installed.packages())) { BiocManager::install("RCy3")}
if (!("clusterProfiler" %in% installed.packages())) { BiocManager::install("clusterProfiler")}
if (!("AnnotationDbi" %in% installed.packages())) { BiocManager::install("AnnotationDbi")}
if (!("org.Hs.eg.db" %in% installed.packages())) { BiocManager::install("org.Hs.eg.db")}
if (!("ComplexHeatmap" %in% installed.packages())) { BiocManager::install("ComplexHeatmap")}
if (!("readr" %in% installed.packages())) { BiocManager::install("readr")}
if (!("gplots" %in% installed.packages())) {BiocManager::install("gplots")}
if (!("survminer" %in% installed.packages())) {install.packages("survminer")}

library(readxl, quietly = TRUE)
library(WGCNA, quietly = TRUE)
library(rstudioapi, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(biomaRt, quietly = TRUE)
library(DESeq2, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(gridExtra, quietly = TRUE)
library(CorLevelPlot, quietly = TRUE)
library(purrr, quietly = TRUE)
library(RCy3, quietly = TRUE)
library(clusterProfiler, quietly = TRUE)
library(AnnotationDbi, quietly = TRUE)
library(org.Hs.eg.db, quietly = TRUE)
library(ComplexHeatmap, quietly = TRUE)
library(readr, quietly = TRUE)
library(gplots, quietly = TRUE)
library(survival)
library(survminer)

#===============================================
# Set working directory
#===============================================
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#===============================================
# Settings
#===============================================
# Specify the TOM type 
TOMType = "unsigned" # the sign of correlation not considered 
pc.filtered = TRUE

# Set the min module size
min.module = 100

# Create the run name
run.name <- paste0(TOMType,"_filtered_power_8_min_",min.module)

# Create an output folder to save results
if (!file.exists(file.path(getwd(), "Output"))) {
  dir.create(file.path(getwd(), "Output"))
}

# Create a folder to save Heatmaps 
if (!file.exists(file.path(getwd(), "Output", "Heatmaps"))) {
  dir.create(file.path(getwd(), "Output", "Heatmaps"))
}

# Create three subfolders within the "Heatmaps" folder
# to store heatmaps with all genes 
if (!file.exists(file.path(getwd(), "Output", "Heatmaps", "All_genes"))) {
  dir.create(file.path(getwd(), "Output", "Heatmaps", "All_genes"))
}
# to store heeatmaps for hub genes alone (kME > 0.9)
if (!file.exists(file.path(getwd(), "Output", "Heatmaps", "Hub_genes_kME>0.9"))) {
  dir.create(file.path(getwd(), "Output", "Heatmaps", "Hub_genes_kME>0.9"))
}
# to store heatmaps for genes with kME > 0.7
if (!file.exists(file.path(getwd(), "Output", "Heatmaps", "kME>0.7"))) {
  dir.create(file.path(getwd(), "Output", "Heatmaps", "kME>0.7"))
}

#===============================================
# Read the trait data file 
#===============================================
traitData <- read.table("metadata.txt", header = T, sep="\t")
# Convert age to numeric 
traitData$Age <- as.numeric(traitData$Age)
# Binarize categorical variables 
traitData$Sex <- ifelse(traitData$Sex == "M", 0, 1)
# Check the number of traits 
ntraits <- ncol(traitData)-1

#===============================================
# Read and process data file
#===============================================
# Access the Ensembl database and retrieve information about genes in the human genome 
ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Read data file
data <- read.table("data.txt", header=TRUE) 
rownames(data) <- data[,1]
data <- data[,-1]

# Remove samples without trait information
data <- data[, colnames(data) %in% traitData$SampleID]
data_rounded <- round(data)  # Round to the nearest integer

# Normalization DeSeq2
group <- as.data.frame(colnames(data_rounded))
# Organize the data into a suitable format for subsequent analysis by DESeq2 - combine count data with sample metadata (SampleID)
dds <- DESeqDataSetFromMatrix(countData = as.matrix(data_rounded), colData = group, design = ~ 1) 
# Perform differential expression analysis (based on the negative binomial generalized linear model)
dds <- DESeq(dds)
# Extract normalized count data
normalized_counts <- counts(dds, normalized = TRUE)

# Calculate percentage of zeros in each row
percentage_zeros <- apply(normalized_counts, 1, function(x) mean(x == 0) * 100)
hist(percentage_zeros)
# Save the histogram
png(file.path(getwd(), "Output", "histogram_percentage_zeros_all_genes.png"))
hist(percentage_zeros)
dev.off()

# Calculate row sums
row_sums <- rowSums(normalized_counts)
# Filter rows: keep rows with less than 90 % zeros and row sum of 20 or more
normalized_counts_filt <- as.data.frame(normalized_counts[percentage_zeros < 90 & row_sums >= 20, ])
# Log transform the data
log_transformed_data <- as.data.frame(log2(normalized_counts_filt+1))
# Density plot for log-transformed expression data (filtered)
png(file.path(getwd(), "Output", "Density_plot_all_genes_no_lowly-expressed_genes.png"))
plot(density(as.matrix(log_transformed_data)))
dev.off()
# Boxplot for log-transformed expression data (filtered) 
png(file.path(getwd(), "Output", "Boxplot_log_transformed_expression_data_no_lowly-expressed_genes.png"))
boxplot(log_transformed_data, las=2, main="Boxplot of Log-transformed Expression Data")
dev.off()

# PCA 
pca_res <- prcomp(t(log_transformed_data))
# plot PCA
png(file.path(getwd(), "Output", "PCA_expression_data_all_genes.png"), width = 4000, height = 3000, res = 300)
pc1_var <- var(pca_res$x[,1])/sum(var(pca_res$x))
pc2_var <- var(pca_res$x[,2])/sum(var(pca_res$x))
plot(pca_res$x[,1], pca_res$x[,2], 
     xlab = paste0("PC1 (", (round(pc1_var, digits = 3)*100), "% variance)"), 
     ylab = paste0("PC2 (", (round(pc2_var, digits = 3)*100), "% variance)"),
     main = "PCA of Expression Data", cex = 0.5)
text(pca_res$x[,1], pca_res$x[,2], labels=colnames(log_transformed_data), pos=4, cex = 0.5)
dev.off()

# Filter out non-protein-coding genes
if (pc.filtered == TRUE) {
  protein_coding_genes <- biomaRt::getBM(attributes = c("ensembl_gene_id"), filters = c("chromosome_name","biotype"),
                                         values = list(c(as.character(1:22), "X", "Y", "MT"),"protein_coding"), mart = ensembl)
  log_transformed_data_filt <- log_transformed_data[rownames(log_transformed_data) %in% protein_coding_genes$ensembl_gene_id,]
  
  # Density plot for protein-coding genes only 
  png(file.path(getwd(), "Output", "Density_plot_protein-coding_genes.png"))
  plot(density(as.matrix(log_transformed_data_filt)))
  dev.off()
  
  # Boxplot for protein-coding genes only 
  png(file.path(getwd(), "Output", "Boxplot_protein-coding_genes.png"))
  boxplot(log_transformed_data_filt, las=2, main="Boxplot of Log-transformed Expression Data")
  dev.off()
  
  # PCA for protein-coding genes only 
  pca_res <- prcomp(t(log_transformed_data_filt))
  png(file.path(getwd(), "Output", "PCA_plot_protein-coding_genes.png"), width = 5000, height = 3500, res = 300)
  pc1_var <- var(pca_res$x[,1])/sum(var(pca_res$x))
  pc2_var <- var(pca_res$x[,2])/sum(var(pca_res$x))
  plot(pca_res$x[,1], pca_res$x[,2], 
       xlab = paste0("PC1 (", (round(pc1_var, digits = 3)*100), "% variance)"), 
       ylab = paste0("PC2 (", (round(pc2_var, digits = 3)*100), "% variance)"),
       main = "PCA of Expression Data", cex = 0.5)
  text(pca_res$x[,1], pca_res$x[,2], labels=colnames(log_transformed_data), pos=4, cex = 0.4)
  dev.off()
} else {
  log_transformed_data_filt <- log_transformed_data
}

# check for outliers
gsg <- goodSamplesGenes(t(log_transformed_data_filt))
gsg$allOK
table(gsg$goodGenes)
table(gsg$goodSamples)

# sample tree
sampleTree <- hclust(dist(t(log_transformed_data_filt)), method="average")
plot(sampleTree)

#===============================================
# Network construction and module detection
#===============================================
# Convert the data into the appropriate format for WGCNA
wgcna_data <- t(log_transformed_data_filt)


# Choose a set of soft thresholding powers
powers <- seq(2,20, 1)

# Analysis of scale-free topology for multiple soft thresholding powers
sft <- pickSoftThreshold(wgcna_data, 
                         powerVector = powers,
                         networkType = "unsigned",
                         verbose = 5)
sft.data <- sft$fitIndices


# Scale-free topology model fit plot
plot_1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) + 
  geom_point() + 
  geom_text(nudge_y = 0.1) + 
  geom_hline(yintercept = 0.8, color = "red") +
  labs(x = 'Power', y = 'Scale Free Topology Model Fit, Signed R^2') +
  theme_classic()

# Main connectivity plot  
plot_2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 100) +
  labs(x = 'Power', y = 'Main Connectivity') +
  theme_classic()

png(file.path(getwd(), "Output", "Pick_soft_threshold_protein-coding_genes.png"))
gridExtra::grid.arrange(plot_1, plot_2, nrow = 2)
dev.off()

# Power 8 selected as the soft threshold 
soft_power <- 8
# Record the start time of the WGCNA 
start_time <- Sys.time()
# Import the 'cor' function from the WGCNA package 
cor <- WGCNA::cor

# Run blockwiseModules function 
bwnet <- blockwiseModules(wgcna_data, power = soft_power,
                          TOMType = TOMType, 
                          maxBlockSize = 17000, 
                          minModuleSize = min.module, 
                          numericLabels = FALSE, 
                          verbose = 3)

# Change name based on settings
saveRDS(bwnet, file.path("Output", paste0(run.name,".RDS")))
# Assign the 'cor' function to the stats package
cor <- stats::cor
# Record the end time of the WGCNA 
end_time <- Sys.time()
# Calculate the duration of the WGCNA
duration <- end_time - start_time
# Print the duration
print(paste("Duration:", duration))

#===============================================
# Investigate modules
#===============================================
# Retrieve eigengenes
module_eigengenes <- bwnet$MEs
table(bwnet$colors)

# Convert labels to colors for plotting
mergedColors <- labels2colors(bwnet$colors)

# Plot the dendrogram and the module colors underneath
png(file.path(getwd(), "Output", "Cluster_Dendrogram_protein-coding_genes.png"))
plotDendroAndColors(bwnet$dendrograms[[1]], mergedColors[bwnet$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Visualize module-trait associations for all modules
heatmap.data <- merge(module_eigengenes, traitData, by.x = 'row.names', by.y="SampleID")
heatmap.data <- heatmap.data %>% tibble::column_to_rownames("Row.names")
sizeGrWindow(12, 68)
png(file.path(getwd(), "Output", paste0("Heatmap_module-trait_all_modules_protein-coding_genes_", run.name, ".png")), width = 6000, height = 3000, res = 300)
CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[(ncol(heatmap.data)-ntraits+1):ncol(heatmap.data)],
             y = names(heatmap.data)[1:(ncol(heatmap.data)-ntraits)],
             col = c("blue1", "skyblue", "white", "pink", "red"))
dev.off()

# Visualize module-trait associations for significant modules only
source("make_heatmap_module_trait_relationship.R")
png(file.path(getwd(), "Output", paste0("Heatmap_module-trait_significant_modules_protein-coding_genes_", run.name, ".png")), width = 6000, height = 3000, res = 300)
significant_modules <- CorLevelPlotSig(heatmap.data,
                                       x = names(heatmap.data)[(ncol(heatmap.data)-ntraits+1):ncol(heatmap.data)],
                                       y = names(heatmap.data)[1:(ncol(heatmap.data)-ntraits)],
                                       col = c("blue1", "skyblue", "white", "pink", "red"))
significant_modules
dev.off()

#===============================================
# Retrieve module names and rearrange trait data
#===============================================

# Retrieve names of all modules 
gene_module_key_trait <- tibble::enframe(bwnet$colors, name = "gene", value = "module") %>%
  dplyr::mutate(module = paste0("ME", module))

# Get the trait data
traitData_df <- read.table("metadata.txt", header = T, sep="\t")
traitData_df$Hospital <- ifelse(grepl("^PANCO", traitData$SampleID), "MUMC", "UMCU")
traitData_df <- traitData_df %>%
  mutate(Subtype = ifelse(ICGCrnaseq < 0, "basal", "classical")) 

#===============================================
# Calculation of eigengene-based connectivity
#===============================================

# Calculate eigengene-based connectivity for identification of hub genes 
kMEs <- cor(wgcna_data, module_eigengenes, method = "pearson")
kMEs <- as.data.frame(kMEs)

#===============================================
# Age distribution across samples
#===============================================
age_dist_data <- subset(traitData, select = c(SampleID, Age))
ggplot(age_dist_data, aes(x = Age)) +
  geom_density(fill = "pink") +
  theme_minimal() +
  labs(title = "Age Distribution",
       x = "Age",
       y = "Density") +
  theme(plot.title = element_text(hjust = 0.5))
ggplot(age_dist_data, aes(x = SampleID, y = Age)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() +
  labs(title = "Age of Patients", x = "Patient", y = "Age") +
  theme(plot.title = element_text(hjust = 0.5))
