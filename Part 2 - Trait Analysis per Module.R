# Script: WGCNA pancreatic cancer - body composition - script per module 
# Date: 2024-07-05
# Author: mkutmon, arina

# Complete WGCNA 
# Identify significant modules 

#===============================================================================
# Provide module name and color name (e.g., "MEmagenta" and "magenta")
#===============================================================================

module.name <- "MEcyan"
color.name <- "cyan"
#===============================================================================
# Make a directory to store the results 
#===============================================================================  
if (!file.exists(file.path(getwd(), "Module Trait Analysis"))) {
  dir.create(file.path(getwd(), "Module Trait Analysis"))
}

dir.create(file.path(getwd(), "Module Trait Analysis", module.name))

#===============================================================================
# Filter genes in the module of interest 
#===============================================================================

module_genes <- gene_module_key_trait$gene[gene_module_key_trait$module == module.name]

#===============================================================================
# Identification of hub genes - kME > 0.9
#===============================================================================

module_kMEs <- kMEs[module_genes, module.name, drop = FALSE]
hub_genes_filt <- module_genes[abs(module_kMEs) > 0.9]
wgcna_data_hub_0.9 <- wgcna_data[, hub_genes_filt]
gene_module_key_trait_hub_0.9 <- gene_module_key_trait[gene_module_key_trait$gene %in% hub_genes_filt, ]

#===============================================================================
# Identification of strongly correlated genes - kME > 0.7
#===============================================================================

st_cor_genes_filt <- module_genes[abs(module_kMEs) > 0.7]
wgcna_data_st_cor_0.7 <- wgcna_data[, st_cor_genes_filt]
gene_module_key_trait_st_cor_0.7 <- gene_module_key_trait[gene_module_key_trait$gene %in% st_cor_genes_filt, ]

#===============================================================================
# Heatmaps with z-scores - all genes
#===============================================================================

# Load source file 
source("make_module_heatmap_for_each_module.R")

png_file <- file.path(getwd(), "Module Trait Analysis", module.name, paste0(module.name, "_z.score_expression_all_genes.png"))
png(filename = png_file, width = 1500, height = 1000)
# Heatmap of z-score expression for all genes 
heatmap_z_score_all <- make_module_heatmap(module_name = module.name)
heatmap_z_score_all
dev.off() 

#===============================================================================
# Heatmaps with z-scores - hub genes (kME > 0.9)
#===============================================================================

# Provide file path to save the heatmap
png_file <- file.path(getwd(), "Module Trait Analysis", module.name, paste0(module.name, "_z.score_expression_hub_0.9.png"))
png(filename = png_file, width = 1500, height = 1000)
# Heatmap of z-score expression for hub genes 
heatmap_z_score_hub <- make_module_heatmap(module_name = module.name, 
                                       expression_mat = wgcna_data_hub_0.9, 
                                       gene_module_key_df = gene_module_key_trait_hub_0.9, 
                                       show_row_names = TRUE)
heatmap_z_score_hub
dev.off()  

#===============================================================================
# Heatmaps with z-scores - strongly-correlated genes (kME > 0.7)
#===============================================================================

# Provide file path to save the heatmap
png_file <- file.path(getwd(), "Module Trait Analysis", module.name, paste0(module.name, "_z.score_expression_st_cor_0.7.png"))
png(filename = png_file, width = 1500, height = 1000)
# Heatmap of z-score expression for hub genes 
heatmap_z_score_st_cor_0.7 <- make_module_heatmap(module_name = module.name, 
                                           expression_mat = wgcna_data_st_cor_0.7, 
                                           gene_module_key_df = gene_module_key_trait_st_cor_0.7, 
                                           show_row_names = FALSE)
heatmap_z_score_st_cor_0.7
dev.off()  # Turn off the PNG device

#===============================================================================
# Heatmaps with abs expression - all genes 
#===============================================================================

source("make_module_heatmap_for_each_module_absolute_expression.R")
png_file <- file.path(getwd(), "Module Trait Analysis", module.name, paste0(module.name, "_absolute_expression_all_genes.png"))
png(filename = png_file, width = 1500, height = 1000)
# Heatmap of absolute expression for all genes 
heatmap_abs_all <- make_module_heatmap(module_name = module.name)
heatmap_abs_all
dev.off() 


#===============================================================================
# Heatmaps with abs expression - hub genes (kME > 0.9)
#===============================================================================

# Provide file path to save the heatmap
png_file <- file.path(getwd(), "Module Trait Analysis", module.name, paste0(module.name, "_absolute_expression_hub_0.9.png"))
png(filename = png_file, width = 1500, height = 1000)
# Heatmap of z-score expression for hub genes 
heatmap_abs_hub <- make_module_heatmap(module_name = module.name, 
                                           expression_mat = wgcna_data_hub_0.9, 
                                           gene_module_key_df = gene_module_key_trait_hub_0.9, 
                                           show_row_names = TRUE)
heatmap_abs_hub
dev.off()  

#===============================================================================
# Heatmaps with abs expression - strongly-correlated genes (kME > 0.7)
#===============================================================================

# Provide file path to save the heatmap
png_file <- file.path(getwd(), "Module Trait Analysis", module.name, paste0(module.name, "_absolute_expression_st_cor_0.7.png"))
png(filename = png_file, width = 1500, height = 1000)
# Heatmap of z-score expression for hub genes 
heatmap_abs_st_cor_0.7 <- make_module_heatmap(module_name = module.name, 
                                                  expression_mat = wgcna_data_st_cor_0.7, 
                                                  gene_module_key_df = gene_module_key_trait_st_cor_0.7, 
                                                  show_row_names = FALSE)
heatmap_abs_st_cor_0.7
dev.off() 

#===============================================================================
# GO enrichment analysis - all genes
#===============================================================================

# Perform GO enrichment analysis
GO_results_module <- enrichGO(gene = module_genes, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP")

png_file <- file.path(getwd(), "Module Trait Analysis", module.name, paste0(module.name, "_GO_barplot_all_genes.png"))
png(filename = png_file, width = 1500, height = 1000)
# Make barplots
GO_barplots_module <- barplot(GO_results_module, showCategory = 20)
GO_barplots_module
dev.off()  

#===============================================================================
# GO enrichment analysis - hub genes (kME > 0.9)
#===============================================================================

# Perform GO enrichment analysis
GO_results_module_hub <- enrichGO(gene = hub_genes_filt, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP")

png_file <- file.path(getwd(), "Module Trait Analysis", module.name, paste0(module.name, "_GO_barplot_hub_genes_0.9.png"))
png(filename = png_file, width = 1500, height = 1000)
# Make barplots
GO_barplots_module_hub <- barplot(GO_results_module_hub, showCategory = 20)
GO_barplots_module_hub
dev.off()  

#===============================================================================
# GO enrichment analysis - hub genes (kME > 0.7)
#===============================================================================

# Perform GO enrichment analysis
GO_results_module_st_cor <- enrichGO(gene = st_cor_genes_filt, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP")

png_file <- file.path(getwd(), "Module Trait Analysis", module.name, paste0(module.name, "_GO_barplot_st_cor_genes_0.7.png"))
png(filename = png_file, width = 1500, height = 1000)
# Make barplots
GO_barplots_module_st_cor <- barplot(GO_results_module_st_cor, showCategory = 20)
GO_barplots_module_st_cor
dev.off()  

#===============================================================================
# PPI in Cytoscape - all genes
#===============================================================================
    
# Extract the list of genes and the name for the module
genes <- module_genes
name <- color.name
print(length(genes))

# create PPI network with gene list (might be a problem for very large GO terms)
query <- format_csv(as.data.frame(genes), col_names=F, escape = "double", eol =",")
RCy3::commandsRun(paste0('string protein query cutoff=0.4 species="Homo sapiens" newNetName="module-',name,'" query="',query,'" limit=0'))
suid <- RCy3::getNetworkSuid()
# update the visual style
style.name <- paste0("viz-",name)
RCy3::createVisualStyle(style.name)
RCy3::setNodeLabelMapping("display name", style.name=style.name)
RCy3::lockNodeDimensions(TRUE, style.name=style.name)
RCy3::setNodeShapeDefault("ELLIPSE", style.name = style.name)
RCy3::setNodeColorDefault(gplots::col2hex(name), style.name = style.name)
setVisualStyle(style.name)
RCy3::exportPNG(file.path(getwd(), "Module Trait Analysis", module.name, paste0(gsub(':', '_', name), "-PPI-network_all_genes.png")), zoom = 500)
RCy3::saveSession(file.path(getwd(), "Module Trait Analysis", module.name, paste0("PPI-networks_all_genes.cys")))

#===============================================================================
# PPI in Cytoscape - hub genes (kME > 0.9)
#===============================================================================
# Extract the list of genes and the name for the module
genes <- hub_genes_filt
name <- color.name
print(length(genes))

# create PPI network with gene list (might be a problem for very large GO terms)
query <- format_csv(as.data.frame(genes), col_names=F, escape = "double", eol =",")
RCy3::commandsRun(paste0('string protein query cutoff=0.4 species="Homo sapiens" newNetName="module-',name,'" query="',query,'" limit=0'))
suid <- RCy3::getNetworkSuid()
# update the visual style
style.name <- paste0("viz-",name)

RCy3::createVisualStyle(style.name)
RCy3::setNodeLabelMapping("display name", style.name=style.name)
RCy3::lockNodeDimensions(TRUE, style.name=style.name)
RCy3::setNodeShapeDefault("ELLIPSE", style.name = style.name)
RCy3::setNodeColorDefault(gplots::col2hex(name), style.name = style.name)
setVisualStyle(style.name)
RCy3::exportPNG(file.path(getwd(), "Module Trait Analysis", module.name, paste0(gsub(':', '_', name), "-PPI-network_hub_genes_0.9.png")), zoom = 500)
RCy3::saveSession(file.path(getwd(), "Module Trait Analysis", module.name, paste0("PPI-networks_hub_genes_0.9.cys")))

#===============================================================================
# PPI in Cytoscape - strongly-correlated genes (kME > 0.7)
#===============================================================================

# Extract the list of genes and the name for the module
genes <- st_cor_genes_filt
name <- color.name
print(length(genes))

# create PPI network with gene list (might be a problem for very large GO terms)
query <- format_csv(as.data.frame(genes), col_names=F, escape = "double", eol =",")
RCy3::commandsRun(paste0('string protein query cutoff=0.4 species="Homo sapiens" newNetName="module-',name,'" query="',query,'" limit=0'))
suid <- RCy3::getNetworkSuid()
# update the visual style
style.name <- paste0("viz-",name)
RCy3::setCurrentNetwork(suid)
RCy3::createVisualStyle(style.name)
RCy3::setNodeLabelMapping("display name", style.name=style.name)
RCy3::lockNodeDimensions(TRUE, style.name=style.name)
RCy3::setNodeShapeDefault("ELLIPSE", style.name = style.name)
RCy3::setNodeColorDefault(gplots::col2hex(name), style.name = style.name)
setVisualStyle(style.name)
RCy3::exportPNG(file.path(getwd(), "Module Trait Analysis", module.name, paste0(gsub(':', '_', name), "-PPI-network_st_cor_genes_0.7.png")), zoom = 500)
RCy3::saveSession(file.path(getwd(), "Module Trait Analysis", module.name, paste0("PPI-networks_st_cor_genes_0.7.cys")))

#===============================================================================
# Kaplan-Meier Curves - Survival Data Analysis 
#===============================================================================

# Read survival data
data_survival <- read_excel("survival-data.xlsx") %>% 
  as.data.frame()

# Restructure survival data frame 
data_survival <- data_survival %>% mutate(Death = 1) %>%
  arrange(Survival)
data_survival$PT_ID <- gsub("\\-", "", data_survival$PT_ID)


# Retrieve expression data of the eigengene for the current module 
expr_data <- module_eigengenes[module.name]

# Divide patients into low- and high-expression groups 

# calculate median value from the eigengene profile 
median_value <- median(expr_data[[module.name]])
expr_data$Expression <- ifelse(expr_data[[module.name]] >= median_value, "HIGH", "LOW")
expr_data$Sample <- rownames(expr_data)
rownames(expr_data) <- NULL

# Merge survival data and expression data 
data_survival <- merge(data_survival, expr_data, by.x = "PT_ID", by.y = "Sample")

# Make Kaplan-Meier curve 
fit <- survfit(Surv(Survival, Death) ~ Expression, 
                    data = data_survival,
                   type = "kaplan-meier")

png_file <- file.path(getwd(), "Module Trait Analysis", module.name, paste0(module.name, "- Kaplan-Meier Plot.png"))
png(filename = png_file, width = 1500, height = 1000)
ggsurvplot(fit, 
           data = data_survival,
           pval = TRUE, 
           risk.table = FALSE,
           title = paste("Survival Analysis of", module.name),
           palette = "Dark1",                 
           ggtheme = theme_minimal() + theme(plot.title = element_text(hjust = 0.5), 
                                             legend.title = element_blank()),
           legend.labs = c("High Expression", "Low Expression"), 
           xlab = "Time (Days)",               
           ylab = "Survival Probability",     
           font.title = c(16, "bold"),         
           font.x = c(14, "bold"),                     
           font.y = c(14, "bold"),                     
)
dev.off()

