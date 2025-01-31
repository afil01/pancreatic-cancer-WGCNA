# Pancreatic Cancer - WGCNA and functional analysis on identified modules 
This repository contains the analysis pipeline for identifying key gene modules and performing functional analysis on RNA-seq data. The main method employed is weighted gene co-expression network analysis (WGCNA) to determine modules significantly associated with disease progression and phenotypic traits. 

# Part 1 - WGCNA 
The first part of the code performs all steps of the standard WGCNA analysis to identify significant modules of co-expressed genes. This includes pre-processing the data (filtering functionally relevant genes), constructing the WGCNA network, identifying gene modules using hierarchical clustering, calculating module-trait relationships, and producing a heatmap to visualise module-trait correlations. 

# Part 2 - Trait Analysis per Module 
The second part of the code focuses on identifying hub genes and strongly correlated genes based on intermodular connectivity, gene ontology (GO) enrichemnt analysis on all genes, strongly correlated genes, and hub genes, visualisation of the enriched pathways using bar plots, construction of PPI networks, visualisation of significant gene modules with heatmaps (on relative and absolute gene expression), and generating Kaplan-Meier survival curves based on high or low eigengene expression. 
