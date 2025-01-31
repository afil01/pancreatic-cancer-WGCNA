make_module_heatmap <- function(module_name,
                                expression_mat = wgcna_data,
                                metadata_df = traitData_df,
                                gene_module_key_df = gene_module_key_trait,
                                module_eigengenes_df = module_eigengenes,
                                show_row_names = FALSE) {
  # Create a summary heatmap of a given module.
  #
  # Args:
  # module_name: a character indicating what module should be plotted, e.g. "MEdarkgreen"
  # expression_mat: The full gene expression matrix. Default is `normalized_counts`.
  # metadata_df: a data frame with new_Sample_ID and time_point
  #              as columns. Default is `metadata`.
  # gene_module_key: a data.frame indicating what genes are a part of what modules. Default is `gene_module_key`.
  # module_eigengenes: a sample x eigengene data.frame with samples as row names. Default is `module_eigengenes`.
  #
  # Returns:
  # A heatmap of expression matrix for a module's genes, with a barplot of the
  # eigengene expression for that module.
  
  # Set up the module eigengene with its SampleID
 module_eigengene <- module_eigengenes_df %>%
    dplyr::select(all_of(module_name)) %>%
    tibble::rownames_to_column("SampleID") 
  
  # Set up column annotation from metadata
  col_annot_df <- metadata_df %>%
    # Only select the treatment and sample ID columns
    dplyr::select(SampleID, Hospital, Subtype) %>%
    # Add on the eigengene expression by joining with sample IDs
    dplyr::inner_join(module_eigengene, by = "SampleID") %>%
    # Arrange by patient and time point
    dplyr::arrange(Subtype) %>%
    # Store sample
    tibble::column_to_rownames("SampleID")
 
  # Set up column annotation from metadata
  
  col_score_df <- metadata_df %>%
    # Only select the treatment and sample ID columns
    dplyr::select(SampleID, Subtype, ICGCrnaseq, SMI.Z.score) %>%
    # Arrange by patient and time point
    dplyr::arrange(Subtype) %>%
    # Store sample
    tibble::column_to_rownames("SampleID")
  
  # Create the ComplexHeatmap column annotation object
  col_annot <- ComplexHeatmap::HeatmapAnnotation(
    # Supply treatment labels
    Subtype = col_annot_df$Subtype,
    Hospital = col_annot_df$Hospital,
    
    # Add annotation barplot
    module_eigengene = ComplexHeatmap::anno_barplot(dplyr::select(col_annot_df, module_name)),
    ICGCrnaseq = ComplexHeatmap::anno_barplot(dplyr::select(col_score_df, ICGCrnaseq)),
    SMI.Z.score = ComplexHeatmap::anno_barplot(dplyr::select(col_score_df, SMI.Z.score)),
    
    
    # Pick colors for each experimental group in subtype
    col = list(Hospital = c("UMCU" = "#f1a340", "MUMC" = "#998ec3"), Subtype = c("classical" = "#f1a340", "basal" = "#998ec3"))
  )
  
  # Get a vector of the Ensembl gene IDs that correspond to this module
  module_genes <- gene_module_key_df %>%
    dplyr::filter(module == module_name) %>%
    dplyr::pull(gene)
  
  # Set up the gene expression data frame
  mod_mat <- expression_mat %>%
    t() %>%
    as.data.frame() %>%
    # Only keep genes from this module
    dplyr::filter(rownames(.) %in% module_genes) %>%
    # Order the samples to match col_annot_df
    dplyr::select(rownames(col_annot_df)) %>%
    # Data needs to be a matrix
    as.matrix()
  
  # Normalize the gene expression values
  mod_mat <- mod_mat %>%
    # Scale can work on matrices, but it does it by column so we will need to
    # transpose first
    t() %>%
    scale() %>%
    # And now we need to transpose back
    t()
  
  # Create a color function based on standardized scale
  color_func <- circlize::colorRamp2(
    c(-2, 0, 2),
    c("#67a9cf", "#f7f7f7", "#ef8a62")
  )
  
  # Plot on a heatmap
  heatmap <- ComplexHeatmap::Heatmap(mod_mat,
                                     name = module_name,
                                     # Supply color function
                                     col = color_func,
                                     # Supply column annotation
                                     bottom_annotation = col_annot,
                                     # We don't want to cluster samples
                                     cluster_columns = TRUE,
                                     # We don't need to show sample or gene labels
                                     show_row_names = show_row_names,
                                     show_column_names = TRUE,
                                     column_names_side = "top",
                                     row_names_side = "left"
  )
  
  # Return heatmap
  return(heatmap)
}