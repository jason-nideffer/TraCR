
#' Calculate Clonotype Frequencies
#'
#' This function takes a seurat object with meta data that already includes
#' TCR annotations and calculates the frequencies of different clonotypes
#' on a per-sample basis.
#'
#' @param seurat_object seurat object to operate on
#' @param clonotype_column Character describing name of column to use for clone calling
#' @param sample_column Character describing name of column to use for sample discrimination
#' @param identity_column Character describing name of column to use for clonotype characterization
#' @return seurat_object with meta data updated to include frequency data
#' @export
calculate_freq <- function(seurat_object, clonotype_column, sample_column, identity_column=NULL) {
  
  # Store rownames and create new metadata object
  meta <- seurat_object@meta.data
  rownames <- row.names(meta)
  
  # Make order_column so order can be preserved and rownames reassigned
  meta$order_column <- 1:nrow(meta)
  
  # Incorporate total cell counts into meta
  cells_per_sample <- meta %>% 
    group_by(pick(all_of(sample_column))) %>% 
    summarise(cells_per_sample=n())
  
  meta <- merge(meta, cells_per_sample, by=c(sample_column))
  
  # Incorporate clonotype frequencies
  clonotype_count_per_sample <- meta %>% 
    group_by(pick(c(sample_column,clonotype_column))) %>% 
    summarise(clonotype_count_per_sample=n())
  
  meta <- merge(meta, clonotype_count_per_sample, 
                by=c(sample_column, clonotype_column), 
                all.x = TRUE)
  
  meta$clonotype_freq_per_sample <-  meta$clonotype_count_per_sample / meta$cells_per_sample
  
  # Calculate Clone Identity Modes
  clonotype_identity_modes <- meta %>%
    group_by(pick(c(clonotype_column,identity_column))) %>%
    summarise(n_clono=n()) %>%
    group_by(pick(clonotype_column))%>%
    slice_max(n=1, order_by = n_clono, na_rm = TRUE, with_ties = FALSE)
  
  clonotype_identity_modes <- select(clonotype_identity_modes, -c(n_clono))
  colnames(clonotype_identity_modes) <- c(clonotype_column, "clonotype_identity_mode")
  
  meta <- merge(meta, clonotype_identity_modes, by=c(clonotype_column))
  
  # Restore metadata order
  meta <- meta[order(meta$order_column), ]
  
  # Remove order column
  drops <- c("order_column")
  meta <- meta[ , !(names(meta) %in% drops)]
  
  # Add back rownames
  row.names(meta) <- rownames
  
  # Reassign metadata
  seurat_object@meta.data <- meta
  
  return (seurat_object)
}
