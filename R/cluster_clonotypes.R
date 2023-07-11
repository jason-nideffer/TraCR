
#' Cluster Clonotypes By Distribution
#'
#' This function plots a UMAP where cells can be colored according to a variable.
#' Multiple clonotypes will also be highlighted.
#'
#' @param seurat_object Seurat object to operate on
#' @param minimum_count Integer specifying minimum number fo clones for a clonotype to be included in analysis
#' @param plot Boolean specifying whether or not to plot a clustered heatmap
#' @return a dataframe with clonotypes as rows and their propotional distribution across celltypes (or the variables contained within 'celltype_column') as columns
#' @export
cluster_clonotypes <- function(seurat_object,
                               minimum_count=5,
                               plot=TRUE) {

  clonotype_column <- seurat_object@misc$clonotype_column
  cell_type_column <- seurat_object@misc$cell_type_column
  sample_column <- seurat_object@misc$sample_column

  meta <- seurat_object@meta.data[is.na(seurat_object@meta.data[clonotype_column])==FALSE,]

  clone_count <- meta %>% group_by(meta[clonotype_column]) %>%
    summarise(clone_count=n())

  # Filter by clone size
  clones_to_include <- clone_count[clone_count$clone_count >= minimum_count, clonotype_column]
  meta <- meta[meta[,clonotype_column] %in% clones_to_include[[clonotype_column]],]

  ## Add Pheno Count
  pheno_clone_count <- meta %>% group_by(meta[clonotype_column], meta[sample_column], meta[cell_type_column]) %>%
    summarise(pheno_clone_count=n())
  pheno_clone_count[is.na(pheno_clone_count$CTaa),]$pheno_clone_count <- NA

  # Add sample_CTaa
  pheno_clone_count$sample_clonotype <- paste(pheno_clone_count[[sample_column]], pheno_clone_count[[clonotype_column]], sep = "_")

  # Make Template
  sample_clonotypes <- unique(pheno_clone_count$sample_clonotype)
  cell_types <- unique(pheno_clone_count[[cell_type_column]])

  template <- expand.grid(cell_types, sample_clonotypes)
  colnames(template) <- c(cell_type_column, "sample_clonotype")

  # Combine template and count data
  pheno_clone_count <- pheno_clone_count[c(cell_type_column, "sample_clonotype", "pheno_clone_count")]

  pheno_clone_count <- left_join(x = template, y = pheno_clone_count, by = c(cell_type_column, "sample_clonotype"))
  pheno_clone_count[is.na(pheno_clone_count)] <- 0

  # Reshape to wide
  pheno_clone_count <- reshape(pheno_clone_count, idvar = c("sample_clonotype"), timevar = cell_type_column, direction = "wide")

  # Format for output
  row.names(pheno_clone_count) <- pheno_clone_count$sample_clonotype
  pheno_clone_count <- subset(pheno_clone_count, select=-c(sample_clonotype))

  colnames(pheno_clone_count) <- gsub("pheno_clone_count.", "", colnames(pheno_clone_count))

  # Divide by total number of clones
  pheno_clone_count <- pheno_clone_count / rowSums(pheno_clone_count)

  if (plot==TRUE) {
    heatmap(t(as.matrix(pheno_clone_count)), scale='column')
  }

  return(pheno_clone_count)

}
