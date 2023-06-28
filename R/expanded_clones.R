
#' Identifying Expanded Clones
#'
#' This function takes a seurat object with meta data that already includes
#' TCR annotations and calculates the expansion of clones in sample.1 compared to sample.2
#'
#' @param seurat_object seurat object to operate on
#' @param clonotype_column Character describing name of column to use for clone calling
#' @param sample_column Character describing name of column to use for sample discrimination
#' @param sample.1 Character describing a sample in 'sample_column' for comparison
#' @param sample.2 Character describing a sample in 'sample_column' for comparison
#' @param minimum_count integer value specifying the minimum number of clones for a clonotype to be included in analysis
#' @param decreasing Boolean value specifying whether or not to sort by decreasing (or increasing) expansion
#' @param identity A value in 'clonotype_identity_mode' to specify clonotypes of interest (must first run the 'calculate_freq' function)
#' @return a dataframe of clones sorted by their expansion in sample.1 compared to sample.2
#' @export
expanded_clones <- function(seurat_object, clonotype_column, sample_column, sample.1, sample.2,
                            minimum_count=0, decreasing=TRUE, identity=NULL) {

  # Extract meta data
  clone_freqs <- seurat_object@meta.data[is.na(seurat_object@meta.data[clonotype_column])==FALSE,
                                         c(clonotype_column, sample_column, "clonotype_count_per_sample","clonotype_freq_per_sample","clonotype_identity_mode")]

  # Get distinct rows (because of duplicates)
  clone_freqs <- distinct(clone_freqs)

  # Remove row names
  row.names(clone_freqs) <- NULL

  # Select clone identities of interest
  if (is.null(identity)==FALSE) {
    clone_freqs <- subset(clone_freqs, clonotype_identity_mode==identity)
  }

  # Select samples of interest
  clone_freqs <- clone_freqs[clone_freqs[,sample_column] %in% c(sample.1, sample.2),]

  # Reshape to wide
  clone_freqs <- reshape(clone_freqs, idvar = c(clonotype_column), timevar = sample_column, direction = "wide")

  # Fill NA with 0
  clone_freqs[is.na(clone_freqs)] <- 0

  # Extract clone counts
  sample.1.clonecounts <- clone_freqs[, paste("clonotype_count_per_sample", sample.1, sep=".")]
  sample.2.clonecounts <- clone_freqs[, paste("clonotype_count_per_sample", sample.2, sep=".")]

  # Extract clone freqs
  sample.1.clonefreqs <- clone_freqs[, paste("clonotype_freq_per_sample", sample.1, sep=".")]
  sample.2.clonefreqs <- clone_freqs[, paste("clonotype_freq_per_sample", sample.2, sep=".")]

  # Initialize output
  output <- clone_freqs[c(clonotype_column)]

  # Calculate change in freq
  output$freq_change <- sample.1.clonefreqs - sample.2.clonefreqs

  # Calculate fold change in freq
  output$freq_FC <- sample.1.clonefreqs / sample.2.clonefreqs

  # Calculate change in count
  output$count_change <- sample.1.clonecounts - sample.2.clonecounts

  # Calculate fold change in count
  output$count_FC <- sample.1.clonecounts / sample.2.clonecounts

  # Include sample counts
  output$sample.1_count <- sample.1.clonecounts
  output$sample.2_count <- sample.2.clonecounts

  # Include total count
  output$total_count <- sample.1.clonecounts + sample.2.clonecounts

  # Sort dataframe
  output <- output[with(output, order(freq_FC, count_change, decreasing=decreasing)),]

  # Exclude clones that are below minimum_count
  output <- output[output$total_count > minimum_count,]

  # Remove row names
  row.names(output) <- NULL

  return(output)
}
