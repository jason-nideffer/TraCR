
#' Highlighting a Single Clonotype
#'
#' This function plots a UMAP where cells can be colored according to a variable.
#' A single clonotype will also be highlighted.
#'
#' @param seurat_object seurat object to operate on
#' @param clonotype Character specifying a clonotype within the 'clonotype_column'
#' @param clonotype_column Character describing name of column to use for clone calling
#' @param group.by Character describing a column to use for the coloring of cells
#' @param highlight_color Color of the highlighted clonotype
#' @param other_colors Vector of colors to use for coloring groups in 'group.by'
#' @param highlight_size Number representing the point size of the highlighted clonotype
#' @param other_size Number representing the point size of the non-highlighted cells
#' @param other_alpha Number (0-1) corresponding to the transparency of non-highlighted cells
#' @return a ggplot object that is plotted and can be stored/modified
#' @export
highlight_clonotype <- function(seurat_object, clonotype, clonotype_column, group.by, highlight_color='black', other_colors, highlight_size=1, other_size=0.3, other_alpha=0.1) {

  clone_data = subset(seurat_object@meta.data, seurat_object@meta.data[clonotype_column]==clonotype)

  clone_data = seurat_object@reductions$umap@cell.embeddings[rownames(clone_data),]
  clone_data = as.data.frame(clone_data)
  clone_data$clonotype <- clonotype

  clone_plot <- geom_point(data=clone_data, aes(x=UMAP_1, y=UMAP_2), size=highlight_size, alpha=1, color=highlight_color)

  umap_data <- seurat_object@reductions$umap@cell.embeddings
  umap_data <- as.data.frame(umap_data)
  umap_data[group.by] <- seurat_object@meta.data[group.by]

  umap_plot <- ggplot2::geom_point(data=umap_data, aes(x=UMAP_1, y=UMAP_2,
                                                       group=umap_data[[group.by]], color=umap_data[[group.by]]),
                                   size=other_size, alpha=other_alpha)

  plot <- ggplot() +
    umap_plot + scale_color_manual(values=other_colors, name=group.by) +
    clone_plot +
    theme(
      # Hide panel borders and remove grid lines
      panel.border = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # Change axis line
      axis.line = element_line(colour = "black"),
      text = element_text(size = 15),
      # Remove gray background on legend
      legend.key=element_rect(fill="white")
    ) +
    guides(colour = guide_legend(override.aes = list(alpha=1, size=2)
    ))

  return(plot)
}

