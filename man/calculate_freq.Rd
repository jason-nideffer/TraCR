% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_freq.R
\name{calculate_freq}
\alias{calculate_freq}
\title{Calculate Clonotype Frequencies}
\usage{
calculate_freq(
  seurat_object,
  clonotype_column,
  sample_column,
  cell_type_column = NULL,
  subject_column = NULL,
  exclude_na = TRUE
)
}
\arguments{
\item{seurat_object}{seurat object to operate on}

\item{clonotype_column}{Character describing name of column to use for clone calling}

\item{sample_column}{Character describing name of column to use for sample discrimination}

\item{cell_type_column}{Character describing name of column to use for cell_type identification (could be and annotation or cluster)}

\item{exclude_na}{Boolean; if TRUE (default), cells without clonotype designation will be excluded from frequency calculations}
}
\value{
seurat_object with meta data updated to include frequency data
}
\description{
This function takes a seurat object with meta data that already includes
TCR annotations and calculates the frequencies of different clonotypes
on a per-sample basis.
}
