% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{load_fetal_organ_ref_data}
\alias{load_fetal_organ_ref_data}
\title{Load Fetal Human Organ Identity Inference Reference Data}
\usage{
load_fetal_organ_ref_data(
  organ_ref_dir = "./",
  do_projection = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{organ_ref_dir}{A character string to specify the path to the folder with
the fetal reference (Yu et al, Cell, 2021) downloaded from Mendeley data
https://doi.org/10.17632/x53tts3zfr.
Default is \code{'../Ref_data_for_projection_to_fetal_atlas/'}}

\item{do_projection}{Logical. Whether to project the query cells to the UMAP
embedding of the reference atlas for visualization. If \code{'TRUE'}, then
the CSS-based UMAP model of the fetal reference will be loaded. Otherwise not.}

\item{verbose}{Logical. Whether to print log message}
}
\value{
A list \code{'fetal_reference'} with five components:
\itemize{
\item \code{'fetal_hvg'}: highly variable genes defined in the fetal reference
scRNA-seq data
\item \code{'fetal_embedding'}: UMAP embedding of fetal reference cells
\item \code{'ref_idx'}: Organ identity of fetal reference cells
\item \code{'css_model'}: Cluster Similarity Spectrum (CSS) model of the fetal
reference. For more detailed explanation of CSS, please refer to He et al,
Genome Biology, 2020 (https://doi.org/10.1186/s13059-020-02147-4)
\item \code{'fetal_umap_res'}: CSS-based UMAP model of the fetal reference
}
}
\description{
This function loads fetal human multi-endodermal-organ single cell meta data as
reference for organ identity inference
}
