#' Print Diagnostic Message
#'
#' @param ... Text to print
#' @param verbose Logical. Whether to display messages
log_message <- function(..., verbose=T){
  if (verbose){
    message(paste0(...))
  }
}



#' Load Fetal Human Organ Identity Inference Reference Data
#' 
#' This function loads fetal human multi-endodermal-organ single cell meta data as 
#' reference for organ identity inference
#' 
#' @import uwot
#' 
#' @param organ_ref_dir A character string to specify the path to the folder with 
#' the fetal reference (Yu et al, Cell, 2021) downloaded from Mendeley data 
#' https://doi.org/10.17632/x53tts3zfr. 
#' Default is \code{'../Ref_data_for_projection_to_fetal_atlas/'}
#' @param do_projection Logical. Whether to project the query cells to the UMAP 
#' embedding of the reference atlas for visualization. If \code{'TRUE'}, then
#' the CSS-based UMAP model of the fetal reference will be loaded. Otherwise not. 
#' @param verbose Logical. Whether to print log message
#'  
#' @return A list \code{'fetal_reference'} with five components:
#' * \code{'fetal_hvg'}: highly variable genes defined in the fetal reference 
#' scRNA-seq data
#' * \code{'fetal_embedding'}: UMAP embedding of fetal reference cells
#' * \code{'ref_idx'}: Organ identity of fetal reference cells
#' * \code{'css_model'}: Cluster Similarity Spectrum (CSS) model of the fetal
#' reference. For more detailed explanation of CSS, please refer to He et al, 
#' Genome Biology, 2020 (https://doi.org/10.1186/s13059-020-02147-4)
#' * \code{'fetal_umap_res'}: CSS-based UMAP model of the fetal reference
#' 
#' @rdname load_fetal_organ_ref_data
#' @export

load_fetal_organ_ref_data <- function(
  organ_ref_dir="./", 
  do_projection=TRUE,
  verbose=TRUE
){
  
  log_message("Reading fetal human multi-endoderm-organ reference data", verbose=verbose)
  
  fetal_hvg <- readLines(paste0(organ_ref_dir, "Res_fetal_atlas_highly_variable_genes.csv"))
  fetal_meta <- read.csv(paste0(organ_ref_dir, "Table_fetal_atlas_meta_data.csv"))
  fetal_embeddings <- fetal_meta[,c("UMAP_X", "UMAP_Y")]
  ref_idx <- fetal_meta$Corrected_organ_group
  css_model <- readRDS(paste0(organ_ref_dir, "Res_fetal_CSS_model.rds"))
  if(do_projection){
    fetal_umap_res <- uwot::load_uwot(paste0(organ_ref_dir, "Res_fetal_all_cell_type_CSS_based_UMAP_res.uwot"))
  }else{
    fetal_umap_res <- NULL
  }
  fetal_reference <- list("fetal_hvg"=fetal_hvg,
                          "fetal_embeddings"=fetal_embeddings,
                          "ref_idx"=ref_idx,
                          "css_model"=css_model,
                          "fetal_umap_res"=fetal_umap_res
  )
  return(fetal_reference)
}



#' Define Default Color Palette 
prettyrainbow <- c("#8966A9","#3FA4D9","#3EB8B4","#B2C224","#F8CB31","#F6A22B","#EC6325","#DC3838")
