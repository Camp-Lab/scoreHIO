#' Estimate the organ fidelity of Human Intestinal Organoids
#' 
#' This function estimates the organ fidelity of cells in pluripotent stem cell
#' derived human intestinal organoids by comparing to a developing human multi-
#' endodermal-organ cell atlas
#' 
#' @import Seurat
#' @import RANN
#' @import uwot
#' 
#' @param que_obj Seurat object in query
#' @param organ_ref_dir A character string to specify the path to the folder with the fetal 
#' reference (Yu et al, Cell, 2021) downloaded from Mendeley data https://doi.org/10.17632/x53tts3zfr.
#' Default is \code{'../Ref_data_for_projection_to_fetal_atlas/'} 
#' @param k Numeric. Specify the number of k-nearest neighboring cells in the reference atlas
#' for each query organoid cell for fidelity estimation.
#' @param group_by A character string to specify the name of the meta.data field of a combined 
#' Seurat object from more than one groups. Only applicable when the \code{'que_obj'} is a 
#' combination of different groups and you would like to separate them. Default is \code{'NULL'} 
#' @param group_cols A vector to specify the color(s) of query group(s) for plotting. Names of
#' the vector specify the group, values specify colors in Hex format.
#' @param organ_cols A vector to specify the colors of organ identity in reference atlas. Names 
#' of the vector specify the group, values specify colors in Hex format.  
#' @param return_object Logical. Whether to return a Seurat object with inferred organ identity
#' for each query cell
#' @param output_dir A character string to specify output directory
#' @param do_projection Logical. Whether to project the query cells to the UMAP embedding of 
#' the reference atlas for visualization
#' @param probability_file_name A character string to specify the file name of the plot showing
#' the projected organ probability
#' @param projection_file_name A character string to specify the file name of the plot projecting
#' query cells to the UMAP embedding of reference atlas 
#' @param verbose Logical. Whether to print log message
#'  
#' @return A bar plot to show the organ projection probability, a UMAP plot to show projected UMAP
#' embedding of the query cell on reference atlas and an updated Seurat object.\cr
#' The updated Seurat has three new components:
#' * Meta.data field \code{'Mapped_fetal_organ'} indicating the inferred organ identity of each query cell;
#' * Dimension reduction object \code{'rss'} representing similarity of each query cell to the fetal
#' similarity spectrum;
#' * Dimension reduction object \code{'rss_umap'} representing projected UMAP embedding of each query cell
#' on the fetal reference atlas
#' 
#'  
#' 
#' @rdname score_fidelity
#' @export

score_fidelity <- function(
  que_obj,
  organ_ref_dir = "../Ref_data_for_projection_to_fetal_atlas/",
  k=20,
  group_by=NULL,
  group_cols=NULL,
  organ_cols=setNames(c("#9D0142", "#5D4EA2", "#FCB163", "#1B7635", "#4ea7b0", "#696969"), 
                      c("Intestine", "Stomach", "Lung", "Esophagus", "Liver", "Pancreas")),
  return_object = TRUE,
  output_dir="./",
  probability_file_name="Plot_barplot_projected_developing_organ_probability.pdf",
  do_projection=TRUE,
  projection_file_name="Plot_UMAP_RSS_projection_to_developing_reference.png",
  verbose=TRUE
){
  
  if(!exists("fetal_reference")){
    fetal_reference <- load_fetal_organ_ref_data(
      organ_ref_dir=organ_ref_dir, 
      do_projection=do_projection, 
      verbose=verbose
    )
  }
  fetal_hvg <- fetal_reference$fetal_hvg
  css_model <- fetal_reference$css_model
  shared_genes <- intersect(fetal_hvg, rownames(que_obj))
  log_message("Number of feature genes = ", length(shared_genes), verbose=verbose)
  que_data <- as.matrix(que_obj@assays$RNA@data[shared_genes,])
  log_message("Calculating RSS of query cells", verbose=verbose)
  rss_list <- lapply(seq(length(css_model$model$profiles)), function(i){
    ref <- css_model$model$profiles[[i]][shared_genes,]
    cor_mat <- cor(que_data, ref, method = "spearman")
    cor_z <- t(scale(t(cor_mat)))
    return(cor_z)
  })
  rss_mat <- do.call('cbind', rss_list)
  que_obj[["rss"]] <- Seurat::CreateDimReducObject(embeddings = rss_mat, key="RSS_", assay=DefaultAssay(que_obj))
  
  css <- css_model$sim2profiles 
  log_message("Calculating Euclidean distance between query cells and reference cells on the similarity spectrum space", verbose=verbose)
  log_message("Get k-nearest neighbors in reference for each query cell (k = ", k, ")", verbose=verbose)
  knn <- RANN::nn2(css, rss_mat, k = k)$nn.idx
  ref_idx <- fetal_reference$ref_idx
  nn_idx <- matrix(ref_idx[as.vector(knn)], nrow=nrow(knn))
  trans_id <- apply(nn_idx, 1, function(vec){
    freq <- table(vec)
    names(which.max(freq))
  })
  que_obj@meta.data$Mapped_fetal_organ <- trans_id
  log_message("Get inferred organ identity of each query cell", verbose=verbose)
  
  log_message("Calculate organ projection probability for the query group(s)", verbose=verbose)
  organ_vec <- que_obj@meta.data$Mapped_fetal_organ
  organ_order <- c("Intestine", "Stomach", "Lung", "Esophagus", "Pancreas")
  if(is.null(group_by)){
    log_message("No grouping within user provided Seurat object", verbose=verbose)
    num1 <- sapply(organ_order, function(organ){ sum(organ_vec==organ) })
    p1 <- num1/sum(num1)
  }else{
    log_message("Group by ", group_by, verbose=verbose)
    group_vec <- que_obj@meta.data[[group_by]]
    num1 <- sapply(sort(unique(group_vec)), function(group){
      sapply(organ_order, function(organ){
        sum(group_vec==group & organ_vec==organ)
      })
    })
    p1 <- t(t(num1)/colSums(num1))
    colnames(p1) <- sort(unique(group_vec))
  }
  
  log_message("Plot organ projection probability for the query group(s)", verbose=verbose)
  pdf(probability_file_name, height=5, width=8)
  par(xpd=TRUE, mar=c(5,5,3,9))
  barplot(p1, col=organ_cols[organ_order], border = NA, beside = T, ylab="Organ projection probability")
  legend("right", title="Organ", inset=c(-0.25,0), legend=names(organ_cols), fill=organ_cols)
  dev.off()
  
  if(do_projection){
    
    fetal_umap_res <- fetal_reference$fetal_umap_res
    fetal_embeddings <- fetal_reference$fetal_embeddings
    log_message("Project query cells to the CSS-based UMAP embedding of fetal atlas data", verbose=verbose)
    que_obj_umap <- uwot::umap_transform(rss_mat, fetal_umap_res)
    rownames(que_obj_umap) <- rownames(rss_mat)
    que_obj[["rss_umap"]] <- Seurat::CreateDimReducObject(embeddings = que_obj_umap, key="RSSUMAP_", assay=DefaultAssay(que_obj))
    
    organ_col_vec <- organ_cols[organ_vec]
    if(is.null(group_by)){
      log_message("Generate projected UMAP plot for query cells", verbose=verbose)
      png(projection_file_name, height=2000, width=2000)
      plot(fetal_embeddings, pch=16, col="#e0e0e0", bty="n", xaxt="n", yaxt="n", xlab="", ylab="", cex=2)
      points(Embeddings(que_obj, reduction = "rss_umap"), pch=21, col="#303030", bg=organ_col_vec, cex=2)
      legend("topleft", pch=16, col=organ_cols, legend = names(organ_cols), bty="n", cex=5)
      dev.off()
      
    }else{
      group_vec <- que_obj@meta.data[[group_by]]
      if(is.null(group_cols)){
        group_cols <- setNames(colorRampPalette(prettyrainbow)(length(unique(group_vec))), sort(unique(group_vec)))
        group_col_vec <- group_cols[group_vec]
      }
      
      log_message("Generate projected UMAP plot for query cells", verbose=verbose)
      png(projection_file_name, height=2000, width=2000*2)
      par(mfrow=c(1,2))
      plot(fetal_embeddings, pch=16, col="#e0e0e0", bty="n", xaxt="n", yaxt="n", xlab="", ylab="", cex=2)
      points(Embeddings(que_obj, reduction = "rss_umap"), pch=21, col="#303030", bg=organ_col_vec, cex=2)
      legend("topleft", pch=21, col="#303030", pt.bg=organ_cols, legend = names(organ_cols), cex=5, title="Organ")
      plot(fetal_embeddings, pch=16, col="#e0e0e0", bty="n", xaxt="n", yaxt="n", xlab="", ylab="", cex=2)
      points(Embeddings(que_obj, reduction = "rss_umap"), pch=21, col="#303030", bg=group_col_vec, cex=2)
      legend("topleft", pch=21, col="#303030", pt.bg=group_cols, legend = names(group_cols), cex=5, title=group_by)
      dev.off()
      
    }
    
  }
  
  if(return_object){
    log_message("Return updated Seurat object", verbose=verbose)
    return(que_obj)
  }
}
