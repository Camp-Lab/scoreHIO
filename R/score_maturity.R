#' Estimate the Stem Cell Maturity of Human Intestinal Organoids
#' 
#' This function provides two approaches to estimate the maturity of stem cells in pluripotent stem cell
#' derived human intestinal organoids. One is to calculate transcriptome similarity (Spearman's correlation
#' coefficient) between query cells and reference intestinal stem cells (ISCs). The other is to decompose 
#' the cellular transcriptome into spheroid and maturing fetal (19 post conceptional week (PCW)) ISC 
#' components using quadratic programming and quantifies the maturity of query cells as the proportion 
#' of 19 PCW fetal ISC components.     
#' 
#' @param que_obj Seurat object in query
#' @param group_by A character string to specify the name of the meta.data field of a combined 
#' Seurat object from more than one groups. Only applicable when the \code{'que_obj'} is a 
#' combination of different groups and you would like to separate them. Default is \code{'NULL'} 
#' @param score_method Specify the method to estimate stem cell maturity. Currently support 
#' \code{'similarity_based'} and \code{'deconvolution_based'}. Default is \code{'similarity_based'}.
#' @param reference_age A character string to specify which age(s) as reference for transcriptome
#' similarity calculation. Currently support one or several values from \code{'thio-4-week'}, 
#' \code{'thio-8-week'}, \code{'fetal-7-week'}, \code{'fetal-8-week'}, \code{'fetal-10-week'}, 
#' \code{'fetal-11-week'}, \code{'fetal-12-week'}, \code{'fetal-14-week'}, \code{'fetal-17-week'}, 
#' \code{'fetal-18-week'}, \code{'fetal-19-week'}, \code{'adult'}. If choose multiple reference
#' ages, only median similarity to each reference for each query group is shown in the output plot.
#' @param add_landmark_data Logical. Whether to add reference ISC cellular data to the query set
#' @param selected_landmark_age A character string to specify which reference ISC cellular data 
#' to be added in the query set
#' @param group_cols A vector to specify the color(s) of query group(s) for plotting. Names of
#' the vector specify the group, values specify colors in Hex format.  
#' @param do_plot Logical. Whether to do plotting for result visualization
#' @param output_dir A character string to specify output directory
#' @param plot_name A character string to specify the file name of the plot
#' @param verbose Logical. Whether to print log message
#'  
#' @return A plot to show the stem cell maturity and a list of numeric values of maturity 
#' 
#' @rdname score_maturity
#' @export

score_maturity <- function(
  que_obj,
  group_by=NULL,
  score_method=c("similarity_based", "deconvolution_based"), 
  reference_age="adult",
  add_landmark_data=TRUE,
  selected_landmark_age=c("thio-4-week", "fetal-7-week", "fetal-8-week", "adult"),
  group_cols=NULL,
  do_plot=TRUE,
  output_dir="./",
  plot_name="Plot_ISC_maturity.pdf",
  verbose=TRUE
){
  score_method <- match.arg(score_method)
  log_message("Choose ", score_method, " method to estimate stem cell maturity", verbose=verbose)
  
  if(is.null(group_by)){
    log_message("No grouping within user provided Seurat object", verbose=verbose)
    que_list <- list(
      "query"=que_obj
    )
  }else{
    log_message("Group by ", group_by, verbose=verbose)
    que_list <- Seurat::SplitObject(que_obj, split.by = group_by)
  }
  
  que_expr_list <- lapply(names(que_list), function(k){
    obj <- que_list[[k]]
    que_expr <- obj@assays$RNA@data
  })
  names(que_expr_list) <- names(que_list)
  
  if(add_landmark_data){
    log_message("Selected landmark data as background: ",paste(selected_landmark_age, collapse = ","), verbose=verbose)
    
    age_vec <- sapply(colnames(duo_sc_expr), function(x){
      strsplit(x,"_")[[1]][3]
    })
    
    for(x in selected_landmark_age){
     idx <- which(age_vec==x)
     que_expr_list[[x]] <- duo_sc_expr[, idx]
     }
  }
  
  if(score_method=="deconvolution_based"){
    feature_gene <- intersect(feature_gene_set$deconvolution_feature_genes, rownames(que_obj))
    log_message("Number of feature genes = ", length(feature_gene), verbose=verbose)
    X <- decov_ref[feature_gene,]
    
    log_message("Run quadratic programming for transcriptome deconvolution", verbose=verbose)
    isc_maturity_list <- lapply(names(que_expr_list), function(k){
      que_expr <- as.matrix(que_expr_list[[k]][feature_gene,])
      identity_que <- matrix(NA, nrow=ncol(que_expr), ncol=4)
      rownames(identity_que) <- colnames(que_expr)
      colnames(identity_que) <- c("frxn_cl12", "frxn_fetLatSC", "Lagrangian", "Error")
      for (j in seq(ncol(que_expr))){
        Y <- as.matrix(que_expr[,j])
        Rinv <- solve(chol(t(X) %*% X));
        C <- cbind(rep(1,2), diag(2))
        b <- c(1,rep(0,2))
        d <- t(Y) %*% X
        QP <- quadprog::solve.QP(Dmat = Rinv, factorized = TRUE, dvec = d, Amat = C, bvec = b, meq = 1)
        Error <-sum(abs(Y-X%*%QP$solution))
        identity_que[j,] <- c(QP$solution[1],QP$solution[2],QP$Lagrangian[1],Error)
      }
      isc_score <- identity_que[,2]
      isc_score[which(isc_score>1)] <- 1
      isc_score[which(isc_score<0)] <- 0
      return(isc_score)
    })
    names(isc_maturity_list) <- names(que_expr_list)
    
  }else if(score_method=="similarity_based"){
    
    log_message("Selected reference age(s) for transcriptome similarity calculation: ",paste(reference_age, collapse = ","), verbose=verbose)
    
    feature_gene <- intersect(feature_gene_set$top_phase_markers, rownames(que_obj))
    log_message("Feature gene number =", length(feature_gene), verbose=verbose)
    X <- stem_cell_expr_by_age[feature_gene, reference_age]
    
    log_message("Calculate transcriptome similarity between query cells and reference", verbose=verbose)
    isc_maturity_list <- lapply(names(que_expr_list), function(k){
      que_expr <- as.matrix(que_expr_list[[k]][feature_gene,])
      scc_res <- cor(que_expr, X, method="spearman")
      if(length(reference_age)==1){
        scc_res <- scc_res[,1]
      }
      return(scc_res)
    })
    names(isc_maturity_list) <- names(que_expr_list)
  }
  if(do_plot){
    log_message("Plot distribution of estimated stem cell maturity", verbose=verbose)
    if(score_method=="similarity_based" & length(reference_age)>1){
      log_message("Multiple reference time points are selected. Only median Spearman's correlation coefficient to each time point will be shown", verbose=verbose)
      if(is.null(group_cols)){
        group_cols <- setNames(colorRampPalette(prettyrainbow)(length(names(isc_maturity_list))), names(isc_maturity_list))
      }
      median_mat <- t(sapply(seq(length(isc_maturity_list)), function(i){
        apply(isc_maturity_list[[i]], 2, median)
      }))
      pdf(paste0(output_dir,plot_name), height=5, width=10)
      par(mar=c(6,5,3,10), xpd=TRUE)
      plot(rep(seq(length(reference_age)), each=length(isc_maturity_list)), as.vector(median_mat), 
           pch=16, col=rep(group_cols, length(reference_age)), bty="n", xaxt="n", xlab="", ylab="Median transcriptome similarity")
      for(i in seq(nrow(median_mat))){
        lines(seq(length(reference_age)), median_mat[i,], col=group_cols[i])
      }
      text(seq(ncol(median_mat)), rep(min(median_mat)-0.05, ncol(median_mat)), labels = colnames(median_mat), srt=30, adj = 1)
      title(xlab="Reference age", line=4.5, cex.lab=1)
      legend("right", inset=c(-0.25,0), legend=names(group_cols), col=group_cols, pch=16, title="Query group")
      dev.off()
    }else{
      if(score_method=="similarity_based" & length(reference_age)==1){
        ylab_title <- paste("Transcriptome similarity to", reference_age)
      }else if(score_method=="deconvolution_based"){
        ylab_title <- "Deconvolution-based stem cell maturity" 
      }
      pdf(paste0(output_dir,plot_name), height=5, width=7)
      par(mar=c(7,5,4,4))
      beanplot::beanplot(isc_maturity_list, what=c(1,1,1,0), ylab=ylab_title, las=2)
      dev.off()
    }
  }
  
  log_message("Return estimated stem cell maturity", verbose=verbose)
  return(isc_maturity_list)
}

