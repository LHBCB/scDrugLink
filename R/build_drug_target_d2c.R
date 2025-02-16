#' @title Build Drug2Cell matrix based on the disease atlas and drug targets.
#' @description The Drug2Cell matrix is constructed by calculating the average gene expression for each drug's target group in each cell and adjusting for background expression by subtracting baseline biases. This is an R implementation of the Drug2Cell Python pipeline designed by Kanemaru et al. (PMID: 37438528).
#' @param dat A Seurat object for the log-normalised, integrated scRNA-seq disease atlas.
#' @param gene_list The list of intersect genes between the atlas and the CMAP L1000 disease tissue specific references.
#' @param drug_target_df The data.frame of each drug's targets.
#' @param n_bins The number of bins determines how finely the background expression is segmented into groups of genes with similar expression patterns.
#' @param ctrl_size The number of genes that are randomly selected from the background (control) group to normalize the expression data for each drug target group.
#' @param seed The random seed for reproducibility, 42 by default.
#' @param out_path The path where the RDS file of the constructed Drug2Cell matrix is saved.
#' @return The data.frame of the Drug2Cell matrix, whose rows are drugs and columns are cell barcodes.
#' @references Kanemaru K. et al. (2023). *Spatially resolved multiomics of human cardiac niches*. Nature.
#' @export


build_drug_target_d2c <- function(dat, gene_list, drug_target_df, n_bins = 25, ctrl_size = 50, seed = 42, out_path="results") {
  
  DefaultAssay(dat) <- "RNA"
  dat_subset <- subset(x = dat, features = gene_list)
  exp_mat <- dat_subset@assays$RNA@data
  exp_mat <- as.matrix(exp_mat)
  
  targets <- list()
  for (drug in drug_target_df$drug_name_lower) {
    targets[[drug]] <- gene_list %in% strsplit(drug_target_df[drug_target_df$drug_name_lower==drug, ]$gene_names, ";")[[1]]
  }
  
  weights <- as.data.frame(targets, row.names = gene_list, check.names = FALSE)
  weights <- sweep(weights, 2, colSums(weights) + 1e-6, `/`)
  
  X <- t(exp_mat)
  scores <- X %*% as.matrix(weights)
  
  # seurat scoring mechanism 
  obs_avg <- colMeans(X)
  n_items <- round(length(obs_avg) / (n_bins - 1))
  obs_rank <- rank(obs_avg, ties.method = "min") - 1
  obs_cut  <- obs_rank %/% n_items
  
  set.seed(seed)
  control_groups <- list()
  for (cut_value in unique(obs_cut)) {
    # Identify the genes belonging to this bin
    mask <- (obs_cut == cut_value)
    r_genes <- which(mask)
    
    # Shuffle genes, then keep only ctrl_size of them
    r_genes <- sample(r_genes)
    if (length(r_genes) > ctrl_size) {
      mask[] <- FALSE
      mask[r_genes[1:ctrl_size]] <- TRUE
    }
    
    # Store the final logical vector (which genes survive in this bin)
    control_groups[[paste0("bin", as.character(cut_value))]] <- mask
  }
  
  control_gene_weights <- as.data.frame(control_groups, row.names = gene_list)
  control_gene_weights <- sweep(
    control_gene_weights,
    2,
    colSums(control_gene_weights) + 1e-6,
    FUN = "/"
  )
  
  control_profiles <- X %*% as.matrix(control_gene_weights)
  
  drug_bins <- list()
  for (drug in colnames(weights)) {
    # targets[[drug]] is assumed to be a logical vector the length of var_names
    bins_used <- unique(obs_cut[targets[[drug]]])
    bins_used <- paste0("bin", bins_used)
    
    # Mark which columns in 'control_gene_weights' correspond to these bins
    drug_bins[[drug]] <- colnames(control_gene_weights) %in% bins_used
  }
  
  drug_weights <- as.data.frame(drug_bins, row.names = colnames(control_gene_weights), check.names = F)
  
  drug_weights <- sweep(
    drug_weights,
    2,
    colSums(drug_weights) + 1e-6,
    FUN = "/"
  )
  
  seurat <- control_profiles %*% as.matrix(drug_weights)
  scores <- scores - seurat
  
  colnames(scores) <- drug_label_df$drug_name[match(colnames(scores), drug_label_df$drug_name_lower)]
  d2c_mat <- t(scores)
  saveRDS(d2c_mat, paste0(out_path, "/drug_target_d2c.rds"))
  return(d2c_mat)
}