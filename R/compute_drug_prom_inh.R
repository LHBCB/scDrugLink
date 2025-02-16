#' @title Estimate the promotion/inhibition effect of each drug on different cell types in the disease atlas.
#' @description Each drug's promotion/inhibition effect on each cell type is calculated by comparing scores derived from control and diseased cells using within-cell-type Wilcoxon rank-sum test, and adjusted for multiple comparisons using the Benjamini-Hochberg procedure.
#' @param dat A Seurat object for the log-normalised, integrated scRNA-seq disease atlas.
#' @param d2c_mat The data.frame of the Drug2Cell matrix constructed by the build_drug_target_d2c function.
#' @param disease The name of the investigated disease.
#' @param out_type The type of returned promotion/inhibition effect scores, by default "weight", as computed in Eq.(7) of the scDrugLink paper.
#' @param out_path The path where the csv files of adjusted p-values and weights are saved.
#' @return The data.frame of promotion/inhibition effect weights (scores) for each pair of cell type (row) and drug (column).
#' @export
#' @import effsize


compute_drug_prom_inh <- function(dat, d2c_mat, disease, out_type="weight", out_path="results"){
  
  Idents(dat) <- "cell_type"
  cell_type_list <- unique(Idents(dat))

  total_num_cells <- ncol(dat)
  
  weight_list <- list()
  prop_weight_list <- list()
  p_adj_list <- list()
  for (cell_type in cell_type_list) {
    print(paste0("Processing cell type: ", cell_type))
    subset_dat <- subset(dat, idents = cell_type)
    subset_cell_labels <- ifelse(subset_dat@meta.data$disease==disease, 1, 0)
    subset_d2c <- d2c_mat[, colnames(subset_dat)]
    
    print(ncol(subset_dat))
    
    diseased_cells <- colnames(subset_d2c)[which(subset_cell_labels==1)]
    healthy_cells <- colnames(subset_d2c)[which(subset_cell_labels==0)]
    
    diseased_d2c <- subset_d2c[, diseased_cells]
    healthy_d2c  <- subset_d2c[, healthy_cells]
    
    p_vals <- c()
    cliff_delta_vals <- c()
    for (drug in rownames(diseased_d2c)) {
      diseased <- diseased_d2c[drug, ]
      healthy <- healthy_d2c[drug, ]
      
      if (any(var(diseased) == 0, var(healthy) == 0)) {
        warning(paste0(drug, ": Zero variance in at least one of the groups; skipping wilcox.test and cliff.delta computation."))
        p_value <- 1
        cliff_delta_value <- 0
      } else {
        wilcox_result <- wilcox.test(diseased, healthy,
                                     alternative = "two.sided",
                                     paired = FALSE,
                                     conf.int = TRUE,
                                     conf.level = 0.95)
        p_value <- wilcox_result$p.value
        cliff_delta <- cliff.delta(diseased, healthy)
        cliff_delta_value <- cliff_delta$estimate
      }
      
      p_vals <- c(p_vals, p_value + 1e-6) # add 1e-6 to prevent log10(0) = inf
      cliff_delta_vals <- c(cliff_delta_vals, cliff_delta_value)
    }
    p_adj <- p.adjust(p_vals, method = "BH")
    
    # weight = effect size x -log10(p_adj)
    # weights can be neg (drug inhibition) or pos (drug promotion), taking into account resistance/sensitivity of cells
    prop_weights <- ncol(subset_dat) / total_num_cells * cliff_delta_vals * (-log10(p_adj)) 
    prop_weight_list[[cell_type]] <- prop_weights
    
    weights <- cliff_delta_vals * (-log10(p_adj)) 
    weight_list[[cell_type]] <- weights
    
    p_adj_list[[cell_type]] <- p_adj
  }
  
  # p_adj
  temp_list <- lapply(p_adj_list, function(x) as.data.frame(t(unlist(x))))
  p_adj_df <- do.call(rbind, temp_list)
  colnames(p_adj_df) <- tolower(rownames(d2c_mat))
  write.csv(p_adj_df, paste0(out_path, "/", disease, "_p_adj_target.csv"), row.names=T)

  # weight
  temp_list <- lapply(weight_list, function(x) as.data.frame(t(unlist(x))))
  weight_df <- do.call(rbind, temp_list)
  colnames(weight_df) <- tolower(rownames(d2c_mat))
  write.csv(weight_df, paste0(out_path, "/", disease, "_prom_inh_weight_target.csv"), row.names=T)
  
  weight_names <- colnames(weight_df)
  weight_sum <- colSums(weight_df)
  weight_sum_df <- data.frame(weight_sum, row.names = weight_names, check.names = F)
  
  # prop_weight
  temp_list <- lapply(prop_weight_list, function(x) as.data.frame(t(unlist(x))))
  prop_weight_df <- do.call(rbind, temp_list)
  colnames(prop_weight_df) <- rownames(d2c_mat)
  
  prop_weight_names <- colnames(prop_weight_df)
  prop_weight_sum <- colSums(prop_weight_df)
  prop_weight_sum_df <- data.frame(prop_weight_sum, row.names = prop_weight_names, check.names = F)
  
  if (out_type == "weight") {
    return(weight_df)
  } else if (out_type == "prop_weight") {
    return(prop_weight_df)
  } else if (out_type == "weight_sum") {
    return(weight_sum_df)
  } else if (out_type == "prop_weight_sum") {
    return(prop_weight_sum_df)
  } else {
    stop("Invalid output type.")
  }
}