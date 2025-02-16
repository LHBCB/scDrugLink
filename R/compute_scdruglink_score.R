#' @title Link the drug promotion/inhibition and drug sensitivity/resistance estimations to generate an overall therapeutic score
#' @description The final drug score is calculated by summing the promotion/inhibition-weighted sensitivity/resistance scores for each cell type, as described by Eq. (7-9) in the scDrugLink paper. This function extends the DrugScore function of Asgard to facilitate the linking mechanism proposed in the scDrugLink paper.
#' @param dat A Seurat object for a log-normalised, integrated scRNA-seq disease atlas
#' @param deg_list The list of intra-cell-type DEGs identified by the get_intra_cell_type_degs function.
#' @param drug_sens_res_pval the p-values and the adjusted p-values (FDRs) for statistical tests that estimate sensitivity/resistance effect for each pair of cell type and drug.
#' @param tissue The name of the disease-specific tissue in CMAP L1000.
#' @param gse70138_gctx_path The path of the CMAP L1000 GSE70138 data.
#' @param gse92742_gctx_path The path of the CMAP L1000 GSE92742 data.
#' @param drug_prom_inh_weight The promotion/inhibition effect weights for each pair of cell type and drug.
#' @param clusters The specific cell types/clusters on which computation of the final drug scores is based, NULL by default for all cell types/clusters.
#' @param case The name of the disease in the atlas. 
#' @param fda_drugs_only Whether only FDA approved drugs are considred for repurposing, TRUE by default. This parameter is more relevant in Asgard, as scDrugLink determines the drugs and fix them from the beginning of the pipeline.
#' @param out_path The path where the csv files of the final drug scores at the cell-type level and the whole-atlas levl are saved.
#' @return The data.frame of the final drug scores at the whole-atlas level.
#' @export
#' @import cmapR

compute_scdruglink_score <- function(dat, deg_list, drug_sens_res_pval, tissue,
					  gse70138_gctx_path, gse92742_gctx_path, drug_prom_inh_weight,
					  clusters = NULL, case = NULL, fda_drugs_only = TRUE, out_path="results") {

  dat_subset <- subset(dat, subset = cell_type %in% names(deg_list))
  cell_metadata <- dat_subset@meta.data

  cell_metadata$cluster <- dat_subset@meta.data$cell_type

  cluster_drugs <- drug_sens_res_pval

	# Subset input data to the set of clusters we are interested in 
  if (length(clusters) > 0) {
    clusters = intersect(clusters, unique(cell_metada$cluster))
      cell_metadata = subset(cell_metadata, cluster %in% clusters)
      cluster_drugs = cluster_drugs[clusters]
      deg_list = deg_list[clusters]
  }

  # Calculate cluster proportions in diseased tissue
  case = c(case)
  if (length(case) > 0) {
      cell_metadata <- subset(cell_metadata, disease %in% case)
  }
  clustering <- cell_metadata$cluster
  cluster_sizes <- table(clustering)
  cluster_sizes <- cluster_sizes[which(cluster_sizes > 3)]
  cluster_prop <- round(100*cluster_sizes/nrow(cell_metadata), 2) 

  library(Asgard)
  # Combine cluster drugs into a single data frame
  drug_list <- data.frame()
  for (i in names(cluster_drugs)) {
    ith_cluster_drugs <- cluster_drugs[[i]]
    drug_names <- ith_cluster_drugs$Drug.name
    ith_cluster_drugs <- ith_cluster_drugs[!duplicated(drug_names), ]

    # Subset to FDA drugs
    if (fda_drugs_only) {
      drug_names <- intersect(drug_names, FDA.drug) # after intersect, drug_names become unique
    }

    if (length(drug_names)>0) {
      ith_cluster_drugs <- subset(ith_cluster_drugs, Drug.name %in% drug_names)
      fdrs <- ith_cluster_drugs$FDR
      p_values <- ith_cluster_drugs$P.value
      
      temp <- data.frame(
        drug = drug_names, 
        cluster = i,
        cluster_prop = cluster_prop[i],
        p_value = p_values,
        fdr = fdrs,
        row.names = NULL
      )
        drug_list <- rbind(drug_list, temp)
    }
  }
  drug_list <- unique(drug_list)
  drug_list$weighted_prop <- drug_list$cluster_prop*(-log10(drug_list$fdr))
  drug_list[is.na(drug_list)] <- 0

  drug_coverage <- tapply(drug_list$weighted_prop, drug_list$drug, sum)
  drugs <- rownames(drug_coverage)
  
  CombineP <- function (p){
    keep <- (p > 0) & (p <= 1)
    invalid <- sum(1L * keep) < 2
    if (invalid) {
      warning("Must have at least two valid p values")
      res <- list(chisq = NA_real_, df = NA_integer_, p = NA_real_, 
                  validp = p[keep])
    }
    else {
      lnp <- log(p[keep])
      chisq <- (-2) * sum(lnp)
      df <- 2 * length(lnp)
      if (length(lnp) != length(p)) {
        warning("Some studies omitted")
      }
      res <- pchisq(chisq,df, lower.tail = FALSE)
    }
    return(res)
  }

  # Combine cluster spesific p-values of drugs
  if(length(unique(names(cluster_drugs)))>1){
      combined_p_values <- tapply(drug_list$p_value, drug_list$drug, CombineP)
  }else{
      combined_p_values <- drug_list$p_value
      names(combined_p_values) <- drug_list$drug
  }

# Cell line information
  cell_lines <- subset(cell_data, primary_site == tissue)$cell_id

  # Load drugs metadata for GSE92742 and subset it to tissue of interest and 
# drugs of interest
  drug_metadata_92742 <- col_meta_GSE92742[, c("sig_id", "pert_iname")]
  row.names(drug_metadata_92742) <- drug_metadata_92742$sig_id
  idx <- which(col_meta_GSE92742$cell_id %in% cell_lines & 
        col_meta_GSE92742$pert_iname %in% drugs)
  sig_ids <- col_meta_GSE92742$sig_id[idx]
  drug_metadata_92742 <- drug_metadata_92742[sig_ids, ]

  # Load drug response for GSE92742
  exprs <- as.data.frame(parse_gctx(gse92742_gctx_path, cid=sig_ids)@mat)
  treatments <- colnames(exprs)
  exprs$gene_id <- row.names(exprs)
  tmp <- merge(exprs, gene_meta, by.x="gene_id", by.y="pr_gene_id")
  drug_responses_92742 <- tmp[, c("pr_gene_symbol", treatments)]

  # Load drugs metadata for GSE70138 and subset it to tissue of interest and 
# drugs of interest
  drug_metadata_70138 <- col_meta_GSE70138[, c("sig_id", "pert_iname")]
  row.names(drug_metadata_70138) <- drug_metadata_70138$sig_id
  idx <- which(col_meta_GSE70138$cell_id %in% cell_lines & 
        col_meta_GSE70138$pert_iname %in% drugs)
  sig_ids <- col_meta_GSE70138$sig_id[idx]
  drug_metadata_70138 <- drug_metadata_70138[sig_ids, ]

  # Load drug response for GSE70138
  exprs <- as.data.frame(parse_gctx(gse70138_gctx_path, cid=sig_ids)@mat)
  treatments <- colnames(exprs)
  exprs$gene_id <- row.names(exprs)
  tmp <- merge(exprs, gene_meta, by.x="gene_id", by.y="pr_gene_id")
  drug_responses_70138 <- tmp[, c("pr_gene_symbol", treatments)]

  drug_responses <- merge(drug_responses_92742, drug_responses_70138, 
            by="pr_gene_symbol")
  row.names(drug_responses) <- drug_responses[, 1]
  drug_responses <- drug_responses[, -1]
  drug_metadata <- rbind(drug_metadata_92742, drug_metadata_70138)

# Find DEGs that are common to all clusters
  common_degs <- list()
  for (i in names(deg_list)) {
    ith_deg_list <- deg_list[[i]]
      ith_deg_list <- subset(ith_deg_list, adj.P.Val < 0.05)
  if (length(ith_deg_list) > 0) {
      common_degs[[i]] <- rownames(ith_deg_list)
  }
  }
  common_degs <- Reduce(intersect, common_degs)

# Combine cluster specific DEG scores into a matrix
deg_scores <- data.frame()
  for (i in names(deg_list)) {
    ith_deg_list <- deg_list[[i]]
    if (nrow(deg_scores) == 0) {
      deg_scores <- data.frame(score = ith_deg_list[common_degs, "score"])
    } else {
      tmp <- data.frame(score = ith_deg_list[common_degs,"score"])
      deg_scores <- cbind(deg_scores, tmp)
    }
  }
  deg_scores <- as.matrix(deg_scores)
  row.names(deg_scores) <- common_degs

  deg_scores_mean <- apply(deg_scores, 1, mean)
  names(deg_scores_mean) <- common_degs

  # Calculate drug score
  counter <- 0 
  drug_scores <- list()
  cluster_drug_scores <- list() # for individual cluster analysis
  for (drug in drugs) {
    # Get response from CMap
    treatments <- subset(drug_metadata, pert_iname == drug)$sig_id
    if (length(treatments) > 1) {
      curr_drug_response <- drug_responses[, treatments]
      mean_response <- apply(curr_drug_response, 1, mean)
    } else {
      curr_drug_response <- drug_responses[, treatments]
      mean_response <- curr_drug_response
    }

    drug_stats <- drug_list[drug_list$drug == drug, ]
    drug_score <- 0
  
    colnames(drug_prom_inh_weight) <- tolower(colnames(drug_prom_inh_weight))
    colnames(drug_prom_inh_weight) <- gsub("\\.", "-", colnames(drug_prom_inh_weight))
    drug_lower <- tolower(drug)

    if (drug_lower %in% colnames(drug_prom_inh_weight)) {
      counter <- counter + 1
      cluster_weight_drug <- drug_prom_inh_weight[, drug_lower, drop = FALSE]
    } else {
      cluster_weight_drug  <- data.frame(rep(1, nrow(drug_prom_inh_weight)))
      colnames(cluster_weight_drug) <- c(drug_lower)
      rownames(cluster_weight_drug) <- rownames(drug_prom_inh_weight)
    }
    #print(counter)

    weighted_cluster_scores <- c() # for individual cluster analysis
    for (i in names(deg_list)) {
      #print(i)
      cluster_prop <- drug_stats[drug_stats$cluster == i, "cluster_prop"]
      fdr <- drug_stats[drug_stats$cluster == i, "fdr"]
      p_value <- drug_stats[drug_stats$cluster == i, "p_value"]

      ith_deg_list <- deg_list[[i]]
      ith_deg_list <- subset(ith_deg_list, adj.P.Val < 0.05)

      treatable_degs <- intersect(row.names(ith_deg_list), names(mean_response))
    
      if (length(treatable_degs > 0)) {
        deg_scores <- ith_deg_list[treatable_degs, "score"]

        treated_degs <- -deg_scores*mean_response[treatable_degs]
        treated_degs <- treated_degs[which(treated_degs > 0)]

        treated_degs_ratio <- length(treated_degs)/length(treatable_degs)
        drug_score <- drug_score + (cluster_prop/100) * (-log10(fdr)) * treated_degs_ratio * exp(cluster_weight_drug[i, ]) # linking

        weighted_cluster_scores <- c(weighted_cluster_scores, (cluster_prop/100) * (-log10(fdr)) * treated_degs_ratio * exp(cluster_weight_drug[i, ])) # for individual cluster analysis
      } else {
        weighted_cluster_scores <- c(weighted_cluster_scores, 0) # for individual cluster analysis
      }
    }
    drug_scores[[drug]] <- drug_score
    cluster_drug_scores[[drug]] <- weighted_cluster_scores
  }

  cluster_drug_scores <- as.data.frame(cluster_drug_scores) 
  rownames(cluster_drug_scores) <- names(deg_list)
  colnames(cluster_drug_scores) <- gsub("\\.", "-", colnames(cluster_drug_scores))
  cluster_drug_scores <- t(cluster_drug_scores[, colnames(drug_prom_inh_weight)])
  write.csv(cluster_drug_scores, file = paste0(out_path, "/", disease, "_", "individual_cell_type_drug_scores.csv"), sep = "\t", quote = FALSE)

  drug_scores <- t(as.data.frame(drug_scores))
  out <- data.frame(
    drug_score = drug_scores,
    p_val = combined_p_values[drugs],
    fdr = p.adjust(combined_p_values[drugs], method = "BH")
    )
  rownames(out) <- gsub("\\.", "-", rownames(out))
  out <- out[colnames(drug_prom_inh_weight), ]
  write.csv(out, file = paste0(out_path, "/", disease, "_", "drug_scores.csv"), sep = "\t", quote = FALSE)
  return(out)
}
