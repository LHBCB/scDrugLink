library(Seurat)
library(Asgard)
library(cmapR)
library(effsize)

# This code is an example of using scDrugLink to repurpose drugs for glioblastoma
# Functions will be factorised and documentations will be prepared upon publication of the paper.

set.seed(42)

# read data
diz <- "GBM"
dat <- readRDS("data/gbm_final.rds")

dat@meta.data$cell_type <- dat@meta.data$finalCelltype_15Sept
dat@meta.data$diz <- dat@meta.data$histology
ctl <- "control"

print(unique(dat@meta.data$cell_type))

# build drug target d2c_mat
drug_label_df <- read.csv("drug_target_sig_label/drug_labels_273.csv")
drug_target_df <- read.csv("drug_target_sig_label/drug_targets_273.csv", row.names = 1)

gene_info <- read.table(file="DrugReference/central-nervous-system_gene_info.txt",sep="\t",header = T,quote = "")
gene_list <- Reduce(intersect, list("seurat"=rownames(dat@assays["RNA"]$RNA), "drug"=gene_info$Gene.Symbol))

compute_drug_target_d2c <- function(dat, gene_list, drug_target_df, out_path, n_bins = 25, ctrl_size = 50, seed = 42) {
  
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

d2c_mat <- compute_drug_target_d2c(dat, gene_list = gene_list, drug_target_df = drug_target_df, out_path = "scdruglink_results")


# compute drug promotion/inhibition weight
get_cell_type_drug_weight <- function(dat, d2c_mat, disease, out_type="weight"){
  
  Idents(dat) <- "cell_type"
  cell_type_list <- unique(Idents(dat))

  total_num_cells <- ncol(dat)
  
  weight_list <- list()
  prop_weight_list <- list()
  p_adj_list <- list()
  for (cell_type in cell_type_list) {
    print(paste0("Processing cell type: ", cell_type))
    subset_dat <- subset(dat, idents = cell_type)
    subset_cell_labels <- ifelse(subset_dat@meta.data$diz==disease, 1, 0)
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
  write.csv(p_adj_df, paste0("scdruglink_results/", disease, "_p_adj_target.csv"), row.names=T)

  # weight
  temp_list <- lapply(weight_list, function(x) as.data.frame(t(unlist(x))))
  weight_df <- do.call(rbind, temp_list)
  colnames(weight_df) <- rownames(d2c_mat)
  
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

cluster_weight_df <- get_cell_type_drug_weight(dat, d2c_mat, disease="GBM", out_type="weight")

# Drug sensitivity/resistance estimated based on the Asgard pipeline
tissue <- "central nervous system"
Case=c(diz)
Control=c(ctl)
DefaultAssay(dat) <- "RNA"

Gene.list <- list()
C_names <- NULL
for(i in unique(dat@meta.data$cell_type)){

  print(i)

  Idents(dat) <- "cell_type"
  c_cells <- subset(dat, cell_type == i)
  Idents(c_cells) <- "diz"

  cells.1 <- WhichCells(c_cells, ident = diz)
  cells.2 <- WhichCells(c_cells, ident = ctl)
  if (length(cells.1) < 3 || length(cells.2) < 3) {
      message("Skipping ", i, " due to insufficient cell count.")
      next
    }

  C_data <- FindMarkers(c_cells, ident.1 = diz, ident.2 = ctl)
  
  C_data_for_drug <- data.frame(row.names=row.names(C_data), score=C_data$avg_log2FC, adj.P.Val=C_data$p_val_adj, P.Value=C_data$p_val)
  Gene.list[[i]] <- C_data_for_drug

  C_names <- c(C_names,i)
}
names(Gene.list) <- C_names

# Load tissue specific drug reference produced by PrepareReference function in ASGARD
print("#############################")
print("Compute drug effect based on perturbation signatures")

gene_info<-read.table(file="DrugReference/central-nervous-system_gene_info.txt",sep="\t",header = T,quote = "")
drug_info<-read.table(file="DrugReference/central-nervous-system_drug_info.txt",sep="\t",header = T,quote = "")
drug.ref.profiles = GetDrugRef(drug.response.path = 'DrugReference/central-nervous-system_rankMatrix.txt',
                               probe.to.genes = gene_info, 
                               drug.info = drug_info)
                          
Drug.ident.res = GetDrug(gene.data = Gene.list, 
                        drug.ref.profiles = drug.ref.profiles, 
                        repurposing.unit = "drug", 
                        connectivity = "negative", 
                        drug.type = "FDA")

# compute drug therapeutic score by linking drug promotion/inhibition and sensitivity/resistance estimations
# this function inherits the DrugScore function of ASGARD
DrugLinkScore <- function(cell_metadata, cluster_degs, cluster_drugs, tissue,
					  gse70138_gctx_path, gse92742_gctx_path, cluster_weight_df,
					  clusters = NULL, case = NULL, fda_drugs_only = TRUE) {

	# Subset input data to the set of clusters we are interested in 
    if (length(clusters) > 0) {
    	clusters = intersect(clusters, unique(cell_metada$cluster))
      	cell_metadata = subset(cell_metadata, cluster %in% clusters)
      	cluster_drugs = cluster_drugs[clusters]
      	cluster_degs = cluster_degs[clusters]
    }

    # Calculate cluster proportions in diseased tissue
    if (length(case) > 0) {
      	cell_metadata <- subset(cell_metadata, diz %in% case)
    }
    clustering <- cell_metadata$cluster
    cluster_sizes <- table(clustering)
    cluster_sizes <- cluster_sizes[which(cluster_sizes > 3)]
    cluster_prop <- round(100*cluster_sizes/nrow(cell_metadata), 2) 

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
    for (i in names(cluster_degs)) {
    	ith_cluster_degs <- cluster_degs[[i]]
      	ith_cluster_degs <- subset(ith_cluster_degs, adj.P.Val < 0.05)
		if (length(ith_cluster_degs) > 0) {
	    	common_degs[[i]] <- rownames(ith_cluster_degs)
		}
    }
	  common_degs <- Reduce(intersect, common_degs)

	# Combine cluster specific DEG scores into a matrix
	deg_scores <- data.frame()
    for (i in names(cluster_degs)) {
    	ith_cluster_degs <- cluster_degs[[i]]
    	if (nrow(deg_scores) == 0) {
    		deg_scores <- data.frame(score = ith_cluster_degs[common_degs, "score"])
    	} else {
    	    tmp <- data.frame(score = ith_cluster_degs[common_degs,"score"])
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
    
      colnames(cluster_weight_df) <- tolower(colnames(cluster_weight_df))
      colnames(cluster_weight_df) <- gsub("\\.", "-", colnames(cluster_weight_df))
      drug_lower <- tolower(drug)

      if (drug_lower %in% colnames(cluster_weight_df)) {
        counter <- counter + 1
        cluster_weight_drug <- cluster_weight_df[, drug_lower, drop = FALSE]
      } else {
        cluster_weight_drug  <- data.frame(rep(1, nrow(cluster_weight_df)))
        colnames(cluster_weight_drug) <- c(drug_lower)
        rownames(cluster_weight_drug) <- rownames(cluster_weight_df)
      }
      print(counter)

      for (i in names(cluster_degs)) {
        #print(i)
        cluster_prop <- drug_stats[drug_stats$cluster == i, "cluster_prop"]
        fdr <- drug_stats[drug_stats$cluster == i, "fdr"]
        p_value <- drug_stats[drug_stats$cluster == i, "p_value"]

        ith_cluster_degs <- cluster_degs[[i]]
        ith_cluster_degs <- subset(ith_cluster_degs, adj.P.Val < 0.05)

        treatable_degs <- intersect(row.names(ith_cluster_degs), names(mean_response))
     
        if (length(treatable_degs > 0)) {
          deg_scores <- ith_cluster_degs[treatable_degs, "score"]

          treated_degs <- -deg_scores*mean_response[treatable_degs]
          treated_degs <- treated_degs[which(treated_degs > 0)]

          treated_degs_ratio <- length(treated_degs)/length(treatable_degs)
          drug_score <- drug_score + (cluster_prop/100) * (-log10(fdr)) * treated_degs_ratio * exp(cluster_weight_drug[i, ]) # linking
        }
      }
      drug_scores[[drug]] <- drug_score
    }
	drug_scores <- t(as.data.frame(drug_scores))
  out <- data.frame(
		drug_score = drug_scores,
		p_val = combined_p_values[drugs],
		fdr = p.adjust(combined_p_values[drugs], method = "BH")
	)
  return(out)
}

# compute drug score
gse92742_gctx_path <- "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
gse70138_gctx_path <- "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx"

SC_subset <- subset(dat, subset = cell_type %in% C_names)
cell_metadata <- SC_subset@meta.data
cell_metadata$cluster <- SC_subset@meta.data$cell_type

Drug.score <- DrugLinkScore(cell_metadata, cluster_degs = Gene.list, 
                        cluster_drugs = Drug.ident.res, tissue = tissue, 
                        case = Case, gse92742_gctx_path = gse92742_gctx_path, 
                        gse70138_gctx_path = gse70138_gctx_path, cluster_weight = cluster_weight_df,
                        fda_drugs_only=TRUE) # clusters can be set to compute drug score for each cluster
rownames(Drug.score) <- gsub("\\.", "-", rownames(Drug.score))
write.csv(Drug.score, file = paste0("scdruglink_results/", diz, "_", "drug_scores_improved.csv"), sep = "\t", quote = FALSE)

# eval
library(pROC)
library(ROCR)
library(PRROC)

eval <- function(drug_scores){
  
  if (diz == "Alzheimer disease") {
    labels <- drug_scores$ad_label_final
  } else if (diz == "GBM") {
    labels <- drug_scores$gbm_label_final
  } else if (diz == "True") {
    labels <- drug_scores$ms_label_final
  }

  pred <- prediction(predictions = drug_scores$drug_score, labels = labels)
  perf <- performance(pred, measure = "tpr", x.measure = "fpr")
  
  auc_perf <- performance(pred, measure = "auc")
  auc_value <- auc_perf@y.values[[1]]
  
  pos_scores <- drug_scores$drug_score[labels == 1]
  neg_scores <- drug_scores$drug_score[labels == 0]
  pr_result <- pr.curve(scores.class0 = pos_scores, scores.class1 = neg_scores, curve = FALSE)
  aupr_value <- pr_result$auc.integral

  return(c(auc_value, aupr_value))
}

eval_diz <- function(diz, Drug.score) {
  print(diz)

  Drug.score.sub <- merge(drug_label_df, Drug.score, by.x = "drug_name_lower", by.y = "row.names", all = FALSE)
  eval_values <- eval(Drug.score.sub)
  print(eval_values)
  return(eval_values)
}

results <- eval_diz(diz, Drug.score)

df <- data.frame(
  disease = diz,
  auc     = results[1],
  aupr    = results[2],
  stringsAsFactors = FALSE
)
write.csv(df, paste0("scdruglink_results/", diz, "_metrics.csv"), row.names = F)
