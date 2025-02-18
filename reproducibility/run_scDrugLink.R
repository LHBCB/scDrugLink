library(Seurat)
library(scDrugLink)

# read the disease atlas
dat <- readRDS("gbm_final.rds")

# standardise the column names in meta data
dat@meta.data$cell_type <- dat@meta.data$finalCelltype_15Sept
dat@meta.data$disease <- dat@meta.data$histology
disease <- "GBM"
control <- "control"

# build Drug2Cell matrix based on drug targets
print('##########################################')
print('Build d2c')

gene_info <- read.table(file = "DrugReference/central-nervous-system_gene_info.txt", sep="\t", header = T, quote = "")
gene_list <- Reduce(intersect, list("seurat" = rownames(dat@assays["RNA"]$RNA), "drug" = gene_info$Gene.Symbol))

dir.create("results", showWarnings = FALSE)
d2c_mat <- build_drug_target_d2c(dat, gene_list = gene_list, drug_target_df = cns_drug_targets, out_path = "results")

# estimate drug promotion/inhibition effects
print('##########################################')
print('Estimate drug promotion/inhibition')
drug_prom_inh_weight <- compute_drug_prom_inh(dat, d2c_mat, disease=disease, out_type="weight", out_path = "results")
print(drug_prom_inh_weight[1:6, 1:6])

# identify intra-cell type degs
print('##########################################')
print('Identify intra-cell type degs')
deg_list <- get_intra_cell_type_degs(dat, disease= disease, control = control)

# estimate drug sensitivity/resistance effects
print('##########################################')
print('Estimate drug sensitivity/resistance')
drug_info <- read.table(file = "DrugReference/central-nervous-system_drug_info.txt", sep="\t", header = T, quote = "")
perturbation_matrix_path <- "DrugReference/central-nervous-system_rankMatrix.txt"

drug_sens_res_pval <- compute_drug_sens_res(disease = disease, 
                                            perturbation_matrix_path = perturbation_matrix_path, 
                                            gene_info = gene_info, 
                                            drug_info = drug_info, 
                                            deg_list = deg_list, 
                                            out_path = "results")
pval_df <- read.csv(paste0("results/" , disease, "_p_adj_perturb_sig.csv"))
print(pval_df[1:6, 1:6])

# compute drug score
print('##########################################')
print('Compute final drug score')

gse92742_gctx_path <- "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
gse70138_gctx_path <- "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx"

drug_score <- compute_scdruglink_score(dat, 
                                      deg_list = deg_list, 
                                      drug_sens_res_pval = drug_sens_res_pval,
                                      tissue = "central nervous system", 
                                      case = disease,
                                      gse92742_gctx_path = gse92742_gctx_path, 
                                      gse70138_gctx_path = gse70138_gctx_path, 
                                      drug_prom_inh_weight = drug_prom_inh_weight,
                                      out_path="results")
print(head(drug_score))

cell_type_drug_score <- read.csv(paste0("results/" , disease, "_individual_cell_type_drug_scores.csv"), row.names = 1)
cell_type_drug_score[1:6, 1:6]
