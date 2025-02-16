#' @title Estimate the sensitivity/resistance effect of each drug on different cell types in the disease atlas.
#' @description This function computes the p-values and adjusted p-values for each drug-cell type comparison via reverse gene expression pattern matching and K-S test (following the Asgard pipeline), which is an essential step for estimating the sensitivity/resistance effect.
#' @param disease The name of the investigated disease.
#' @param perturbation_matrix_path The path of the CMAP L1000 disease tissue-specific pertuabtion signature matrix.
#' @param gene_info The meta data for genes in the perturbation matrix.
#' @param drug_info The meta data for drugs in the perturbation matrix.
#' @param deg_list The list of intra-cell-type DEGs identified by the get_intra_cell_type_degs function.
#' @param out_path The path where the csv file of the p-values and the adjusted p-values (FDRs) for each pair of drug and cell type is saved.
#' @return The data.frame of the p-values and the adjusted p-values (FDRs) for statistical tests that estimate sensitivity/resistance effect for each pair of cell type (row) and drug (column).
#' @references He B. et al. (2023). *ASGARD is A Single-cell Guided Pipeline to Aid Repurposing of Drugs*. Nature Communications.
#' @export
#' @import Asgard

compute_drug_sens_res <- function(disease, perturbation_matrix_path, gene_info, drug_info, deg_list, out_path="results") {

  library(Asgard)
  drug_ref_profiles = GetDrugRef(drug.response.path = perturbation_matrix_path,
                               probe.to.genes = gene_info, 
                               drug.info = drug_info)
                          
  drug_sens_res = GetDrug(gene.data = deg_list, 
                          drug.ref.profiles = drug_ref_profiles, 
                          repurposing.unit = "drug", 
                          connectivity = "negative", 
                          drug.type = "FDA")
  drug_ident_res_df <- drug_sens_res
  for (ci in names(drug_ident_res_df)) {
    df <- drug_ident_res_df[[ci]]
    df <- df[!duplicated(df$Drug.name), ]
    drug_ident_res_df[[ci]] <- df
  }
  write.csv(drug_ident_res_df, file = paste0(out_path, "/", disease, "_", "p_adj_perturb_sig.csv"), row.names = F)
  return(drug_sens_res)
}