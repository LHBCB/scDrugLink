library(scDrugPrio)

lit_ppi <- readRDS("lit_ppin_gene_symbol.rds")
dt <- read.csv("data/drug_targets_273.csv")

###############################
## Make drug target matrix
###############################
print("###############################")
print("Make drug target matrix")
library(dplyr)
library(tidyr)
dt_long <- dt %>% separate_rows(gene_names, sep = ";")

drug_target_matrix <- create_drug_target_matrix(drugID = dt_long$drug_name, target = dt_long$gene_names)
print(drug_target_matrix[1:5,1:6])

############################################################################################################################
## Selection of the largest connected component (LCC) of the protein-protein interaction network (PPIN)
############################################################################################################################
print("###############################")
print("Selection of the largest connected component (LCC) of the protein-protein interaction network (PPIN)")
lit_ppi <- as.matrix(lit_ppi[,3:4])
print(head(lit_ppi))
print(dim(lit_ppi))
# n unique proteins in PPIN
print(length(unique(as.vector(lit_ppi))))

# Select LCC
ppin <- ppin_formatting(lit_ppi)
print(dim(ppin))
# n unique proteins in LCC of PPIN
print(length(unique(as.vector(ppin))))

# Select unique drug target combinations found in LCC of PPIN
ppin <- ppin_formatting(lit_ppi)
dim(ppin)
length(unique(as.vector(ppin)))

sum(!is.na(drug_target_matrix))# number of drug targets in drug_target_matrix
sum(drug_target_matrix[!is.na(drug_target_matrix)] %in% unique(as.vector(ppin)))# number of drug targets in the drug_target_matrix that were found in PPIN
sum(colSums(!is.na(drug_target_matrix)) == 0)# are there drugs with no drug targets in PPIN?
sum(duplicated(t(drug_target_matrix))) # are there drugs with identical drug targets?

## Preparation of drug target matrix for network distance calculation.
drug_target_matrix <- prepare_drug_target_matrix_for_network_distance_calculation(ppin, drug_target_matrix, file_name = "in_lit_ppin", out_dir = "scdrugprio_results")
same_drugs <- read.table(file = "scdrugprio_results/SAME_DRUG_TARGETS_in_lit_ppin.txt", sep="\t", header = T, stringsAsFactors = F)
head(same_drugs)

############################################################################################################################
## Compute DEGs
############################################################################################################################
print("###############################")
print("Compute DEGs")
library(Seurat)

diz <- "GBM"
SC.integrated <- readRDS("gbm_final.rds")

if (diz == "GBM") {
  SC.integrated@meta.data$celltype <- SC.integrated@meta.data$finalCelltype_15Sept
  SC.integrated@meta.data$type <- SC.integrated@meta.data$histology
  ctl <- "control"
} else if (diz == "AD") {
  SC.integrated@meta.data$celltype <- SC.integrated@meta.data$author_cell_type
  SC.integrated@meta.data$type <- SC.integrated@meta.data$disease
  ctl <- "TLE"
} else if (diz == "MS") {
  SC.integrated@meta.data$celltype <- SC.integrated@meta.data$labels
  SC.integrated@meta.data$type <- SC.integrated@meta.data$MS
  ctl <- "IIH"
}
SC.integrated@meta.data$sample <- SC.integrated@meta.data$type
print(SC.integrated)
tissue <- "central nervous system"
#Case sample names
Case=c(diz)
#Control sample names
Control=c(ctl)
#Get differential genes from Seurat (Wilcoxon Rank Sum test)
DefaultAssay(SC.integrated) <- "RNA"

Gene.list <- list()
C_names <- NULL

for(i in unique(SC.integrated@meta.data$celltype)){
  
  print(i)
  
  Idents(SC.integrated) <- "celltype"
  c_cells <- subset(SC.integrated, celltype == i)
  Idents(c_cells) <- "type"
  
  cells.1 <- WhichCells(c_cells, ident = diz)
  cells.2 <- WhichCells(c_cells, ident = ctl)
  if (length(cells.1) < 3 || length(cells.2) < 3) {
    # Print a message and skip the current iteration
    message("Skipping ", i, " due to insufficient cell count.")
    next
  }
  
  C_data <- FindMarkers(c_cells, ident.1 = diz, ident.2 = ctl)
  
  C_data_for_drug <- data.frame(row.names=row.names(C_data), score=C_data$avg_log2FC, adj.P.Val=C_data$p_val_adj, P.Value=C_data$p_val) ##for Seurat version > 4.0, please use avg_log2FC instead of avg_logFC
  Gene.list[[i]] <- C_data_for_drug
  
  if(exists("C_data_for_drug")){
    write.table(C_data_for_drug, file = paste("scdrugprio_results/Cluster_", i, ".txt", sep=""), sep="\t",col.names=NA,row.names=T)
    rm(C_data_for_drug)
  }
  C_names <- c(C_names,i)
}
head(read.table(file = "scdrugprio_results/Cluster_Mg_1_1 .txt", sep="\t", header = T, stringsAsFactors = F))

# DEG summary for all clusters
outdir = "scdrugprio_results"
files <- list.files(path = outdir, pattern = "Cluster_")

degs.clusters.h.vs.s <- matrix(NA, nrow = SC.integrated@assays$RNA@counts@Dim[1], ncol = length(files))

for(i in 1:ncol(degs.clusters.h.vs.s)){
  temp <- read.table(file = paste(outdir, "/", files[i], sep=""), sep="\t", header = T)
  temp <- temp[temp$adj.P.Val < 0.05,]
  temp <- as.character(temp[order(as.numeric(temp$adj.P.Val), decreasing = F),1]) # order by significance
  if(length(temp)>0){
    degs.clusters.h.vs.s[1:length(temp),i] <- temp
  }
  rm(temp)
}
colnames(degs.clusters.h.vs.s) <- unlist(strsplit(files, split = ".txt"))

# remove empty rows
degs.clusters.h.vs.s <- degs.clusters.h.vs.s[rowSums(!is.na(degs.clusters.h.vs.s))>0,]
head(degs.clusters.h.vs.s)
write.table(degs.clusters.h.vs.s, file = paste(outdir, "/SUMMARY_sig_adj_SEURAT_DEGs_log(1.5)_all_clusters.txt", sep=""),sep="\t", col.names = T, row.names = F)

############################################################################################################################
## NicheNet ligand activity analysis
############################################################################################################################
print("###############################")
print("NicheNet ligand activity analysis")

cell_IDs <- SC.integrated$finalCelltype_15Sept
unique_cl <- unique(SC.integrated$finalCelltype_15Sept)
cell_IDs <- foreach(i = unique_cl) %do% {
  temp <- cell_IDs[as.character(cell_IDs) == i]
  return(names(temp))
}
names(cell_IDs) <- paste("Cluster_", unique_cl, sep="")

# background gene calculation
bg <- background_genes_NicheNet(data = as.matrix(SC.integrated@assays$RNA@counts), cell_IDs = cell_IDs)

# make sure DEGs and background genes are annotated in human entrez symbols
degs_transl <- degs.clusters.h.vs.s
head(degs_transl)
degs_transl <- degs_transl[rowSums(!is.na(degs_transl))>0, colSums(!is.na(degs_transl))>0]
bg <- bg[rowSums(!is.na(bg))>0,]
head(bg)

dir.create(path = "scdrugprio_results/NicheNet", showWarnings = F)

ligand_target <- as.matrix(readRDS("scdrugprio_data/ligand_target_matrix.rds"))
lr_network = as.matrix(readRDS("scdrugprio_data/lr_network.rds"))

all_ligand_activity <- NicheNet_ligand_activity_analysis(degs = degs_transl, background_genes = bg, 
                                                         out_dir = "scdrugprio_results/NicheNet",
                                                         ligand_target = ligand_target,
                                                         lr_network = lr_network,
                                                         cores = 2)
print(head(all_ligand_activity))

############################################################################################################################
## Intercellular centrality
############################################################################################################################
print("###############################")
print("Intercellular centrality")

# Calculate centrality for cell types based on ligand activity outcomes
intercellular_centrality <- NicheNet_cell_type_centrality(all_ligand_activity, out_dir = "scdrugprio_results/NicheNet/")
print(intercellular_centrality)

############################################################################################################################
## Intracellular disease models
############################################################################################################################
print("###############################")
print("Intracellular disease models")

degs.clusters.h.vs.s <- degs.clusters.h.vs.s[rowSums(!is.na(degs.clusters.h.vs.s))>0, colSums(!is.na(degs.clusters.h.vs.s))>0]

# run intracellular centrality analysis
dir.create("scdrugprio_results/intracellular_centrality", showWarnings = F)

intra_cent <- intracellular_drug_centrality(ppin = ppin,
                                            drug_target_matrix = drug_target_matrix,
                                            degs = degs.clusters.h.vs.s,
                                            file_name = "intracellular_centrality_drugs_lit_PPIN",
                                            centrality_alg = "eigenvector centralities",
                                            out_dir = "scdrugprio_results/intracellular_centrality")
print(head(intra_cent))

############################################################################################################################
## Network distance calculation
############################################################################################################################
print("###############################")
print("Network distance calculation")

dir.create("scdrugprio_results/network_distances")
for(i in 1:ncol(degs.clusters.h.vs.s)){ # For every cell types DEGs
  # calculate average closest distances between all drugs and the degs
  # if possible on your machine this should be parallelized (cores > 1), though the function will be memory heavy
  average_closest_distance_network_drug_screening(ppin = ppin,
                                                  drug_target_matrix = drug_target_matrix,
                                                  disease_genes = degs.clusters.h.vs.s[,i], # using symbols rather than IDs
                                                  file_name = colnames(degs.clusters.h.vs.s)[i],
                                                  cores = 25,
                                                  out_dir = "scdrugprio_results/network_distances")
}

#data <- read.table(file = paste("scdrugprio_results/network_distances/drug-disease_closest_distances_vs_random_bin_adjusted__Cluster_0.txt", sep=""), sep ="\t", header = T, stringsAsFactors = F)
#print(head(data))

lf  <- list.files(path = "scdrugprio_results/network_distances/", pattern = "drug-disease_closest_distances_vs_random_bin_adjusted__")
separate_unique_drug_target_combinations_into_individual_drugs(files = lf, 
                                                               remove_part_of_file_name_for_output = "drug-disease_closest_distances_vs_random_bin_adjusted__",
                                                               same_drugs = same_drugs, 
                                                               disease_specific_drugs = RA_drugs, 
                                                               output_file_name_add_on = "INDIVIDUAL_DRUGS_", 
                                                               in_dir = "scdrugprio_results/network_distances/", 
                                                               out_dir = "scdrugprio_results/network_distances/")

############################################################################################################################
## Evaluate pharmacological actions on targeted DEGs fold change
############################################################################################################################
print("###############################")
print("Evaluate pharmacological actions on targeted DEGs fold change")

# drugID = dt_long$drug_name, target = dt_long$gene_names

# translate DEG files
lf <- list.files(path = "scdrugprio_results", pattern = "Cluster_")
dir.create(path = "scdrugprio_results/translated_DEGs", showWarnings = F)
for(i in 1:length(lf)){
  temp <- read.table(file = paste("scdrugprio_results/",lf[i],sep=""), sep = "\t", header = T, stringsAsFactors = F)
  temp <- temp[temp$adj.P.Val < 0.05,] # watch out for the different names given by Seurat
  #temp[,1] <- translation_mouse_human[match(temp[,1], translation_mouse_human[,3]),1]
  if(any(!is.na(temp[,1]))){
    temp <- temp[!is.na(temp[,1]),]
    write.table(temp, file = paste("scdrugprio_results/translated_DEGs/",lf[i],sep=""), sep="\t", col.names = T, row.names = F)
  }
  rm(temp)
}

# get paths to translated DEG files
deg_files <- list.files("scdrugprio_results/translated_DEGs", pattern = "Cluster_", full.names = T)
print(deg_files)
#> [1] "scdrugprio_results/translated_DEGs/Cluster_0.txt"
#> [2] "scdrugprio_results/translated_DEGs/Cluster_1.txt"
# get paths to average closest netork distances
drug_dists <- list.files("scdrugprio_results/network_distances", pattern = "drug-disease_closest_distance", full.names = T)
print(drug_dists)
#> [1] "scdrugprio_results/network_distances/drug-disease_closest_distances_vs_random_bin_adjusted__Cluster_0.txt"
#> [2] "scdrugprio_results/network_distances/drug-disease_closest_distances_vs_random_bin_adjusted__Cluster_1.txt"
# save names
save_name <- unlist(strsplit(list.files("scdrugprio_results/translated_DEGs", pattern = "Cluster_"), split = ".txt"))
# make output folder
dir.create("scdrugprio_results/FC_criteria_checking", showWarnings = F)

fc_evaluation <- pharma_effect_on_fold_change(drug_dist_files = drug_dists, deg_files = deg_files, pharma_effect = dt_long[,c(1, 2, 3)], save_names = save_name, out_dir = "scdrugprio_results/FC_criteria_checking")

############################################################################################################################
## Drug candidate selection based on network distance calculations
############################################################################################################################
print("###############################")
print("Drug candidate selection based on network distance calculations")

# Apply drug selection network criteria
#fc_evaluation <- fc_evaluation[as.numeric(fc_evaluation[,3]) < 1,]
#fc_evaluation <- fc_evaluation[as.numeric(fc_evaluation[,7]) < 0.05,]
#fc_evaluation <- fc_evaluation[order(fc_evaluation[,1]),]
#fc_evaluation <- cbind(Drug_name = dt_long[match(fc_evaluation[,1], dt_long[,1]),2], fc_evaluation)
#write.table(fc_evaluation, file = "scdrugprio_results/FC_criteria_checking/SUMMARY_only_drugs_passing_network_criteria.txt", sep="\t", col.names = T, row.names = F)

############################################################################################################################
## Final drug candidate selection and ranking
############################################################################################################################
print("###############################")
print("Final drug candidate selection and ranking")

dir.create("scdrugprio_results/Final_ranking")
drug_rank <- final_drug_prioritization(fc_evaluation = fc_evaluation, 
                                       #pos_DrugID = 2,
                                       #pos_clusterID = 19, 
                                       #keep = c(1:4,18), 
                                       inter_cent = intercellular_centrality[5,],  #
                                       intra_cent = intra_cent[,ncol(intra_cent)], #
                                       file_name = "FINAL_drug_ranking", 
                                       out_dir = "scdrugprio_results/Final_ranking")
head(drug_rank[,-c(5:7)])

write.csv(drug_rank, "scdrugprio_results/final_drug_scores.csv")
