library('Seurat')
library('Asgard')

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
tissue <- "central nervous system"

#Case sample names
Case=c(diz)
#Control sample names
Control=c(ctl)
#Get differential genes from Seurat (Wilcoxon Rank Sum test)
DefaultAssay(SC.integrated) <- "RNA"

set.seed(42)

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

  C_names <- c(C_names,i)
}

names(Gene.list) <- C_names

############################################
#Load tissue specific drug reference produced by PrepareReference function as mentioned above. Please select proper tissue accroding to the disease.
my_gene_info<-read.table(file="DrugReference/central-nervous-system_gene_info.txt",sep="\t",header = T,quote = "")
my_drug_info<-read.table(file="DrugReference/central-nervous-system_drug_info.txt",sep="\t",header = T,quote = "")
drug.ref.profiles = GetDrugRef(drug.response.path = 'DrugReference/central-nervous-system_rankMatrix.txt',
                               probe.to.genes = my_gene_info, 
                               drug.info = my_drug_info)

#Repurpose mono-drugs for every cell type                               
Drug.ident.res = GetDrug(gene.data = Gene.list, 
                        drug.ref.profiles = drug.ref.profiles, 
                        repurposing.unit = "drug", 
                        connectivity = "negative", 
                        drug.type = "FDA")

############################################
# Change the following two lines with the paths on your computer
gse92742_gctx_path <- "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
gse70138_gctx_path <- "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx"

SC_subset <- subset(SC.integrated, subset = celltype %in% C_names)
cell_metadata <- SC_subset@meta.data
cell_metadata$cluster <- SC_subset@meta.data$celltype

#cell_metadata <- SC.integrated@meta.data
#cell_metadata$cluster <- SC.integrated@meta.data$celltype

Drug.score <- DrugScore(cell_metadata, cluster_degs = Gene.list, 
                        cluster_drugs = Drug.ident.res, tissue = tissue, 
                        case = Case, gse92742_gctx_path = gse92742_gctx_path, 
                        gse70138_gctx_path = gse70138_gctx_path,
                        fda_drugs_only=TRUE) # clusters can be set to compute drug score for each cluster
rownames(Drug.score) <- gsub("\\.", "-", rownames(Drug.score))
write.csv(Drug.score, "asgard_results/results.csv")
