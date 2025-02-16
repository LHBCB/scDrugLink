# scDrugLink: Single-Cell Drug Repurposing for CNS Diseases via Computationally Linking Drug Targets and Perturbation Signatures
`scDrugLink` is an R package designed to integrate disease single-cell transcriptomics data with drug target information (for drug promotion/inhibition effect estimation) and drug perturbation signatures (for drug sensitivity/resistance effect estimation) to compute robust drug therapeutic scores. It is applicable to a variety of diseases and tissues, such as glioblastoma (GBM), multiple sclerosis (MS), and Alzheimer's disease (AD) in central nervous system (CNS). Drug repurposing by therapeutic score computation and ranking can be performed both at the cell-type level and the whole-atlas level.

![scDrugLink pipeline](scDrugLink_pipeline.jpg)

## Installation
`scDrugLink` requires several R packages: `Asgard`, `cmapR`, `effsize`, and `Seurat`. Please first install `devtools` and `BiocManager` if it is not already done.
```
install.packages('devtools')
devtools::install_github("lanagarmire/Asgard")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("cmapR")

install.packages("effsize")
install.packages("Seurat")
```
Then install scDrugLink
```
devtools::install_github("lhbcb/scDrugLink")
```
Load scDrugLink
```
library(scDrugLink)
```

## Tutorial
### 1. Load disease scRNA-seq data and prepare drug targets and perturbation signatures
Here, we take drug repurposing for GBM as an example to demonstrate the steps required to reproduce the results for the disease in the scDrugLink paper. The disease scRNA-seq data can be downloaded from GliomaAtlas (https://gbmvisium.snu.ac.kr/seuratObjects/); two GBM patient samples (SNU24 and SNU25) and one control sample (SNU57normal) are used as the GBM dataset for our study. 
```
library(Seurat)
seurat_obj <- readRDS("gbm.rds")
disease <- "GBM" # the disease condition
control <- "control" # the control condition
# standardise the column names in meta data
seurat_obj@meta.data$cell_type <- seurat_obj@meta.data$finalCelltype_15Sept
seurat_obj@meta.data$disease <- seurat_obj@meta.data$histology
# remove the "Doublet" cell type and cell types with < 3 diseased or control cells
seurat_obj_subset <- subset(x = seurat_obj, subset = cell_type != "Doublet") 
meta_data <- seurat_obj_subset@meta.data
counts <- table(meta_data$cell_type, meta_data$disease)
remain_cell_types <- rownames(counts)[apply(counts, 1, function(x) all(x >= 3))]
dat <- subset(seurat_obj_subset, cell_type %in% remain_cell_types)
```
The drug targets can be sourced as follows: first download the complete dataset of drugs from DrugBank’s releases (available in XML format at https://go.drugbank.com/releases); then, use an XML parser (such as the xml2 package in R) to extract the text data, and apply keyword or pattern matching techniques to identify the gene targets for each drug. In our study, we have prepared the targets for 273 drugs effective in the CNS tissue, directly accessible in the `scDrugLink` package.
```
head(cns_drug_targets)
#    drug_name
#1      Biotin
#2  Calcitriol
#3 Calcifediol
#4   Icosapent
#5   Menadione
#6 Pravastatin
#                                                                                                                                                        gene_names
#1                                                                                                   PCCB;HLCS;MCCC2;ACACB;MCCC1;PC;PCCA;ACACA;CYP1B1;SLC5A6;SLC5A8
#2                                                                                                                                     VDR;HOXA10;CYP24A1;CYP3A4;GC
#3                                                                                                                                              VDR;CYP27B1;CYP24A1
#4                                                                                     PTGS2;PTGS1;PPARG;PPARD;FFAR1;SLC8A1;FABP7;FADS1;ACSL4;TRPV1;PPARA;ALOX5;ALB
#5 GGCX;VKORC1;VKORC1L1;F2;F9;PROC;PROS1;NQO2;NQO1;BGLAP;CYP1A2;CYP2A6;CYP1B1;CYP2B6;CYP2C8;CYP2C9;CYP2C19;CYP2D6;CYP2E1;CYP3A4;CYP3A5;CYP3A7;XDH;AOX1;MTHFR;CYP1A1
#6                                                    HMGCR;HDAC2;SLCO1B1;SLCO2B1;ABCB1;SLCO1A2;SLC22A6;SLC22A8;ABCC2;SLC22A11;ABCG2;SLC22A7;SLC16A1;ABCB11;SLCO1B3
```
The drug perturbation signatures can be obtained following Asgard's step (https://github.com/lanagarmire/Asgard): download the Connectivity Map (CMAP) L1000 perturbational profiles GSE70138 and GSE92742 from GEO; then generate tissue specific drug references using the `PrepareReference` function.
```
library(Asgard)
PrepareReference(cell.info="GSE70138_Broad_LINCS_cell_info_2017-04-28.txt",
                 gene.info="GSE70138_Broad_LINCS_gene_info_2017-03-06.txt",
                 GSE70138.sig.info = "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt",
                 GSE92742.sig.info = "GSE92742_Broad_LINCS_sig_info.txt",
                 GSE70138.gctx = "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx",
                 GSE92742.gctx = "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx",
                 Output.Dir = "DrugReference/"
                )
```
The "central nervous system" tissue will be utlised for GBM drug repurposig in the subsequent steps.
### 2. Build Drug2Cell matrix based on drug targets
The Drug2Cell matrix is constructed by calculating the average gene expression for each drug's target group in each cell and adjusting for background expression by subtracting baseline biases. This is an R implementation of the Drug2Cell Python pipeline designed by Kanemaru et al. (PMID: 37438528).
```
gene_info <- read.table(file = "DrugReference/central-nervous-system_gene_info.txt", sep="\t", header = T, quote = "")
gene_list <- Reduce(intersect,
                    list("seurat" = rownames(dat@assays["RNA"]$RNA),
                         "drug" = gene_info$Gene.Symbol))

dir.create("results", showWarnings = FALSE)
d2c_mat <- build_drug_target_d2c(dat,
                                 gene_list = gene_list,
                                 drug_target_df = cns_drug_targets,
                                 out_path = "results")
```
### 3. Estimate drug promotion/inhibition effects
Each drug's promotion/inhibition effect on each cell type is calculated by comparing scores derived from control and diseased cells using within-cell-type Wilcoxon rank-sum test, and adjusted for multiple comparisons using the Benjamini-Hochberg procedure.
```
drug_prom_inh_weight <- compute_drug_prom_inh(dat,
                                              d2c_mat = d2c_mat,
                                              disease = disease,
                                              out_type = "weight",
                                              out_path = "results")
```
### 4. Identify intra-cell-type differentially expressed genes (DEGs)
This step processes each cell type individually, performs statistical testing to find markers, and stores the results (log fold change and adjusted p-values) for each cell type in a list.
```
deg_list <- get_intra_cell_type_degs(dat, disease= disease, control = control)
```
### 5. Estimate drug sensitivity/resistance effects
This step computes the p-values and adjusted p-values for each drug-cell type comparison via reverse gene expression pattern matching and K-S test (following the Asgard pipeline), which is an essential step for estimating the sensitivity/resistance effect.
```
drug_info <- read.table(file = "DrugReference/central-nervous-system_drug_info.txt", sep = "\t", header = T, quote = "")
perturbation_matrix_path <- "DrugReference/central-nervous-system_rankMatrix.txt"
drug_sens_res_pval <- compute_drug_sens_res(disease = disease, 
                                            perturbation_matrix_path = perturbation_matrix_path, 
                                            gene_info = gene_info, 
                                            drug_info = drug_info, 
                                            deg_list = deg_list, 
                                            out_path = "results")
```
### 6. Compute final drug therapeutic score
The final drug score is calculated by summing the promotion/inhibition-weighted sensitivity/resistance scores for each cell type, as described by Eq. (7-9) in the scDrugLink paper.
```
gse92742_gctx_path <- "../scDD/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
gse70138_gctx_path <- "../scDD/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx"
drug_score <- compute_scdruglink_score(dat, 
                                      deg_list = deg_list, 
                                      drug_sens_res_pval = drug_sens_res_pval,
                                      tissue = "central nervous system", 
                                      case = disease,
                                      gse92742_gctx_path = gse92742_gctx_path, 
                                      gse70138_gctx_path = gse70138_gctx_path, 
                                      drug_prom_inh_weight = drug_prom_inh_weight,
                                      out_path="results")
head(drug_score)
#             drug_score       p_val        fdr
#biotin      0.000000000 0.773365987 1.00000000
#calcitriol  0.013820636 0.001709220 0.01435745
#calcifediol 0.010978312 0.228458279 0.58153016
#icosapent   0.001226290 0.794338210 1.00000000
#menadione   0.000000000 0.055883363 0.22353345
#pravastatin 0.001372659 0.002719219 0.01928024
```
To repurpose drugs for GBM, rank the drugs by their therapeutic scores, with higher scores indicating a greater likelihood of having an effect on the disease.

### Citation
> Huang, L., Lu, X., Chen, D. scDrugLink: Single-Cell Drug Repurposing for CNS Diseases via Computationally Linking Drug Targets and Perturbation Signatures. 2025.
### Acknowledgement
We thank He et al. for their `Asgard` package (https://github.com/lanagarmire/Asgard), which helps prepare drug perturbation signatures and forms the basis for estimating drug sensitivity/resistance effects.
