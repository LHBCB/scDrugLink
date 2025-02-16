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

## Quick start
