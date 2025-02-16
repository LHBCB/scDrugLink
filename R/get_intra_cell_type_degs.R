#' @title Identify within cell-type differentially expressed genes (DEGs) between diseased and control cells.
#' @description This function processes each cell type individually, performs statistical testing to find markers, and stores the results (log fold change and adjusted p-values) for each cell type in a list.
#' @param dat A Seurat object for a log-normalised, integrated scRNA-seq disease atlas.
#' @param disease The disease term in the atlas for identifying diseased cells.
#' @param control The control term in the atlas for identifying control cells.
#' @return A list of DEGs for each cell type.
#' @export
#' @import Seurat


get_intra_cell_type_degs <- function(dat, disease, control) {
  DefaultAssay(dat) <- "RNA"

  # identify DEGs
  deg_list <- list()
  cell_type_list <- NULL
  for(i in unique(dat@meta.data$cell_type)){

    print(paste0("Identify DEGs for cell_type/cluster: ", i))

    Idents(dat) <- "cell_type"
    ct_cells <- subset(dat, cell_type == i)
    Idents(ct_cells) <- "disease"

    cells_1 <- WhichCells(ct_cells, ident = disease)
    cells_2 <- WhichCells(ct_cells, ident = control)
    if (length(cells_1) < 3 || length(cells_2) < 3) {
        message("Skipping cell_type/cluster: ", i, " due to insufficient cell count.")
        next
      }

    ct_markers <- FindMarkers(ct_cells, ident.1 = disease, ident.2 = control)
    
    ct_markers_for_drug <- data.frame(row.names=row.names(ct_markers), score=ct_markers$avg_log2FC, adj.P.Val=ct_markers$p_val_adj, P.Value=ct_markers$p_val)
    deg_list[[i]] <- ct_markers_for_drug

    cell_type_list <- c(cell_type_list,i)
  }
  names(deg_list) <- cell_type_list
  return(deg_list)
}