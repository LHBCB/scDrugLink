library(DrugReSC)

set.seed(42)

drug_label_df <- read.csv("data/drug_labels_273.csv")

dat <- readRDS("gbm_final.rds")
cell_labels <- ifelse(dat@meta.data$histology=="GBM", 1, 0) # label

drug_signature_list <- readRDS("drug_sig_list_273.rds")

drug_score <- DrugReSC(sc_data, drug_signature_list, cell_labels, improtance_method = "rf", ncore = 8)
drug_score_df <- data.frame(drug_name = names(drug_signature_list), drug_score = drug_score, stringsAsFactors = FALSE)

write.csv(df, "drugresc_results/results.csv", row.names = F)
