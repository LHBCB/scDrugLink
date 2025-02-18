library(pROC)
library(ROCR)
library(PRROC)
library(scDrugLink)

eval <- function(drug_score){
  
  if (disease == "AD") {
    labels <- drug_score$ad_label_final
  } else if (disease == "GBM") {
    labels <- drug_score$gbm_label_final
  } else if (disease == "MS") {
    labels <- drug_score$ms_label_final
  }

  pred <- prediction(predictions = drug_score$drug_score, labels = labels)
  perf <- performance(pred, measure = "tpr", x.measure = "fpr")
  
  auc_perf <- performance(pred, measure = "auc")
  auc_value <- auc_perf@y.values[[1]]
  
  pos_scores <- drug_score$drug_score[labels == 1]
  neg_scores <- drug_score$drug_score[labels == 0]
  pr_result <- pr.curve(scores.class0 = pos_scores, scores.class1 = neg_scores, curve = FALSE)
  aupr_value <- pr_result$auc.integral

  return(c(auc_value, aupr_value))
}

eval_disease <- function(disease, drug_score) {
  drug_score_sub <- merge(cns_drug_labels, drug_score, by.x = "drug_name_lower", by.y = "row.names", all = FALSE)
  eval_values <- eval(drug_score_sub)
  return(eval_values)
}

results <- eval_disease(disease, drug_score)

df <- data.frame(
  disease = disease,
  auc     = results[1],
  aupr    = results[2],
  stringsAsFactors = FALSE
)
write.csv(df, paste0("results/", disease, "_metrics.csv"), row.names = F)
