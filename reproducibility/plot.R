
drug_label_df <- read.csv("data/drug_labels_273.csv")

##################################################################################################################
# plot functions
##################################################################################################################
library(ggplot2)
library(pROC)

plot_roc_auc <- function(diz, roc_list){
  
  if (diz == "GBM") {
    full_diz <- "Glioblastoma"
  } else if (diz == "MS") {
    full_diz <- "Multiple Sclerosis"
  } else if (diz == "AD") {
    full_diz <- "Alzheimer's Disease"
  }
  
  auc1 <- auc(roc_list$scDrugLink)
  auc2 <- auc(roc_list$ASGARD)
  auc3 <- auc(roc_list$DrugReSC)
  auc4 <- auc(roc_list$scDrugPrio)
  
  g <- ggroc(roc_list, linewidth = 1.2) +
    geom_abline(slope = 1, intercept = 1,
                linetype = "dashed", color = "gray", linewidth = 1) +
    labs(
      title = NULL,
      x = "Specificity",
      y = "Sensitivity"
    ) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.line = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 18), 
      axis.title = element_text(color = "black", size = 20),
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 1),
      axis.text = element_text(color = "black", size = 20),
      legend.position = "none"
    ) +
    scale_color_manual(values = c(
      "ASGARD" = "#9BC985",  # Teal
      "DrugReSC" = "#F7D58B",  # Orange
      "scDrugPrio" = "#B595BF",  # Purple
      "scDrugLink" = "#797BB7"   # Reddish pink
    ))
  
  g <- g + annotate("text", x = 0.6, y = 0.33, label = paste("ASGARD AUC =", sprintf("%.4f", auc2)), color = "#9BC985", size = 6, fontface = "bold", hjust = 0) +
    annotate("text",x = 0.6, y = 0.27, label = paste("DrugReSC AUC =", sprintf("%.4f", auc3)),color = "#F7D58B",size = 6,fontface = "bold", hjust = 0) + 
    annotate("text",x = 0.6, y = 0.21,label = paste("scDrugPrio AUC =", sprintf("%.4f", auc4)),color = "#B595BF",size = 6,fontface = "bold", hjust = 0) +
    annotate("text", x = 0.6, y = 0.15,label = paste("scDrugLink AUC =", sprintf("%.4f", auc1)),color = "#797BB7",size = 6, fontface = "bold", hjust = 0)
  g$layers <- rev(g$layers)
  ggsave(paste0("Figure 2/", diz, "_AUROC.png"), plot = g, width = 6, height = 6, dpi = 300, units = "in")
  return("Plot done")
}

plot_pr_auc <- function(diz, pr_df_all, aupr_list, baseline_precision){
  
  if (diz == "GBM") {
    full_diz <- "Glioblastoma"
  } else if (diz == "MS") {
    full_diz <- "Multiple Sclerosis"
  } else if (diz == "AD") {
    full_diz <- "Alzheimer's Disease"
  }
  
  g <- ggplot(pr_df_all, aes(x = recall, y = precision, colour = method)) +
    geom_line(size = 1.2) +                     # Draw lines
    geom_hline(
      yintercept = baseline_precision,
      linetype   = "dashed",
      color = "gray", 
      linewidth = 1
    ) +
    theme_minimal() +           # Use a minimal theme
    labs(
      title = NULL,
      x     = "Recall",
      y     = "Precision",
      color = "Method"
    ) + 
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.line = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 18), 
      axis.title = element_text(color = "black", size = 20),
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 1),
      axis.text = element_text(color = "black", size = 20),
      legend.position = "none"
    ) + 
    scale_color_manual(values = c(
      "ASGARD" = "#9BC985",  # Teal
      "DrugReSC" = "#F7D58B",  # Orange
      "scDrugPrio" = "#B595BF",  # Purple
      "scDrugLink" = "#797BB7"   # Reddish pink
    ))
  
  g <- g + annotate("text", x = 0.35, y = 0.83, label = paste("ASGARD AUPR =", sprintf("%.4f", aupr_list[2])), color = "#9BC985", size = 6, fontface = "bold", hjust = 0) +
    annotate("text",x = 0.35, y = 0.77, label = paste("DrugReSC AUPR =", sprintf("%.4f", aupr_list[3])),color = "#F7D58B",size = 6,fontface = "bold", hjust = 0) + 
    annotate("text",x = 0.35, y = 0.71,label = paste("scDrugPrio AUPR =", sprintf("%.4f", aupr_list[4])),color = "#B595BF",size = 6,fontface = "bold", hjust = 0) +
    annotate("text", x = 0.35, y = 0.65,label = paste("scDrugLink AUPR =", sprintf("%.4f", aupr_list[1])),color = "#797BB7",size = 6, fontface = "bold", hjust = 0)
  g$layers <- rev(g$layers)
  ggsave(paste0("Figure 2/", diz, "_PRROC.png"), plot = g, width = 6, height = 6, dpi = 300, units = "in")
  return("Plot done")
}

make_prroc_plot <- function(diz, drug_scores_labels, pr_df_all, aupr_list) {
  if (diz == "AD") {
    labels <- drug_scores_labels$ad_label_final
  } else if (diz == "GBM") {
    labels <- drug_scores_labels$gbm_label_final
  } else if (diz == "MS") {
    labels <- drug_scores_labels$ms_label_final
  }
  baseline_precision <- mean(labels == 1)
  plot_pr_auc(diz, pr_df_all, aupr_list, baseline_precision)
}
##################################################################################################################
# Figure 2: AUC & AUPR comparison of scDrugLink, ASGARD, DrugReSC, scDrugPrio on diseases: GBM, MS, AD
##################################################################################################################
### GBM
diz <- "GBM"
## scDrugLink
drug_scores_df <- read.csv("asgard_scdruglink/GBM_final_scores_df_improved.csv")
drug_scores_df$drug_score <- drug_scores_df$Drug.therapeutic.score
drug_scores_labels <- merge(drug_scores_df[, c("drug_name", "drug_score")], drug_label_df, by.x = "drug_name", by.y = "drug_name", all = FALSE)
eval(drug_scores_labels, diz) # 0.7241731 0.5484190
roc_obj1 <- make_roc(drug_scores_labels, diz)
pr_df1 <- make_aupr_df(drug_scores_labels, diz, "scDrugLink")
aupr1 <- make_aupr(drug_scores_labels, diz, "scDrugLink")
## ASGARD
drug_scores_df <- read.csv("asgard_scdruglink/GBM_final_scores_df.csv")
drug_scores_df$drug_score <- drug_scores_df$Drug.therapeutic.score
drug_scores_labels <- merge(drug_scores_df[, c("drug_name", "drug_score")], drug_label_df, by.x = "drug_name", by.y = "drug_name", all = FALSE)
eval(drug_scores_labels, diz) # 0.7089726 0.5023980
roc_obj2 <- make_roc(drug_scores_labels, diz)
pr_df2 <- make_aupr_df(drug_scores_labels, diz, "ASGARD")
aupr2 <- make_aupr(drug_scores_labels, diz, "ASGARD")
## DrugReSC
drug_scores_df <- readRDS("drugresc_scdruglink/glioblastoma_drug_scores.rds")
drug_scores_labels <- merge(drug_scores_df[, c("drug_name", "drug_score")], drug_label_df, by.x = "drug_name", by.y = "drug_name", all = FALSE)
eval(drug_scores_labels, diz) # 0.5082336 0.2770647
roc_obj3 <- make_roc(drug_scores_labels, diz)
pr_df3 <- make_aupr_df(drug_scores_labels, diz, "DrugReSC")
aupr3 <- make_aupr(drug_scores_labels, diz, "DrugReSC")
## scDrugPrio
drug_scores_df <- read.csv("scdrugprio_scdruglink/gbm/final_drug_scores.csv")
drug_scores_df$drug_score <- drug_scores_df$combined_centrality_score
drug_scores_df$drug_name <- drug_scores_df$Drug
drug_scores_labels <- merge(drug_scores_df[, c("drug_name", "drug_score")], drug_label_df, by.x = "drug_name", by.y = "drug_name", all = FALSE)
eval(drug_scores_labels, diz) # 0.6897959 0.4514044
roc_obj4 <- make_roc(drug_scores_labels, diz)
pr_df4 <- make_aupr_df(drug_scores_labels, diz, "scDrugPrio")
aupr4 <- make_aupr(drug_scores_labels, diz, "scDrugPrio")

roc_list <- list(
  "ASGARD" = roc_obj2,
  "DrugReSC" = roc_obj3,
  "scDrugPrio" = roc_obj4,
  "scDrugLink" = roc_obj1
)
plot_roc_auc(diz, roc_list)

pr_df_all <- rbind(pr_df1, pr_df2, pr_df3, pr_df4)
aupr_list <- c(aupr1, aupr2, aupr3, aupr4)
make_prroc_plot(diz, drug_scores_labels, pr_df_all, aupr_list)

### MS
diz <- "MS"
## scDrugLink
drug_scores_df <- read.csv("asgard_scdruglink/True_final_scores_df_improved.csv")
drug_scores_df$drug_score <- drug_scores_df$Drug.therapeutic.score
drug_scores_labels <- merge(drug_scores_df[, c("drug_name", "drug_score")], drug_label_df, by.x = "drug_name", by.y = "drug_name", all = FALSE)
eval(drug_scores_labels, diz) # 0.6285957 0.3412467
roc_obj1 <- make_roc(drug_scores_labels, diz)
pr_df1 <- make_aupr_df(drug_scores_labels, diz, "scDrugLink")
aupr1 <- make_aupr(drug_scores_labels, diz, "scDrugLink")
## ASGARD
drug_scores_df <- read.csv("asgard_scdruglink/True_final_scores_df.csv")
drug_scores_df$drug_score <- drug_scores_df$Drug.therapeutic.score
drug_scores_labels <- merge(drug_scores_df[, c("drug_name", "drug_score")], drug_label_df, by.x = "drug_name", by.y = "drug_name", all = FALSE)
eval(drug_scores_labels, diz) # 0.6146768 0.3071688
roc_obj2 <- make_roc(drug_scores_labels, diz)
pr_df2 <- make_aupr_df(drug_scores_labels, diz, "ASGARD")
aupr2 <- make_aupr(drug_scores_labels, diz, "ASGARD")
## DrugReSC
drug_scores_df <- readRDS("drugresc_scdruglink/multiple sclerosis_drug_scores.rds")
drug_scores_labels <- merge(drug_scores_df[, c("drug_name", "drug_score")], drug_label_df, by.x = "drug_name", by.y = "drug_name", all = FALSE)
eval(drug_scores_labels, diz) # 0.5631766 0.2666879
roc_obj3 <- make_roc(drug_scores_labels, diz)
pr_df3 <- make_aupr_df(drug_scores_labels, diz, "DrugReSC")
aupr3 <- make_aupr(drug_scores_labels, diz, "DrugReSC")
## scDrugPrio
drug_scores_df <- read.csv("scdrugprio_scdruglink/ms/final_drug_scores.csv")
drug_scores_df$drug_score <- drug_scores_df$combined_centrality_score
drug_scores_df$drug_name <- drug_scores_df$Drug
drug_scores_labels <- merge(drug_scores_df[, c("drug_name", "drug_score")], drug_label_df, by.x = "drug_name", by.y = "drug_name", all = FALSE)
eval(drug_scores_labels, diz) # 0.5594649 0.2574637
roc_obj4 <- make_roc(drug_scores_labels, diz)
pr_df4 <- make_aupr_df(drug_scores_labels, diz, "scDrugPrio")
aupr4 <- make_aupr(drug_scores_labels, diz, "scDrugPrio")

roc_list <- list(
  "ASGARD" = roc_obj2,
  "DrugReSC" = roc_obj3,
  "scDrugPrio" = roc_obj4,
  "scDrugLink" = roc_obj1
)
plot_roc_auc(diz, roc_list)

pr_df_all <- rbind(pr_df1, pr_df2, pr_df3, pr_df4)
aupr_list <- c(aupr1, aupr2, aupr3, aupr4)
make_prroc_plot(diz, drug_scores_labels, pr_df_all, aupr_list)

### AD
diz <- "AD"
## scDrugLink
drug_scores_df <- read.csv("asgard_scdruglink/Alzheimer disease_final_scores_df_improved.csv")
drug_scores_df$drug_score <- drug_scores_df$Drug.therapeutic.score
drug_scores_labels <- merge(drug_scores_df[, c("drug_name", "drug_score")], drug_label_df, by.x = "drug_name", by.y = "drug_name", all = FALSE)
eval(drug_scores_labels, diz) # 0.6374676 0.3887734
roc_obj1 <- make_roc(drug_scores_labels, diz)
pr_df1 <- make_aupr_df(drug_scores_labels, diz, "scDrugLink")
aupr1 <- make_aupr(drug_scores_labels, diz, "scDrugLink")
## ASGARD
drug_scores_df <- read.csv("asgard_scdruglink/Alzheimer disease_final_scores_df.csv")
drug_scores_df$drug_score <- drug_scores_df$Drug.therapeutic.score
drug_scores_labels <- merge(drug_scores_df[, c("drug_name", "drug_score")], drug_label_df, by.x = "drug_name", by.y = "drug_name", all = FALSE)
eval(drug_scores_labels, diz) # 0.6279469 0.3759152
roc_obj2 <- make_roc(drug_scores_labels, diz)
pr_df2 <- make_aupr_df(drug_scores_labels, diz, "ASGARD")
aupr2 <- make_aupr(drug_scores_labels, diz, "ASGARD")
## DrugReSC
drug_scores_df <- readRDS("drugresc_scdruglink/Alzheimer disease_drug_scores.rds")
drug_scores_labels <- merge(drug_scores_df[, c("drug_name", "drug_score")], drug_label_df, by.x = "drug_name", by.y = "drug_name", all = FALSE)
eval(drug_scores_labels, diz) # 0.5154145 0.3063060
roc_obj3 <- make_roc(drug_scores_labels, diz)
pr_df3 <- make_aupr_df(drug_scores_labels, diz, "DrugReSC")
aupr3 <- make_aupr(drug_scores_labels, diz, "DrugReSC")
## scDrugPrio
drug_scores_df <- read.csv("scdrugprio_scdruglink/ad/final_drug_scores.csv")
drug_scores_df$drug_score <- drug_scores_df$combined_centrality_score
drug_scores_df$drug_name <- drug_scores_df$Drug
drug_scores_labels <- merge(drug_scores_df[, c("drug_name", "drug_score")], drug_label_df, by.x = "drug_name", by.y = "drug_name", all = FALSE)
eval(drug_scores_labels, diz) # 0.5934262 0.3819951
roc_obj4 <- make_roc(drug_scores_labels, diz)
pr_df4 <- make_aupr_df(drug_scores_labels, diz, "scDrugPrio")
aupr4 <- make_aupr(drug_scores_labels, diz, "scDrugPrio")

roc_list <- list(
  "ASGARD" = roc_obj2,
  "DrugReSC" = roc_obj3,
  "scDrugPrio" = roc_obj4,
  "scDrugLink" = roc_obj1
)
plot_roc_auc(diz, roc_list)

pr_df_all <- rbind(pr_df1, pr_df2, pr_df3, pr_df4)
aupr_list <- c(aupr1, aupr2, aupr3, aupr4)
make_prroc_plot(diz, drug_scores_labels, pr_df_all, aupr_list)

##################################################################################################################
# Figure 3: GBM comparison by cluster: scDrugLink, ASGARD, scDrugPrio, with table on average+-std of AUC & AUPR
##################################################################################################################

eval_cluster <- function(drug_scores_labels, cluster, diz){
  drug_score_list <- drug_scores_labels[, cluster]
  if (diz == "AD") {
    labels <- drug_scores_labels$ad_label_final
  } else if (diz == "GBM") {
    labels <- drug_scores_labels$gbm_label_final
  } else if (diz == "MS") {
    labels <- drug_scores_labels$ms_label_final
  }
  pred <- prediction(predictions = drug_score_list, labels = labels)
  perf <- performance(pred, measure = "tpr", x.measure = "fpr")
  #plot(perf, col = "red", lwd = 2, main = "ROC Curve",
  #     xlab = "False Positive Rate", ylab = "True Positive Rate")
  auc_perf <- performance(pred, measure = "auc")
  auc_value <- auc_perf@y.values[[1]]
  #legend("bottomright", legend = paste("AUC =", round(auc_value, 3)), 
  #       col = "red", lty = 1, bty = "n")
  pos_scores <- drug_score_list[labels == 1]
  neg_scores <- drug_score_list[labels == 0]
  pr_result <- pr.curve(scores.class0 = pos_scores, scores.class1 = neg_scores, curve = FALSE)
  aupr_value <- pr_result$auc.integral
  return(c(auc_value, aupr_value))
}

library(dplyr)
library(stringr)
eval_method_by_cluster <- function(drug_scores_labels, diz, method) {
  aucs <- c()
  auprs <- c()
  clusters <- c()
  for (cluster in colnames(drug_scores_df)) {
    #print(cluster)
    clusters <- c(clusters, cluster)
    eval_values <- eval_cluster(drug_scores_labels, cluster, diz)
    #print(eval_values)
    aucs <- c(aucs, eval_values[1])
    auprs <- c(auprs, eval_values[2])
  }
  #clusters <- c(clusters, "Average")
  #aucs <- c(aucs, mean(aucs))
  #auprs <- c(auprs, mean(auprs))
  metric_df <- as.data.frame(list(cell_type = clusters, method = rep(method, length(clusters)), AUC = aucs, AUPR = auprs))
  #return(c(mean(aucs), sd(aucs), mean(auprs), sd(auprs)))
  metric_df <- metric_df %>%
    mutate(
      cell_type = cell_type %>%
        str_replace_all("\\s+", "") %>%  # remove all spaces
        str_remove("^Cluster_") %>%      # remove prefix
        str_remove("\\.$")               # remove trailing period
    )
  return(metric_df)
}

plot_auc_bar <- function(df_all, bar_width) {
  df_auc <- df_all[, c("cell_type", "method", "AUC")]
  df_auc$AUC <- as.numeric(df_all$AUC)
  df_avg <- df_auc %>% group_by(method) %>% summarise(AUC = mean(AUC)) %>% mutate(cell_type = "Average")
  df_auc <- bind_rows(df_auc, df_avg)
  df_auc <- df_auc %>% mutate(method = factor(method, levels = c("ASGARD", "scDrugPrio", "scDrugLink")))
  
  df_notavg <- df_auc %>% filter(cell_type != "Average")
  df_sd <- df_notavg %>% group_by(method) %>% summarize(sd_value = sd(AUC, na.rm = TRUE))
  df_avg <- df_auc %>% filter(cell_type == "Average") %>% left_join(df_sd, by = "method") %>% rename(sd_AUC = sd_value)
  
  nonavg_types   <- setdiff(unique(df_auc$cell_type), "Average")
  df_auc <- df_auc %>% mutate(cell_type = factor(cell_type, levels = c(sort(nonavg_types), "Average")))
  
  g <- ggplot(df_auc, aes(x = cell_type, y = AUC, fill = method, group = method)) +
    geom_bar(stat = "identity", position = position_dodge(width = bar_width), width = bar_width) +
    labs(
      title = NULL,
      x = NULL,
      y = "AUC",
      fill = "Method"
    ) +
    scale_fill_manual(name = "Methods:", values = c(
      "scDrugLink" = "#797BB7",  # Reddish pink
      "ASGARD" = "#9BC985",  # Teal
      "scDrugPrio" = "#B595BF"  # Purple
    )) +
    geom_errorbar(
      data = df_avg,  # only the rows for Average
      aes(
        x    = cell_type,
        ymin = AUC - sd_AUC,
        ymax = AUC + sd_AUC
      ),
      width    = max(c(bar_width-0.3, 0.2)),
      position = position_dodge(width = bar_width),
      linewidth = 1
    ) +
    scale_x_discrete(labels = function(x) {
      sapply(x, function(lbl) {
        if (lbl == "Average") {
          bquote(bold(.(lbl)))
        } else {
          bquote(.(lbl))
        }
      })
    }) +
    geom_vline(xintercept = length(nonavg_types) + 0.5, linetype = "dashed", color = "gray", linewidth = 1.2) +
    theme_minimal() +
    theme(
      #plot.title = element_text(hjust = 0.5, size = 18), 
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 1),
      axis.text = element_text(color = "black", size = 14),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.title = element_text(color = "black", size = 18),
      axis.line = element_blank(),
      legend.position = "none",
      #legend.position = "bottom",
      #legend.text     = element_text(size = 18),
      #legend.key.width  = unit(1.5, "cm"),
      #legend.title = element_text(size = 18, margin = margin(r = 20))
    ) +
    scale_y_continuous(breaks=seq(0.25, 0.75, by=0.25), limits = c(0., 0.75), expand = c(0, 0))
  
  file_name <- paste0("Figure 3/", diz, "_AUC.png")
  ggsave(file_name, plot = g, width = 20, height = 4, dpi = 300, units = "in")
}

plot_aupr_bar <- function(df_all, bar_width) {
  df_AUPR <- df_all[, c("cell_type", "method", "AUPR")]
  df_AUPR$AUPR <- as.numeric(df_all$AUPR)
  df_avg <- df_AUPR %>% group_by(method) %>% summarise(AUPR = mean(AUPR)) %>% mutate(cell_type = "Average")
  df_AUPR <- bind_rows(df_AUPR, df_avg)
  df_AUPR <- df_AUPR %>% mutate(method = factor(method, levels = c("ASGARD", "scDrugPrio", "scDrugLink")))
  
  df_notavg <- df_AUPR %>% filter(cell_type != "Average")
  df_sd <- df_notavg %>% group_by(method) %>% summarize(sd_value = sd(AUPR, na.rm = TRUE))
  df_avg <- df_AUPR %>% filter(cell_type == "Average") %>% left_join(df_sd, by = "method") %>% rename(sd_AUPR = sd_value)
  
  nonavg_types   <- setdiff(unique(df_AUPR$cell_type), "Average")
  df_AUPR <- df_AUPR %>% mutate(cell_type = factor(cell_type, levels = c(sort(nonavg_types), "Average")))
  
  g <- ggplot(df_AUPR, aes(x = cell_type, y = AUPR, fill = method, group = method)) +
    geom_bar(stat = "identity", position = position_dodge(width = bar_width), width = bar_width) +
    labs(
      title = NULL,
      x = NULL,
      y = "AUPR",
      fill = "Method"
    ) +
    scale_fill_manual(name = "Methods:", values = c(
      "scDrugLink" = "#797BB7",  # Reddish pink
      "ASGARD" = "#9BC985",  # Teal
      "scDrugPrio" = "#B595BF"  # Purple
    )) +
    geom_errorbar(
      data = df_avg,  # only the rows for Average
      aes(
        x    = cell_type,
        ymin = AUPR - sd_AUPR,
        ymax = AUPR + sd_AUPR
      ),
      width    = max(c(bar_width-0.3, 0.2)),
      position = position_dodge(width = bar_width),
      linewidth = 1
    ) +
    scale_x_discrete(labels = function(x) {
      sapply(x, function(lbl) {
        if (lbl == "Average") {
          bquote(bold(.(lbl)))
        } else {
          bquote(.(lbl))
        }
      })
    }) +
    geom_vline(xintercept = length(nonavg_types) + 0.5, linetype = "dashed", color = "gray", linewidth = 1.2) +
    theme_minimal() +
    theme(
      #plot.title = element_text(hjust = 0.5, size = 18), 
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 1),
      axis.text = element_text(color = "black", size = 14),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.title = element_text(color = "black", size = 18),
      axis.line = element_blank(),
      legend.position = "none",
      #legend.position = "bottom",
      #legend.text     = element_text(size = 18),
      #legend.key.width  = unit(1.5, "cm"),
      #legend.title = element_text(size = 18, margin = margin(r = 20))
    ) +
    scale_y_continuous(breaks=seq(0.25, 0.75, by=0.25), limits = c(0., 0.75), expand = c(0, 0))
  
  file_name <- paste0("Figure 3/", diz, "_AUPR.png")
  ggsave(file_name, plot = g, width = 20, height = 4, dpi = 300, units = "in")
}

plot_aupr_bar_with_legend <- function(df_all, bar_width) {
  df_AUPR <- df_all[, c("cell_type", "method", "AUPR")]
  df_AUPR$AUPR <- as.numeric(df_all$AUPR)
  df_avg <- df_AUPR %>% group_by(method) %>% summarise(AUPR = mean(AUPR)) %>% mutate(cell_type = "Average")
  df_AUPR <- bind_rows(df_AUPR, df_avg)
  df_AUPR <- df_AUPR %>% mutate(method = factor(method, levels = c("ASGARD", "scDrugPrio", "scDrugLink")))
  
  df_notavg <- df_AUPR %>% filter(cell_type != "Average")
  df_sd <- df_notavg %>% group_by(method) %>% summarize(sd_value = sd(AUPR, na.rm = TRUE))
  df_avg <- df_AUPR %>% filter(cell_type == "Average") %>% left_join(df_sd, by = "method") %>% rename(sd_AUPR = sd_value)
  
  nonavg_types   <- setdiff(unique(df_AUPR$cell_type), "Average")
  df_AUPR <- df_AUPR %>% mutate(cell_type = factor(cell_type, levels = c(sort(nonavg_types), "Average")))
  
  g <- ggplot(df_AUPR, aes(x = cell_type, y = AUPR, fill = method, group = method)) +
    geom_bar(stat = "identity", position = position_dodge(width = bar_width), width = bar_width) +
    labs(
      title = NULL,
      x = NULL,
      y = "AUPR",
      fill = "Method"
    ) +
    scale_fill_manual(name = "Methods:", values = c(
      "scDrugLink" = "#797BB7",  # Reddish pink
      "ASGARD" = "#9BC985",  # Teal
      "scDrugPrio" = "#B595BF"  # Purple
    )) +
    geom_errorbar(
      data = df_avg,  # only the rows for Average
      aes(
        x    = cell_type,
        ymin = AUPR - sd_AUPR,
        ymax = AUPR + sd_AUPR
      ),
      width    = max(c(bar_width-0.3, 0.2)),
      position = position_dodge(width = bar_width),
      linewidth = 1
    ) +
    scale_x_discrete(labels = function(x) {
      sapply(x, function(lbl) {
        if (lbl == "Average") {
          bquote(bold(.(lbl)))
        } else {
          bquote(.(lbl))
        }
      })
    }) +
    geom_vline(xintercept = length(nonavg_types) + 0.5, linetype = "dashed", color = "gray", linewidth = 1.2) +
    theme_minimal() +
    theme(
      #plot.title = element_text(hjust = 0.5, size = 18), 
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 1),
      axis.text = element_text(color = "black", size = 14),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.title = element_text(color = "black", size = 18),
      axis.line = element_blank(),
      #legend.position = "none",
      legend.position = "bottom",
      legend.text     = element_text(size = 18),
      legend.key.width  = unit(1.5, "cm"),
      legend.title = element_text(size = 18, margin = margin(r = 20))
    ) +
    scale_y_continuous(breaks=seq(0.25, 0.75, by=0.25), limits = c(0., 0.75), expand = c(0, 0))
  
  file_name <- paste0("Figure 3/", diz, "_AUPR_with_legend.png")
  ggsave(file_name, plot = g, width = 20, height = 4, dpi = 300, units = "in")
}

diz <- "GBM"
### scDrugLink
drug_scores_df <- read.csv("asgard_scdruglink_clusters/gbm/weighted_drug_cluster_scores.csv", row.names = 1)
drug_scores_df <- t(as.matrix(drug_scores_df))
rownames(drug_scores_df) <- gsub("\\.", "-", rownames(drug_scores_df))
drug_scores_labels <- merge(drug_scores_df, drug_label_df, by.x = "row.names", by.y = "drug_name_lower", all = FALSE)
df1 <- eval_method_by_cluster(drug_scores_labels, "GBM", "scDrugLink") # 0.66445913 0.05678218 0.41912587 0.05588805

### ASGARD
drug_scores_df <- read.csv("asgard_scdruglink_clusters/gbm/drug_cluster_scores.csv", row.names = 1)
drug_scores_df <- t(as.matrix(drug_scores_df))
rownames(drug_scores_df) <- gsub("\\.", "-", rownames(drug_scores_df))
drug_scores_labels <- merge(drug_scores_df, drug_label_df, by.x = "row.names", by.y = "drug_name_lower", all = FALSE)
df2 <- eval_method_by_cluster(drug_scores_labels, "GBM", "ASGARD") # 0.65535086 0.05373208 0.38954649 0.05102038

### scDrugPrio
drug_scores_df <- read.csv("scdrugprio_clusters/gbm_cluster_scores.csv", row.names = 1)
drug_scores_df <- as.matrix(drug_scores_df)
rownames(drug_scores_df) <- gsub("\\.", "-", rownames(drug_scores_df))
drug_scores_labels <- merge(drug_scores_df, drug_label_df, by.x = "row.names", by.y = "drug_name", all = FALSE)
df3 <- eval_method_by_cluster(drug_scores_labels, "GBM", "scDrugPrio") # 0.55667461 0.06083130 0.32655631 0.06145681

df_all <- rbind(df1, df2, df3)

unique_cell_types <- unique(df_all$cell_type)
sub_df <- df_all[df_all$method == "scDrugPrio", ]
# assign 0 to NA cell types in scDrugPrio
for (cell_type in unique_cell_types) {
  if (!(cell_type %in% sub_df$cell_type)) {
    df_all[nrow(df_all)+1, ] <- c(cell_type, "scDrugPrio", 0, 0)
  }
}

plot_auc_bar(df_all, 0.7)
plot_aupr_bar(df_all, 0.7)

##################################################################################################################
# Figure 4: MS comparison by cluster: scDrugLink, ASGARD, scDrugPrio, with table on average+-std of AUC & AUPR
##################################################################################################################
diz <- "MS"
### scDrugLink
drug_scores_df <- read.csv("asgard_scdruglink_clusters/ms/weighted_drug_cluster_scores.csv", row.names = 1)
drug_scores_df <- t(as.matrix(drug_scores_df))
rownames(drug_scores_df) <- gsub("\\.", "-", rownames(drug_scores_df))
drug_scores_labels <- merge(drug_scores_df, drug_label_df, by.x = "row.names", by.y = "drug_name_lower", all = FALSE)
df1 <- eval_method_by_cluster(drug_scores_labels, "MS", "scDrugLink") # 0.56806372 0.05187308 0.27916325 0.04768145

### ASGARD
drug_scores_df <- read.csv("asgard_scdruglink_clusters/ms/drug_cluster_scores.csv", row.names = 1)
drug_scores_df <- t(as.matrix(drug_scores_df))
rownames(drug_scores_df) <- gsub("\\.", "-", rownames(drug_scores_df))
drug_scores_labels <- merge(drug_scores_df, drug_label_df, by.x = "row.names", by.y = "drug_name_lower", all = FALSE)
df2 <- eval_method_by_cluster(drug_scores_labels, "MS", "ASGARD") # 0.56613311 0.05031224 0.26990665 0.03508120

### scDrugPrio
drug_scores_df <- read.csv("scdrugprio_clusters/ms_cluster_scores.csv", row.names = 1)
drug_scores_df <- as.matrix(drug_scores_df)
rownames(drug_scores_df) <- gsub("\\.", "-", rownames(drug_scores_df))
drug_scores_labels <- merge(drug_scores_df, drug_label_df, by.x = "row.names", by.y = "drug_name", all = FALSE)
df3 <- eval_method_by_cluster(drug_scores_labels, "MS", "scDrugPrio") # 0.51021616 0.02354195 0.23097959 0.02110718

df_all <- rbind(df1, df2, df3)

unique_cell_types <- unique(df_all$cell_type)
sub_df <- df_all[df_all$method == "scDrugPrio", ]
# assign 0 to NA cell types in scDrugPrio
for (cell_type in unique_cell_types) {
  if (!(cell_type %in% sub_df$cell_type)) {
    df_all[nrow(df_all)+1, ] <- c(cell_type, "scDrugPrio", 0, 0)
  }
}

plot_auc_bar(df_all, 0.7)
plot_aupr_bar(df_all, 0.7)

##################################################################################################################
# Figure 5: AD comparison by cluster: scDrugLink, ASGARD, scDrugPrio, with table on average+-std of AUC & AUPR
##################################################################################################################
diz <- "AD"
### scDrugLink
drug_scores_df <- read.csv("asgard_scdruglink_clusters/ad/weighted_drug_cluster_scores.csv", row.names = 1)
drug_scores_df <- t(as.matrix(drug_scores_df))
rownames(drug_scores_df) <- gsub("\\.", "-", rownames(drug_scores_df))
drug_scores_labels <- merge(drug_scores_df, drug_label_df, by.x = "row.names", by.y = "drug_name_lower", all = FALSE)
df1 <- eval_method_by_cluster(drug_scores_labels, "AD", "scDrugLink") # 0.60997814 0.03095573 0.36920953 0.02155681

### ASGARD
drug_scores_df <- read.csv("asgard_scdruglink_clusters/ad/drug_cluster_scores.csv", row.names = 1)
drug_scores_df <- t(as.matrix(drug_scores_df))
rownames(drug_scores_df) <- gsub("\\.", "-", rownames(drug_scores_df))
drug_scores_labels <- merge(drug_scores_df, drug_label_df, by.x = "row.names", by.y = "drug_name_lower", all = FALSE)
df2 <- eval_method_by_cluster(drug_scores_labels, "AD", "ASGARD") # 0.60965431 0.03082130 0.37336473 0.02286845

### scDrugPrio
drug_scores_df <- read.csv("scdrugprio_clusters/ad_cluster_scores.csv", row.names = 1)
drug_scores_df <- as.matrix(drug_scores_df)
rownames(drug_scores_df) <- gsub("\\.", "-", rownames(drug_scores_df))
drug_scores_labels <- merge(drug_scores_df, drug_label_df, by.x = "row.names", by.y = "drug_name", all = FALSE)
df3 <- eval_method_by_cluster(drug_scores_labels, "AD", "scDrugPrio") # 0.53348041 0.04074165 0.32072096 0.03812069

df_all <- rbind(df1, df2, df3)
df_all$cell_type <- gsub("\\.cluster\\.", "", df_all$cell_type)
df_all$cell_type <- gsub("cluster", "", df_all$cell_type)

unique_cell_types <- unique(df_all$cell_type)
sub_df <- df_all[df_all$method == "scDrugPrio", ]
# assign 0 to NA cell types in scDrugPrio
for (cell_type in unique_cell_types) {
  if (!(cell_type %in% sub_df$cell_type)) {
    df_all[nrow(df_all)+1, ] <- c(cell_type, "scDrugPrio", 0, 0)
  }
}

plot_auc_bar(df_all, bar_width = 0.4)
plot_aupr_bar(df_all, 0.4)
plot_aupr_bar_with_legend(df_all, 0.4)
