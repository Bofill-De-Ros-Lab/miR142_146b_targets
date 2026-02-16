setwd("~/OneDrive - Aarhus universitet/My_Lab/1_Projects/Collaborations/C_Fillat/")
require(data.table)
#targets_142 <- read_excel("TargetScan7.2__miR-142-5p.predicted_targets.xlsx", sheet = "TargetScan7")
targets_142 <- read_excel("extended_TargetScan8.0__miR-142-5p.Human.predicted_targets.xlsx", sheet = "TargetScanHuman_8.0")

#targets_146 <- read_excel("TargetScan7.2__miR-146-5p.predicted_targets.xlsx", sheet = "TargetScan7")
targets_146 <- read_excel("extended_TargetScan8.0__miR-146-5p.Human.predicted_targets.xlsx", sheet = "TargetScanHuman_8.0")

dataset <- fread(file = "Gene_counts_filtered_normalized_THP-1.csv", header = T)
geneid <- dataset[,c(1,17)]
colnames(geneid) <- c("gene_name","V1")
#gene_id <- sub("\\..*", "", geneid$V1)

#library(biomaRt)
#mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://useast.ensembl.org") 

geneid$T142 <- geneid$gene_name %in% targets_142$`Target gene`
geneid$T146 <- geneid$gene_name %in% targets_146$`Target gene`
geneid$T142 <- gsub("TRUE","T142",geneid$T142 )
geneid$T146 <- gsub("TRUE","T146",geneid$T146 )
geneid$targets <- paste0(geneid$T142,"_",geneid$T146)


setwd("~/OneDrive - Aarhus universitet/My_Lab/1_Projects/Collaborations/C_Fillat/THP-1_results/")
THP_neg_luc <- fread(file = "topTable1_vs_2.csv", header = T)
THP_neg_luc <- merge(THP_neg_luc,geneid)
# plot(THP_neg_luc$logFC,-log10(THP_neg_luc$adj.P.Val))
# abline(h=1.30103)

THP_luc_L142x4 <- fread(file = "topTable2_vs_4.csv", header = T)
THP_luc_L142x4 <- merge(THP_luc_L142x4,geneid)
# plot(THP_luc_L142x4$logFC,-log10(THP_luc_L142x4$adj.P.Val))
# abline(h=1.30103)

# THP_neg_L142x4 <- fread(file = "topTable1_vs_4.csv", header = T)
# THP_neg_L142x4 <- merge(THP_neg_L142x4,geneid)
# plot(THP_neg_L142x4$logFC,-log10(THP_neg_L142x4$adj.P.Val))
# abline(h=1.30103)

THP_luc_L142x2 <- fread(file = "topTable2_vs_3.csv", header = T)
THP_luc_L142x2 <- merge(THP_luc_L142x2,geneid)

THP_luc_L142_146 <- fread(file = "topTable2_vs_5.csv", header = T)
THP_luc_L142_146 <- merge(THP_luc_L142_146,geneid)


library(ggplot2)
to_plot <- THP_neg_luc
to_plot <- THP_luc_L142x2
to_plot <- THP_luc_L142x44
to_plot <- THP_luc_L142_146

to_plot$both <- "FALSE"
i <- 1676
for (i in 1:nrow(to_plot)) {
  test_142 <- to_plot$T142[i] == "T142" 
  test_146 <- to_plot$T146[i] == "T146" 
  test_both <- to_plot$targets[i] == "T142_T146" 
  test_pval <- to_plot$adj.P.Val[i]<0.05
  
  if(test_142==T & test_pval==T & test_both==F){
    to_plot$T142[i] <- -log10(to_plot$adj.P.Val[i])
  }
  if(test_146==T & test_pval==T & test_both==F){
    to_plot$T146[i] <- -log10(to_plot$adj.P.Val[i])
  } 
  if(test_both==T & test_pval==T){
    to_plot$both[i] <- -log10(to_plot$adj.P.Val[i])
  }
  print(i)
  rm(test_142,test_146,test_both, test_pval)
}


# Plot (P.Value OR adj.P.Val)
ggplot(to_plot, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = targets), alpha = 0.7) +
  geom_hline(yintercept = 1.30103, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("FALSE_FALSE" = "grey", "T142_FALSE" = "red","FALSE_T146" = "blue","T142_T146" = "purple" )) +
  labs(title = "Volcano Plot",
       x = "log2 Fold Change",
       y = "-log10(adj. P-value)",
       color = "Predicted targets") +
  theme_minimal()

to_plot <- to_plot[order(to_plot$both, to_plot$T142, to_plot$T146), ]
to_plot$mlog10pval <- -log10(to_plot$adj.P.Val)
to_plot$T142 <- gsub("T142","",to_plot$T142)
to_plot$T146 <- gsub("T146","",to_plot$T146)
to_plot$targets <- gsub("FALSE_T146","",to_plot$targets)
to_plot$targets <- gsub("T142_FALSE","",to_plot$targets)
to_plot$targets <- gsub("T142_T146","",to_plot$targets)

write.table(to_plot, "THP_luc_L142_146.tsv", append = F, row.names = F, sep = "\t")

colnames(targets_142)[1] <- "gene_name"
colnames(targets_146)[1] <- "gene_name"
best_142_4X <- merge(targets_142,THP_luc_L142x4,by="gene_name")
best_142_4X <- best_142_4X[,c(1,3,5:12,23,24)]

best_146 <- merge(targets_146,THP_luc_L142_146,by="gene_name")
best_146 <- best_146[,c(1,3,5:12,23,24)]
best_146 <- best_146[which(best_146$`Conserved sites total`>0 & best_146$logFC>0),]

best_142_2X <- merge(targets_142,THP_luc_L142x2,by="gene_name")
best_142_2X <- best_142_2X[,c(1,3,5:12,23,24)]

best_double <- merge(targets_142,THP_luc_L142_146,by="gene_name")
best_double <- best_double[,c(1,3,5:12,23,24)]

best_142_4X <- best_142_4X[which(best_142_4X$`Conserved sites total`>0 & best_142_4X$logFC>0),]
best_142_2X <- best_142_2X[which(best_142_2X$`Conserved sites total`>0 & best_142_2X$logFC>0),]
best_double <- best_double[which(best_double$`Conserved sites total`>0 & best_double$logFC>0),]
best_142_4X <- best_142_4X[best_142_4X$gene_name %in% best_142_2X$gene_name,]
best_142_4X <- best_142_4X[best_142_4X$gene_name %in% best_double$gene_name,]


best_targets <- best_142_4X
best_targets$sites <- as.numeric(best_targets$`Conserved sites total`) + as.numeric(best_targets$`Poorly conserved sites total`)

ggplot(best_targets, aes(x = AveExpr, y = logFC)) +
  geom_point(aes(color = sites), alpha = 0.7) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "black") +
  labs(title = "Volcano Plot",
       x = "Average Expression",
       y = "log2 Fold Change",
       color = "Predicted targets") +
  theme_minimal()

best_targets <- best_targets[order(-best_targets$sites, -best_targets$AveExpr), ]
write.table(best_targets, paste0("THP_best_142","_all",".tsv"), append = F, row.names = F, sep = "\t")
write.table(best_146, paste0("THP_best_146","_all",".tsv"), append = F, row.names = F, sep = "\t")
