setwd("~/OneDrive - Aarhus universitet/My_Lab/1_Projects/Collaborations/C_Fillat/")
require(data.table)
dataset <- fread(file = "Gene_counts_filtered_normalized_HL-60.csv", header = T)

dataset <- dataset[which(dataset$H_C_1>0&dataset$H_C_2>0&dataset$H_C_4>0),]
dataset <- dataset[which(dataset$H_1_2>0&dataset$H_1_3>0&dataset$H_1_4>0),]
dataset <- dataset[which(dataset$H_2_1>0&dataset$H_2_2>0&dataset$H_2_4>0),]
dataset <- dataset[which(dataset$H_3_2>0&dataset$H_3_3>0&dataset$H_3_4>0),]
dataset <- dataset[which(dataset$H_4_1>0&dataset$H_4_2>0&dataset$H_4_3>0),]

median <- (median(dataset$H_C_1)+median(dataset$H_C_2)+median(dataset$H_C_4)
           +median(dataset$H_1_2)+median(dataset$H_1_3)+median(dataset$H_1_4)
           +median(dataset$H_2_1)+median(dataset$H_2_2)+median(dataset$H_2_4)
           +median(dataset$H_3_2)+median(dataset$H_3_3)+median(dataset$H_3_4)
           +median(dataset$H_4_1)+median(dataset$H_4_2)+median(dataset$H_4_3))/15

dataset$H_C_1 <- dataset$H_C_1-(median(dataset$H_C_1)-median)
dataset$H_C_2 <- dataset$H_C_2-(median(dataset$H_C_2)-median)
dataset$H_C_4 <- dataset$H_C_4-(median(dataset$H_C_4)-median)

dataset$H_1_2 <- dataset$H_1_2-(median(dataset$H_1_2)-median)
dataset$H_1_3 <- dataset$H_1_3-(median(dataset$H_1_3)-median)
dataset$H_1_4 <- dataset$H_1_4-(median(dataset$H_1_4)-median)

dataset$H_2_1 <- dataset$H_2_1-(median(dataset$H_2_1)-median)
dataset$H_2_2 <- dataset$H_2_2-(median(dataset$H_2_2)-median)
dataset$H_2_4 <- dataset$H_2_4-(median(dataset$H_2_4)-median)

dataset$H_3_2 <- dataset$H_3_2-(median(dataset$H_3_2)-median)
dataset$H_3_3 <- dataset$H_3_3-(median(dataset$H_3_3)-median)
dataset$H_3_4 <- dataset$H_3_4-(median(dataset$H_3_4)-median)

dataset$H_4_1 <- dataset$H_4_1-(median(dataset$H_4_1)-median)
dataset$H_4_2 <- dataset$H_4_2-(median(dataset$H_4_2)-median)
dataset$H_4_3 <- dataset$H_4_3-(median(dataset$H_4_3)-median)

dataset$neg_cont <- (dataset$H_C_1+dataset$H_C_2+dataset$H_C_4)/3
dataset$luc_cont <- (dataset$H_1_2+dataset$H_1_3+dataset$H_1_4)/3
dataset$m142T2X <- (dataset$H_2_1+dataset$H_2_2+dataset$H_2_4)/3
dataset$m142T4X <- (dataset$H_3_2+dataset$H_3_3+dataset$H_3_4)/3
dataset$m142146 <- (dataset$H_4_1+dataset$H_4_2+dataset$H_4_3)/3

dataset <- dataset[which(dataset$neg_cont>0&dataset$luc_cont>0&dataset$m142T2X>0&dataset$m142T4X>0&dataset$m142146>0),]

plot(log10(dataset$neg_cont),log10(dataset$luc_cont))
cor(dataset$neg_cont,dataset$luc_cont)

plot(log10(dataset$luc_cont),log10(dataset$m142T2X))
#cor(dataset$luc_cont,dataset$m142T2X)
cor(dataset$neg_cont,dataset$m142T2X)

plot(log10(dataset$luc_cont),log10(dataset$m142T4X))
#cor(dataset$luc_cont,dataset$m142T4X)
cor(dataset$neg_cont,dataset$m142T4X)

plot(log10(dataset$luc_cont),log10(dataset$m142146))
#cor(dataset$luc_cont,dataset$m142146)
cor(dataset$neg_cont,dataset$m142146)

dataset$gene_name <- toupper(dataset$gene_name)

#targets_142 <- fread(file = "TargetScan7.2__miR-142-5p.predicted_targets.txt", header = T)
#library(readxl)
targets_142 <- read_excel("TargetScan7.2__miR-142-5p.predicted_targets.xlsx", sheet = "TargetScan7")
targets_146 <- read_excel("TargetScan7.2__miR-146-5p.predicted_targets.xlsx", sheet = "TargetScan7")
colnames(targets_142)[1] <- "gene_name"
colnames(targets_146)[1] <- "gene_name"

summary(dataset$luc_cont)
#cuttoff <- 598.9
#dataset <- dataset[which(dataset$luc_cont>=cuttoff),]
dataset$log2FC_142X2 <- log2(dataset$m142T2X/dataset$luc_cont)
dataset$log2FC_142X4 <- log2(dataset$m142T4X/dataset$luc_cont)
dataset$log2FC_both <- log2(dataset$m142146/dataset$luc_cont)

dataset$log2FC_selected <- dataset$log2FC_both

dataset_142 <- dataset[dataset$gene_name %in% targets_142$gene_name,]
dataset_146 <- dataset[dataset$gene_name %in% targets_146$gene_name,]

dataset_142 <- merge(dataset_142,targets_142)
dataset_146 <- merge(dataset_146,targets_146)

#not in function
`%notin%` <- Negate(`%in%`)
dataset_NT <- dataset[dataset$gene_name %notin% dataset_142$gene_name,]
dataset_NT <- dataset_NT[dataset_NT$gene_name %notin% dataset_146$gene_name,]

breaks <-  seq(-10, 10, by=0.1)

## Baseline
par(mfrow=c(1,1))

#cumulative tables negetive
set_NT <-  as.numeric(dataset_NT$log2FC_selected) 
set_NT_cumulative = c(0, cumsum(table(cut(set_NT, breaks, right=TRUE))/nrow(dataset_NT))) 

plot(breaks, set_NT_cumulative, type="n", xlab = "log2(FC)", ylab = "Cumulative fraction", ylim=c(0,1), xlim=c(-3,3))
lines(breaks, set_NT_cumulative,col="black")

selected_miR <- dataset_142
selected_miR <- dataset_146
selected_miR2 <- selected_miR[which(selected_miR$`Conserved sites total`>=2),]
set_sponge2 <-  as.numeric(selected_miR2$log2FC_selected) 
set_sponge2_cumulative = c(0, cumsum(table(cut(set_sponge2, breaks, right=TRUE))/nrow(selected_miR2))) 
lines(breaks, set_sponge2_cumulative,col="blue")

selected_miR <- selected_miR[which(selected_miR$`Conserved sites total`==1),]
set_sponge <-  as.numeric(selected_miR$log2FC_selected) 
set_sponge_cumulative = c(0, cumsum(table(cut(set_sponge, breaks, right=TRUE))/nrow(selected_miR))) 
lines(breaks, set_sponge_cumulative,col="red")




wilcox.test(dataset_NT$log2FC_selected,selected_miR$log2FC_selected)
wilcox.test(dataset_NT$log2FC_selected,selected_miR2$log2FC_selected)


write.table(rbind(breaks, set_NT_cumulative, set_sponge_cumulative, set_sponge2_cumulative), "cumulative-log2FC_142-146-2X-dataset_146.txt", sep="\t", append = FALSE)


