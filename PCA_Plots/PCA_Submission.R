## PCA plots to show similarities between samples

rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
workd<-getwd()

library(readxl)
library(dplyr)
library(Hmisc)
library(ggplot2)
library(limma)  # for plotMDS

# Read in data

baseline = paste0(workd,"/pub_2024-07-07_aging_baseline_cyt_cytof_combined_wide.xlsx")
timepoint1 = paste0(workd,"/pub_2024-07-07_aging_fold.change_cyt_cytof_combined_wide.xlsx")
foldchange = paste0(workd,"/pub_2024-07-07_aging_timepoint1_cyt_cytof_combined_wide.xlsx")
metadata = paste0(workd,"/md.xlsx")

baseline <- read_xlsx(baseline)
baseline_cytokine<-baseline[,1:33] #subset cytokine columns
baseline_props<-baseline[,c(1,34:49)] #subset cytof props columns

timepoint1 <- read_xlsx(timepoint1)
timepoint1_cytokine<-timepoint1[,1:33] #subset cytokine columns
timepoint1_props<-timepoint1[,c(1,34:49)] #subset cytof props columns

foldchange <- read_xlsx(foldchange)
foldchange_cytokine<-foldchange[,1:33] #subset cytokine columns
foldchange_props<-foldchange[,c(1,34:49)] #subset cytof props columns

md <- read_xlsx(metadata)

#### All Samples Baseline ####

expr_mean_sample_tbl_baseline <- data.frame(baseline) # make the dataframe
expr_mean_sample_baseline <- t(expr_mean_sample_tbl_baseline[, -1]) # remove the pt_id column & transpose
colnames(expr_mean_sample_baseline) <- expr_mean_sample_tbl_baseline$pt_id # add pt_ids as column names
data_baseline <- na.omit(expr_mean_sample_baseline)

#### All samples Timepoint 1 ####
expr_mean_sample_tbl_timepoint1 <- data.frame(timepoint1) # make the dataframe
expr_mean_sample_timepoint1 <- t(expr_mean_sample_tbl_timepoint1[, -1]) # remove the pt_id column & transpose
colnames(expr_mean_sample_timepoint1) <- expr_mean_sample_tbl_timepoint1$pt_id # add pt_ids as column names
data_timepoint1 <- na.omit(expr_mean_sample_timepoint1)


#### All samples FC ####
expr_mean_sample_tbl_fc <- data.frame(foldchange) # make the dataframe
expr_mean_sample_fc <- t(expr_mean_sample_tbl_fc[, -1]) # remove the pt_id column & transpose
colnames(expr_mean_sample_fc) <- expr_mean_sample_tbl_fc$pt_id # add pt_ids as column names
data_fc <- na.omit(expr_mean_sample_fc)

#### normalizing data across both modalities and capping outlier data  #####

rng <- colQuantiles(t(expr_mean_sample_baseline), probs = c(0, 0.95), na.rm = T)
data01_0 <- t((expr_mean_sample_baseline - rng[, 1]) / (rng[, 2] - rng[, 1]))
data01_0[data01_0 < 0] <- 0; data01_0[data01_0 > 1] <- 1
data01_0 <- t(data01_0)

rng <- colQuantiles(t(expr_mean_sample_timepoint1), probs = c(0, 0.95), na.rm = T)
data01_1 <- t((expr_mean_sample_timepoint1 - rng[, 1]) / (rng[, 2] - rng[, 1]))
data01_1[data01_1 < 0] <- 0; data01_1[data01_1 > 1] <- 1
data01_1 <- na.omit(t(data01_1))

rng <- colQuantiles(t(expr_mean_sample_fc), probs = c(0, 0.95), na.rm = T)
data01_fc <- t((expr_mean_sample_fc - rng[, 1]) / (rng[, 2] - rng[, 1]))
data01_fc[data01_fc < 0] <- 0; data01_fc[data01_fc > 1] <- 1
data01_fc <- na.omit(t(data01_fc))

#### select data to plot ####

data01 <- t(data01_0); ggtitle <- "baseline"
data01 <- t(data01_1); ggtitle <- "timepoint1"
data01 <- t(data01_fc); ggtitle <- "fold_change"

#### Running PCA ####

pca_res <- prcomp(data01, scale. = T)
ggdf<- data.frame(pca_res$x, pt_id = rownames(data01)) #scores
agedlevels=c("Yes","No")
ggdf$Aged <- factor(md$Aged[match(ggdf$pt_id,md$pt_id)], levels=agedlevels)
responselevels=c("Yes","No","NA")
ggdf$response <- factor(md$Response[match(ggdf$pt_id,md$pt_id)], levels=responselevels)
ggdf <- ggdf[ggdf$response %nin% c("NA"),]

htmp <- Heatmap(pca_res$rotation, #loadings
                name='loadings',
                cluster_rows = F,
                cluster_columns = F)

#### sorting by aged and save plot ####

means <- ggdf %>% group_by(Aged) %>% summarise_all(mean,na.rm=T)
diff <- means[1,] - means[2,]
diff <- unlist(diff)

summ <- summary(pca_res)
top_pc <- rev(sort(abs(diff[2:48])))

ggp <- ggplot(ggdf, aes(x = !!sym(names(top_pc[1])), y = !!sym(names(top_pc[2])), color = Aged, shape = response)) +
  ggtitle(ggtitle) +
  geom_point(size = 1.5) +
  scale_shape_manual(values = c(19, 17))+
  xlab(paste0(names(top_pc[1]), " (", summ$importance[,names(top_pc[1])][2], ")"))+
  ylab(paste0(names(top_pc[2]), " (", summ$importance[,names(top_pc[2])][2], ")"))+
  theme(plot.background = element_rect(fill="white"),
        panel.background = element_rect(fill="white"),
        panel.grid = element_blank(),
        axis.line = element_line(color="black", linewidth=0.25),
        axis.ticks = element_line(linewidth=0.25),
        axis.text = element_text(color="black"),
        axis.title = element_text(color="black"),
        legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill="white"),
        legend.text = element_text(color="black"))

#### sorting by response and save plot ####

means <- ggdf %>% group_by(response) %>% summarise_all(mean,na.rm=T)
diff <- means[1,] - means[2,]
diff <- unlist(diff)

top_pc <- rev(sort(abs(diff[2:48])))

ggp2 <- ggplot(ggdf, aes(x = !!sym(names(top_pc[1])), y = !!sym(names(top_pc[2])), color = response, shape = Aged)) +
  ggtitle(ggtitle) +
  geom_point(size = 1.5) +
  scale_shape_manual(values = c(19, 17))+
  xlab(paste0(names(top_pc[1]), " (", summ$importance[,names(top_pc[1])][2], ")"))+
  ylab(paste0(names(top_pc[2]), " (", summ$importance[,names(top_pc[2])][2], ")"))+
  theme(plot.background = element_rect(fill="white"),
        panel.background = element_rect(fill="white"),
        panel.grid = element_blank(),
        axis.line = element_line(color="black", linewidth=0.25),
        axis.ticks = element_line(linewidth=0.25),
        axis.text = element_text(color="black"),
        axis.title = element_text(color="black"),
        legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill="white"),
        legend.text = element_text(color="black"))

#### print out plots ####

pdf(paste0('pca_',ggtitle,'.pdf'),width=5,height=3.5);ggp;ggp2;dev.off()
pdf(paste0('pca_',ggtitle,'_loadings.pdf'),width=8,height=8);htmp;dev.off()


