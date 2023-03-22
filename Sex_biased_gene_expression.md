# Analysis of masculinization of gene expression (Kallisto, DeSeq2)
```R
library(edgeR)
library(tximport)
library('edgeR')
library('rhdf5')
library('readxl')
library('ggplot2')
library(grid)
require('gridExtra')
library("org.Xl.eg.db")
library(PCAtools)
library("HTSFilter")
library(tidyverse)
library(purrr)
library(writexl)

setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_Supergene/2022_KO_tad_RNAseq/2022_EdgeR_and_DeSeq2/2023_Kallisto_DeSeq2_done")
dir <- "/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_Supergene/2022_KO_tad_RNAseq/2022_EdgeR_and_DeSeq2/2023_Kallisto_DeSeq2_done"
list.files(dir)


f_files<- list.files(".", pattern = "Kallisto_DeSeq2_unfiltered.csv", full.names = T);f_files
# import into a list
myfiles = lapply(f_files, read.delim, sep = ",")

# rename the columns so they are sensible
colnames(myfiles[[1]]) <- c("gene","MF_ccdc_baseMean","MF_ccdc_logFC","MF_ccdc_lfcSE","MF_ccdc_stat","MF_ccdc_pvalue","MF_ccdc_padj")
colnames(myfiles[[2]]) <- c("gene","MF_dmrt1L_baseMean","MF_dmrt1L_logFC","MF_dmrt1L_lfcSE","MF_dmrt1L_stat","MF_dmrt1L_pvalue","MF_dmrt1L_padj")
colnames(myfiles[[3]]) <- c("gene","MF_dmrt1S_baseMean","MF_dmrt1S_logFC","MF_dmrt1S_lfcSE","MF_dmrt1S_stat","MF_dmrt1S_pvalue","MF_dmrt1S_padj")
colnames(myfiles[[4]]) <- c("gene","wtko_ccdc_baseMean","wtko_ccdc_logFC","wtko_ccdc_lfcSE","wtko_ccdc_stat","wtko_ccdc_pvalue","wtko_ccdc_padj")
colnames(myfiles[[5]]) <- c("gene","wtko_dmw_baseMean","wtko_dmw_logFC","wtko_dmw_lfcSE","wtko_dmw_stat","wtko_dmw_pvalue","wtko_dmw_padj")
colnames(myfiles[[6]]) <- c("gene","wtko_scan_baseMean","wtko_scan_logFC","wtko_scan_lfcSE","wtko_scan_stat","wtko_scan_pvalue","wtko_scan_padj")

library(plyr)
alldata<-join_all(myfiles, by = "gene", type = "full", match = "all")
library(ggplot2)
library(GGally)

# get rid of outliers
# https://www.r-bloggers.com/2020/01/how-to-remove-outliers-in-r/
boxplot(alldata$MF_ccdc_logFC, plot=FALSE)$out
# these outliers are the first quartile - 1.5 the interquartile range
# and the third quartile plus 1.5 the interquartile range
#MF_ccdc
outliers <- boxplot(alldata$MF_ccdc_logFC, plot=FALSE)$out
MF_ccdc_logFC_trim<-alldata[,c(1,3)]
if(any(outliers)) {
  MF_ccdc_logFC_trim<- MF_ccdc_logFC_trim[-which(MF_ccdc_logFC_trim$MF_ccdc_logFC %in% outliers),]
}
#MF_dmrt1L
outliers <- boxplot(alldata$MF_dmrt1L_logFC, plot=FALSE)$out
MF_dmrt1L_logFC_trim<-alldata[,c(1,9)]
if(any(outliers)) {
  MF_dmrt1L_logFC_trim<- MF_dmrt1L_logFC_trim[-which(MF_dmrt1L_logFC_trim$MF_dmrt1L_logFC %in% outliers),]
}
#MF_dmrt1S
outliers <- boxplot(alldata$MF_dmrt1S_logFC, plot=FALSE)$out
MF_dmrt1S_logFC_trim<-alldata[,c(1,15)]
if(any(outliers)) {
  MF_dmrt1S_logFC_trim<- MF_dmrt1S_logFC_trim[-which(MF_dmrt1S_logFC_trim$MF_dmrt1S_logFC %in% outliers),]
}

#wtko_dmw
outliers <- boxplot(alldata$wtko_dmw_logFC, plot=FALSE)$out
wtko_dmw_logFC_trim<- alldata[,c(1,27)]
if(any(outliers)) {
  wtko_dmw_logFC_trim<- wtko_dmw_logFC_trim[-which(wtko_dmw_logFC_trim$wtko_dmw_logFC %in% outliers),]
}
#wtko_scan
outliers <- boxplot(alldata$wtko_scan_logFC, plot=FALSE)$out
wtko_scan_logFC_trim<- alldata[,c(1,33)]
if(any(outliers)) {
  wtko_scan_logFC_trim<- wtko_scan_logFC_trim[-which(wtko_scan_logFC_trim$wtko_scan_logFC %in% outliers),]
}
#wtko_ccdc
outliers <- boxplot(alldata$wtko_ccdc_logFC, plot=FALSE)$out
wtko_ccdc_logFC_trim<- alldata[,c(1,21)]
if(any(outliers)) {
  wtko_ccdc_logFC_trim<- wtko_ccdc_logFC_trim[-which(wtko_ccdc_logFC_trim$wtko_ccdc_logFC %in% outliers),]
}

# combine no outlier files
# make a list of df
df_list <- list(MF_ccdc_logFC_trim, MF_dmrt1L_logFC_trim, MF_dmrt1S_logFC_trim, wtko_dmw_logFC_trim, wtko_scan_logFC_trim, wtko_ccdc_logFC_trim)
#merge all data frames in list
alldata_no_outliers <- df_list %>% reduce(full_join, by='gene')

colnames(alldata_no_outliers) <- c("gene","MF_1","MF_2",
                                   "MF_3","dmw","scan",
                                   "ccdc")
library(ggplot2)

my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
  #  geom_smooth(method=loess, fill="red", color="red", ...) +
    geom_smooth(method=lm, fill="blue", color="blue", ...)
  p
}

my_fn2 <- function(data, mapping, method="p", use="pairwise", ...){
  
  # grab data
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  # calculate correlation
  corr <- cor(x, y, method=method, use=use)
  
  # calculate colour based on correlation value
  # Here I have set a correlation of minus one to blue, 
  # zero to white, and one to red 
  # Change this to suit: possibly extend to add as an argument of `my_fn`
  colFn <- colorRampPalette(c("blue", "white", "red"), interpolate ='spline')
  fill <- colFn(100)[findInterval(corr, seq(-1, 1, length=100))]
  
  ggally_cor(data = data, mapping = mapping, ...) + 
    theme_void() +
    theme(panel.background = element_rect(fill=fill))
}

my_custom_smooth <- function(data, mapping, ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(color = I("blue")) + 
    geom_smooth(method = "lm", fill="blue", color="blue", ...)
  
  lmModel <- eval(substitute(lm(y ~ x, data = data), mapping))
  fs <- summary(lmModel)$fstatistic
  pValue <- pf(fs[1], fs[2], fs[3], lower.tail = FALSE)
  
  if (pValue < 0.05) {
    p <- p + theme(
      panel.border = element_rect(
        color = "red", 
        size = 3,
        linetype = "solid",
        fill = "transparent"
      )
    )
  }
  
  p
}
p_ <- GGally::print_if_interactive
g<-ggpairs(alldata_no_outliers[,c(2:7)], 
        #upper = list(continuous = "density", combo = "box_no_facet"),
        #upper = list(continuous = wrap(ggally_cor, size = 2)), 
        # upper = list(continuous = my_fn2),
        upper = "blank",
        lower = list(continuous = my_fn)) +
        #lower = list(continuous = my_custom_smooth)) +
  #theme_bw() +
  theme(strip.background = element_rect(
      color="white", fill="white", size=1.5, linetype="solid")) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank());g

box1_2 <- ggally_text("\nr = 0.285\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box1_3 <- ggally_text("\nr = 0.134\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box1_4 <- ggally_text("\nr = 0.468*\np = 0.013*\n",geom_text = ggplot2::aes(size = 6), color = I("red"))
box1_5 <- ggally_text("\nr = 0.091\np = 0.340\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box1_6 <- ggally_text("\nr = 0.347*\np = 0.547\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box2_3 <- ggally_text("\nr = -0.127\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box2_4 <- ggally_text("\nr = 0.205\np = 0.391\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box2_5 <- ggally_text("\nr = -0.250\np = 0.924\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box2_6 <- ggally_text("\nr = -0.068\np = 0.458\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box3_4 <- ggally_text("\nr = 0.613*\np = 0.001*\n",geom_text = ggplot2::aes(size = 6), color = I("red"))
box3_5 <- ggally_text("\nr = 0.119\np = 0.390\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box3_6 <- ggally_text("\nr = -0.204\np = 0.420\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box4_5 <- ggally_text("\nr = -0.050\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box4_6 <- ggally_text("\nr = 0.158\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box5_6 <- ggally_text("\nr = -0.021\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
g[1, 2] <- box1_2
g[1, 3] <- box1_3
g[1, 4] <- box1_4
g[1, 5] <- box1_5
g[1, 6] <- box1_6
g[2, 3] <- box2_3
g[2, 4] <- box2_4
g[2, 5] <- box2_5
g[2, 6] <- box2_6
g[3, 4] <- box3_4
g[3, 5] <- box3_5
g[3, 6] <- box3_6
g[4, 5] <- box4_5
g[4, 6] <- box4_6
g[5, 6] <- box5_6
# small function to display plots only if it's interactive

p_(g)


ggsave(file="Kallisto_DeSeq2_sexrelated_pairwise_unfiltered_new.pdf", g, width=10, height=4)

```