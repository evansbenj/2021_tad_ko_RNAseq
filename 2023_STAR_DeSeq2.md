# STAR DeSeq2 plus permutations
```
## RNA-seq analysis with DESeq2
## Adapted from Stephen Turner, @genetics_blog
library(DESeq2)
library(tidyverse)
library(purrr)

# https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#input-data
# if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("gplots")
# BiocManager::install("DESeq2")
library(apeglm)
# RNA-seq data from PRJNA315516
# https://www.ncbi.nlm.nih.gov//bioproject/PRJNA315516.
# 3 control samples ("ctl"), 3 samples grown under drought condition ("dro")

# Import & pre-process ----------------------------------------------------

# Import data from featureCounts
## Previously ran at command line something like this:
setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_Supergene/2022_KO_tad_RNAseq/2022_EdgeR_and_DeSeq2/STAR_done")
dir <- "/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_Supergene/2022_KO_tad_RNAseq/2022_EdgeR_and_DeSeq2/STAR_done"
list.files(dir)


# load the count data from STAR
f_files<- list.files(".", pattern = "_counts.out", full.names = T);f_files

read_in_feature_counts<- function(file){
  cnt<- read_tsv(file, col_names =T, comment = "#")
  cnt<- cnt %>% dplyr::select(-Chr, -Start, -End, -Strand, -Length)
  return(cnt)
}

raw_counts<- map(f_files, read_in_feature_counts)
raw_counts_df<- purrr::reduce(raw_counts, inner_join) 

dim(raw_counts_df)
# [1] 44478    67
# View(raw_counts_df)

# get rid of any rows that have incomplete data
counts <- raw_counts_df[complete.cases(raw_counts_df), ]
dim(counts)
# [1]  44478    67 # 2023_STAR

colnames(counts) <- c("geneID","ccdc_11_wt_M","ccdc_12_wt_M","ccdc_13_wt_M",
                      "ccdc_14_ko_F","ccdc_2_wt_M","ccdc_25_wt_M","ccdc_3_wt_F",
                      "ccdc_30_ko_F","ccdc_32_ko_F","ccdc_34_wt_F","ccdc_35_ko_F",
                      "ccdc_36_ko_F","ccdc_42_ko_F","ccdc_9_wt_M","dmrt1L_11_ko_F",
                      "dmrt1L_13_wt_M","dmrt1L_15_wt_M","dmrt1L_17_ko_M",
                      "dmrt1L_19_ko_M","dmrt1L_24_ko_F","dmrt1L_25_ko_F",
                      "dmrt1L_26_ko_M","dmrt1L_27_ko_F","dmrt1L_30_wt_F",
                      "dmrt1L_34_wt_M","dmrt1L_35_ko_M","dmrt1L_41_ko_M",
                      "dmrt1L_43_wt_M","dmrt1L_50_wt_F","dmrt1L_55_wt_F",
                      "dmrt1L_59_ko_M","dmrt1L_6_ko_F","dmrt1L_7_ko_F",
                      "dmrt1L_8_wt_M","dmrt1S_10_ko_F","dmrt1S_14_wt_M",
                      "dmrt1S_15_ko_F","dmrt1S_18_wt_M","dmrt1S_20_wt_M",
                      "dmrt1S_23_ko_F","dmrt1S_24_wt_F","dmrt1S_29_wt_F",
                      "dmrt1S_30_ko_F","dmrt1S_4_wt_F","dmrt1S_5_ko_M",
                      "dmw_12_wt_F","dmw_14_ko_F","dmw_15_wt_F","dmw_16_ko_F",
                      "dmw_17_wt_F","dmw_20_wt_F","dmw_22_wt_F","dmw_26_ko_F",
                      "dmw_28_ko_F","dmw_29_ko_F","dmw_35_ko_F","dmw_9_wt_F",
                      "scanw_14_wt_F","scanw_15_wt_F","scanw_19_ko_F",
                      "scanw_27_ko_F","scanw_28_ko_F","scanw_30_ko_F",
                      "scanw_31_ko_F","scanw_6_wt_F","scanw_9_wt_F")

# this is useful info about batch effects:
# https://support.bioconductor.org/p/96627/

# make design matrix
# Here we want to test for differential expression between KO and
# WT, while adjusting for differences between batches. In statistical
# terms, this "may be" an additive linear model with batch as the blocking factor:

# countz
countz <- counts[,-c(1)]
gene_names <- counts$geneID
row.names(countz) <- gene_names
#View(countz)
# samples
samples <- read.table(file.path(dir, "samples.txt"), header = T)
samples
# sexez
sexez <- factor(samples$sex)
sexez <- relevel(sexez, ref="F")
# batch
batchez <- factor(samples$batch)
row.names(samples) <- colnames(countz)



# MF dmrt1S ----
colnames(countz)
new_counts <- as.data.frame(countz[,-c(1:35,37,40,43,45,46:66) ])
row.names(new_counts) <- gene_names
new_samples <-as.data.frame(samples[-c(1:35,37,40,43,45,46:66), ]);new_samples
new_sexez <- factor(samples$sex[-c(1:35,37,40,43,45,46:66)])
new_sexez <- relevel(new_sexez, ref="F")
new_batchez <- factor(samples$batch[-c(1:35,37,40,43,45,46:66)])
# relevel
new_samples$sex <-factor(new_samples$sex, levels = c("F","M"))
new_samples$sex <- relevel(new_samples$sex, ref="F")
new_samples$batch <-factor(new_samples$batch, levels = c("dmrt1S"))
dds <- DESeqDataSetFromMatrix(countData = new_counts,
                              colData = new_samples,
                              design= ~sex)
# save the unfiltered logFC to a dataframe 
dds$sex <- relevel(dds$sex, ref="F") # this makes the expression levels relative to M
dds <- DESeq(dds)
res <- results(dds)
MF_dmrt1S_unfiltered <- res;MF_dmrt1S_unfiltered
# Only sex related
sex_related_MF_dmrt1S <- res[c('XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274'),]
write.csv(sex_related_MF_dmrt1S, file="Sex_related_MF_dmrt1S_Kallisto_DeSeq2_unfiltered.csv", row.names = T)
# Write counts of sex related to a file
sex_related_MF_dmrt1S_counts <- new_counts[c('XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274'),]
write.csv(sex_related_MF_dmrt1S_counts, file="Sex_related_MF_dmrt1S_STAR_edgeR_counts_unfiltered.csv", row.names = T)

# Now do analysis of differential expression; 
# first remove transcripts where the average count per sample is 2 or less:
keep <- rowSums(counts(dds)) >= 2* length(colnames(dds))
dds <- dds[keep,]
# relevel
dds$sex <- relevel(dds$sex, ref="F") # this makes the expression levels relative to M
# now do the analysis
dds <- DESeq(dds)
res <- results(dds)
summary(res)
resOrdered <- res[order(res$pvalue),]
p<-resOrdered[1:1053,];p
write.csv(p, file="MF_STAR_dmrt1Sonly_DE_DeSeq2.csv", row.names = T)


# MF dmrt1L ----
colnames(countz)
new_counts <- as.data.frame(countz[,-c(1:15,18:23,26:27,31:33,35:66) ])
row.names(new_counts) <- gene_names
new_samples <-as.data.frame(samples[-c(1:15,18:23,26:27,31:33,35:66), ]);new_samples
new_sexez <- factor(samples$sex[-c(1:15,18:23,26:27,31:33,35:66)])
new_sexez <- relevel(new_sexez, ref="F")
new_batchez <- factor(samples$batch[-c(1:15,18:23,26:27,31:33,35:66)])
# relevel
new_samples$sex <-factor(new_samples$sex, levels = c("F","M"))
new_samples$sex <- relevel(new_samples$sex, ref="F")
new_samples$batch <-factor(new_samples$batch, levels = c("dmrt1S"))
dds <- DESeqDataSetFromMatrix(countData = new_counts,
                              colData = new_samples,
                              design= ~sex)
# save the unfiltered logFC to a dataframe 
dds$sex <- relevel(dds$sex, ref="F") # this makes the expression levels relative to M
dds <- DESeq(dds)
res <- results(dds)
MF_dmrt1L_unfiltered <- res;MF_dmrt1L_unfiltered
# Only sex related
sex_related_MF_dmrt1L <- res[c('XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274'),]
write.csv(sex_related_MF_dmrt1L, file="Sex_related_MF_dmrt1L_Kallisto_DeSeq2_unfiltered.csv", row.names = T)

# Now do analysis of differential expression; 
# first remove transcripts where the average count per sample is 2 or less:
keep <- rowSums(counts(dds)) >= 2* length(colnames(dds))
dds <- dds[keep,]
# relevel
dds$sex <- relevel(dds$sex, ref="F") # this makes the expression levels relative to M
# now do the analysis
dds <- DESeq(dds)
res <- results(dds)
summary(res)
resOrdered <- res[order(res$pvalue),]
p<-resOrdered[1:47,];p
write.csv(p, file="MF_STAR_dmrt1Lonly_DE_DeSeq2.csv", row.names = T)


# MF ccdc ----
colnames(countz)
new_counts <- as.data.frame(countz[,-c(4,8:9,11:13,15:66) ])
row.names(new_counts) <- gene_names
new_samples <-as.data.frame(samples[-c(4,8:9,11:13,15:66), ]);new_samples
new_sexez <- factor(samples$sex[-c(4,8:9,11:13,15:66)])
new_sexez <- relevel(new_sexez, ref="F")
new_batchez <- factor(samples$batch[-c(4,8:9,11:13,15:66)])
# relevel
new_samples$sex <-factor(new_samples$sex, levels = c("F","M"))
new_samples$sex <- relevel(new_samples$sex, ref="F")
new_samples$batch <-factor(new_samples$batch, levels = c("ccdc"))
dds <- DESeqDataSetFromMatrix(countData = new_counts,
                              colData = new_samples,
                              design= ~sex)
# save the unfiltered logFC to a dataframe 
dds$sex <- relevel(dds$sex, ref="F") # this makes the expression levels relative to M
dds <- DESeq(dds)
res <- results(dds)
MF_ccdc_unfiltered <- res;MF_ccdc_unfiltered
# Only sex related
sex_related_MF_ccdc <- res[c('XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274'),]
write.csv(sex_related_MF_ccdc, file="Sex_related_MF_ccdc_Kallisto_DeSeq2_unfiltered.csv", row.names = T)

# Now do analysis of differential expression; 
# first remove transcripts where the average count per sample is 2 or less:
keep <- rowSums(counts(dds)) >= 2* length(colnames(dds))
dds <- dds[keep,]
# relevel
dds$sex <- relevel(dds$sex, ref="F") # this makes the expression levels relative to M
# now do the analysis
dds <- DESeq(dds)
res <- results(dds)
summary(res)
resOrdered <- res[order(res$pvalue),]
p<-resOrdered[1:7,];p
write.csv(p, file="MF_STAR_ccdconly_DE_DeSeq2.csv", row.names = T)


# wtko dmw ----
colnames(countz)
new_counts <- as.data.frame(countz[,c(46:57) ])
row.names(new_counts) <- gene_names
new_samples <-as.data.frame(samples[c(46:57), ]);new_samples
new_sexez <- factor(samples$sex[c(46:57)])
new_sexez <- relevel(new_sexez, ref="F")
new_batchez <- factor(samples$batch[c(46:57)])
# relevel
new_samples$sex <-factor(new_samples$sex, levels = c("F","M"))
new_samples$sex <- relevel(new_samples$sex, ref="F")
new_samples$genotype <-factor(new_samples$genotype, levels = c("wt","ko"))
new_samples$genotype <- relevel(new_samples$genotype, ref="wt")
new_samples$batch <-factor(new_samples$batch, levels = c("dmw"))
dds <- DESeqDataSetFromMatrix(countData = new_counts,
                              colData = new_samples,
                              design= ~genotype)
# save the unfiltered logFC to a dataframe 
dds$genotype <- relevel(dds$genotype, ref="wt") # this makes the expression levels relative to M
dds <- DESeq(dds)
res <- results(dds)
wtko_dmw_unfiltered <- res;wtko_dmw_unfiltered
# Only sex related
sex_related_wtko_dmw <- res[c('XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274'),]
write.csv(sex_related_wtko_dmw, file="Sex_related_wtko_dmw_Kallisto_DeSeq2_unfiltered.csv", row.names = T)

# Now do analysis of differential expression; 
# first remove transcripts where the average count per sample is 2 or less:
keep <- rowSums(counts(dds)) >= 2* length(colnames(dds))
dds <- dds[keep,]
# relevel
dds$genotype <- relevel(dds$genotype, ref="wt") # this makes the expression levels relative to M
# now do the analysis
dds <- DESeq(dds)
res <- results(dds)
summary(res)
resOrdered <- res[order(res$pvalue),]
p<-resOrdered[1:18,];p
write.csv(p, file="wt_ko_STAR_dmw_DE_DeSeq2.csv", row.names = T)


# wtko scan ----
colnames(countz)
new_counts <- as.data.frame(countz[,c(58:66) ])
row.names(new_counts) <- gene_names
new_samples <-as.data.frame(samples[c(58:66), ]);new_samples
# relevel
new_samples$genotype <-factor(new_samples$genotype, levels = c("wt","ko"))
new_samples$genotype <- relevel(new_samples$genotype, ref="wt")
dds <- DESeqDataSetFromMatrix(countData = new_counts,
                              colData = new_samples,
                              design= ~genotype)
# save the unfiltered logFC to a dataframe 
dds$genotype <- relevel(dds$genotype, ref="wt") # this makes the expression levels relative to M
dds <- DESeq(dds)
res <- results(dds)
wtko_scan_unfiltered <- res;wtko_scan_unfiltered
# Only sex related
sex_related_wtko_scan <- res[c('XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274'),]
write.csv(sex_related_wtko_scan, file="Sex_related_wtko_scan_Kallisto_DeSeq2_unfiltered.csv", row.names = T)

# Now do analysis of differential expression; 
# first remove transcripts where the average count per sample is 2 or less:
keep <- rowSums(counts(dds)) >= 2* length(colnames(dds))
dds <- dds[keep,]
# relevel
dds$genotype <- relevel(dds$genotype, ref="wt") # this makes the expression levels relative to M
# now do the analysis
dds <- DESeq(dds)
res <- results(dds)
summary(res)
resOrdered <- res[order(res$pvalue),]
p<-resOrdered[1:32,];p
write.csv(p, file="wt_ko_STAR_scan_DE_DeSeq2.csv", row.names = T)

# wtko ccdc ----
colnames(countz)
new_counts <- as.data.frame(countz[,c(4,7:13) ])
row.names(new_counts) <- gene_names
new_samples <-as.data.frame(samples[c(4,7:13), ]);new_samples
# relevel
new_samples$genotype <-factor(new_samples$genotype, levels = c("wt","ko"))
new_samples$genotype <- relevel(new_samples$genotype, ref="wt")
dds <- DESeqDataSetFromMatrix(countData = new_counts,
                              colData = new_samples,
                              design= ~genotype)
# save the unfiltered logFC to a dataframe 
dds$genotype <- relevel(dds$genotype, ref="wt") # this makes the expression levels relative to M
dds <- DESeq(dds)
res <- results(dds)
wtko_ccdc_unfiltered <- res;wtko_ccdc_unfiltered
# Only sex related
sex_related_wtko_ccdc <- res[c('XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274'),]
write.csv(sex_related_wtko_ccdc, file="Sex_related_wtko_ccdc_Kallisto_DeSeq2_unfiltered.csv", row.names = T)

# Now do analysis of differential expression; 
# first remove transcripts where the average count per sample is 2 or less:
keep <- rowSums(counts(dds)) >= 2* length(colnames(dds))
dds <- dds[keep,]
# relevel
dds$genotype <- relevel(dds$genotype, ref="wt") # this makes the expression levels relative to M
# now do the analysis
dds <- DESeq(dds)
res <- results(dds)
summary(res)
resOrdered <- res[order(res$pvalue),]
p<-resOrdered[1:263,];p
write.csv(p, file="wt_ko_STAR_ccdc_DE_DeSeq2.csv", row.names = T)

# dmw permutations ----

# MF_ccdc vs dmw
correlations <- c()
# get rownames of sexrelated transcripts
SL_rownames <- c('XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274')

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(counts)[1], 74, replace = F);indexes
  rownames <- counts$geneID[indexes]
  # remove outliers from MF
  MF_ccdc_trim <- MF_ccdc_unfiltered[rownames,]
  outliers <- boxplot(MF_ccdc_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_ccdc_trim<- MF_ccdc_trim[-which(MF_ccdc_trim$log2FoldChange %in% outliers),]
  }  
  # remove outliers from wtko
  wtko_dmw_trim <- wtko_dmw_unfiltered[rownames,]
  outliers <- boxplot(wtko_dmw_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    wtko_dmw_trim<- wtko_dmw_trim[-which(wtko_dmw_trim$log2FoldChange %in% outliers),]
  }
  correlations[x] <- cor(MF_ccdc_trim[rownames,'log2FoldChange'],
                         wtko_dmw_trim[rownames,'log2FoldChange'], 
                         method = "pearson", use="pairwise")
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_ccdc_trim <- sex_related_MF_ccdc[SL_rownames,]
outliers <- boxplot(sex_related_MF_ccdc_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_ccdc_trim <- sex_related_MF_ccdc_trim[-which(sex_related_MF_ccdc_trim$log2FoldChange %in% outliers),]
}  
# remove outliers from wtko
sex_related_wtko_dmw_trim <- sex_related_wtko_dmw[SL_rownames,]
outliers <- boxplot(sex_related_wtko_dmw_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_wtko_dmw_trim <- sex_related_wtko_dmw_trim[-which(sex_related_wtko_dmw_trim$log2FoldChange %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_ccdc_trim[SL_rownames,'log2FoldChange'],
                          sex_related_wtko_dmw_trim[SL_rownames,'log2FoldChange'], 
                          method = "pearson", use="pairwise")
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1] 0.3556444


# MF_dmrt1L vs dmw
correlations <- c()
# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 74, replace = F);indexes
  rownames <- counts$geneID[indexes]
  # remove outliers from MF
  MF_dmrt1L_trim <- MF_dmrt1L_unfiltered[rownames,]
  outliers <- boxplot(MF_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_dmrt1L_trim<- MF_dmrt1L_trim[-which(MF_dmrt1L_trim$log2FoldChange %in% outliers),]
  }  
  # remove outliers from wtko
  wtko_dmw_trim <- wtko_dmw_unfiltered[rownames,]
  outliers <- boxplot(wtko_dmw_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    wtko_dmw_trim<- wtko_dmw_trim[-which(wtko_dmw_trim$log2FoldChange %in% outliers),]
  }
  correlations[x] <- cor(MF_dmrt1L_trim[rownames,'log2FoldChange'],
                         wtko_dmw_trim[rownames,'log2FoldChange'], 
                         method = "pearson", use="pairwise")
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_dmrt1L_trim <- sex_related_MF_dmrt1L[SL_rownames,]
outliers <- boxplot(sex_related_MF_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_dmrt1L_trim <- sex_related_MF_dmrt1L_trim[-which(sex_related_MF_dmrt1L_trim$log2FoldChange %in% outliers),]
}  
# remove outliers from wtko
sex_related_wtko_dmw_trim <- sex_related_wtko_dmw[SL_rownames,]
outliers <- boxplot(sex_related_wtko_dmw_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_wtko_dmw_trim <- sex_related_wtko_dmw_trim[-which(sex_related_wtko_dmw_trim$log2FoldChange %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_dmrt1L_trim[SL_rownames,'log2FoldChange'],
                          sex_related_wtko_dmw_trim[SL_rownames,'log2FoldChange'], 
                          method = "pearson", use="pairwise")
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1] 0.5874126

# MF_dmrt1S vs dmw
correlations <- c()
# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 74, replace = F);indexes
  rownames <- counts$geneID[indexes]
  # remove outliers from MF
  MF_dmrt1S_trim <- MF_dmrt1S_unfiltered[rownames,]
  outliers <- boxplot(MF_dmrt1S_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_dmrt1S_trim<- MF_dmrt1S_trim[-which(MF_dmrt1S_trim$log2FoldChange %in% outliers),]
  }  
  # remove outliers from wtko
  wtko_dmw_trim <- wtko_dmw_unfiltered[rownames,]
  outliers <- boxplot(wtko_dmw_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    wtko_dmw_trim<- wtko_dmw_trim[-which(wtko_dmw_trim$log2FoldChange %in% outliers),]
  }
  correlations[x] <- cor(MF_dmrt1S_trim[rownames,'log2FoldChange'],
                         wtko_dmw_trim[rownames,'log2FoldChange'], 
                         method = "pearson", use="pairwise")
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_dmrt1S_trim <- sex_related_MF_dmrt1S[SL_rownames,]
outliers <- boxplot(sex_related_MF_dmrt1S_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_dmrt1S_trim <- sex_related_MF_dmrt1S_trim[-which(sex_related_MF_dmrt1S_trim$log2FoldChange %in% outliers),]
}  
# remove outliers from wtko
sex_related_wtko_dmw_trim <- sex_related_wtko_dmw[SL_rownames,]
outliers <- boxplot(sex_related_wtko_dmw_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_wtko_dmw_trim <- sex_related_wtko_dmw_trim[-which(sex_related_wtko_dmw_trim$log2FoldChange %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_dmrt1S_trim[SL_rownames,'log2FoldChange'],
                          sex_related_wtko_dmw_trim[SL_rownames,'log2FoldChange'], 
                          method = "pearson", use="pairwise")
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1] 0.02497502


# scan permutations ----

# MF_ccdc vs scan
correlations <- c()
# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 74, replace = F);indexes
  rownames <- counts$geneID[indexes]
  # remove outliers from MF
  MF_ccdc_trim <- MF_ccdc_unfiltered[rownames,]
  outliers <- boxplot(MF_ccdc_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_ccdc_trim<- MF_ccdc_trim[-which(MF_ccdc_trim$log2FoldChange %in% outliers),]
  }  
  # remove outliers from wtko
  wtko_scan_trim <- wtko_scan_unfiltered[rownames,]
  outliers <- boxplot(wtko_scan_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    wtko_scan_trim<- wtko_scan_trim[-which(wtko_scan_trim$log2FoldChange %in% outliers),]
  }
  correlations[x] <- cor(MF_ccdc_trim[rownames,'log2FoldChange'],
                         wtko_scan_trim[rownames,'log2FoldChange'], 
                         method = "pearson", use="pairwise")
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_ccdc_trim <- sex_related_MF_ccdc[SL_rownames,]
outliers <- boxplot(sex_related_MF_ccdc_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_ccdc_trim <- sex_related_MF_ccdc_trim[-which(sex_related_MF_ccdc_trim$log2FoldChange %in% outliers),]
}  
# remove outliers from wtko
sex_related_wtko_scan_trim <- sex_related_wtko_scan[SL_rownames,]
outliers <- boxplot(sex_related_wtko_scan_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_wtko_scan_trim <- sex_related_wtko_scan_trim[-which(sex_related_wtko_scan_trim$log2FoldChange %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_ccdc_trim[SL_rownames,'log2FoldChange'],
                          sex_related_wtko_scan_trim[SL_rownames,'log2FoldChange'], 
                          method = "pearson", use="pairwise")
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1] 0.2577423


# MF_dmrt1L vs scan
correlations <- c()
# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 74, replace = F);indexes
  rownames <- counts$geneID[indexes]
  # remove outliers from MF
  MF_dmrt1L_trim <- MF_dmrt1L_unfiltered[rownames,]
  outliers <- boxplot(MF_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_dmrt1L_trim<- MF_dmrt1L_trim[-which(MF_dmrt1L_trim$log2FoldChange %in% outliers),]
  }  
  # remove outliers from wtko
  wtko_scan_trim <- wtko_scan_unfiltered[rownames,]
  outliers <- boxplot(wtko_scan_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    wtko_scan_trim<- wtko_scan_trim[-which(wtko_scan_trim$log2FoldChange %in% outliers),]
  }
  correlations[x] <- cor(MF_dmrt1L_trim[rownames,'log2FoldChange'],
                         wtko_scan_trim[rownames,'log2FoldChange'], 
                         method = "pearson", use="pairwise")
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_dmrt1L_trim <- sex_related_MF_dmrt1L[SL_rownames,]
outliers <- boxplot(sex_related_MF_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_dmrt1L_trim <- sex_related_MF_dmrt1L_trim[-which(sex_related_MF_dmrt1L_trim$log2FoldChange %in% outliers),]
}  
# remove outliers from wtko
sex_related_wtko_scan_trim <- sex_related_wtko_scan[SL_rownames,]
outliers <- boxplot(sex_related_wtko_scan_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_wtko_scan_trim <- sex_related_wtko_scan_trim[-which(sex_related_wtko_scan_trim$log2FoldChange %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_dmrt1L_trim[SL_rownames,'log2FoldChange'],
                          sex_related_wtko_scan_trim[SL_rownames,'log2FoldChange'], 
                          method = "pearson", use="pairwise")
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1] 0.2177822

# MF_dmrt1S vs scan
correlations <- c()
# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 74, replace = F);indexes
  rownames <- counts$geneID[indexes]
  # remove outliers from MF
  MF_dmrt1S_trim <- MF_dmrt1S_unfiltered[rownames,]
  outliers <- boxplot(MF_dmrt1S_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_dmrt1S_trim<- MF_dmrt1S_trim[-which(MF_dmrt1S_trim$log2FoldChange %in% outliers),]
  }  
  # remove outliers from wtko
  wtko_scan_trim <- wtko_scan_unfiltered[rownames,]
  outliers <- boxplot(wtko_scan_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    wtko_scan_trim<- wtko_scan_trim[-which(wtko_scan_trim$log2FoldChange %in% outliers),]
  }
  correlations[x] <- cor(MF_dmrt1S_trim[rownames,'log2FoldChange'],
                         wtko_scan_trim[rownames,'log2FoldChange'], 
                         method = "pearson", use="pairwise")
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_dmrt1S_trim <- sex_related_MF_dmrt1S[SL_rownames,]
outliers <- boxplot(sex_related_MF_dmrt1S_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_dmrt1S_trim <- sex_related_MF_dmrt1S_trim[-which(sex_related_MF_dmrt1S_trim$log2FoldChange %in% outliers),]
}  
# remove outliers from wtko
sex_related_wtko_scan_trim <- sex_related_wtko_scan[SL_rownames,]
outliers <- boxplot(sex_related_wtko_scan_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_wtko_scan_trim <- sex_related_wtko_scan_trim[-which(sex_related_wtko_scan_trim$log2FoldChange %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_dmrt1S_trim[SL_rownames,'log2FoldChange'],
                          sex_related_wtko_scan_trim[SL_rownames,'log2FoldChange'], 
                          method = "pearson", use="pairwise")
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1] 0.3116883


# ccdc permutations ----

# MF_ccdc vs ccdc
correlations <- c()
# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 74, replace = F);indexes
  rownames <- counts$geneID[indexes]
  # remove outliers from MF
  MF_ccdc_trim <- MF_ccdc_unfiltered[rownames,]
  outliers <- boxplot(MF_ccdc_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_ccdc_trim<- MF_ccdc_trim[-which(MF_ccdc_trim$log2FoldChange %in% outliers),]
  }  
  # remove outliers from wtko
  wtko_ccdc_trim <- wtko_ccdc_unfiltered[rownames,]
  outliers <- boxplot(wtko_ccdc_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    wtko_ccdc_trim<- wtko_ccdc_trim[-which(wtko_ccdc_trim$log2FoldChange %in% outliers),]
  }
  correlations[x] <- cor(MF_ccdc_trim[rownames,'log2FoldChange'],
                         wtko_ccdc_trim[rownames,'log2FoldChange'], 
                         method = "pearson", use="pairwise")
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_ccdc_trim <- sex_related_MF_ccdc[SL_rownames,]
outliers <- boxplot(sex_related_MF_ccdc_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_ccdc_trim <- sex_related_MF_ccdc_trim[-which(sex_related_MF_ccdc_trim$log2FoldChange %in% outliers),]
}  
# remove outliers from wtko
sex_related_wtko_ccdc_trim <- sex_related_wtko_ccdc[SL_rownames,]
outliers <- boxplot(sex_related_wtko_ccdc_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_wtko_ccdc_trim <- sex_related_wtko_ccdc_trim[-which(sex_related_wtko_ccdc_trim$log2FoldChange %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_ccdc_trim[SL_rownames,'log2FoldChange'],
                          sex_related_wtko_ccdc_trim[SL_rownames,'log2FoldChange'], 
                          method = "pearson", use="pairwise")
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1] 0.962038


# MF_dmrt1L vs ccdc
correlations <- c()
# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 74, replace = F);indexes
  rownames <- counts$geneID[indexes]
  # remove outliers from MF
  MF_dmrt1L_trim <- MF_dmrt1L_unfiltered[rownames,]
  outliers <- boxplot(MF_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_dmrt1L_trim<- MF_dmrt1L_trim[-which(MF_dmrt1L_trim$log2FoldChange %in% outliers),]
  }  
  # remove outliers from wtko
  wtko_ccdc_trim <- wtko_ccdc_unfiltered[rownames,]
  outliers <- boxplot(wtko_ccdc_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    wtko_ccdc_trim<- wtko_ccdc_trim[-which(wtko_ccdc_trim$log2FoldChange %in% outliers),]
  }
  correlations[x] <- cor(MF_dmrt1L_trim[rownames,'log2FoldChange'],
                         wtko_ccdc_trim[rownames,'log2FoldChange'], 
                         method = "pearson", use="pairwise")
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_dmrt1L_trim <- sex_related_MF_dmrt1L[SL_rownames,]
outliers <- boxplot(sex_related_MF_dmrt1L_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_dmrt1L_trim <- sex_related_MF_dmrt1L_trim[-which(sex_related_MF_dmrt1L_trim$log2FoldChange %in% outliers),]
}  
# remove outliers from wtko
sex_related_wtko_ccdc_trim <- sex_related_wtko_ccdc[SL_rownames,]
outliers <- boxplot(sex_related_wtko_ccdc_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_wtko_ccdc_trim <- sex_related_wtko_ccdc_trim[-which(sex_related_wtko_ccdc_trim$log2FoldChange %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_dmrt1L_trim[SL_rownames,'log2FoldChange'],
                          sex_related_wtko_ccdc_trim[SL_rownames,'log2FoldChange'], 
                          method = "pearson", use="pairwise")
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1] 0.6083916

# MF_dmrt1S vs ccdc
correlations <- c()
# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 74, replace = F);indexes
  rownames <- counts$geneID[indexes]
  # remove outliers from MF
  MF_dmrt1S_trim <- MF_dmrt1S_unfiltered[rownames,]
  outliers <- boxplot(MF_dmrt1S_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_dmrt1S_trim<- MF_dmrt1S_trim[-which(MF_dmrt1S_trim$log2FoldChange %in% outliers),]
  }  
  # remove outliers from wtko
  wtko_ccdc_trim <- wtko_ccdc_unfiltered[rownames,]
  outliers <- boxplot(wtko_ccdc_trim$log2FoldChange, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    wtko_ccdc_trim<- wtko_ccdc_trim[-which(wtko_ccdc_trim$log2FoldChange %in% outliers),]
  }
  correlations[x] <- cor(MF_dmrt1S_trim[rownames,'log2FoldChange'],
                         wtko_ccdc_trim[rownames,'log2FoldChange'], 
                         method = "pearson", use="pairwise")
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_dmrt1S_trim <- sex_related_MF_dmrt1S[SL_rownames,]
outliers <- boxplot(sex_related_MF_dmrt1S_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_dmrt1S_trim <- sex_related_MF_dmrt1S_trim[-which(sex_related_MF_dmrt1S_trim$log2FoldChange %in% outliers),]
}  
# remove outliers from wtko
sex_related_wtko_ccdc_trim <- sex_related_wtko_ccdc[SL_rownames,]
outliers <- boxplot(sex_related_wtko_ccdc_trim$log2FoldChange, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_wtko_ccdc_trim <- sex_related_wtko_ccdc_trim[-which(sex_related_wtko_ccdc_trim$log2FoldChange %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_dmrt1S_trim[SL_rownames,'log2FoldChange'],
                          sex_related_wtko_ccdc_trim[SL_rownames,'log2FoldChange'], 
                          method = "pearson", use="pairwise")
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1] 0.8811189

```
