# Analysis of differential expression
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


setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_Supergene/2022_KO_tad_RNAseq/2022_EdgeR_and_DeSeq2/STAR_done")
dir <- "/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_Supergene/2022_KO_tad_RNAseq/2022_EdgeR_and_DeSeq2/STAR_done"
list.files(dir)


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
# samples
samples <- read.table(file.path(dir, "samples.txt"), header = T)
samples
# sexez
sexez <- factor(samples$sex)
sexez <- relevel(sexez, ref="F")
# batch
batchez <- factor(samples$batch)

# MF dmrt1S ----
colnames(countz)
new_counts <- as.data.frame(countz[,-c(1:35,37,40,43,45,46:66) ])
row.names(new_counts) <- gene_names
new_samples <-as.data.frame(samples[-c(1:35,37,40,43,45,46:66), ]);new_samples
new_sexez <- factor(samples$sex[-c(1:35,37,40,43,45,46:66)])
new_sexez <- relevel(new_sexez, ref="F")
# batch
new_batchez <- factor(samples$batch[-c(1:35,37,40,43,45,46:66)])
# Create DGEList object - this is a data structure that is used for 
# the analysis of differential expression
d0 <- DGEList(new_counts, group = new_sexez, remove.zeros = TRUE)
dim(d0$counts) #each row is a transcript - here is the number before filtering
# [1] 34261     6 # 2023 STAR ccdc dmrt1L no dmrt1S

# save the unfiltered logFC to a dataframe 
d0 <- calcNormFactors(d0, method="TMM")
design <- model.matrix(~ 0 + new_sexez, data=d0$samples) # last coefficient = difference between sexes)
d0 <- estimateDisp(d0, design, robust=TRUE)
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
MF_dmrt1S_unfiltered <- exacttest$table;MF_dmrt1S_unfiltered
# Write sex_related to a file
sex_related_MF_dmrt1S <- data.frame(exacttest$table[c('XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274'),])
write.csv(sex_related_MF_dmrt1S, file="Sex_related_MF_dmrt1S_STAR_edgeR_unfiltered.csv", row.names = T)
# Write counts of sex related to a file
sex_related_MF_dmrt1S_counts <- new_counts[c('XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274'), ]
write.csv(sex_related_MF_dmrt1S_counts, file="Sex_related_MF_dmrt1S_STAR_edgeR_counts_unfiltered.csv", row.names = T)

# Now do analysis of differential expression; 
# here we remove transcripts where the average count per sample is 2 or less:
d0$counts <- d0$counts[rowSums(d0$counts)> 2* ncol(d0$counts),] 
# Now we have far fewer transcripts:
dim(d0$counts)
# [1] 28918     6 # 2023 STAR ccdc dmrt1L no dmrt1S
# many rows with low expression were eliminated
# TMM normalization is applied to this dataset to account for compositional difference between
# the libraries.
d0 <- calcNormFactors(d0, method="TMM")
# check the normalization factors
d0$samples
# plot by sex
#plotMDS(d0,labels=c(rep("M",3),rep("F",3)),
#        col=c(rep("green",3),rep("blue",3)))
# design matrix: this is used for the model of differential expression
design <- model.matrix(~ 0 + new_sexez, data=d0$samples) # last coefficient = difference between sexes)
design
# estimate dispersion
d0 <- estimateDisp(d0, design, robust=TRUE)
#d0$common.dispersion
#d0$tagwise.dispersion
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
summary(decideTests(object = exacttest, p.value = 0.1))
# M-F
# Down     209
# NotSig 28443
# Up       266
dmrt1S_DE <- as.data.frame(topTags(exacttest, n=475))
write.csv(dmrt1S_DE, file="MF_STAR_dmrt1S_DE_edgeR.csv", row.names = T)

# pairwise scatterplot
# This is how you get the normalized counts in edgeR (https://support.bioconductor.org/p/44300/#44301)
effective.lib.size <- d0$samples$lib.size * d0$samples$norm.factors
normalized_counts <- log2( t(t(d0$counts+0.5) / (effective.lib.size+0.5)) )

# get Rsquare value for all pairwise comparisons
rsquare <- data.frame(matrix(ncol = ncol(normalized_counts), nrow = ncol(normalized_counts)))
for(i in 1:(ncol(normalized_counts)-1)) {       # for-loop over columns
  for(j in (i+1):ncol(normalized_counts)) { 
    print(paste(i," ",j))
    x <- cor.test(normalized_counts[ , i], 
                  normalized_counts[ , j], 
                  method = 'spearman')
    rsquare[i,j] <- x$estimate
  }
}
colnames(rsquare) <- colnames(normalized_counts)
rownames(rsquare) <- colnames(normalized_counts)

View(rsquare)
write_xlsx(rsquare, "./MF_rsquare_dmrt1S_STAR_edgeR.xls")

# negative logFC indicates higher expression in wt because wt is the
# denominator (reference)

#new_counts[c('XBXL10_1g7999','XBXL10_1g10668','XBXL10_1g34625'),]
#exacttest$table[c('XBXL10_1g7999','XBXL10_1g10668','XBXL10_1g34625'),]
#female related
# View(exacttest$table[c('XBXL10_1g34625','XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g29076','XBXL10_1g30057','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g33473','XBXL10_1g35876','XBXL10_1g37293','XBXL10_1g5748','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8966'),])
# male related
# View(exacttest$table[c('XBXL10_1g2070','XBXL10_1g4848','XBXL10_1g11002','XBXL10_1g30252','XBXL10_1g32546','XBXL10_1g605','XBXL10_1g3639','XBXL10_1g37486'),])
# sox genes
# View(exacttest$table[c('XBXL10_1g39526','XBXL10_1g42722','XBXL10_1g35158','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g42158','XBXL10_1g39443','XBXL10_1g42662'),])
# testis differentiation
# View(exacttest$table[c('XBXL10_1g41173','XBXL10_1g43880','XBXL10_1g19698','XBXL10_1g22028','XBXL10_1g815','XBXL10_1g3800','XBXL10_1g8007','XBXL10_1g2154','XBXL10_1g4928','XBXL10_1g27310','XBXL10_1g29128','XBXL10_1g40425','XBXL10_1g43291','XBXL10_1g8118','XBXL10_1g10760','XBXL10_1g8117','XBXL10_1g10758','XBXL10_1g1634','XBXL10_1g4460'),])
# steroidogenic genes
# View(exacttest$table[c('XBXL10_1g22534','XBXL10_1g25047','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g30377','XBXL10_1g13205','XBXL10_1g6054','XBXL10_1g29226','XBXL10_1g34871'),])



# MF dmrt1L ----
colnames(countz)
new_counts <- as.data.frame(countz[,-c(1:15,18:23,26:27,31:33,35:66) ])
row.names(new_counts) <- gene_names
new_samples <-as.data.frame(samples[-c(1:15,18:23,26:27,31:33,35:66), ]);new_samples
new_sexez <- factor(samples$sex[-c(1:15,18:23,26:27,31:33,35:66)])
new_sexez <- relevel(new_sexez, ref="F")
# batch
new_batchez <- factor(samples$batch[-c(1:15,18:23,26:27,31:33,35:66)])
# Create DGEList object - this is a data structure that is used for 
# the analysis of differential expression
d0 <- DGEList(new_counts, group = new_sexez, remove.zeros = TRUE)
dim(d0$counts) #each row is a transcript - here is the number before filtering
# [1] 35000     8 # 2023 STAR ccdc dmrt1L no dmrt1S

# save the unfiltered logFC to a dataframe 
d0 <- calcNormFactors(d0, method="TMM")
design <- model.matrix(~ 0 + new_sexez, data=d0$samples) # last coefficient = difference between sexes)
d0 <- estimateDisp(d0, design, robust=TRUE)
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
MF_dmrt1L_unfiltered <- exacttest$table;MF_dmrt1L_unfiltered
# Write sex_related to a file
sex_related_MF_dmrt1L <- data.frame(exacttest$table[c('XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274'),])
write.csv(sex_related_MF_dmrt1L, file="Sex_related_MF_dmrt1L_STAR_edgeR_unfiltered.csv", row.names = T)
# Write counts of sex related to a file
sex_related_MF_dmrt1L_counts <- new_counts[c('XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274'), ]
write.csv(sex_related_MF_dmrt1L_counts, file="Sex_related_MF_dmrt1L_STAR_edgeR_counts_unfiltered.csv", row.names = T)

# Now do analysis of differential expression; 
# here we remove transcripts where the average count per sample is 2 or less:
d0$counts <- d0$counts[rowSums(d0$counts)> 2* ncol(d0$counts),] 
# Now we have far fewer transcripts:
dim(d0$counts)
# [1] 28721     8 # 2023 STAR ccdc dmrt1L no dmrt1S
# many rows with low expression were eliminated
# TMM normalization is applied to this dataset to account for compositional difference between
# the libraries.
d0 <- calcNormFactors(d0, method="TMM")
# check the normalization factors
d0$samples
# plot by sex
#plotMDS(d0,labels=c(rep("M",2),"F",rep("M",2),rep("F",2),"M"),
#        col=c(rep("blue",2),"red",rep("blue",2),rep("red",2),"blue"))
# design matrix: this is used for the model of differential expression
design <- model.matrix(~ 0 + new_sexez, data=d0$samples) # last coefficient = difference between sexes)
design
# estimate dispersion
d0 <- estimateDisp(d0, design, robust=TRUE)
#d0$common.dispersion
#d0$tagwise.dispersion
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
summary(decideTests(object = exacttest, p.value = 0.1))
# M-F
# Down      26
# NotSig 28689
# Up         6
topTags(exacttest, n=32)
# Coefficient:  new_sexezm 
# logFC   logCPM        F       PValue        FDR
# XBXL10_1g8729  -8.261778  0.05457581 4.300636e-18 1.235186e-13  # this is dmw
# XBXL10_1g31107 -7.306227 -0.24264218 1.651868e-07 2.372166e-03
# XBXL10_1g28440 -6.013350 -1.39124555 2.916002e-07 2.791683e-03
# XBXL10_1g28511 -9.237984  1.04595625 4.961647e-07 3.562587e-03
# XBXL10_1g30091 -7.778468 -1.65769032 2.503454e-06 1.438034e-02
# XBXL10_1g44212 -3.423333  2.81516140 3.149509e-06 1.507617e-02
# XBXL10_1g3968   2.446234  2.66603804 5.028324e-06 2.063121e-02
# XBXL10_1g30717 -2.881452  1.06159385 6.238864e-06 2.119850e-02
# XBXL10_1g8948  -5.190441  0.08807785 6.642753e-06 2.119850e-02
# XBXL10_1g22487  3.917060 -0.11597090 8.623697e-06 2.318810e-02
# XBXL10_1g28385 -2.788075  2.63915640 8.880928e-06 2.318810e-02
# XBXL10_1g33085 -5.537669 -1.26488551 9.711802e-06 2.324439e-02
# XBXL10_1g32966 -4.777311  0.91243053 1.067252e-05 2.357889e-02
# XBXL10_1g25989 -5.657178 -1.78668504 1.168029e-05 2.396212e-02
# XBXL10_1g26938 -7.245238 -2.09741196 1.432764e-05 2.591665e-02
# XBXL10_1g11202 -9.159648 -0.43501781 1.458572e-05 2.591665e-02
# XBXL10_1g6540  -4.020248  2.90823349 1.534010e-05 2.591665e-02
# XBXL10_1g23366 -8.709820 -0.84574163 2.216987e-05 3.359738e-02
# XBXL10_1g15042 -2.328457  3.39615396 2.222591e-05 3.359738e-02
# XBXL10_1g37344 -2.947692  8.78789893 2.626404e-05 3.771648e-02
# XBXL10_1g7963  -3.141757  1.65351506 3.476186e-05 4.754264e-02
# XBXL10_1g31213 -5.471000  1.02750585 3.752027e-05 4.898271e-02
# XBXL10_1g24521 -7.557744 -1.84599713 4.059185e-05 5.068863e-02
# XBXL10_1g7804  -2.459094  1.29238280 4.537734e-05 5.430345e-02
# XBXL10_1g32518  2.429673  0.82211029 5.201557e-05 5.773739e-02
# XBXL10_1g17575  4.625349  1.69631694 5.360835e-05 5.773739e-02
# XBXL10_1g40319 -1.643922  4.25911016 5.427769e-05 5.773739e-02
# XBXL10_1g17518 -6.968335 -2.30556537 5.714516e-05 5.861665e-02
# XBXL10_1g5635   7.879352 -0.95229338 6.793244e-05 6.727888e-02
# XBXL10_1g8098  -6.946950 -2.32856224 8.528415e-05 8.164820e-02
# XBXL10_1g17574  3.632461  1.24137345 9.512396e-05 8.563408e-02
# XBXL10_1g44213 -4.075446  2.09771627 9.541069e-05 8.563408e-02
# XBXL10_1g43578  5.353484 -0.34327546 1.313532e-04 1.127739e-01
dmrt1L_DE <- as.data.frame(topTags(exacttest, n=32))
write.csv(dmrt1L_DE, file="MF_STAR_dmrt1L_DE_edgeR.csv", row.names = T)

# Let's try to get the logFC for all the genes in MF dmrt1L 
# that are sigDE in dmrt1S
is.data.frame(exacttest$table)
# first read in the file from the dmrt1S
dmrt1_DE_list <- read.csv(file.path(dir, "MF_STAR_dmrt1S_DE_edgeR.csv"), header = T)
dmrt1_DE_list_names <- dmrt1_DE_list$X
write.csv(exacttest$table[dmrt1_DE_list_names,], file="MF_STAR_edgeR_dmrt1Lexpression_of_dmrt1S_DEs.csv", row.names = T)


# pairwise scatterplot
# This is how you get the normalized counts in edgeR (https://support.bioconductor.org/p/44300/#44301)
effective.lib.size <- d0$samples$lib.size * d0$samples$norm.factors
normalized_counts <- log2( t(t(d0$counts+0.5) / (effective.lib.size+0.5)) )

# get Rsquare value for all pairwise comparisons
rsquare <- data.frame(matrix(ncol = ncol(normalized_counts), nrow = ncol(normalized_counts)))
for(i in 1:(ncol(normalized_counts)-1)) {       # for-loop over columns
  for(j in (i+1):ncol(normalized_counts)) { 
    print(paste(i," ",j))
    x <- cor.test(normalized_counts[ , i], 
                  normalized_counts[ , j], 
                  method = 'spearman')
    rsquare[i,j] <- x$estimate
  }
}
colnames(rsquare) <- colnames(normalized_counts)
rownames(rsquare) <- colnames(normalized_counts)

View(rsquare)
library(writexl)
write_xlsx(rsquare, "./MF_rsquare_dmrt1L_STAR_edgeR.xls")
# all are above 0.8.

# MF ccdc ----
colnames(countz)
new_counts <- as.data.frame(countz[,-c(4,8:9,11:13,15:66) ])
row.names(new_counts) <- gene_names
new_samples <-as.data.frame(samples[-c(4,8:9,11:13,15:66), ]);new_samples
new_sexez <- factor(samples$sex[-c(4,8:9,11:13,15:66)]);new_sexez
new_sexez <- relevel(new_sexez, ref="F")
# batch
new_batchez <- factor(samples$batch[-c(4,8:9,11:13,15:66)])
# Create DGEList object - this is a data structure that is used for 
# the analysis of differential expression
d0 <- DGEList(new_counts, group = new_sexez, remove.zeros = TRUE)
dim(d0$counts) #each row is a transcript - here is the number before filtering
# [1] 34467     8 # 2023 STAR ccdc only

# save the unfiltered logFC to a dataframe 
d0 <- calcNormFactors(d0, method="TMM")
design <- model.matrix(~ 0 + new_sexez, data=d0$samples) # last coefficient = difference between sexes)
d0 <- estimateDisp(d0, design, robust=TRUE)
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
MF_ccdc_unfiltered <- exacttest$table;MF_ccdc_unfiltered
# Write sex_related to a file
sex_related_MF_ccdc <- data.frame(exacttest$table[c('XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274'),])
write.csv(sex_related_MF_ccdc, file="Sex_related_MF_ccdc_STAR_edgeR_unfiltered.csv", row.names = T)
# Write counts of sex related to a file
sex_related_MF_ccdc_counts <- new_counts[c('XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274'), ]
write.csv(sex_related_MF_ccdc_counts, file="Sex_related_MF_dmrt1L_STAR_edgeR_counts_unfiltered.csv", row.names = T)

# Now do analysis of differential expression; 
# here we remove transcripts where the average count per sample is 2 or less:
d0$counts <- d0$counts[rowSums(d0$counts)> 2* ncol(d0$counts),] 
# Now we have far fewer transcripts:
dim(d0$counts)
# [1] 29073     8 # 2023 STAR ccdc only
# many rows with low expression were eliminated
# TMM normalization is applied to this dataset to account for compositional difference between
# the libraries.
d0 <- calcNormFactors(d0, method="TMM")
# check the normalization factors
d0$samples
# plot by sex
#plotMDS(d0,labels=c(rep("M",5),rep("F",2),"M"),
#        col=c(rep("blue",5),rep("red",2),"blue"))
# design matrix: this is used for the model of differential expression
design <- model.matrix(~ 0 + new_sexez, data=d0$samples) # last coefficient = difference between sexes)
design
# estimate dispersion
d0 <- estimateDisp(d0, design, robust=TRUE)
#d0$common.dispersion
#d0$tagwise.dispersion
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
summary(decideTests(object = exacttest, p.value = 0.1))
# M-F
# Down       2
# NotSig 29070
#  Up         1
topTags(exacttest, n=3)
# Coefficient:  new_sexezm 
# logFC   logCPM        F       PValue        FDR
# XBXL10_1g3050  -7.520457 -0.6587820 1.307588e-06 0.01585958
# XBXL10_1g43040 -8.439321 -1.6729722 1.477786e-06 0.01585958
# XBXL10_1g3473   9.507057  0.7988971 1.636526e-06 0.01585958
ccdc_DE <- as.data.frame(topTags(exacttest, n=3))
write.csv(ccdc_DE, file="MF_STAR_ccdc_DE_edgeR.csv", row.names = T)


# pairwise scatterplot
# This is how you get the normalized counts in edgeR (https://support.bioconductor.org/p/44300/#44301)
effective.lib.size <- d0$samples$lib.size * d0$samples$norm.factors
normalized_counts <- log2( t(t(d0$counts+0.5) / (effective.lib.size+0.5)) )

# get Rsquare value for all pairwise comparisons
rsquare <- data.frame(matrix(ncol = ncol(normalized_counts), nrow = ncol(normalized_counts)))
for(i in 1:(ncol(normalized_counts)-1)) {       # for-loop over columns
  for(j in (i+1):ncol(normalized_counts)) { 
    print(paste(i," ",j))
    x <- cor.test(normalized_counts[ , i], 
                  normalized_counts[ , j], 
                  method = 'spearman')
    rsquare[i,j] <- x$estimate
     }
}
colnames(rsquare) <- colnames(normalized_counts)
rownames(rsquare) <- colnames(normalized_counts)

View(rsquare)
library(writexl)
write_xlsx(rsquare, "./MF_rsquare_ccdc_STAR_edgeR.xls")
## all comparisons are above 0.8



# wtko dmw ----
colnames(countz)
new_counts <- as.data.frame(countz[,c(46:57) ])
row.names(new_counts) <- gene_names
new_samples <-as.data.frame(samples[c(46:57), ]);new_samples
new_genotypez <- factor(samples$genotype[c(46:57)]);new_genotypez
new_genotypez <- relevel(new_genotypez, ref="wt")
# Create DGEList object - this is a data structure that is used for 
# the analysis of differential expression
d0 <- DGEList(new_counts, group = new_genotypez, remove.zeros = TRUE)
dim(d0$counts) #each row is a transcript - here is the number before filtering
# [1] 36327    12 # 2023 STAR dmw only

# save the unfiltered logFC to a dataframe 
d0 <- calcNormFactors(d0, method="TMM")
design <- model.matrix(~ 0 + new_genotypez, data=d0$samples) # last coefficient = difference between sexes)
d0 <- estimateDisp(d0, design, robust=TRUE)
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
wtko_dmw_unfiltered <- exacttest$table;wtko_dmw_unfiltered
# Write sex_related to a file
sex_related_wtko_dmw <- data.frame(exacttest$table[c('XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274'),])
write.csv(sex_related_wtko_dmw, file="Sex_related_wtko_dmw_STAR_edgeR_unfiltered.csv", row.names = T)
# Write counts of sex related to a file
sex_related_wtko_dmw_counts <- new_counts[c('XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274'), ]
write.csv(sex_related_wtko_dmw_counts, file="Sex_related_wtko_dmw_STAR_edgeR_counts_unfiltered.csv", row.names = T)

# Now do analysis of differential expression; 
# here we remove transcripts where the average count per sample is 2 or less:
d0$counts <- d0$counts[rowSums(d0$counts)> 2* ncol(d0$counts),] 
# Now we have far fewer transcripts:
dim(d0$counts)
# [1] 29051    12 # 2023 STAR dmw only
# many rows with low expression were eliminated
# TMM normalization is applied to this dataset to account for compositional difference between
# the libraries.
d0 <- calcNormFactors(d0, method="TMM")
# check the normalization factors
d0$samples
# plot by sex
#plotMDS(d0,labels=c("wt","ko","wt","ko","wt","wt","wt","ko","ko","ko","ko","wt"),
#        col=c("blue","green","blue","green","blue","blue","blue","green","green","green","green","blue"))
# design matrix: this is used for the model of differential expression
design <- model.matrix(~ 0 + new_genotypez, data=d0$samples) # last coefficient = difference between sexes)
design
# estimate dispersion
d0 <- estimateDisp(d0, design, robust=TRUE)
#d0$common.dispersion
#d0$tagwise.dispersion
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
summary(decideTests(object = exacttest, p.value = 0.1))
#        ko-wt
# Down       0
# NotSig 29043
# Up         8
topTags(exacttest, n=8)
# Comparison of groups:  ko-wt 
# logFC    logCPM       PValue        FDR
# XBXL10_1g15465 5.447262 -1.995802 5.548449e-07 0.01045141
# XBXL10_1g7424  3.277253  2.152194 7.195213e-07 0.01045141
# XBXL10_1g3709  4.585148 -1.764584 4.462500e-06 0.04321336
# XBXL10_1g4777  1.827103  1.842327 9.498220e-06 0.06898320
# XBXL10_1g26938 6.427278 -2.296535 1.401094e-05 0.08140634
# XBXL10_1g19716 1.860419  1.673926 2.168770e-05 0.09371137
# XBXL10_1g10773 1.655619  1.913439 2.258028e-05 0.09371137
# XBXL10_1g40801 4.648535 -2.594841 2.639426e-05 0.09584746
dmw_DE <- as.data.frame(topTags(exacttest, n=8))
write.csv(dmw_DE, file="wtko_STAR_dmw_DE_edgeR.csv", row.names = T)
# how does expression of dmw look?
new_counts[c('XBXL10_1g8729'),]
#how about dmrt1L (XBXL10_1g2070)
new_counts[c('XBXL10_1g2070'),]
exacttest$table['XBXL10_1g2070',]
#how about dmrt1S (XBXL10_1g4848)
new_counts[c('XBXL10_1g4848'),]
exacttest$table['XBXL10_1g4848',]


# negative logFC indicates higher expression in wt because wt is the
# denominator (reference)
#female related
new_counts[c('XBXL10_1g7999','XBXL10_1g10668','XBXL10_1g34625'),]
exacttest$table[c('XBXL10_1g7999','XBXL10_1g10668','XBXL10_1g34625'),]
View(exacttest$table[c('XBXL10_1g34625','XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g29076','XBXL10_1g30057','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g33473','XBXL10_1g35876','XBXL10_1g37293','XBXL10_1g5748','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8966'),])
# logFC of female related genes is mostly positive, suggesting
# higher expression in ko. This could be due to the female pathway
# being associated with inhibition of male related genes (that
# are thus not inhibited in the ko)

# male related
View(exacttest$table[c('XBXL10_1g2070','XBXL10_1g4848','XBXL10_1g11002','XBXL10_1g30252','XBXL10_1g32546','XBXL10_1g605','XBXL10_1g3639','XBXL10_1g37486'),])
# logFC of male related genes is ~50/50 positive and negative

# sox genes
View(exacttest$table[c('XBXL10_1g39526','XBXL10_1g42722','XBXL10_1g35158','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g42158','XBXL10_1g39443','XBXL10_1g42662'),])
# logFC of sox genes is mostly positive

# testis differentiation
View(exacttest$table[c('XBXL10_1g41173','XBXL10_1g43880','XBXL10_1g19698','XBXL10_1g22028','XBXL10_1g815','XBXL10_1g3800','XBXL10_1g8007','XBXL10_1g2154','XBXL10_1g4928','XBXL10_1g27310','XBXL10_1g29128','XBXL10_1g40425','XBXL10_1g43291','XBXL10_1g8118','XBXL10_1g10760','XBXL10_1g8117','XBXL10_1g10758','XBXL10_1g1634','XBXL10_1g4460'),])
# logFC of male related genes is ~60/40 positive and negative
# this is somewhat consistent with dmw inhibiting male related genes

# steroidogenic genes
View(exacttest$table[c('XBXL10_1g22534','XBXL10_1g25047','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g30377','XBXL10_1g13205','XBXL10_1g6054','XBXL10_1g29226','XBXL10_1g34871'),])
# logFC of steriodogenic genes is ~50/50 positive and negative


# Now let's get the logFC for all the genes in dmw wt_ko
# that are sigDE in M vs F dmrt1S
# first read in the file from the dmrt1S
dmrt1_DE_list <- read.csv(file.path(dir, "MF_STAR_dmrt1S_DE_edgeR.csv"), header = T)
dmrt1_DE_list_names <- dmrt1_DE_list$X
write.csv(exacttest$table[dmrt1_DE_list_names,], file="wt_ko_STAR_edgeR_dmw_expression_of_dmrt1S_MF_DEs.csv", row.names = T)

# pairwise scatterplot
# This is how you get the normalized counts in edgeR (https://support.bioconductor.org/p/44300/#44301)
effective.lib.size <- d0$samples$lib.size * d0$samples$norm.factors
normalized_counts <- log2( t(t(d0$counts+0.5) / (effective.lib.size+0.5)) )

# get Rsquare value for all pairwise comparisons
rsquare <- data.frame(matrix(ncol = ncol(normalized_counts), nrow = ncol(normalized_counts)))
for(i in 1:(ncol(normalized_counts)-1)) {       # for-loop over columns
  for(j in (i+1):ncol(normalized_counts)) { 
    print(paste(i," ",j))
    x <- cor.test(normalized_counts[ , i], 
                  normalized_counts[ , j], 
                  method = 'spearman')
    rsquare[i,j] <- x$estimate
  }
}
colnames(rsquare) <- colnames(normalized_counts)
rownames(rsquare) <- colnames(normalized_counts)

View(rsquare)
library(writexl)
write_xlsx(rsquare, "./wt_ko_rsquare_dmw_STAR_edgeR.xls")
# all comparisons are above 0.8


# wtko scan ----
colnames(countz)
new_counts <- as.data.frame(countz[,c(58:66) ])
row.names(new_counts) <- gene_names
new_samples <-as.data.frame(samples[c(58:66), ]);new_samples
new_genotypez <- factor(samples$genotype[c(58:66)]);new_genotypez
new_genotypez <- relevel(new_genotypez, ref="wt")
# Create DGEList object - this is a data structure that is used for 
# the analysis of differential expression
d0 <- DGEList(new_counts, group = new_genotypez, remove.zeros = TRUE)
dim(d0$counts) #each row is a transcript - here is the number before filtering
# [1] 35697     9 # 2023 STAR scan only

# save the unfiltered logFC to a dataframe 
d0 <- calcNormFactors(d0, method="TMM")
design <- model.matrix(~ 0 + new_genotypez, data=d0$samples) # last coefficient = difference between sexes)
d0 <- estimateDisp(d0, design, robust=TRUE)
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
wtko_scan_unfiltered <- exacttest$table;wtko_scan_unfiltered
# Write sex_related to a file
sex_related_wtko_scan <- data.frame(exacttest$table[c('XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274'),])
write.csv(sex_related_wtko_scan, file="Sex_related_wtko_scan_STAR_edgeR_unfiltered.csv", row.names = T)
# Write counts of sex related to a file
sex_related_wtko_scan_counts <- new_counts[c('XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274'), ]
write.csv(sex_related_wtko_scan_counts, file="Sex_related_wtko_scan_STAR_edgeR_counts_unfiltered.csv", row.names = T)

# Now do analysis of differential expression; 
# here we remove transcripts where the average count per sample is 2 or less:
d0$counts <- d0$counts[rowSums(d0$counts)> 2* ncol(d0$counts),] 
# Now we have far fewer transcripts:
dim(d0$counts)
# [1] 29479     9 # 2023 STAR scan
# many rows with low expression were eliminated
# TMM normalization is applied to this dataset to account for compositional difference between
# the libraries.
d0 <- calcNormFactors(d0, method="TMM")
# check the normalization factors
d0$samples
# plot by sex
#plotMDS(d0,labels=c("wt","ko","wt","ko","wt","wt","wt","ko","ko","ko","ko","wt"),
#        col=c("blue","green","blue","green","blue","blue","blue","green","green","green","green","blue"))
# design matrix: this is used for the model of differential expression
design <- model.matrix(~ 0 + new_genotypez, data=d0$samples) # last coefficient = difference between sexes)
#design
# estimate dispersion
d0 <- estimateDisp(d0, design, robust=TRUE)
#d0$common.dispersion
#d0$tagwise.dispersion
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
summary(decideTests(object = exacttest, p.value = 0.1))
#        ko-wt
# Down      11
# NotSig 29462
# Up         6
topTags(exacttest, n=17)
scan_DE <- as.data.frame(topTags(exacttest, n=17))
write.csv(scan_DE, file="wtko_STAR_scan_DE_edgeR.csv", row.names = T)

# pairwise scatterplot
# This is how you get the normalized counts in edgeR (https://support.bioconductor.org/p/44300/#44301)
effective.lib.size <- d0$samples$lib.size * d0$samples$norm.factors
normalized_counts <- log2( t(t(d0$counts+0.5) / (effective.lib.size+0.5)) )

# get Rsquare value for all pairwise comparisons
rsquare <- data.frame(matrix(ncol = ncol(normalized_counts), nrow = ncol(normalized_counts)))
for(i in 1:(ncol(normalized_counts)-1)) {       # for-loop over columns
  for(j in (i+1):ncol(normalized_counts)) { 
    print(paste(i," ",j))
    x <- cor.test(normalized_counts[ , i], 
                  normalized_counts[ , j], 
                  method = 'spearman')
    rsquare[i,j] <- x$estimate
  }
}
colnames(rsquare) <- colnames(normalized_counts)
rownames(rsquare) <- colnames(normalized_counts)

View(rsquare)
write_xlsx(rsquare, "./wt_ko_rsquare_scan_STAR_edgeR.xls")


# wtko ccdc ----
colnames(countz)
new_counts <- as.data.frame(countz[,c(4,7:13) ])
row.names(new_counts) <- gene_names
new_samples <-as.data.frame(samples[c(4,7:13), ]);new_samples
new_genotypez <- factor(samples$genotype[c(4,7:13)]);new_genotypez
new_genotypez <- relevel(new_genotypez, ref="wt")
# Create DGEList object - this is a data structure that is used for 
# the analysis of differential expression
d0 <- DGEList(new_counts, group = new_genotypez, remove.zeros = TRUE)
dim(d0$counts) #each row is a transcript - here is the number before filtering
# [1] 34358     8 # 2023 STAR ccdc only

# save the unfiltered logFC to a dataframe 
d0 <- calcNormFactors(d0, method="TMM")
design <- model.matrix(~ 0 + new_genotypez, data=d0$samples) # last coefficient = difference between sexes)
d0 <- estimateDisp(d0, design, robust=TRUE)
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
wtko_ccdc_unfiltered <- exacttest$table;wtko_ccdc_unfiltered
# Write sex_related to a file
sex_related_wtko_ccdc <- data.frame(exacttest$table[c('XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274'),])
write.csv(sex_related_wtko_ccdc, file="Sex_related_wtko_ccdc_STAR_edgeR_unfiltered.csv", row.names = T)
# Write counts of sex related to a file
sex_related_wtko_ccdc_counts <- new_counts[c('XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274'), ]
write.csv(sex_related_wtko_ccdc_counts, file="Sex_related_wtko_ccdc_STAR_edgeR_counts_unfiltered.csv", row.names = T)

# Now do analysis of differential expression; 
# here we remove transcripts where the average count per sample is 2 or less:
d0$counts <- d0$counts[rowSums(d0$counts)> 2* ncol(d0$counts),] 
# Now we have far fewer transcripts:
dim(d0$counts)
# [1] 29101     8 # 2023 STAR ccdc only
# many rows with low expression were eliminated
# TMM normalization is applied to this dataset to account for compositional difference between
# the libraries.
d0 <- calcNormFactors(d0, method="TMM")
# check the normalization factors
d0$samples
# plot by sex
#plotMDS(d0,labels=c("wt","ko","wt","ko","wt","wt","wt","ko","ko","ko","ko","wt"),
#        col=c("blue","green","blue","green","blue","blue","blue","green","green","green","green","blue"))
# design matrix: this is used for the model of differential expression
design <- model.matrix(~ 0 + new_genotypez, data=d0$samples) # last coefficient = difference between sexes)
#design
# estimate dispersion
d0 <- estimateDisp(d0, design, robust=TRUE)
#d0$common.dispersion
#d0$tagwise.dispersion
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
summary(decideTests(object = exacttest, p.value = 0.1))
#        ko-wt
# Down       8
# NotSig 29082
# Up        11
topTags(exacttest, n=19)
ccdc_DE <- as.data.frame(topTags(exacttest, n=19))
write.csv(ccdc_DE, file="wtko_STAR_ccdc_DE_edgeR.csv", row.names = T)

# pairwise scatterplot
# This is how you get the normalized counts in edgeR (https://support.bioconductor.org/p/44300/#44301)
effective.lib.size <- d0$samples$lib.size * d0$samples$norm.factors
normalized_counts <- log2( t(t(d0$counts+0.5) / (effective.lib.size+0.5)) )

# get Rsquare value for all pairwise comparisons
rsquare <- data.frame(matrix(ncol = ncol(normalized_counts), nrow = ncol(normalized_counts)))
for(i in 1:(ncol(normalized_counts)-1)) {       # for-loop over columns
  for(j in (i+1):ncol(normalized_counts)) { 
    print(paste(i," ",j))
    x <- cor.test(normalized_counts[ , i], 
                  normalized_counts[ , j], 
                  method = 'spearman')
    rsquare[i,j] <- x$estimate
  }
}
colnames(rsquare) <- colnames(normalized_counts)
rownames(rsquare) <- colnames(normalized_counts)

View(rsquare)
write_xlsx(rsquare, "./wt_ko_rsquare_ccdc_STAR_edgeR.xls")

# permutations ----

# OK now do some permutations to assess significance of the correlations
# between the log2FC of each wtko and each MF

# These permutations will randomly select 74 logFC from each MF and wtko
# 1000 times and calculate the correlation.  Then this will be compared to
# the observed

# get rownames of sexrelated transcripts
SL_rownames <- c('XBXL10_1g10089','XBXL10_1g10668','XBXL10_1g10675','XBXL10_1g10758','XBXL10_1g10760','XBXL10_1g11002','XBXL10_1g13205','XBXL10_1g13810','XBXL10_1g15286','XBXL10_1g15724','XBXL10_1g1634','XBXL10_1g19698','XBXL10_1g2070','XBXL10_1g2154','XBXL10_1g22028','XBXL10_1g22534','XBXL10_1g22535','XBXL10_1g23152','XBXL10_1g24241','XBXL10_1g24554','XBXL10_1g25046','XBXL10_1g25047','XBXL10_1g25243','XBXL10_1g26060','XBXL10_1g26280','XBXL10_1g27265','XBXL10_1g27310','XBXL10_1g29076','XBXL10_1g29128','XBXL10_1g29226','XBXL10_1g30057','XBXL10_1g30252','XBXL10_1g30377','XBXL10_1g31301','XBXL10_1g3211','XBXL10_1g32392','XBXL10_1g32546','XBXL10_1g33473','XBXL10_1g34625','XBXL10_1g34871','XBXL10_1g35158','XBXL10_1g35876','XBXL10_1g3639','XBXL10_1g37293','XBXL10_1g37486','XBXL10_1g37811','XBXL10_1g3800','XBXL10_1g38013','XBXL10_1g38893','XBXL10_1g39443','XBXL10_1g39526','XBXL10_1g40425','XBXL10_1g41173','XBXL10_1g42158','XBXL10_1g42662','XBXL10_1g42722','XBXL10_1g43291','XBXL10_1g43880','XBXL10_1g4460','XBXL10_1g4848','XBXL10_1g4928','XBXL10_1g5748','XBXL10_1g605','XBXL10_1g6054','XBXL10_1g6566','XBXL10_1g7278','XBXL10_1g7999','XBXL10_1g8007','XBXL10_1g8117','XBXL10_1g8118','XBXL10_1g815','XBXL10_1g8430','XBXL10_1g8966','XBXL10_1g9274')


# the two functions below were developed by Ian Dworkin and colleagues
# https://github.com/DworkinLab/Trypoxylus_RNAseq/blob/master/analysis_scripts/re_analysis_scripts_jan_2018.Rmd

# This function calculates the Euclidean Distances (the L2 norm) 
# so we can use unit vectors in the estimation of the vector 
# correlation.
# These functions follow the logic of Kuruvilla et al 2002, 
# and were adapted from Pitchers et al 2013
PD <- function(x) { 
  sqrt(t(x)%*%x)}

#this function gives the vector correlation and angle, 
# and vector magnitude ratio, alpha, between two vectors
# alpha of 1 means that the length of vector 2 is the same 
# as the length as vector 1. 
# lower than one means that it is smaller, greater than 
# 1 means that vector 2 is larger
ang.vec.alph <- function(vec1, vec2) {
  vec1 <- vec1 - mean(vec1)
  vec2 <- vec2 - mean(vec2)
  vec.cor <- abs((t(vec1) %*% vec2)/(PD(vec1)*PD(vec2)))
  vec.angle <- acos(vec.cor)*(180/pi)
  vec.alpha <- PD(vec1)/PD(vec2) 
  vec.ED <- PD(vec2-vec1) #Subtract vector one from vector two and then calculate PD for Euclidean Distance. 
  return(c(vector.cor=vec.cor, vec.angle=vec.angle, vec.alpha=vec.alpha, vector.ED=vec.ED))} 



# dmw permutations ----
# MF_ccdc vs dmw
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(counts)[1], 74, replace = F);indexes
  rownames <- counts$geneID[indexes]
  # remove outliers from MF
  MF_ccdc_trim <- MF_ccdc_unfiltered[rownames,]
  outliers <- boxplot(MF_ccdc_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_ccdc_trim<- MF_ccdc_trim[-which(MF_ccdc_trim$logFC %in% outliers),]
  }  
  # remove outliers from wtko
  wtko_dmw_trim <- wtko_dmw_unfiltered[rownames,]
  outliers <- boxplot(wtko_dmw_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    wtko_dmw_trim<- wtko_dmw_trim[-which(wtko_dmw_trim$logFC %in% outliers),]
  }
  correlations[x] <- cor(MF_ccdc_trim[rownames,'logFC'],
                         wtko_dmw_trim[rownames,'logFC'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(wtko_dmw_trim[,'logFC'], # ko:wt first
             MF_ccdc_trim[,'logFC'], # reference M:F second
             by = 'row.names', 
             incomparables = NA)
  b <- a[complete.cases(a), ];b
  magnitudes[x] <-ang.vec.alph(b$x,b$y)[3] # this is the ratio of the magnitudes of each vector
  # the reference (wildtype M:F is the denominator)
  # if this is >1 then the ko:wt has a bigger effect
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_ccdc_trim <- sex_related_MF_ccdc[SL_rownames,]
outliers <- boxplot(sex_related_MF_ccdc_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_ccdc_trim <- sex_related_MF_ccdc_trim[-which(sex_related_MF_ccdc_trim$logFC %in% outliers),]
}  
# remove outliers from wtko
sex_related_wtko_dmw_trim <- sex_related_wtko_dmw[SL_rownames,]
outliers <- boxplot(sex_related_wtko_dmw_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_wtko_dmw_trim <- sex_related_wtko_dmw_trim[-which(sex_related_wtko_dmw_trim$logFC %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_ccdc_trim[SL_rownames,'logFC'],
                          sex_related_wtko_dmw_trim[SL_rownames,'logFC'], 
                          method = "pearson", use="pairwise")
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1] 0.1168831

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_wtko_dmw_trim[SL_rownames,'logFC'],
           sex_related_MF_ccdc_trim[SL_rownames,'logFC'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.2027972


# MF_dmrt1L vs dmw
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(counts)[1], 74, replace = F);indexes
  rownames <- counts$geneID[indexes]
  # remove outliers from MF
  MF_dmrt1L_trim <- MF_dmrt1L_unfiltered[rownames,]
  outliers <- boxplot(MF_dmrt1L_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_dmrt1L_trim<- MF_dmrt1L_trim[-which(MF_dmrt1L_trim$logFC %in% outliers),]
  }  
  # remove outliers from wtko
  wtko_dmw_trim <- wtko_dmw_unfiltered[rownames,]
  outliers <- boxplot(wtko_dmw_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    wtko_dmw_trim<- wtko_dmw_trim[-which(wtko_dmw_trim$logFC %in% outliers),]
  }
  correlations[x] <- cor(MF_dmrt1L_trim[rownames,'logFC'],
                         wtko_dmw_trim[rownames,'logFC'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(wtko_dmw_trim[,'logFC'], # ko:wt first
             MF_dmrt1L_trim[,'logFC'], # reference M:F second
             by = 'row.names', 
             incomparables = NA)
  b <- a[complete.cases(a), ];b
  magnitudes[x] <-ang.vec.alph(b$x,b$y)[3] # this is the ratio of the magnitudes of each vector
  # the reference (wildtype M:F is the denominator)
  # if this is >1 then the ko:wt has a bigger effect
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_dmrt1L_trim <- sex_related_MF_dmrt1L[SL_rownames,]
outliers <- boxplot(sex_related_MF_dmrt1L_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_dmrt1L_trim <- sex_related_MF_dmrt1L_trim[-which(sex_related_MF_dmrt1L_trim$logFC %in% outliers),]
}  
# remove outliers from wtko
sex_related_wtko_dmw_trim <- sex_related_wtko_dmw[SL_rownames,]
outliers <- boxplot(sex_related_wtko_dmw_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_wtko_dmw_trim <- sex_related_wtko_dmw_trim[-which(sex_related_wtko_dmw_trim$logFC %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_dmrt1L_trim[SL_rownames,'logFC'],
                          sex_related_wtko_dmw_trim[SL_rownames,'logFC'], 
                          method = "pearson", use="pairwise")
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1] 0.5294705

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_wtko_dmw_trim[SL_rownames,'logFC'],
           sex_related_MF_dmrt1L_trim[SL_rownames,'logFC'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.1618382


# MF_dmrt1S vs dmw
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(counts)[1], 74, replace = F);indexes
  rownames <- counts$geneID[indexes]
  # remove outliers from MF
  MF_dmrt1S_trim <- MF_dmrt1S_unfiltered[rownames,]
  outliers <- boxplot(MF_dmrt1S_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_dmrt1S_trim<- MF_dmrt1S_trim[-which(MF_dmrt1S_trim$logFC %in% outliers),]
  }  
  # remove outliers from wtko
  wtko_dmw_trim <- wtko_dmw_unfiltered[rownames,]
  outliers <- boxplot(wtko_dmw_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    wtko_dmw_trim<- wtko_dmw_trim[-which(wtko_dmw_trim$logFC %in% outliers),]
  }
  correlations[x] <- cor(MF_dmrt1S_trim[rownames,'logFC'],
                         wtko_dmw_trim[rownames,'logFC'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(wtko_dmw_trim[,'logFC'], # ko:wt first
             MF_dmrt1S_trim[,'logFC'], # reference M:F second
             by = 'row.names', 
             incomparables = NA)
  b <- a[complete.cases(a), ];b
  magnitudes[x] <-ang.vec.alph(b$x,b$y)[3] # this is the ratio of the magnitudes of each vector
  # the reference (wildtype M:F is the denominator)
  # if this is >1 then the ko:wt has a bigger effect
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_dmrt1S_trim <- sex_related_MF_dmrt1S[SL_rownames,]
outliers <- boxplot(sex_related_MF_dmrt1S_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_dmrt1S_trim <- sex_related_MF_dmrt1S_trim[-which(sex_related_MF_dmrt1S_trim$logFC %in% outliers),]
}  
# remove outliers from wtko
sex_related_wtko_dmw_trim <- sex_related_wtko_dmw[SL_rownames,]
outliers <- boxplot(sex_related_wtko_dmw_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_wtko_dmw_trim <- sex_related_wtko_dmw_trim[-which(sex_related_wtko_dmw_trim$logFC %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_dmrt1S_trim[SL_rownames,'logFC'],
                          sex_related_wtko_dmw_trim[SL_rownames,'logFC'], 
                          method = "pearson", use="pairwise")
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1] 0.02597403

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_wtko_dmw_trim[SL_rownames,'logFC'],
           sex_related_MF_dmrt1S_trim[SL_rownames,'logFC'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.7352647


# scan permutations ----

# MF_ccdc vs scan
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(counts)[1], 74, replace = F);indexes
  rownames <- counts$geneID[indexes]
  # remove outliers from MF
  MF_ccdc_trim <- MF_ccdc_unfiltered[rownames,]
  outliers <- boxplot(MF_ccdc_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_ccdc_trim<- MF_ccdc_trim[-which(MF_ccdc_trim$logFC %in% outliers),]
  }  
  # remove outliers from wtko
  wtko_scan_trim <- wtko_scan_unfiltered[rownames,]
  outliers <- boxplot(wtko_scan_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    wtko_scan_trim<- wtko_scan_trim[-which(wtko_scan_trim$logFC %in% outliers),]
  }
  correlations[x] <- cor(MF_ccdc_trim[rownames,'logFC'],
                         wtko_scan_trim[rownames,'logFC'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(wtko_scan_trim[,'logFC'], # ko:wt first
             MF_ccdc_trim[,'logFC'], # reference M:F second
             by = 'row.names', 
             incomparables = NA)
  b <- a[complete.cases(a), ];b
  magnitudes[x] <-ang.vec.alph(b$x,b$y)[3] # this is the ratio of the magnitudes of each vector
  # the reference (wildtype M:F is the denominator)
  # if this is >1 then the ko:wt has a bigger effect
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_ccdc_trim <- sex_related_MF_ccdc[SL_rownames,]
outliers <- boxplot(sex_related_MF_ccdc_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_ccdc_trim <- sex_related_MF_ccdc_trim[-which(sex_related_MF_ccdc_trim$logFC %in% outliers),]
}  
# remove outliers from wtko
sex_related_wtko_scan_trim <- sex_related_wtko_scan[SL_rownames,]
outliers <- boxplot(sex_related_wtko_scan_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_wtko_scan_trim <- sex_related_wtko_scan_trim[-which(sex_related_wtko_scan_trim$logFC %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_ccdc_trim[SL_rownames,'logFC'],
                          sex_related_wtko_scan_trim[SL_rownames,'logFC'], 
                          method = "pearson", use="pairwise")
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1] 0.2827173

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_wtko_scan_trim[SL_rownames,'logFC'],
           sex_related_MF_ccdc_trim[SL_rownames,'logFC'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.2387612



# MF_dmrt1L vs scan
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(counts)[1], 74, replace = F);indexes
  rownames <- counts$geneID[indexes]
  # remove outliers from MF
  MF_dmrt1L_trim <- MF_dmrt1L_unfiltered[rownames,]
  outliers <- boxplot(MF_dmrt1L_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_dmrt1L_trim<- MF_dmrt1L_trim[-which(MF_dmrt1L_trim$logFC %in% outliers),]
  }  
  # remove outliers from wtko
  wtko_scan_trim <- wtko_scan_unfiltered[rownames,]
  outliers <- boxplot(wtko_scan_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    wtko_scan_trim<- wtko_scan_trim[-which(wtko_scan_trim$logFC %in% outliers),]
  }
  correlations[x] <- cor(MF_dmrt1L_trim[rownames,'logFC'],
                         wtko_scan_trim[rownames,'logFC'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(wtko_scan_trim[,'logFC'], # ko:wt first
             MF_dmrt1L_trim[,'logFC'], # reference M:F second
             by = 'row.names', 
             incomparables = NA)
  b <- a[complete.cases(a), ];b
  magnitudes[x] <-ang.vec.alph(b$x,b$y)[3] # this is the ratio of the magnitudes of each vector
  # the reference (wildtype M:F is the denominator)
  # if this is >1 then the ko:wt has a bigger effect
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_dmrt1L_trim <- sex_related_MF_dmrt1L[SL_rownames,]
outliers <- boxplot(sex_related_MF_dmrt1L_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_dmrt1L_trim <- sex_related_MF_dmrt1L_trim[-which(sex_related_MF_dmrt1L_trim$logFC %in% outliers),]
}  
# remove outliers from wtko
sex_related_wtko_scan_trim <- sex_related_wtko_scan[SL_rownames,]
outliers <- boxplot(sex_related_wtko_scan_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_wtko_scan_trim <- sex_related_wtko_scan_trim[-which(sex_related_wtko_scan_trim$logFC %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_dmrt1L_trim[SL_rownames,'logFC'],
                          sex_related_wtko_scan_trim[SL_rownames,'logFC'], 
                          method = "pearson", use="pairwise")
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1] 0.2147852

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_wtko_scan_trim[SL_rownames,'logFC'],
           sex_related_MF_dmrt1L_trim[SL_rownames,'logFC'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.1498501


# MF_dmrt1S vs scan
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(counts)[1], 74, replace = F);indexes
  rownames <- counts$geneID[indexes]
  # remove outliers from MF
  MF_dmrt1S_trim <- MF_dmrt1S_unfiltered[rownames,]
  outliers <- boxplot(MF_dmrt1S_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_dmrt1S_trim<- MF_dmrt1S_trim[-which(MF_dmrt1S_trim$logFC %in% outliers),]
  }  
  # remove outliers from wtko
  wtko_scan_trim <- wtko_scan_unfiltered[rownames,]
  outliers <- boxplot(wtko_scan_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    wtko_scan_trim<- wtko_scan_trim[-which(wtko_scan_trim$logFC %in% outliers),]
  }
  correlations[x] <- cor(MF_dmrt1S_trim[rownames,'logFC'],
                         wtko_scan_trim[rownames,'logFC'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(wtko_scan_trim[,'logFC'], # ko:wt first
             MF_dmrt1S_trim[,'logFC'], # reference M:F second
             by = 'row.names', 
             incomparables = NA)
  b <- a[complete.cases(a), ];b
  magnitudes[x] <-ang.vec.alph(b$x,b$y)[3] # this is the ratio of the magnitudes of each vector
  # the reference (wildtype M:F is the denominator)
  # if this is >1 then the ko:wt has a bigger effect
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_dmrt1S_trim <- sex_related_MF_dmrt1S[SL_rownames,]
outliers <- boxplot(sex_related_MF_dmrt1S_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_dmrt1S_trim <- sex_related_MF_dmrt1S_trim[-which(sex_related_MF_dmrt1S_trim$logFC %in% outliers),]
}  
# remove outliers from wtko
sex_related_wtko_scan_trim <- sex_related_wtko_scan[SL_rownames,]
outliers <- boxplot(sex_related_wtko_scan_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_wtko_scan_trim <- sex_related_wtko_scan_trim[-which(sex_related_wtko_scan_trim$logFC %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_dmrt1S_trim[SL_rownames,'logFC'],
                          sex_related_wtko_scan_trim[SL_rownames,'logFC'], 
                          method = "pearson", use="pairwise")
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1] 0.2277722

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_wtko_scan_trim[SL_rownames,'logFC'],
           sex_related_MF_dmrt1S_trim[SL_rownames,'logFC'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.6003996

# ccdc permutations ----

# MF_ccdc vs ccdc
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(counts)[1], 74, replace = F);indexes
  rownames <- counts$geneID[indexes]
  # remove outliers from MF
  MF_ccdc_trim <- MF_ccdc_unfiltered[rownames,]
  outliers <- boxplot(MF_ccdc_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_ccdc_trim<- MF_ccdc_trim[-which(MF_ccdc_trim$logFC %in% outliers),]
  }  
  # remove outliers from wtko
  wtko_ccdc_trim <- wtko_ccdc_unfiltered[rownames,]
  outliers <- boxplot(wtko_ccdc_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    wtko_ccdc_trim<- wtko_ccdc_trim[-which(wtko_ccdc_trim$logFC %in% outliers),]
  }
  correlations[x] <- cor(MF_ccdc_trim[rownames,'logFC'],
                         wtko_ccdc_trim[rownames,'logFC'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(wtko_ccdc_trim[,'logFC'], # ko:wt first
             MF_ccdc_trim[,'logFC'], # reference M:F second
             by = 'row.names', 
             incomparables = NA)
  b <- a[complete.cases(a), ];b
  magnitudes[x] <-ang.vec.alph(b$x,b$y)[3] # this is the ratio of the magnitudes of each vector
  # the reference (wildtype M:F is the denominator)
  # if this is >1 then the ko:wt has a bigger effect
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_ccdc_trim <- sex_related_MF_ccdc[SL_rownames,]
outliers <- boxplot(sex_related_MF_ccdc_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_ccdc_trim <- sex_related_MF_ccdc_trim[-which(sex_related_MF_ccdc_trim$logFC %in% outliers),]
}  
# remove outliers from wtko
sex_related_wtko_ccdc_trim <- sex_related_wtko_ccdc[SL_rownames,]
outliers <- boxplot(sex_related_wtko_ccdc_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_wtko_ccdc_trim <- sex_related_wtko_ccdc_trim[-which(sex_related_wtko_ccdc_trim$logFC %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_ccdc_trim[SL_rownames,'logFC'],
                          sex_related_wtko_ccdc_trim[SL_rownames,'logFC'], 
                          method = "pearson", use="pairwise")
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1] 0.972028

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_wtko_ccdc_trim[SL_rownames,'logFC'],
           sex_related_MF_ccdc_trim[SL_rownames,'logFC'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.1278721



# MF_dmrt1L vs ccdc
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(counts)[1], 74, replace = F);indexes
  rownames <- counts$geneID[indexes]
  # remove outliers from MF
  MF_dmrt1L_trim <- MF_dmrt1L_unfiltered[rownames,]
  outliers <- boxplot(MF_dmrt1L_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_dmrt1L_trim<- MF_dmrt1L_trim[-which(MF_dmrt1L_trim$logFC %in% outliers),]
  }  
  # remove outliers from wtko
  wtko_ccdc_trim <- wtko_ccdc_unfiltered[rownames,]
  outliers <- boxplot(wtko_ccdc_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    wtko_ccdc_trim<- wtko_ccdc_trim[-which(wtko_ccdc_trim$logFC %in% outliers),]
  }
  correlations[x] <- cor(MF_dmrt1L_trim[rownames,'logFC'],
                         wtko_ccdc_trim[rownames,'logFC'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(wtko_ccdc_trim[,'logFC'], # ko:wt first
             MF_dmrt1L_trim[,'logFC'], # reference M:F second
             by = 'row.names', 
             incomparables = NA)
  b <- a[complete.cases(a), ];b
  magnitudes[x] <-ang.vec.alph(b$x,b$y)[3] # this is the ratio of the magnitudes of each vector
  # the reference (wildtype M:F is the denominator)
  # if this is >1 then the ko:wt has a bigger effect
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_dmrt1L_trim <- sex_related_MF_dmrt1L[SL_rownames,]
outliers <- boxplot(sex_related_MF_dmrt1L_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_dmrt1L_trim <- sex_related_MF_dmrt1L_trim[-which(sex_related_MF_dmrt1L_trim$logFC %in% outliers),]
}  
# remove outliers from wtko
sex_related_wtko_ccdc_trim <- sex_related_wtko_ccdc[SL_rownames,]
outliers <- boxplot(sex_related_wtko_ccdc_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_wtko_ccdc_trim <- sex_related_wtko_ccdc_trim[-which(sex_related_wtko_ccdc_trim$logFC %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_dmrt1L_trim[SL_rownames,'logFC'],
                          sex_related_wtko_ccdc_trim[SL_rownames,'logFC'], 
                          method = "pearson", use="pairwise")
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1] 0.7852148

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_wtko_ccdc_trim[SL_rownames,'logFC'],
           sex_related_MF_dmrt1L_trim[SL_rownames,'logFC'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.06593407



# MF_dmrt1S vs ccdc
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(counts)[1], 74, replace = F);indexes
  rownames <- counts$geneID[indexes]
  # remove outliers from MF
  MF_dmrt1S_trim <- MF_dmrt1S_unfiltered[rownames,]
  outliers <- boxplot(MF_dmrt1S_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    MF_dmrt1S_trim<- MF_dmrt1S_trim[-which(MF_dmrt1S_trim$logFC %in% outliers),]
  }  
  # remove outliers from wtko
  wtko_ccdc_trim <- wtko_ccdc_unfiltered[rownames,]
  outliers <- boxplot(wtko_ccdc_trim$logFC, plot=FALSE)$out;outliers
  # check if there are any outliers
  if(any(outliers)) {
    wtko_ccdc_trim<- wtko_ccdc_trim[-which(wtko_ccdc_trim$logFC %in% outliers),]
  }
  correlations[x] <- cor(MF_dmrt1S_trim[rownames,'logFC'],
                         wtko_ccdc_trim[rownames,'logFC'], 
                         method = "pearson", use="pairwise")
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(wtko_ccdc_trim[,'logFC'], # ko:wt first
             MF_dmrt1S_trim[,'logFC'], # reference M:F second
             by = 'row.names', 
             incomparables = NA)
  b <- a[complete.cases(a), ];b
  magnitudes[x] <-ang.vec.alph(b$x,b$y)[3] # this is the ratio of the magnitudes of each vector
  # the reference (wildtype M:F is the denominator)
  # if this is >1 then the ko:wt has a bigger effect
}
# now figure out where the observed ranks within the correlations vector
# remove outliers from MF
sex_related_MF_dmrt1S_trim <- sex_related_MF_dmrt1S[SL_rownames,]
outliers <- boxplot(sex_related_MF_dmrt1S_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_MF_dmrt1S_trim <- sex_related_MF_dmrt1S_trim[-which(sex_related_MF_dmrt1S_trim$logFC %in% outliers),]
}  
# remove outliers from wtko
sex_related_wtko_ccdc_trim <- sex_related_wtko_ccdc[SL_rownames,]
outliers <- boxplot(sex_related_wtko_ccdc_trim$logFC, plot=FALSE)$out;outliers
# check if there are any outliers
if(any(outliers)) {
  sex_related_wtko_ccdc_trim <- sex_related_wtko_ccdc_trim[-which(sex_related_wtko_ccdc_trim$logFC %in% outliers),]
}
correlations[1001] <- cor(sex_related_MF_dmrt1S_trim[SL_rownames,'logFC'],
                          sex_related_wtko_ccdc_trim[SL_rownames,'logFC'], 
                          method = "pearson", use="pairwise")
print("pvalue: "); 1-rank(correlations)[1001]/1001
# [1] "pvalue: "
# [1] 0.8911089


# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_wtko_ccdc_trim[SL_rownames,'logFC'],
           sex_related_MF_dmrt1S_trim[SL_rownames,'logFC'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.4145854
```

# Plotting

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

setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_Supergene/2022_KO_tad_RNAseq/2022_EdgeR_and_DeSeq2/STAR_done")
dir <- "/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_Supergene/2022_KO_tad_RNAseq/2022_EdgeR_and_DeSeq2/STAR_done"
list.files(dir)


# import into a list
temp = list.files(pattern="STAR_edgeR_unfiltered.csv");temp
myfiles = lapply(temp, read.delim, sep = ",")

# rename the columns so they are sensible
colnames(myfiles[[1]]) <- c("gene","MF_ccdc_logFC","MF_ccdc_logCPM","MF_ccdc_PValue")
colnames(myfiles[[2]]) <- c("gene","MF_dmrt1L_logFC","MF_dmrt1L_logCPM","MF_dmrt1L_PValue")
colnames(myfiles[[3]]) <- c("gene","MF_dmrt1S_logFC","MF_dmrt1S_logCPM","MF_dmrt1S_PValue")
colnames(myfiles[[4]]) <- c("gene","wtko_ccdc_logFC","wtko_ccdc_logCPM","wtko_ccdc_PValue")
colnames(myfiles[[5]]) <- c("gene","wtko_dmw_logFC","wtko_dmw_logCPM","wtko_dmw_PValue")
colnames(myfiles[[6]]) <- c("gene","wtko_scan_logFC","wtko_scan_logCPM","wtko_scan_PValue")

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
MF_ccdc_logFC_trim<-alldata[,c(1,2)]
MF_ccdc_logFC_trim<- MF_ccdc_logFC_trim[-which(MF_ccdc_logFC_trim$MF_ccdc_logFC %in% outliers),]
#MF_dmrt1L
outliers <- boxplot(alldata$MF_dmrt1L_logFC, plot=FALSE)$out
MF_dmrt1L_logFC_trim<-alldata[,c(1,5)]
MF_dmrt1L_logFC_trim<- MF_dmrt1L_logFC_trim[-which(MF_dmrt1L_logFC_trim$MF_dmrt1L_logFC %in% outliers),]
#MF_dmrt1S
outliers <- boxplot(alldata$MF_dmrt1S_logFC, plot=FALSE)$out
MF_dmrt1S_logFC_trim<-alldata[,c(1,8)]
MF_dmrt1S_logFC_trim<- MF_dmrt1S_logFC_trim[-which(MF_dmrt1S_logFC_trim$MF_dmrt1S_logFC %in% outliers),]
#wtko_dmw
outliers <- boxplot(alldata$wtko_dmw_logFC, plot=FALSE)$out
wtko_dmw_logFC_trim<- alldata[,c(1,14)]
wtko_dmw_logFC_trim<- wtko_dmw_logFC_trim[-which(wtko_dmw_logFC_trim$wtko_dmw_logFC %in% outliers),]
#wtko_scan
outliers <- boxplot(alldata$wtko_scan_logFC, plot=FALSE)$out
wtko_scan_logFC_trim<- alldata[,c(1,17)]
wtko_scan_logFC_trim<- wtko_scan_logFC_trim[-which(wtko_scan_logFC_trim$wtko_scan_logFC %in% outliers),]
#wtko_ccdc
outliers <- boxplot(alldata$wtko_ccdc_logFC, plot=FALSE)$out
wtko_ccdc_logFC_trim<- alldata[,c(1,11)]
wtko_ccdc_logFC_trim<- wtko_ccdc_logFC_trim[-which(wtko_ccdc_logFC_trim$wtko_ccdc_logFC %in% outliers),]


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
        upper = list(continuous = my_fn2),
        lower = list(continuous = my_fn)) +
        #lower = list(continuous = my_custom_smooth)) +
  #theme_bw() +
  theme(strip.background = element_rect(
      color="white", fill="white", size=1.5, linetype="solid")) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
box1_2 <- ggally_text("\nr = 0.159\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box1_3 <- ggally_text("\nr = -0.098\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box1_4 <- ggally_text("\nr = 0.278*\np = 0.116\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box1_5 <- ggally_text("\nr = 0.123\np = 0.283\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box1_6 <- ggally_text("\nr = 0.212\np = 0.972\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box2_3 <- ggally_text("\nr = 0.002\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box2_4 <- ggally_text("\nr = 0.111\np = 0.162\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box2_5 <- ggally_text("\nr = 0.138\np = 0.215\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box2_6 <- ggally_text("\nr = -0.222\np = 0.785\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box3_4 <- ggally_text("\nr = 0.331*\np = 0.026\n",geom_text = ggplot2::aes(size = 6), color = I("red"))
box3_5 <- ggally_text("\nr = 0.085\np = 0.228\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box3_6 <- ggally_text("\nr = -0.400*\np = 0.891\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box4_5 <- ggally_text("\nr = 0.143\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box4_6 <- ggally_text("\nr = 0.025\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box5_6 <- ggally_text("\nr = -0.167\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
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
ggsave(file="STAR_edgeR_sexrelated_pairwise_unfiltered.pdf", g, width=10, height=4)


```

# Venn diagrams (needs to be updated)
```
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_dmrt1_knockouts/2023_KO_tad_RNAseq/2022_EdgeR_and_DeSeq2/STAR')
library(ggvenn)

# load contig names
MF_ccdc_DeSeq2 = read.table("MF_STAR_ccdconly_DE_DeSeq2.csv", sep = ",", header=T);
colnames(MF_ccdc_DeSeq2)[1] <- "gene"
MF_dmrt1L_DeSeq2 = read.table("MF_STAR_dmrt1Lonly_DE_DeSeq2.csv", sep = ",", header=T);
colnames(MF_dmrt1L_DeSeq2)[1] <- "gene"
MF_dmrt1S_DeSeq2 = read.table("MF_STAR_dmrt1Sonly_DE_DeSeq2.csv", sep = ",", header=T);
colnames(MF_dmrt1S_DeSeq2)[1] <- "gene"
wtkofemsonly_dmrt1L_DeSeq2 = read.table("wtkofemsonly_dmrt1L_STAR_DE_DeSeq2.csv", sep = ",", header=T);
colnames(MF_dmrt1S_DeSeq2)[1] <- "gene"
wtkomalesonly_dmrt1L_DeSeq2 = read.table("wtkomalesonly_dmrt1L_STAR_DE_DeSeq2.csv", sep = ",", header=T);
colnames(wtkomalesonly_dmrt1L_DeSeq2)[1] <- "gene"
wtkofemsonly_dmrt1S_DeSeq2 = read.table("wtko_STAR_dmrt1Sfemsonly_DE_DeSeq2.csv", sep = ",", header=T);
colnames(wtkofemsonly_dmrt1S_DeSeq2)[1] <- "gene"
wtkomalesonly_dmrt1S_DeSeq2 = read.table("wtko_STAR_dmrt1Smalesonly_DE_DeSeq2.csv", sep = ",", header=T);
colnames(wtkomalesonly_dmrt1S_DeSeq2)[1] <- "gene"


# Make list for comparison
x <- list(MF1=MF_ccdc_DeSeq2$gene, 
          MF2=MF_dmrt1L_DeSeq2$gene, 
          MF3=MF_dmrt1S_DeSeq2$gene,
          dmrt1L_F = wtkofemsonly_dmrt1L_DeSeq2$gene)

# plot
#pdf("dmw_wtko_and3_MF_DeSeq2.pdf",w=8, h=4, version="1.4", bg="transparent")
plot1 <- ggvenn(x, 
       show_elements = F,
       label_sep = "\n", 
       text_size = 3,
       fill_color = c("black","grey70", "grey80", "grey90"),
       set_name_size = 3,
       show_percentage = F) + 
        coord_cartesian(clip="off") +
        theme(plot.margin = margin(20,20,20,20))
#dev.off()

# Make list for comparison
x <- list(MF1=MF_ccdc_DeSeq2$gene, 
          MF2=MF_dmrt1L_DeSeq2$gene, 
          MF3=MF_dmrt1S_DeSeq2$gene,
          dmrt1L_M = wtkomalesonly_dmrt1L_DeSeq2$gene)

# plot
#pdf("dmw_wtko_and3_MF_DeSeq2.pdf",w=8, h=4, version="1.4", bg="transparent")
plot2 <- ggvenn(x, 
                show_elements = F,
                label_sep = "\n", 
                text_size = 3,
                fill_color = c("black","grey70", "grey80", "grey90"),
                set_name_size = 3,
                show_percentage = F) + 
                coord_cartesian(clip="off") +
                theme(plot.margin = margin(20,20,20,20))
#dev.off()


# Make list for comparison
x <- list(MF1=MF_ccdc_DeSeq2$gene, 
          MF2=MF_dmrt1L_DeSeq2$gene, 
          MF3=MF_dmrt1S_DeSeq2$gene,
          dmrt1S_F = wtkofemsonly_dmrt1S_DeSeq2$gene)

# plot
#pdf("dmw_wtko_and3_MF_DeSeq2.pdf",w=8, h=4, version="1.4", bg="transparent")
plot3 <- ggvenn(x, 
                show_elements = F,
                label_sep = "\n", 
                text_size = 3,
                fill_color = c("black","grey70", "grey80", "grey90"),
                set_name_size = 3,
                show_percentage = F) + 
                coord_cartesian(clip="off") +
                theme(plot.margin = margin(20,20,20,20))
#dev.off()


# Make list for comparison
x <- list(MF1=MF_ccdc_DeSeq2$gene, 
          MF2=MF_dmrt1L_DeSeq2$gene, 
          MF3=MF_dmrt1S_DeSeq2$gene,
          dmrt1S_M = wtkomalesonly_dmrt1S_DeSeq2$gene)

# plot
#pdf("dmw_wtko_and3_MF_DeSeq2.pdf",w=8, h=4, version="1.4", bg="transparent")
plot4 <- ggvenn(x, 
                show_elements = F,
                label_sep = "\n", 
                text_size = 3,
                fill_color = c("black","grey70", "grey80", "grey90"),
                set_name_size = 3,
                show_percentage = F) + 
                coord_cartesian(clip="off") +
                theme(plot.margin = margin(20,20,20,20))
#dev.off()


# Make list for comparison
x <- list(dmrt1L_F= wtkofemsonly_dmrt1L_DeSeq2$gene,
          dmrt1L_M = wtkomalesonly_dmrt1L_DeSeq2$gene,
          dmrt1S_F= wtkofemsonly_dmrt1S_DeSeq2$gene,
          dmrt1S_M= wtkomalesonly_dmrt1S_DeSeq2$gene)

# plot
#pdf("DMRT1L_and_DMRT1S_STAR_DeSeq2_Venn.pdf",w=8, h=4, version="1.4", bg="transparent")
plot5 <- ggvenn(x, 
       show_elements = F,
       label_sep = "\n", 
       text_size = 3,
       fill_color = c("black","grey70", "grey80", "grey90"),
       set_name_size = 3,
       show_percentage = F) + 
        coord_cartesian(clip="off")+
        theme(plot.margin = margin(20,20,20,20))
#dev.off()




library(gridExtra)
#merge all three plots within one grid (and visualize this)
g <- grid.arrange(arrangeGrob(plot1, plot2, plot3, plot4,ncol=2,nrow=2), plot5, 
             ncol=2,
             top = "Kallisto and DeSeq2")

ggsave(file="DMRT1L_and_DMRT1S_STAR_DeSeq2_Venn.pdf", w=10, h=6, g)


```
