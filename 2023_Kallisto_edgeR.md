# Kallisto_EdgeR plus permutations

```R
library(edgeR)
library(tximport)
library('rhdf5')
library('readxl')
library('ggplot2')
library(grid)
require('gridExtra')
library("org.Xl.eg.db")
library(PCAtools)
library("HTSFilter")
library(tidyverse)
library(writexl)
setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_Supergene/2022_KO_tad_RNAseq/2022_EdgeR_and_DeSeq2/2022_Kallisto_EdgeR_done")
dir <- "/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_Supergene/2022_KO_tad_RNAseq/2022_EdgeR_and_DeSeq2/2022_Kallisto_EdgeR_done"
list.files(dir)

# load MF count data (from Kalisto)
# read the count data from Kallisto that was combined from each sample
# into a single file into a dataframe called "counts"
MFcounts <- read.table("MF_.isoform.TMM.EXPR.matrix", header=T, row.names = 1)

# samples
MFsamples <- read.table(file.path(dir, "MF_samples.txt"), header = F);samples

# MF dmrt1L ----
colnames(MFcounts)
new_MFcounts <- MFcounts[,c(3:5,15:19) ]
new_MFsamples <-as.data.frame(MFsamples[c(3:5,15:19), ]);new_MFsamples
new_MFsexez <- factor(c(rep("f",3),rep("m",5)))
new_MFsexez <- relevel(new_MFsexez, ref="f")
# Create DGEList object - this is a data structure that is used for 
# the analysis of differential expression
d0 <- DGEList(new_MFcounts, group = new_MFsexez, remove.zeros = TRUE)
dim(d0$counts) #each row is a transcript - here is the number before filtering
# [1] 39217     8



# save the unfiltered logFC to a dataframe 
d0 <- calcNormFactors(d0, method="TMM")
design <- model.matrix(~ 0 + new_MFsexez, data=d0$samples) # last coefficient = difference between sexes)
d0 <- estimateDisp(d0, design, robust=TRUE)
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
MF_dmrt1L_unfiltered <- exacttest$table;MF_dmrt1L_unfiltered
# Write sex_related to a file
#exacttest$table[grepl("*\\|XBXL10_1g2070\\|", rownames(exacttest$table)), ]
sex_related_MF_dmrt1L <- exacttest$table[grepl("gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g34871\\|", rownames(exacttest$table)), ]
write.csv(sex_related_MF_dmrt1L, file="Sex_related_MF_dmrt1L_Kallisto_edgeR_unfiltered.csv", row.names = T)
# Write counts of sex related to a file
sex_related_dmrt1L_counts <- new_MFcounts[grepl("gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g34871\\|", rownames(new_MFcounts)), ]
write.csv(sex_related_dmrt1L_counts, file="Sex_related_MF_dmrt1L_Kallisto_edgeR_counts_unfiltered.csv", row.names = T)

# Now do analysis of differential expression; 
# first remove transcripts where the average count per sample is 2 or less:
d0$counts <- d0$counts[rowSums(d0$counts)> 2* ncol(d0$counts),] 
dim(d0$counts)
# 10782     8
# Now we have far fewer transcripts:
# TMM normalization is applied to this dataset to account for 
# compositional difference between the libraries.
d0 <- calcNormFactors(d0, method="TMM")
# check the normalization factors
d0$samples
# plot by sex
#plotMDS(d0,labels=c(rep("F",3),rep("M",5)),col=c(rep("green",3),rep("blue",5)))
# design matrix: this is used for the model of differential expression
design <- model.matrix(~ 0 + new_MFsexez, data=d0$samples) # last coefficient = difference between sexes)
# estimate dispersion
d0 <- estimateDisp(d0, design, robust=TRUE)
#d0$common.dispersion
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
#exacttest
summary(decideTests(object = exacttest, p.value = 0.1)) 
#         m-f
# Down       3
# NotSig 10779
# Up         0
topTags(exacttest, n=3)

# logFC   logCPM       PValue         FDR
# gnl|XBXL10_1g8664|XBmRNA16317|     -4.568098 2.070722 5.284365e-07 0.005697602
# gnl|XBXL10_1g21751|XBmRNA40341|    -2.488928 2.504201 1.404523e-05 0.059567088
# gnl|XBXL10_1g44212|XBmRNA83108|    -3.277658 4.107902 1.657404e-05 0.059567088
# gnl|XBXL10_1g37344|XBmRNA70491|    -2.895938 7.819552 6.190799e-05 0.148644488
MF_dmrt1L_DE <- as.data.frame(topTags(exacttest, n=3))
write.csv(MF_dmrt1L_DE, file="MF_Kallisto_dmrt1L_edgeR.csv", row.names = T)

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
# all are above 0.8 
library(writexl)
write_xlsx(rsquare, "./MF_rsquare_dmrt1L_Kallisto_edgeR.xls")



# dmrt1S ----
colnames(MFcounts)
new_MFcounts <- MFcounts[,c(6:8,20:22) ]
new_MFsamples <-as.data.frame(MFsamples[c(6:8,20:22), ]);new_MFsamples
new_MFsexez <- factor(c(rep("f",3),rep("m",3)))
new_MFsexez <- relevel(new_MFsexez, ref="f")
# Create DGEList object - this is a data structure that is used for 
# the analysis of differential expression
d0 <- DGEList(new_MFcounts, group = new_MFsexez, remove.zeros = TRUE)
dim(d0$counts) #each row is a transcript - here is the number before filtering
# [1] 38394     6

# save the unfiltered logFC to a dataframe 
d0 <- calcNormFactors(d0, method="TMM")
design <- model.matrix(~ 0 + new_MFsexez, data=d0$samples) # last coefficient = difference between sexes)
d0 <- estimateDisp(d0, design, robust=TRUE)
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
MF_dmrt1S_unfiltered <- exacttest$table;MF_dmrt1S_unfiltered
# Write sex_related to a file
sex_related_MF_dmrt1S <- exacttest$table[grepl("gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g34871\\|", rownames(exacttest$table)), ]
write.csv(sex_related_MF_dmrt1S, file="Sex_related_MF_dmrt1S_Kallisto_edgeR_unfiltered.csv", row.names = T)
# Write counts of sex related to a file
sex_related_dmrt1S_counts <- new_MFcounts[grepl("gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g34871\\|", rownames(new_MFcounts)), ]
write.csv(sex_related_dmrt1S_counts, file="Sex_related_MF_dmrt1S_Kallisto_edgeR_counts_unfiltered.csv", row.names = T)

# Now do analysis of differential expression; 
# here we remove transcripts where the average count per sample is 2 or less:
d0$counts <- d0$counts[rowSums(d0$counts)> 2* ncol(d0$counts),] 
# Now we have far fewer transcripts:
dim(d0$counts)
# 10899     6
# TMM normalization is applied to this dataset to account for 
# compositional difference between the libraries.
d0 <- calcNormFactors(d0, method="TMM")
# check the normalization factors
d0$samples
# plot by sex
#plotMDS(d0,labels=c(rep("F",3),rep("M",3)),col=c(rep("green",3),rep("blue",3)))
# design matrix: this is used for the model of differential expression
design <- model.matrix(~ new_MFsexez, data=d0$samples) # last coefficient = difference between sexes)
# estimate dispersion
d0 <- estimateDisp(d0, design, robust=TRUE)
#d0$common.dispersion
#plotBCV(d0)
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
#exacttest
summary(decideTests(object = exacttest, p.value = 0.1)) 
# m-f
# Down      75
# NotSig 10661
# Up       163

topTags(exacttest, n=238)
MF_dmrt1S_DE <- as.data.frame(topTags(exacttest, n=238))
write.csv(MF_dmrt1S_DE, file="MF_Kallisto_dmrt1S_edgeR.csv", row.names = T)

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
# just one below 0.8 
library(writexl)
write_xlsx(rsquare, "./MF_rsquare_dmrt1S_Kallisto_edgeR.xls")

# MF ccdc ---- 
colnames(MFcounts)
new_MFcounts <- MFcounts[,c(1:2,9:14) ]
new_MFsamples <-as.data.frame(MFsamples[c(1:2,9:14), ]);new_MFsamples
new_MFsexez <- factor(c(rep("f",2),rep("m",6)))
new_MFsexez <- relevel(new_MFsexez, ref="f")
# Create DGEList object - this is a data structure that is used for 
# the analysis of differential expression
d0 <- DGEList(new_MFcounts, group = new_MFsexez, remove.zeros = TRUE)
dim(d0$counts) #each row is a transcript - here is the number before filtering
# 39029     8

# save the unfiltered logFC to a dataframe 
d0 <- calcNormFactors(d0, method="TMM")
design <- model.matrix(~ 0 + new_MFsexez, data=d0$samples) # last coefficient = difference between sexes)
d0 <- estimateDisp(d0, design, robust=TRUE)
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
MF_ccdc_unfiltered <- exacttest$table;MF_ccdc_unfiltered
# Write sex_related to a file
sex_related_MF_ccdc <- exacttest$table[grepl("gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g34871\\|", rownames(exacttest$table)), ]
write.csv(sex_related_MF_ccdc, file="Sex_related_MF_ccdc_Kallisto_edgeR_unfiltered.csv", row.names = T)
# Write counts of sex related to a file
sex_related_ccdc_counts <- new_MFcounts[grepl("gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g34871\\|", rownames(new_MFcounts)), ]
write.csv(sex_related_ccdc_counts, file="Sex_related_MF_ccdc_Kallisto_edgeR_counts_unfiltered.csv", row.names = T)

# Now do analysis of differential expression; 
# here we remove transcripts where the average count per sample is 2 or less:
d0$counts <- d0$counts[rowSums(d0$counts)> 2* ncol(d0$counts),] 
# Now we have far fewer transcripts:
dim(d0$counts)
# 11023     8
# TMM normalization is applied to this dataset to account for compositional difference between
# the libraries.
d0 <- calcNormFactors(d0, method="TMM")
# check the normalization factors
d0$samples
# plot by sex
#plotMDS(d0,labels=c(rep("F",2),rep("M",6)),col=c(rep("green",2),rep("blue",6)))
# design matrix: this is used for the model of differential expression
design <- model.matrix(~ new_MFsexez, data=d0$samples) # last coefficient = difference between sexes)
#design
# estimate dispersion
d0 <- estimateDisp(d0, design, robust=TRUE)
#d0$common.dispersion
#plotBCV(d0)
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
#exacttest
#summary(decideTests(object = exacttest, lfc = 1))
summary(decideTests(object = exacttest, p.value = 0.1)) 
#          m-f
# Down       1
# NotSig 11022
# Up         0

topTags(exacttest, n=1)
# logFC   logCPM       PValue          FDR
# gnl|XBXL10_1g23832|XBmRNA44524| -5.005899 2.086205 4.242890e-10 4.676938e-06
# gnl|XBXL10_1g33805|XBmRNA63843|  5.723055 5.683792 2.219891e-05 1.223493e-01

MF_ccdc_DE <- as.data.frame(topTags(exacttest, n=1))
write.csv(MF_ccdc_DE, file="MF_Kallisto_ccdc_edgeR.csv", row.names = T)

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
# five below 0.8 
write_xlsx(rsquare, "./MF_rsquare_ccdc_Kallisto_edgeR.xls")


# dmw ----
dmw_counts <- read.table("2022_dmw_only.isoform.TMM.EXPR.matrix", header=T, row.names = 1)
dim(dmw_counts)
# [1]  44441    12 # 2022_Kalisto 
# get rid of any rows that have incomplete data
dmw_counts <- dmw_counts[complete.cases(dmw_counts), ]
dim(dmw_counts)
# [1]  44441    12 # 2022_Kalisto

# this is useful info about batch effects:
# https://support.bioconductor.org/p/96627/

# make design matrix
# Here we want to test for differential expression between KO and
# WT, while adjusting for differences between batches. In statistical
# terms, this "may be" an additive linear model with batch as the blocking factor:

# samples
dmw_samples <- read.table(file.path(dir, "dmw_samples.txt"), header = F)
dmw_samples

# rename the col names - it is crucial that the order of the samples in the "samples" file
# match the order of the columns in the "counts" dataframe
colnames(dmw_counts) <- dmw_samples$V1

# batch (for dmw we are not able to control for batch effects because the second batch had no ko individuals)
batch <- factor(c("first","first","first","first","first","first",
                  "second","second","second","first","first","second"))
batch <- relevel(batch, ref="first")

# genotypes
genotypez <- factor(c("dmwKO","dmwKO","dmwKO","dmwKO","dmwKO","dmwKO",
                      "WT","WT","WT","WT","WT","WT"))
genotypez <- relevel(genotypez, ref="WT") # this makes the expression levels relative to WT, and not KO

# sexez
sexez <- factor(c("F","F","F","F","F","F","F","F","F","F","F","F"))
sexez <- relevel(sexez, ref="F")

# Create DGEList object - this is a data structure that is used for 
# the analysis of differential expression
d0 <- DGEList(dmw_counts, group = factor(genotypez), remove.zeros = TRUE)

# save the unfiltered logFC to a dataframe 
d0 <- calcNormFactors(d0, method="TMM")
design<-model.matrix(~0+genotypez, data=d0$samples) # last coefficient = difference between genotypez
d0 <- estimateDisp(d0, design, robust=TRUE)
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
wtko_dmw_unfiltered <- exacttest$table;wtko_dmw_unfiltered
# Write sex_related to a file
sex_related_wtko_dmw <- exacttest$table[grepl("gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g34871\\|", rownames(exacttest$table)), ]
write.csv(sex_related_wtko_dmw, file="Sex_related_wtko_dmw_Kallisto_edgeR_unfiltered.csv", row.names = T)
# Write counts of sex related to a file
sex_related_wtko_dmw_counts <- dmw_counts[grepl("gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g34871\\|", rownames(dmw_counts)), ]
write.csv(sex_related_wtko_dmw_counts, file="Sex_related_wtko_dmw_Kallisto_edgeR_counts_unfiltered.csv", row.names = T)

# Now do analysis of differential expression; 
# here we remove transcripts where the average count per sample is 2 or less:
d0$counts <- d0$counts[rowSums(d0$counts)> 2* ncol(d0$counts),] 
# Now we have far fewer transcripts:
dim(d0$counts)
# [1] 10632    12 # 2022 kallisto
# many rows with low expression were eliminated
# TMM normalization is applied to this dataset to account for 
# compositional difference between the libraries.
d0 <- calcNormFactors(d0, method="TMM")
# check the normalization factors
d0$samples

# plot by genotype
#plotMDS(d0,labels=NULL,col=c(rep("green",6),rep("blue",6)))
#dmw_15_wt and dmw_28_ko seem weird

# plot by batch
#plotMDS(d0,labels=NULL,col=c(rep("green",6),rep("blue",3),rep("green",2),"blue"))
#dmw_15_wt and dmw_28_ko seem weird

# pca
# project.pca <- prcomp(t(d0$counts))
# plot(project.pca$x, type="n")
# points(project.pca$x, col=genotypez, pch=16, cex=1)
# points(project.pca$x, col=sexez, pch=16, cex=1)
# points(project.pca$x, col=batch, pch=16, cex=1) +
# theme(legend.position="top")
# var_explained <- project.pca$sdev^2/sum(project.pca $sdev^2)
# var_explained

# project.pca$x %>% 
#  as.data.frame %>%
  # rownames_to_column("continent_country") %>%
  # separate(continent_country,c("continent")) %>%
#  ggplot(aes(x=PC1,y=PC2)) + geom_point(aes(color=genotypez),size=4) +
#  theme_bw(base_size=12) + 
#  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
#       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
#  theme(legend.position="top")

# pca mostly overlaps...

# design matrix: this is used for the model of differential expression
design<-model.matrix(~0+genotypez, data=d0$samples) # last coefficient = difference between genotypez
# estimate dispersion
d0 <- estimateDisp(d0, design, robust=TRUE)
#d0$common.dispersion
#d0$tagwise.dispersion
#plotBCV(d0)
# exact test
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
summary(decideTests(object = exacttest, p.value = 0.1))  # no DEs
# dmwKO-WT
# Down          0
# NotSig    10632
# Up            0
topTags(exacttest, n=1)
wtko_dmw_DE <- as.data.frame(topTags(exacttest, n=1))
wtko_dmw_DE <-wtko_dmw_DE[-c(1),];wtko_dmw_DE
# Just output the headers
write.csv(wtko_dmw_DE, file="wtko_Kallisto_dmw_edgeR.csv", row.names = T)

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
# two below 0.8 
library(writexl)
write_xlsx(rsquare, "./MF_rsquare_dmw_Kallisto_edgeR.xls")



# scan ----
scan_counts <- read.table("scanw.isoform.TMM.EXPR.matrix", header=T, row.names = 1)
dim(scan_counts)
# [1]  44441    9 # 2022_Kalisto 
# get rid of any rows that have incomplete data
scan_counts <- scan_counts[complete.cases(scan_counts), ]
dim(scan_counts)
# [1]  44441    9 # 2022_Kalisto
samples <- read.table(file.path(dir, "scan_samples.txt"), header = F)
samples
# rename the col names - it is crucial that the order of the samples in the "samples" file
# match the order of the columns in the "counts" dataframe
colnames(scan_counts) <- samples$V1
# genotypes
genotypez <- factor(c("wt","wt","ko","ko","ko",
                      "ko","ko","wt","wt"))
genotypez <- relevel(genotypez, ref="wt") # this makes the expression levels relative to WT, and not KO
# sexez
sexez <- factor(c("f","f","f","f","f","f","f","f","f"))
sexez <- relevel(sexez, ref="f")
# Create DGEList object - this is a data structure that is used for 
# the analysis of differential expression
d0 <- DGEList(scan_counts, group = factor(genotypez), remove.zeros = TRUE)
dim(d0$counts) #each row is a transcript - here is the number before filtering
# [1] 39323     9 # 2022 kallisto ccdc

# save the unfiltered logFC to a dataframe 
d0 <- calcNormFactors(d0, method="TMM")
design<-model.matrix(~0+genotypez, data=d0$samples) # last coefficient = difference between genotypez
d0 <- estimateDisp(d0, design, robust=TRUE)
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
wtko_scan_unfiltered <- exacttest$table;wtko_scan_unfiltered
# Write sex_related to a file
sex_related_wtko_scan <- exacttest$table[grepl("gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g34871\\|", rownames(exacttest$table)), ]
write.csv(sex_related_wtko_scan, file="Sex_related_wtko_scan_Kallisto_edgeR_unfiltered.csv", row.names = T)
# Write counts of sex related to a file
sex_related_wtko_scan_counts <- scan_counts[grepl("gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g34871\\|", rownames(scan_counts)), ]
write.csv(sex_related_wtko_scan_counts, file="Sex_related_wtko_scan_Kallisto_edgeR_counts_unfiltered.csv", row.names = T)

# Now do analysis of differential expression; 
# here we remove transcripts where the average count per sample is 2 or less:
d0$counts <- d0$counts[rowSums(d0$counts)> 2* ncol(d0$counts),] 
# Now we have far fewer transcripts:
dim(d0$counts)
# [1] 13034     9 # 2022 kallisto ccdc
# many rows with low expression were eliminated
# TMM normalization is applied to this dataset to account for 
# compositional difference between the libraries.
d0 <- calcNormFactors(d0, method="TMM")
# check the normalization factors
d0$samples

# plot by genotype
# plotMDS(d0,labels=NULL,col=c(rep("green",2),rep("blue",5),rep("green",2)))
#dmw_15_wt and dmw_28_ko seem weird

# plot by batch
# plotMDS(d0,labels=NULL,col=c(rep("green",6),rep("blue",3),rep("green",2),"blue"))
#dmw_15_wt and dmw_28_ko seem weird

# pca
#project.pca <- prcomp(t(d0$counts))
#plot(project.pca$x, type="n")
#points(project.pca$x, col=genotypez, pch=16, cex=1)
#points(project.pca$x, col=sexez, pch=16, cex=1)
#points(project.pca$x, col=batch, pch=16, cex=1) +
#theme(legend.position="top")
#var_explained <- project.pca$sdev^2/sum(project.pca $sdev^2)
#var_explained

#project.pca$x %>% 
#  as.data.frame %>%
# rownames_to_column("continent_country") %>%
# separate(continent_country,c("continent")) %>%
#  ggplot(aes(x=PC1,y=PC2)) + geom_point(aes(color=genotypez),size=4) +
#  theme_bw(base_size=12) + 
#  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
#       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
#  theme(legend.position="top")

# design matrix: this is used for the model of differential expression
design<-model.matrix(~0+genotypez, data=d0$samples) # last coefficient = difference between genotypez
# estimate dispersion
d0 <- estimateDisp(d0, design, robust=TRUE)
#d0$common.dispersion
#d0$tagwise.dispersion
#plotBCV(d0)

# exact test
exacttest <- exactTest(d0, dispersion = "auto") 
summary(decideTests(object = exacttest, p.value = 0.1)) 
# ko-wt
# Down       3
# NotSig 13031
# Up         0

topTags(exacttest, n=3)
# logFC   logCPM       PValue        FDR
# gnl|XBXL10_1g43753|XBmRNA82360| -3.243891 3.995481 5.462670e-06 0.07120044
# gnl|XBXL10_1g26238|XBmRNA49054| -5.602802 2.862834 1.157393e-05 0.07542732
# gnl|XBXL10_1g569|XBmRNA596|     -2.799005 4.521137 1.982904e-05 0.08615057

wtko_scan_DE <- as.data.frame(topTags(exacttest, n=3))
write.csv(wtko_scan_DE, file="wtko_Kallisto_scan_edgeR.csv", row.names = T)

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
# none below 0.8 
write_xlsx(rsquare, "./MF_rsquare_scan_Kallisto_edgeR.xls")


# ccdc ----
ccdc_counts <- read.table("ccdc.isoform.TMM.EXPR.matrix", header=T, row.names = 1)

dim(ccdc_counts)
# [1]  44441    12 # 2022_Kalisto 
# get rid of any rows that have incomplete data
ccdc_counts <- ccdc_counts[complete.cases(ccdc_counts), ]
dim(ccdc_counts)
# [1]  44441    12 # 2022_Kalisto
# samples
ccdc_samples <- read.table(file.path("ccdc_samples.txt"), header = F)
ccdc_samples

# rename the col names - it is crucial that the order of the samples in the "samples" file
# match the order of the columns in the "counts" dataframe
colnames(ccdc_counts) <- ccdc_samples$V1

# genotypes
genotypez <- factor(c("wt","ko","wt","ko","ko","wt",
                      "ko","ko","wt","ko","wt","wt",
                      "wt","wt"))
genotypez <- relevel(genotypez, ref="wt") # this makes the expression levels relative to WT, and not KO

# sexez
sexez <- factor(c("m","f","m","f","f","f","f","f","f","f","m","m","m","m"))
sexez <- relevel(sexez, ref="f")

# subset the females only
ccdc_counts_onlyfemales <- ccdc_counts[,c(2,4:10)]
genotypez_onlyfemales <- genotypez[c(2,4:10)]
sexez_onlyfemales <- rep("f",8)

# Create DGEList object - this is a data structure that is used for 
# the analysis of differential expression
d0 <- DGEList(ccdc_counts_onlyfemales, group = factor(genotypez_onlyfemales), remove.zeros = TRUE)
dim(d0$counts) #each row is a transcript - here is the number before filtering
# [1] 38682     8 # 2022 kallisto ccdc

# save the unfiltered logFC to a dataframe 
d0 <- calcNormFactors(d0, method="TMM")
design<-model.matrix(~0+genotypez_onlyfemales, data=d0$samples) # last coefficient = difference between genotypez
d0 <- estimateDisp(d0, design, robust=TRUE)
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
wtko_ccdc_unfiltered <- exacttest$table;wtko_ccdc_unfiltered
# Write sex_related to a file
sex_related_wtko_ccdc <- exacttest$table[grepl("gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g34871\\|", rownames(exacttest$table)), ]
write.csv(sex_related_wtko_ccdc, file="Sex_related_wtko_ccdc_Kallisto_edgeR_unfiltered.csv", row.names = T)
# Write counts of sex related to a file
sex_related_wtko_ccdc_counts <- ccdc_counts_onlyfemales[grepl("gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g34871\\|", rownames(scan_counts)), ]
write.csv(sex_related_wtko_ccdc_counts, file="Sex_related_wtko_ccdc_Kallisto_edgeR_counts_unfiltered.csv", row.names = T)

# Now do analysis of differential expression; 
# here we remove transcripts where the average count per sample is 2 or less:
d0$counts <- d0$counts[rowSums(d0$counts)> 2* ncol(d0$counts),] 
# Now we have far fewer transcripts:
dim(d0$counts)
# [1] 11702     8 # 2022 kallisto ccdc
# TMM normalization is applied to this dataset to account for 
# compositional difference between the libraries.
d0 <- calcNormFactors(d0, method="TMM")
# check the normalization factors
d0$samples

# plot by genotype
#plotMDS(d0,labels=NULL,col=c(rep("green",3),"blue",rep("green",2),"blue","green"))
#dmw_15_wt and dmw_28_ko seem weird

# plot by batch
# plotMDS(d0,labels=NULL,col=c(rep("green",6),rep("blue",3),rep("green",2),"blue"))
#dmw_15_wt and dmw_28_ko seem weird

# pca
#project.pca <- prcomp(t(d0$counts))
#plot(project.pca$x, type="n")
#points(project.pca$x, col=genotypez, pch=16, cex=1)
#points(project.pca$x, col=sexez, pch=16, cex=1)
#points(project.pca$x, col=batch, pch=16, cex=1) +
#theme(legend.position="top")
#var_explained <- project.pca$sdev^2/sum(project.pca $sdev^2)
#var_explained

#project.pca$x %>% 
#  as.data.frame %>%
# rownames_to_column("continent_country") %>%
# separate(continent_country,c("continent")) %>%
#  ggplot(aes(x=PC1,y=PC2)) + geom_point(aes(color=genotypez),size=4) +
#  theme_bw(base_size=12) + 
#  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
#       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
#  theme(legend.position="top")

# design matrix: this is used for the model of differential expression
design<-model.matrix(~0+genotypez_onlyfemales, data=d0$samples) # last coefficient = difference between genotypez
# estimate dispersion
d0 <- estimateDisp(d0, design, robust=TRUE)
#plotBCV(d0)

# exact test
exacttest <- exactTest(d0, dispersion = "auto") 
exacttest
summary(decideTests(object = exacttest, p.value = 0.1)) 
#        ko-wt
# Down      18
# NotSig 11682
# Up         2
topTags(exacttest, n=20)
wtko_ccdc_DE <- as.data.frame(topTags(exacttest, n=20))
write.csv(wtko_ccdc_DE, file="wtko_Kallisto_ccdc_edgeR.csv", row.names = T)

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
# none below 0.8 
write_xlsx(rsquare, "./MF_rsquare_ccdc_Kallisto_edgeR.xls")

# permutations ----

# OK now do some permutations to assess significance of the correlations
# between the log2FC of each wtko and each MF

# These permutations will randomly select 74 logFC from each MF and wtko
# 1000 times and calculate the correlation.  Then this will be compared to
# the observed

# get rownames of sexrelated transcripts
SL_rownames <- row.names(MFcounts[grepl("gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g34871\\|", rownames(MFcounts)), ])

# dmw permutations ----

# MF_ccdc vs dmw
correlations <- c()
# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 74, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
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
# [1] 0.02097902


# MF_dmrt1L vs dmw
correlations <- c()
# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 74, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
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
# [1] 0.4885115

# MF_dmrt1S vs dmw
correlations <- c()
# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 74, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
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
# [1] 0.1238761


# scan permutations ----

# MF_ccdc vs scan
correlations <- c()
# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 74, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
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
# [1] 0.7562438


# MF_dmrt1L vs scan
correlations <- c()
# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 74, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
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
# [1] 0.4365634

# MF_dmrt1S vs scan
correlations <- c()
# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 74, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
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
# [1] 0.3486513


# ccdc permutations ----

# MF_ccdc vs ccdc
correlations <- c()
# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 74, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
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
# [1] 0.7792208


# MF_dmrt1L vs ccdc
correlations <- c()
# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 74, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
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

# MF_dmrt1S vs ccdc
correlations <- c()
# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 74, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
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
# [1] 0.3426573

```
