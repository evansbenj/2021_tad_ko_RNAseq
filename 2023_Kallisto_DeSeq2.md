# kallisto DeSeq2 with permutations
```R
## RNA-seq analysis with DESeq2
## Adapted from Stephen Turner, @genetics_blog
library(DESeq2)
# https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#input-data
# if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("gplots")
# BiocManager::install("DESeq2")
library(apeglm)
# RNA-seq data from PRJNA315516
# https://www.ncbi.nlm.nih.gov//bioproject/PRJNA315516.
# 3 control samples ("ctl"), 3 samples grown under drought condition ("dro")
library(writexl)
# Import & pre-process ----------------------------------------------------

# Import data from featureCounts
## Previously ran at command line something like this:
setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_Supergene/2022_KO_tad_RNAseq/2022_EdgeR_and_DeSeq2/2023_Kallisto_DeSeq2_done")
dir <- "/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_Supergene/2022_KO_tad_RNAseq/2022_EdgeR_and_DeSeq2/2023_Kallisto_DeSeq2_done"
list.files(dir)


# MF ccdc ----
counts <- read.table("MF_.isoform.TMM.EXPR.matrix", header=T, row.names = 1)
colnames(counts)
new_counts <- counts[,c(1:2,9:14)]
colnames(new_counts)
coldata <- read.table("MF_sample_genotype.txt", header=T, row.names = 1)
new_coldata <- coldata[c(1:2,9:14),]; new_coldata
new_coldata$batch <- as.factor(new_coldata$batch)
new_coldata$sex <- as.factor(new_coldata$sex)

dds <- DESeqDataSetFromMatrix(countData = round(new_counts),
                              colData = new_coldata,
                              design= ~ sex)
# save the unfiltered log2FoldChange to a dataframe 
dds <- DESeq(dds)
dds$sex <- relevel(dds$sex, ref="F") # this makes the expression levels relative to F
res <- results(dds)
MF_ccdc_unfiltered <- res;MF_ccdc_unfiltered
# Only sex related
sex_related_MF_ccdc <- res[grepl("gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g34871\\|", rownames(res)), ]
write.csv(sex_related_MF_ccdc, file="Sex_related_MF_ccdc_Kallisto_DeSeq2_unfiltered.csv", row.names = T)

# Now do analysis of differential expression; 
# first remove transcripts where the average count per sample is 2 or less:
keep <- rowSums(counts(dds)) >= 2* length(colnames(dds))
dds <- dds[keep,]
dim(dds)
#[1] 11182     8
# relevel
dds$sex <- relevel(dds$sex, ref="F") # this makes the expression levels relative to F
# now do the analysis
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$pvalue),]
summary(res)
p<-resOrdered[1,];p
write.csv(p, file="MF_Kallisto_ccdc_DE_DeSeq2.csv", row.names = T)
# get Rsquare value for all pairwise comparisons
# normalized counts for DeSeq2:
# https://bioinformatics.stackexchange.com/questions/193/how-can-i-extract-normalized-read-count-values-from-deseq2-results
norm <- as.data.frame(counts(dds, normalized=T))
rsquare <- data.frame(matrix(ncol = ncol(norm), 
                             nrow = ncol(norm)))
for(i in 1:(ncol(norm)-1)) {       # for-loop over columns
  for(j in (i+1):ncol(norm)) { 
    print(paste(i," ",j))
    x <- cor.test(norm[ , i], 
                  norm[ , j], 
                  method = 'spearman')
    rsquare[i,j] <- x$estimate
  }
}
colnames(rsquare) <- colnames(norm)
rownames(rsquare) <- colnames(norm)
View(rsquare)
write_xlsx(rsquare, "./MF_Kallisto_ccdc_rsquare.xls")


# MF dmrt1L ----
colnames(counts)
new_counts <- counts[,c(3:5,15:19)]
colnames(new_counts)
coldata <- read.table("MF_sample_genotype.txt", header=T, row.names = 1)
new_coldata <- coldata[c(3:5,15:19),]; new_coldata
new_coldata$batch <- as.factor(new_coldata$batch)
new_coldata$sex <- as.factor(new_coldata$sex)
dds <- DESeqDataSetFromMatrix(countData = round(new_counts),
                              colData = new_coldata,
                              design= ~ sex)
# save the unfiltered log2FoldChange to a dataframe 
dds <- DESeq(dds)
dds$sex <- relevel(dds$sex, ref="F") # this makes the expression levels relative to F
res <- results(dds)
MF_dmrt1L_unfiltered <- res;MF_dmrt1L_unfiltered
# Only sex related
sex_related_MF_dmrt1L <- res[grepl("gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g34871\\|", rownames(res)), ]
write.csv(sex_related_MF_dmrt1L, file="Sex_related_MF_dmrt1L_Kallisto_DeSeq2_unfiltered.csv", row.names = T)

# Now do analysis of differential expression; 
# first remove transcripts where the average count per sample is 2 or less:
keep <- rowSums(counts(dds)) >= 2* length(colnames(dds))
dds <- dds[keep,]
dim(dds)
#[1] 10955     8
# relevel
dds$sex <- relevel(dds$sex, ref="F") # this makes the expression levels relative to F
# now do the analysis
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$pvalue),]
summary(res)
p<-resOrdered[1:4,];p
write.csv(p, file="MF_Kallisto_dmrt1L_DE_DeSeq2.csv", row.names = T)
# get Rsquare value for all pairwise comparisons
# normalized counts for DeSeq2:
# https://bioinformatics.stackexchange.com/questions/193/how-can-i-extract-normalized-read-count-values-from-deseq2-results
norm <- as.data.frame(counts(dds, normalized=T))
rsquare <- data.frame(matrix(ncol = ncol(norm), 
                             nrow = ncol(norm)))
for(i in 1:(ncol(norm)-1)) {       # for-loop over columns
  for(j in (i+1):ncol(norm)) { 
    print(paste(i," ",j))
    x <- cor.test(norm[ , i], 
                  norm[ , j], 
                  method = 'spearman')
    rsquare[i,j] <- x$estimate
  }
}
colnames(rsquare) <- colnames(norm)
rownames(rsquare) <- colnames(norm)
View(rsquare)
write_xlsx(rsquare, "./MF_Kallisto_dmrt1L_rsquare.xls")


# MF dmrt1S ----
colnames(counts)
new_counts <- counts[,c(6:8,20:22)]
colnames(new_counts)
coldata <- read.table("MF_sample_genotype.txt", header=T, row.names = 1)
new_coldata <- coldata[c(6:8,20:22),]; new_coldata
new_coldata$batch <- as.factor(new_coldata$batch)
new_coldata$sex <- as.factor(new_coldata$sex)
dds <- DESeqDataSetFromMatrix(countData = round(new_counts),
                              colData = new_coldata,
                              design= ~ sex)
# save the unfiltered log2FoldChange to a dataframe 
dds <- DESeq(dds)
dds$sex <- relevel(dds$sex, ref="F") # this makes the expression levels relative to F
res <- results(dds)
MF_dmrt1S_unfiltered <- res;MF_dmrt1S_unfiltered
# Only sex related
sex_related_MF_dmrt1S <- res[grepl("gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g34871\\|", rownames(res)), ]
write.csv(sex_related_MF_dmrt1S, file="Sex_related_MF_dmrt1S_Kallisto_DeSeq2_unfiltered.csv", row.names = T)

# Now do analysis of differential expression; 
# first remove transcripts where the average count per sample is 2 or less:
keep <- rowSums(counts(dds)) >= 2* length(colnames(dds))
dds <- dds[keep,]
dim(dds)
#[1] 11108     6
# relevel
dds$sex <- relevel(dds$sex, ref="F") # this makes the expression levels relative to F

# now do the analysis
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$pvalue),]
summary(res)
p<-resOrdered[1:237,];p
write.csv(p, file="MF_Kallisto_dmrt1S_DE_DeSeq2.csv", row.names = T)
# get Rsquare value for all pairwise comparisons
# normalized counts for DeSeq2:
# https://bioinformatics.stackexchange.com/questions/193/how-can-i-extract-normalized-read-count-values-from-deseq2-results
norm <- as.data.frame(counts(dds, normalized=T))
rsquare <- data.frame(matrix(ncol = ncol(norm), 
                             nrow = ncol(norm)))
for(i in 1:(ncol(norm)-1)) {       # for-loop over columns
  for(j in (i+1):ncol(norm)) { 
    print(paste(i," ",j))
    x <- cor.test(norm[ , i], 
                  norm[ , j], 
                  method = 'spearman')
    rsquare[i,j] <- x$estimate
  }
}
colnames(rsquare) <- colnames(norm)
rownames(rsquare) <- colnames(norm)
View(rsquare)
write_xlsx(rsquare, "./MF_Kallisto_dmrt1S_rsquare.xls")

# wtko dmw ----
# load count data (from Kalisto)
# read the count data from Kallisto that was combined from each sample
# into a single file into a dataframe called "dds"
dmw_counts <- read.table("dmw_only.isoform.counts.matrix", header=T, row.names = 1)
coldata <- read.table("dmw_sample_genotype.txt", header=T, row.names = 1)
coldata$genotype <- as.factor(coldata$genotype)
dds <- DESeqDataSetFromMatrix(countData = round(dmw_counts),
                              colData = coldata,
                              design= ~ genotype)

# save the unfiltered log2FoldChange to a dataframe 
dds$genotype <- relevel(dds$genotype, ref="wt") # this makes the expression levels relative to WT, and not KO
dds <- DESeq(dds)
res <- results(dds)
wtko_dmw_unfiltered <- res;wtko_dmw_unfiltered
# Only sex related
sex_related_wtko_dmw <- res[grepl("gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g34871\\|", rownames(res)), ]
write.csv(sex_related_wtko_dmw, file="Sex_related_wtko_dmw_Kallisto_DeSeq2_unfiltered.csv", row.names = T)


# Now do analysis of differential expression; 
# here we remove transcripts where the average count per sample is 2 or less:
keep <- rowSums(counts(dds)) >= 2* length(colnames(dds))
# here we remove transcripts where the average count per sample is 1or more
#keep <- rowSums(counts(dds)) >=length(colnames(dds))
dds <- dds[keep,]
dim(dds)
# [1] 32249    12
# relevel
dds$genotype <- relevel(dds$genotype, ref="wt") # this makes the expression levels relative to WT, and not KO
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$pvalue),]
summary(res)
p<-resOrdered[1:33,];p
write.csv(p, file="wtko_Kallisto_dmw_DE_DeSeq2.csv", row.names = T)

# get Rsquare value for all pairwise comparisons
# normalized counts for DeSeq2:
# https://bioinformatics.stackexchange.com/questions/193/how-can-i-extract-normalized-read-count-values-from-deseq2-results
norm <- as.data.frame(counts(dds, normalized=T))
rsquare <- data.frame(matrix(ncol = ncol(norm), 
                             nrow = ncol(norm)))
for(i in 1:(ncol(norm)-1)) {       # for-loop over columns
  for(j in (i+1):ncol(norm)) { 
    print(paste(i," ",j))
    x <- cor.test(norm[ , i], 
                  norm[ , j], 
                  method = 'spearman')
    rsquare[i,j] <- x$estimate
  }
}
colnames(rsquare) <- colnames(norm)
rownames(rsquare) <- colnames(norm)
View(rsquare)
write_xlsx(rsquare, "./wtko_Kallisto_dmw_rsquare.xls")

# wtko scan ----
scan_counts <- read.table("scanw.isoform.counts.matrix", header=T, row.names = 1)
coldata <- read.table("scanw_sample_genotype.txt", header=T, row.names = 1)
colnames(counts) <- row.names(coldata)
coldata$genotype <- as.factor(coldata$genotype)
colnames(scan_counts) <- row.names(coldata)

dds <- DESeqDataSetFromMatrix(countData = round(scan_counts),
                              colData = coldata,
                              design= ~ genotype)
# save the unfiltered log2FoldChange to a dataframe 
dds$genotype <- relevel(dds$genotype, ref="wt") # this makes the expression levels relative to WT, and not KO
dds <- DESeq(dds)
res <- results(dds)
wtko_scan_unfiltered <- res;wtko_scan_unfiltered
# Only sex related
sex_related_wtko_scan <- res[grepl("gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g34871\\|", rownames(res)), ]
write.csv(sex_related_wtko_scan, file="Sex_related_wtko_scan_Kallisto_DeSeq2_unfiltered.csv", row.names = T)


# Now do analysis of differential expression; 
# here we remove transcripts where the average count per sample is 2 or less:
keep <- rowSums(counts(dds)) >= 2* length(colnames(dds))
# here we remove transcripts where the average count per sample is 1or more
# keep <- rowSums(counts(dds)) >=length(colnames(dds))
dds <- dds[keep,]
dim(dds)
#[1] 32469     9
# relevel
dds$genotype <- factor(dds$genotype, levels = c("wt","knockout"))
dds$genotype <- relevel(dds$genotype, ref="wt") # this makes the expression levels relative to WT, and not KO

dds <- DESeq(dds)
res <- results(dds)
res
resOrdered <- res[order(res$pvalue),]
summary(res)
p<-resOrdered[1:34,];p
write.csv(p, file="wtko_Kallisto_scan_DE_DeSeq2.csv", row.names = T)

# get Rsquare value for all pairwise comparisons
# normalized counts for DeSeq2:
# https://bioinformatics.stackexchange.com/questions/193/how-can-i-extract-normalized-read-count-values-from-deseq2-results
norm <- as.data.frame(counts(dds, normalized=T))
rsquare <- data.frame(matrix(ncol = ncol(norm), 
                             nrow = ncol(norm)))
for(i in 1:(ncol(norm)-1)) {       # for-loop over columns
  for(j in (i+1):ncol(norm)) { 
    print(paste(i," ",j))
    x <- cor.test(norm[ , i], 
                  norm[ , j], 
                  method = 'spearman')
    rsquare[i,j] <- x$estimate
  }
}
colnames(rsquare) <- colnames(norm)
rownames(rsquare) <- colnames(norm)
View(rsquare)
write_xlsx(rsquare, "./wtko_Kallisto_scan_rsquare.xls")


# wtko ccdc ----
ccdc_counts <- read.table("ccdc.isoform.counts.matrix", header=T, row.names = 1)
coldata <- read.table("ccdc_sample_genotype.txt", header=T, row.names = 1)
colnames(ccdc_counts) <- row.names(coldata)
coldata$genotype <- as.factor(coldata$genotype)

# subset the data to keep only females
ccdc_counts_fem_only <- ccdc_counts[,c(2,4:10)]
coldata_fem_only <- coldata[c(2,4:10),]
dds <- DESeqDataSetFromMatrix(countData = round(ccdc_counts_fem_only),
                              colData = coldata_fem_only,
                              design= ~ genotype)
# save the unfiltered log2FoldChange to a dataframe 
dds$genotype <- relevel(dds$genotype, ref="wt") # this makes the expression levels relative to WT, and not KO
dds <- DESeq(dds)
res <- results(dds)
wtko_ccdc_unfiltered <- res;wtko_ccdc_unfiltered
# Only sex related
sex_related_wtko_ccdc <- res[grepl("gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g34871\\|", rownames(res)), ]
write.csv(sex_related_wtko_ccdc, file="Sex_related_wtko_ccdc_Kallisto_DeSeq2_unfiltered.csv", row.names = T)

# Now do analysis of differential expression; 
# here we remove transcripts where the average count per sample is 2 or less:
keep <- rowSums(counts(dds)) >= 2* length(colnames(dds))
dds <- dds[keep,]
dim(dds)
# [1] 32326     8
# relevel
dds$genotype <- factor(dds$genotype, levels = c("wt","knockout"))
dds$genotype <- relevel(dds$genotype, ref="wt") # this makes the expression levels relative to WT, and not KO
dds <- DESeq(dds)
res <- results(dds)
res
resOrdered <- res[order(res$pvalue),]
summary(res)
p<-resOrdered[1:128,];p
write.csv(p, file="wtko_Kallisto_ccdc_DE_DeSeq2.csv", row.names = T)

# get Rsquare value for all pairwise comparisons
# normalized counts for DeSeq2:
# https://bioinformatics.stackexchange.com/questions/193/how-can-i-extract-normalized-read-count-values-from-deseq2-results
norm <- as.data.frame(counts(dds, normalized=T))
rsquare <- data.frame(matrix(ncol = ncol(norm), 
                             nrow = ncol(norm)))
for(i in 1:(ncol(norm)-1)) {       # for-loop over columns
  for(j in (i+1):ncol(norm)) { 
    print(paste(i," ",j))
    x <- cor.test(norm[ , i], 
                  norm[ , j], 
                  method = 'spearman')
    rsquare[i,j] <- x$estimate
  }
}
colnames(rsquare) <- colnames(norm)
rownames(rsquare) <- colnames(norm)
View(rsquare)
write_xlsx(rsquare, "./wtko_Kallisto_ccdc_rsquare.xls")

# permutations ----
MFcounts <- read.table("MF_.isoform.TMM.EXPR.matrix", header=T, row.names = 1)

# get rownames of sexrelated transcripts
SL_rownames <- row.names(MFcounts[grepl("gnl\\|XBXL10_1g34625\\||gnl\\|XBXL10_1g37811\\||gnl\\|XBXL10_1g10668\\||gnl\\|XBXL10_1g7999\\||gnl\\|XBXL10_1g24241\\||gnl\\|XBXL10_1g26060\\||gnl\\|XBXL10_1g24554\\||gnl\\|XBXL10_1g26280\\||gnl\\|XBXL10_1g27265\\||gnl\\|XBXL10_1g29076\\||gnl\\|XBXL10_1g30057\\||gnl\\|XBXL10_1g32392\\||gnl\\|XBXL10_1g31301\\||gnl\\|XBXL10_1g33473\\||gnl\\|XBXL10_1g3211\\||gnl\\|XBXL10_1g5748\\||gnl\\|XBXL10_1g35876\\||gnl\\|XBXL10_1g37293\\||gnl\\|XBXL10_1g6566\\||gnl\\|XBXL10_1g8966\\||gnl\\|XBXL10_1g7278\\||gnl\\|XBXL10_1g10089\\||gnl\\|XBXL10_1g23152\\||gnl\\|XBXL10_1g25243\\||gnl\\|XBXL10_1g2070\\||gnl\\|XBXL10_1g4848\\||gnl\\|XBXL10_1g8430\\||gnl\\|XBXL10_1g11002\\||gnl\\|XBXL10_1g30252\\||gnl\\|XBXL10_1g32546\\||gnl\\|XBXL10_1g605\\||gnl\\|XBXL10_1g3639\\||gnl\\|XBXL10_1g37486\\||gnl\\|XBXL10_1g39526\\||gnl\\|XBXL10_1g42722\\||gnl\\|XBXL10_1g35158\\||gnl\\|XBXL10_1g38013\\||gnl\\|XBXL10_1g38893\\||gnl\\|XBXL10_1g42158\\||gnl\\|XBXL10_1g39443\\||gnl\\|XBXL10_1g42662\\||gnl\\|XBXL10_1g41173\\||gnl\\|XBXL10_1g43880\\||gnl\\|XBXL10_1g19698\\||gnl\\|XBXL10_1g22028\\||gnl\\|XBXL10_1g815\\||gnl\\|XBXL10_1g3800\\||gnl\\|XBXL10_1g8007\\||gnl\\|XBXL10_1g10675\\||gnl\\|XBXL10_1g2154\\||gnl\\|XBXL10_1g4928\\||gnl\\|XBXL10_1g27310\\||gnl\\|XBXL10_1g29128\\||gnl\\|XBXL10_1g40425\\||gnl\\|XBXL10_1g43291\\||gnl\\|XBXL10_1g8118\\||gnl\\|XBXL10_1g10760\\||gnl\\|XBXL10_1g8117\\||gnl\\|XBXL10_1g10758\\||gnl\\|XBXL10_1g1634\\||gnl\\|XBXL10_1g4460\\||gnl\\|XBXL10_1g22534\\||gnl\\|XBXL10_1g25047\\||gnl\\|XBXL10_1g22535\\||gnl\\|XBXL10_1g25046\\||gnl\\|XBXL10_1g13810\\||gnl\\|XBXL10_1g15286\\||gnl\\|XBXL10_1g30377\\||gnl\\|XBXL10_1g13205\\||gnl\\|XBXL10_1g15724\\||gnl\\|XBXL10_1g6054\\||gnl\\|XBXL10_1g9274\\||gnl\\|XBXL10_1g29226\\||gnl\\|XBXL10_1g34871\\|", rownames(MFcounts)), ])

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
  indexes <- sample.int(dim(MFcounts)[1], 74, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
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
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(wtko_dmw_trim[,'log2FoldChange'], # ko:wt first
             MF_ccdc_trim[,'log2FoldChange'], # reference M:F second
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
# [1] 0.01298701
# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_wtko_dmw_trim[SL_rownames,'log2FoldChange'],
           sex_related_MF_ccdc_trim[SL_rownames,'log2FoldChange'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.01198801


# MF_dmrt1L vs dmw
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 74, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
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
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(wtko_dmw_trim[,'log2FoldChange'], # ko:wt first
             MF_dmrt1L_trim[,'log2FoldChange'], # reference M:F second
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
# [1] 0.3906094
# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_wtko_dmw_trim[SL_rownames,'log2FoldChange'],
           sex_related_MF_dmrt1L_trim[SL_rownames,'log2FoldChange'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.1918082



# MF_dmrt1S vs dmw
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 74, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
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
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(wtko_dmw_trim[,'log2FoldChange'], # ko:wt first
             MF_dmrt1S_trim[,'log2FoldChange'], # reference M:F second
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
# [1] 0.000999001

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_wtko_dmw_trim[SL_rownames,'log2FoldChange'],
           sex_related_MF_dmrt1S_trim[SL_rownames,'log2FoldChange'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.3236763



# scan permutations ----

# MF_ccdc vs scan
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 74, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
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
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(wtko_scan_trim[,'log2FoldChange'], # ko:wt first
             MF_ccdc_trim[,'log2FoldChange'], # reference M:F second
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
# [1] 0.3396603

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_wtko_scan_trim[SL_rownames,'log2FoldChange'],
           sex_related_MF_ccdc_trim[SL_rownames,'log2FoldChange'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.06293706



# MF_dmrt1L vs scan
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 74, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
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
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(wtko_scan_trim[,'log2FoldChange'], # ko:wt first
             MF_dmrt1L_trim[,'log2FoldChange'], # reference M:F second
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
# [1] 0.9240759

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_wtko_scan_trim[SL_rownames,'log2FoldChange'],
           sex_related_MF_dmrt1L_trim[SL_rownames,'log2FoldChange'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.09490509


# MF_dmrt1S vs scan
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 74, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
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
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(wtko_scan_trim[,'log2FoldChange'], # ko:wt first
             MF_dmrt1S_trim[,'log2FoldChange'], # reference M:F second
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
# [1] 0.3896104
# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_wtko_scan_trim[SL_rownames,'log2FoldChange'],
           sex_related_MF_dmrt1S_trim[SL_rownames,'log2FoldChange'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.4545455



# ccdc permutations ----

# MF_ccdc vs ccdc
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 74, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
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
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(wtko_ccdc_trim[,'log2FoldChange'], # ko:wt first
             MF_ccdc_trim[,'log2FoldChange'], # reference M:F second
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
# [1] 0.5474525
# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_wtko_ccdc_trim[SL_rownames,'log2FoldChange'],
           sex_related_MF_ccdc_trim[SL_rownames,'log2FoldChange'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.06393606


# MF_dmrt1L vs ccdc
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 74, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
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
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(wtko_ccdc_trim[,'log2FoldChange'], # ko:wt first
             MF_dmrt1L_trim[,'log2FoldChange'], # reference M:F second
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
# [1] 0.4575425

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_wtko_ccdc_trim[SL_rownames,'log2FoldChange'],
           sex_related_MF_dmrt1L_trim[SL_rownames,'log2FoldChange'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.2177822




# MF_dmrt1S vs ccdc
correlations <- c()
magnitudes <- c()

# Use a for loop
for (x in 1:1000) {
  indexes <- sample.int(dim(MFcounts)[1], 74, replace = F);indexes
  rownames <- row.names(MFcounts[indexes,])
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
  # calculate and add the ratio of vector lengths to a vector
  a <- merge(wtko_ccdc_trim[,'log2FoldChange'], # ko:wt first
             MF_dmrt1S_trim[,'log2FoldChange'], # reference M:F second
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
# [1] 0.4195804

# now figure out where the observed magnitude ratio is within the permutation magntiude vector
a <- merge(sex_related_wtko_ccdc_trim[SL_rownames,'log2FoldChange'],
           sex_related_MF_dmrt1S_trim[SL_rownames,'log2FoldChange'],
           by = 'row.names', 
           incomparables = NA)
b <- a[complete.cases(a), ];b
magnitudes[1001] <- ang.vec.alph(b$x,b$y)[3]

print("Magnitude pvalue: "); 1-rank(magnitudes)[1001]/1001
# [1] 0.6273726

```
# plotting
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
