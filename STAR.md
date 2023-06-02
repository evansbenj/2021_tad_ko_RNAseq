# Path
```
/home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq
```

# Index genome

```
#!/bin/sh
#SBATCH --job-name=STAR_index
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=256gb
#SBATCH --output=STAR_index.%J.out
#SBATCH --error=STAR_index.%J.err
#SBATCH --account=def-ben

module load StdEnv/2020 star/2.7.9a

STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir /home/ben/projects/rrg-ben/ben/2020_XL_v9.2_refgenome/ \
--genomeFastaFiles /home/ben/projects/rrg-ben/ben/2020_XL_v9.2_refgenome/XENLA_9.2_genome.fa \
--sjdbGTFfile /home/ben/projects/rrg-ben/ben/2020_XL_v9.2_refgenome/XLv9.2_xenbase_annotations.gff \
--sjdbOverhang 99 \
--limitGenomeGenerateRAM=124544990592
```

# Map reads
```
#!/bin/sh
#SBATCH --job-name=STAR_map
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=64gb
#SBATCH --output=STAR_map.%J.out
#SBATCH --error=STAR_map.%J.err
#SBATCH --account=def-ben

module load StdEnv/2020 star/2.7.9a

STAR --genomeDir /home/ben/projects/rrg-ben/ben/2020_XL_v9.2_refgenome/ \
--runThreadN 6 \
--readFilesIn ${1} ${2} \
--outFileNamePrefix ${3} \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--readFilesCommand zcat
```

# Count
First use samtools sort to sort all the bam files

featureCounts package  
http://subread.sourceforge.net/  
https://subread.sourceforge.net/SubreadUsersGuide.pdf

```
#!/bin/sh
#SBATCH --job-name=STAR_count
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --time=3:00:00
#SBATCH --mem=4Gb
#SBATCH --output=STAR_count.%J.out
#SBATCH --error=STAR_count.%J.err
#SBATCH --account=def-ben

# sbatch 2022_STAR_count.sh inputbam output_counts

module load StdEnv/2020 gcc/9.3.0 star/2.7.9a samtools subread/2.0.3

# must use -s 0 because the data are unstranded
# must use -p because the data are paired
# use --countReadPairs to count read pairs instead of reads
# use -C to prevent counting of chimeric reads
# -T is the number of threads
featureCounts -T 4 -s 0 -p --countReadPairs -C -g ID \
  -a /home/ben/projects/rrg-ben/ben/2020_XL_v9.2_refgenome/XLv9.2_xenbase_annotations.gff \
  -o ${2} \
  ${1}
  
```

# Combine and analyze

I downloded the counts for each individual and combined and analyzed them using this R script. 

```R
# if (!requireNamespace('edgeR', quietly = T)) {
#   install.packages('edgeR')
# }
# if (!requireNamespace('rhdf5')) {
#   BiocManager::install("rhdf5")
# }
# if (!requireNamespace('tximportData')) {
#   BiocManager::install('tximportData')
# }
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install('PCAtools')
# BiocManager::install("org.Xl.eg.db")
# BiocManager::install("limma")
# BiocManager::install("edgeR")
# BiocManager::install("tximport")
# BiocManager::install("HTSFilter")
# BiocManager::install("rhdf5")
library(limma)
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


setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_Supergene/2022_KO_tad_RNAseq/2022_EdgeR_and_DeSeq2/STAR")
dir <- "/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_Supergene/2022_KO_tad_RNAseq/2022_EdgeR_and_DeSeq2/STAR"
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
View(countz)
# samples
samples <- read.table(file.path(dir, "samples.txt"), header = T)
samples


# sexez
sexez <- factor(samples$sex)
sexez <- relevel(sexez, ref="F")

# batch
batchez <- factor(samples$batch)

# Trim the data to leave only the samples under consideration
# MvsFwt
new_counts <- as.data.frame(countz[,-c(4,8:9,11:13,15,18:23,26:27,31:33,35,37,40,43,45,46:66) ])
row.names(new_counts) <- gene_names
new_samples <-as.data.frame(samples[-c(4,8:9,11:13,15,18:23,26:27,31:33,35,37,40,43,45,46:66), ]);new_samples

new_sexez <- factor(samples$sex[-c(4,8:9,11:13,15,18:23,26:27,31:33,35,37,40,43,45,46:66)])
new_sexez <- relevel(new_sexez, ref="F")
# batch
new_batchez <- factor(samples$batch[-c(4,8:9,11:13,15,18:23,26:27,31:33,35,37,40,43,45,46:66)])

# Create DGEList object - this is a data structure that is used for 
# the analysis of differential expression
d0 <- DGEList(new_counts, group = new_sexez, remove.zeros = TRUE)

# get rid of rows where there is less than an average of two read per individual
# in the analysis
dim(d0$counts) #each row is a transcript - here is the number before filtering
# [1] 37380    22 # 2023 STAR MF
d0$counts <- d0$counts[rowSums(d0$counts)> 2* ncol(d0$counts),] 
# Now we have far fewer transcripts:
dim(d0$counts)
# [1] 29387    22 # 2023 STAR MF

# TMM normalization is applied to this dataset to account for compositional difference between
# the libraries.
d0 <- calcNormFactors(d0, method="TMM")
# check the normalization factors
d0$samples

# plot by genotype
plotMDS(d0,labels=new_samples$sex,col=c(rep("green",8),rep("blue",8),rep("red",6)))

# explorative analysis
design <- model.matrix(~new_batchez + new_batchez:new_sexez)
logFC <- predFC(d0,design,prior.count=1,dispersion=0.05)
cor(logFC[,4:6])
#                                new_batchezccdc:new_sexezM new_batchezdmrt1L:new_sexezM new_batchezdmrt1S:new_sexezM
# new_batchezccdc:new_sexezM                   1.00000000    0.07570747   -0.10293765
# new_batchezdmrt1L:new_sexezM                 0.07570747    1.00000000   -0.08498052
# new_batchezdmrt1S:new_sexezM                -0.10293765    -0.08498052   1.00000000

# the correlation between ccdc and dmrt1L is positive but the
# correlations between dmrt1S and both of these is negative!!!

# design matrix: this is used for the model of differential expression
design <- model.matrix(~ new_batchez + new_sexez, data=d0$samples) # last coefficient = difference between sexes)
design

# estimate dispersion
d0 <- estimateDisp(d0, design, robust=TRUE)
d0$common.dispersion
#d0$tagwise.dispersion
plotBCV(d0)

fit <- glmQLFit(d0, design, robust=TRUE)
plotQLDisp(fit)

# check if there is a substantial effect of batch
qlf <- glmQLFTest(fit, coef=2:3)
topTags(qlf)
summary(decideTests(qlf))
# there are thousands of DE genes between batches so we need to 
# take this into account (obviously as well based on PCA plots)

# now do test for effect of sex.  Be default the test is for the last
# coefficient in the design matrix, which in this case is the treatment effect
qlf <- glmQLFTest(fit)
summary(decideTests(qlf), p.value = 0.1)
#        new_sexezM
# Down            1
# NotSig      29386
# Up              0
topTags(qlf, n=4)
#                   logFC     logCPM        F       PValue         FDR
# XBXL10_1g8729  -8.142117 -0.3168759 94.61785 4.456324e-10 1.30958e-05 (this is dmw!)
MF_ccdc_dmrt1L_dmrt1S_DE <- as.data.frame(topTags(exacttest, n=2))
write.csv(MF_ccdc_dmrt1L_dmrt1S_DE, file="MF_STAR_ccdcdmrt1Ldmrt1S_DE_edgeR.csv", row.names = T)



# because the dmrt1S batch was negatively correlated; lets remove
# it and try again. This could be because the dmrt1S tads were at a 
# slightly different developmental stage
colnames(countz)

new_counts <- as.data.frame(countz[,-c(4,8:9,11:13,15,18:23,26:27,31:33,35:66) ])
row.names(new_counts) <- gene_names
new_samples <-as.data.frame(samples[-c(4,8:9,11:13,15,18:23,26:27,31:33,35:66), ]);new_samples

new_sexez <- factor(samples$sex[-c(4,8:9,11:13,15,18:23,26:27,31:33,35:66)])
new_sexez <- relevel(new_sexez, ref="F")
# batch
new_batchez <- factor(samples$batch[-c(4,8:9,11:13,15,18:23,26:27,31:33,35:66)])

# Create DGEList object - this is a data structure that is used for 
# the analysis of differential expression
d0 <- DGEList(new_counts, group = new_sexez, remove.zeros = TRUE)
dim(d0$counts) #each row is a transcript - here is the number before filtering
# [1] 36632    16 # 2023 STAR ccdc dmrt1L no dmrt1S
# here we remove transcripts where the average count per sample is 2 or less:
d0$counts <- d0$counts[rowSums(d0$counts)> 2* ncol(d0$counts),] 
# Now we have far fewer transcripts:
dim(d0$counts)
# [1] 29324    16 # 2023 STAR ccdc dmrt1L no dmrt1S
# many rows with low expression were eliminated
# TMM normalization is applied to this dataset to account for compositional difference between
# the libraries.
d0 <- calcNormFactors(d0, method="TMM")
# check the normalization factors
d0$samples
# plot by sex
plotMDS(d0,labels=c(rep("M",5),rep("F",2),rep("M",3),"F",rep("M",2),rep("F",2),"M"),
        col=c(rep("green",8),rep("blue",8)))
# design matrix: this is used for the model of differential expression
design <- model.matrix(~ new_batchez + new_sexez, data=d0$samples) # last coefficient = difference between sexes)
design

# estimate dispersion
d0 <- estimateDisp(d0, design, robust=TRUE)
d0$common.dispersion
#d0$tagwise.dispersion
plotBCV(d0)

fit <- glmQLFit(d0, design, robust=TRUE)
plotQLDisp(fit)

# check if there is a substantial effect of batch
qlf <- glmQLFTest(fit, coef=2:3)
topTags(qlf)
summary(decideTests(qlf))
# there are thousands of DE genes between batches so we need to 
# take this into account (obviously as well based on PCA plots)

# now do test for effect of sex.  Be default the test is for the last
# coefficient in the design matrix, which in this case is the treatment effect
qlf <- glmQLFTest(fit)
summary(decideTests(qlf, p.value = 0.10))
#         sexezm
# Down            1
# NotSig      11054
# Up              0
topTags(qlf, n=2)
# Coefficient:  new_sexezm 
# logFC   logCPM        F       PValue        FDR
# XBXL10_1g8729  -7.384947 -0.8953262 69.36840 5.379099e-08 0.001577367 # this is dmw





# Try dmrt1S alone
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
plotMDS(d0,labels=c(rep("M",3),rep("F",3)),
        col=c(rep("green",3),rep("blue",3)))
# design matrix: this is used for the model of differential expression
design <- model.matrix(~ 0 + new_sexez, data=d0$samples) # last coefficient = difference between sexes)
design
# estimate dispersion
d0 <- estimateDisp(d0, design, robust=TRUE)
d0$common.dispersion
#d0$tagwise.dispersion
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
summary(decideTests(object = exacttest, p.value = 0.1))
# M-F
# Down     209
# NotSig 28443
# Up       266
dmrt1S_DE <- as.data.frame(topTags(exacttest, n=475))
library(rJava)
library(xlsx)
#write.xlsx(dmrt1S_DE, "./MF_STAR_dmrt1S_DE_edgeR.xls", row.names=TRUE)
write.csv(dmrt1S_DE, file="MF_STAR_dmrt1S_DE_edgeR.csv", row.names = T)
# Coefficient:  new_sexezm 
# logFC   logCPM        F       PValue        FDR
# XBXL10_1g8729  -7.384947 -0.8953262 69.36840 5.379099e-08 0.001577367 # this is dmw


# Try dmrt1L alone
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
plotMDS(d0,labels=c(rep("M",2),"F",rep("M",2),rep("F",2),"M"),
        col=c(rep("blue",2),"red",rep("blue",2),rep("red",2),"blue"))
# design matrix: this is used for the model of differential expression
design <- model.matrix(~ 0 + new_sexez, data=d0$samples) # last coefficient = difference between sexes)
design
# estimate dispersion
d0 <- estimateDisp(d0, design, robust=TRUE)
d0$common.dispersion
#d0$tagwise.dispersion
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
summary(decideTests(object = exacttest, p.value = 0.1))
# M-F
# Down      26
# NotSig 28689
# Up         6
topTags(exacttest, n=33)
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
dmrt1L_DE <- as.data.frame(topTags(exacttest, n=33))
write.csv(dmrt1L_DE, file="MF_STAR_dmrt1L_DE_edgeR.csv", row.names = T)

# Let's try to get the logFC for all the genes in MF dmrt1L 
# that are sigDE in dmrt1S
is.data.frame(exacttest$table)
# first read in the file from the dmrt1S
dmrt1_DE_list <- read.csv(file.path(dir, "MF_STAR_dmrt1S_DE_edgeR.csv"), header = T)
dmrt1_DE_list_names <- dmrt1_DE_list$X
write.csv(exacttest$table[dmrt1_DE_list_names,], file="MF_STAR_edgeR_dmrt1Lexpression_of_dmrt1S_DEs.csv", row.names = T)


# Try ccdc alone
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
plotMDS(d0,labels=c(rep("M",5),rep("F",2),"M"),
        col=c(rep("blue",5),rep("red",2),"blue"))
# design matrix: this is used for the model of differential expression
design <- model.matrix(~ 0 + new_sexez, data=d0$samples) # last coefficient = difference between sexes)
design
# estimate dispersion
d0 <- estimateDisp(d0, design, robust=TRUE)
d0$common.dispersion
#d0$tagwise.dispersion
exacttest <- exactTest(d0, dispersion = "auto") # no differentially expressed genes
summary(decideTests(object = exacttest, p.value = 0.1))
# M-F
# Down       2
# NotSig 29070
#  Up         1
topTags(exacttest, n=4)
# Coefficient:  new_sexezm 
# logFC   logCPM        F       PValue        FDR
# XBXL10_1g3050  -7.520457 -0.6587820 1.307588e-06 0.01585958
# XBXL10_1g43040 -8.439321 -1.6729722 1.477786e-06 0.01585958
# XBXL10_1g3473   9.507057  0.7988971 1.636526e-06 0.01585958
ccdc_DE <- as.data.frame(topTags(exacttest, n=33))
write.csv(ccdc_DE, file="MF_STAR_ccdc_DE_edgeR.csv", row.names = T)









# OK now try dmw wt vs ko
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
plotMDS(d0,labels=c("wt","ko","wt","ko","wt","wt","wt","ko","ko","ko","ko","wt"),
        col=c("blue","green","blue","green","blue","blue","blue","green","green","green","green","blue"))
# design matrix: this is used for the model of differential expression
design <- model.matrix(~ 0 + new_genotypez, data=d0$samples) # last coefficient = difference between sexes)
design
# estimate dispersion
d0 <- estimateDisp(d0, design, robust=TRUE)
d0$common.dispersion
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
```
