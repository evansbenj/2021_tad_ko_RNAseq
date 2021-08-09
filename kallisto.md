# Generate counts
```
#!/bin/sh
#SBATCH --job-name=kallisto
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=32gb
#SBATCH --output=kallisto.%J.out
#SBATCH --error=kallisto.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this
# sbatch ./2021_kallisto_withBoot.sh ../raw_data/dmrt1

module load kallisto/0.46.1

#  Always use for-loop, prefix glob, check if exists file.
for file in $1/*R1_trim_001.fastq.gz ; do         # Use ./* ... NEVER bare *
dmw_35_S36_L001_R2_trim_001.fastq.gz
    if [ -e "$file" ] ; then   # Check whether file exists.
      kallisto quant -b 100 -i /home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNA
seq/XL_transcriptome/xlaevisMRNA.idx -o ${file::-26}_kallisto_boot_out ${file::-
26}_L001_R1_trim_001.fastq.gz ${file::-26}_L001_R2_trim_001.fastq.gz ${file::-26
}_L002_R1_trim_001.fastq.gz ${file::-26}_L002_R2_trim_001.fastq.gz
  fi
done 
```

# Combine counts from multiple individuals
```
#!/bin/sh
#SBATCH --job-name=kallisto
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=32gb
#SBATCH --output=kallisto.%J.out
#SBATCH --error=kallisto.%J.err
#SBATCH --account=def-ben

module load r/4.0.5
module load trinity/2.11.0

# run by passing an argument like this
# sbatch 2021_Trinity_combine_kallisto.sh

/home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/ben_scripts/trinityrnaseq-v
2.12.0/util/abundance_estimates_to_matrix.pl --est_method kallisto --out_prefix 
dmw_ccdc_dmrt1L  --gene_trans_map none --name_sample_by_basedir /home/ben/projec
ts/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_data/dmw/dmw_14_S29_kallisto_boot_out/a
bundance.tsv /home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_data/dmw/d
mw_16_S30_kallisto_boot_out/abundance.tsv /home/ben/projects/rrg-ben/ben/2021_XL
_ko_tad_RNAseq/raw_data/dmw/dmw_17_S31_kallisto_boot_out/abundance.tsv /home/ben
/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_data/dmw/dmw_20_S32_kallisto_boo
t_out/abundance.tsv /home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_dat
a/dmw/dmw_26_S33_kallisto_boot_out/abundance.tsv /home/ben/projects/rrg-ben/ben/
2021_XL_ko_tad_RNAseq/raw_data/dmw/dmw_28_S34_kallisto_boot_out/abundance.tsv /h
ome/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_data/dmw/dmw_29_S35_kalli
sto_boot_out/abundance.tsv /home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/
raw_data/dmw/dmw_35_S36_kallisto_boot_out/abundance.tsv /home/ben/projects/rrg-b
en/ben/2021_XL_ko_tad_RNAseq/raw_data/ccdc/ccdc_12_S3_kallisto_boot_out/abundanc
e.tsv /home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_data/ccdc/ccdc_14
_S4_kallisto_boot_out/abundance.tsv /home/ben/projects/rrg-ben/ben/2021_XL_ko_ta
d_RNAseq/raw_data/ccdc/ccdc_25_S5_kallisto_boot_out/abundance.tsv /home/ben/proj
ects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_data/ccdc/ccdc_30_S6_kallisto_boot_ou
t/abundance.tsv /home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_data/cc
dc/ccdc_32_S7_kallisto_boot_out/abundance.tsv /home/ben/projects/rrg-ben/ben/202
1_XL_ko_tad_RNAseq/raw_data/ccdc/ccdc_34_S8_kallisto_boot_out/abundance.tsv /hom
e/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_data/ccdc/ccdc_35_S9_kallis
to_boot_out/abundance.tsv /home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/r
aw_data/ccdc/ccdc_36_S10_kallisto_boot_out/abundance.tsv /home/ben/projects/rrg-
ben/ben/2021_XL_ko_tad_RNAseq/raw_data/ccdc/ccdc_3_S1_kallisto_boot_out/abundanc
e.tsv /home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_data/ccdc/ccdc_42
_S11_kallisto_boot_out/abundance.tsv /home/ben/projects/rrg-ben/ben/2021_XL_ko_t
ad_RNAseq/raw_data/ccdc/ccdc_9_S2_kallisto_boot_out/abundance.tsv /home/ben/proj
ects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_data/dmrt1L/dmrt1L_11_S15_kallisto_bo
ot_out/abundance.tsv /home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_da
ta/dmrt1L/dmrt1L_17_S16_kallisto_boot_out/abundance.tsv /home/ben/projects/rrg-b
en/ben/2021_XL_ko_tad_RNAseq/raw_data/dmrt1L/dmrt1L_19_S17_kallisto_boot_out/abu
ndance.tsv /home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_data/dmrt1L/
dmrt1L_24_S18_kallisto_boot_out/abundance.tsv /home/ben/projects/rrg-ben/ben/202
1_XL_ko_tad_RNAseq/raw_data/dmrt1L/dmrt1L_25_S19_kallisto_boot_out/abundance.tsv
 /home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_data/dmrt1L/dmrt1L_26_
S20_kallisto_boot_out/abundance.tsv /home/ben/projects/rrg-ben/ben/2021_XL_ko_ta
d_RNAseq/raw_data/dmrt1L/dmrt1L_27_S21_kallisto_boot_out/abundance.tsv /home/ben
/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_data/dmrt1L/dmrt1L_30_S22_kallis
to_boot_out/abundance.tsv /home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/r
aw_data/dmrt1L/dmrt1L_35_S23_kallisto_boot_out/abundance.tsv /home/ben/projects/
rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_data/dmrt1L/dmrt1L_41_S24_kallisto_boot_ou
t/abundance.tsv /home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_data/dm
rt1L/dmrt1L_43_S25_kallisto_boot_out/abundance.tsv /home/ben/projects/rrg-ben/be
n/2021_XL_ko_tad_RNAseq/raw_data/dmrt1L/dmrt1L_50_S26_kallisto_boot_out/abundanc
e.tsv /home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_data/dmrt1L/dmrt1
L_55_S27_kallisto_boot_out/abundance.tsv /home/ben/projects/rrg-ben/ben/2021_XL_
ko_tad_RNAseq/raw_data/dmrt1L/dmrt1L_59_S28_kallisto_boot_out/abundance.tsv /hom
e/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_data/dmrt1L/dmrt1L_6_S12_ka
llisto_boot_out/abundance.tsv /home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAs
eq/raw_data/dmrt1L/dmrt1L_7_S13_kallisto_boot_out/abundance.tsv /home/ben/projec
ts/rrg-ben/ben/2021_XL_ko_tad_RNAseq/raw_data/dmrt1L/dmrt1L_8_S14_kallisto_boot_
out/abundance.tsv
```
# Analysis of differential expression (for dmw only)
I assembled the dmw transcriptome from 12 individuals - 6 wt female and 6 dmw ko female.  I used this script to analyze differentially expressed genes:
```
if (!requireNamespace('edgeR', quietly = T)) {
  install.packages('edgeR')
}

if (!requireNamespace('rhdf5')) {
  BiocManager::install("rhdf5")
}
if (!requireNamespace('tximportData')) {
  BiocManager::install('tximportData')
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('PCAtools')

BiocManager::install("org.Xl.eg.db")

BiocManager::install("limma")





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

setwd("/Users/Shared/Previously Relocated Items/Security/projects/2021_KO_tad_RNAseq/Kallisto_EdgeR")
dir <- "/Users/Shared/Previously Relocated Items/Security/projects/2021_KO_tad_RNAseq/Kallisto_EdgeR"
list.files(dir)

# load count data (from Kalisto)
# read the count data from Kallisto that was combined from each sample
# into a single file into a dataframe called "counts"
counts <- read.table("dmw_only.isoform.counts.matrix", header=T, row.names = 1)
dim(counts)
# [1] 704350     12
# get rid of any rows that have incomplete data
counts <- counts[complete.cases(counts), ]
dim(counts)
# [1] 704350     12

# this is useful info about batch effects:
# https://support.bioconductor.org/p/96627/

# make design matrix
# Here we want to test for differential expression between KO and
# WT, while adjusting for differences between batches. In statistical
# terms, this "may be" an additive linear model with batch as the blocking factor:

# samples
samples <- read.table(file.path(dir, "dmw_samples.txt"), header = F)
samples

# rename the col names
colnames(counts) <- samples$V1

# batch
batch <- factor(c("first","first","first","first","first","first",
                  "second","second","second","first","first","second"))
batch <- relevel(batch, ref="first")

# genotypes
genotypez <- factor(c("dmwKO","dmwKO","dmwKO","dmwKO","dmwKO","dmwKO",
                      "WT","WT","WT","WT","WT","WT"))
genotypez <- relevel(genotypez, ref="WT")

# sexez
sexez <- factor(c("F","F","F","F","F","F","F","F","F","F","F","F"))
sexez <- relevel(sexez, ref="F")

grouping <- factor(paste(batch,genotypez, sep=".")) # ignore sex because all dmw samples are female

# Create DGEList object
d0 <- DGEList(counts, group = grouping)

# get rid of rows where there is less than an average of one read per individual
# in the analysis
# this is an important step - lowly expressed genes tend to have very variable
# expression and are unreliable.  THis can also result from assembly errors
# first get rid of scientific notation
dim(d0$counts)
# 704350     12
d0$counts <- d0$counts[rowSums(d0$counts)> ncol(d0$counts),]
dim(d0$counts)
# 397895     12
# many rows with low expression were eliminated

# design matrix
design<-model.matrix(~0+batch+genotypez, data=d0$samples) # last coefficient = difference between genotypez
design

# TMM normalization is applied to this dataset to account for compositional difference between
# the libraries.
d0 <- calcNormFactors(d0, method="TMM")
# check the normalization factors
d0$samples

plotMDS(d0)
plotMDS(d0,labels=NULL,col=c(rep("green",6),rep("blue",6)))
#dmw_15 and dmw_28 seem weird

# pca
project.pca <- prcomp(t(d0$counts))
plot(project.pca$x, type="n")
points(project.pca$x, col=genotypez, pch=16, cex=1)
#points(project.pca$x, col=sexez, pch=16, cex=1)
points(project.pca$x, col=batch, pch=16, cex=1) +
theme(legend.position="top")
var_explained <- project.pca$sdev^2/sum(project.pca $sdev^2)

# OK let's compare wt F to dmw_ko_F
# don't include dmw_wt_m because we don't have any males from that batch


# pca
project.pca <- prcomp(t(d0$counts))
var_explained <- project.pca$sdev^2/sum(project.pca $sdev^2)

library(tidyverse)
project.pca$x %>% 
  as.data.frame %>%
  # rownames_to_column("continent_country") %>%
  # separate(continent_country,c("continent")) %>%
  ggplot(aes(x=PC1,y=PC2)) + geom_point(aes(color=genotypez),size=4) +
  theme_bw(base_size=12) + 
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme(legend.position="top")

# pca mostly overlaps...

# estimate dispersion
d0 <- estimateDisp(d0, design, robust=TRUE)
d0$common.dispersion
plotBCV(d0)

# fit the model
fit <- glmQLFit(d0, design, robust=TRUE)
plotQLDisp(fit)
qlf <- glmQLFTest(fit, coef=3)
topTags(qlf)
#TRINITY_DN15246_c0_g1_i10  11.745462 3.9729331 75.11365 6.771376e-08 0.01586036
#TRINITY_DN5299_c0_g1_i4    11.292092 2.9022767 61.46581 1.120121e-07 0.01586036
#TRINITY_DN2367_c0_g1_i1    -5.490647 2.0221252 58.63102 1.195820e-07 0.01586036
#TRINITY_DN209_c5_g1_i3      9.037564 2.7173802 56.41456 1.633396e-07 0.01624800
#TRINITY_DN1063_c0_g1_i6    11.181356 2.8607074 56.16741 2.487316e-07 0.01979381
#TRINITY_DN3375_c0_g1_i6   -10.581700 0.1669082 51.94127 3.153063e-07 0.02090980
#TRINITY_DN492_c0_g2_i2     10.540963 2.6369223 50.76775 4.942737e-07 0.02809558
#TRINITY_DN3927_c0_g1_i6     3.423790 5.9308253 44.69616 1.002688e-06 0.04479863
#TRINITY_DN5058_c0_g1_i1     9.507158 1.6073057 44.63367 1.013302e-06 0.04479863
#TRINITY_DN1524_c0_g1_i5    -9.697684 0.2788082 42.80529 1.385260e-06 0.05511878

# assembly is on graham here:
# /home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/trinity_assembly_all_batches_dmw
summary(decideTests(qlf))
# genotypezwt
#       genotypezdmwKO
# Down                2
# NotSig         397886
# Up                  7




FDR <- p.adjust(qlf$table$PValue, method="BH")
sum(FDR < 0.05)
# [1] 9
```

Relative to wt, 7 are upregulated in dmw ko and 2 are dowregulated.  I used grep to extract the seqs from the assembly and I blasted these against XL in xenbase and against everythign in NCBI. I put notes after each entry.
```
>TRINITY_DN15246_c0_g1_i10 len=2766 path=[0:0-452 3:453-1606 14:1607-1722 15:1723-1741 17:1742-1792 18:1793-1842 20:1843-1867 21:1868-1960 23:1961-1992 24:1993-2120 26:2121-2153 27:2154-2196 29:2197-2262 30:2263-2359 32:2360-2424 33:2425-2440 34:2441-2581 35:2582-2765] FAD-dependent oxidoreductase domain containing 1 Required for the assembly of the mitochondrial membrane respiratory chain NADH dehydrogenase (Complex I). Involved in mid-late stages of complex I assembly.
TATCAACGCAGAGTACGGGAGACTGATGACGTAATGTGAAAACGATGTATTAAACAGCGAGTAGTGGCGCACGATTCCTAGTTGTCCGCAGCTATGTATCTGGCTGGGATGTATGTCAGTTCATTGAAGTTCCCATTCCTTGGTATTGGAGTGTGGAAGGGGGTTCGGCTGTGGAGACAGAGGTCCTTGGGTACCAGTGCCTGCGCCCTGAAACAAGATGATTTCATAAAAGAACTGGACCAGAATTTTGTGCGGCTCCAGAAGAAGCTGATGGATAGTCTACCTTCCAGCGACTGGAGTCCATTTACACCCACTGGAGACCTACCACCCGAGCGAGCTGATGTAGTCATAGTTGGGGGTGGAGTGATGGGATGGTCTATTGCGTATTGGCTAAAACAGAAAGAAAACCGCAGGGGAGCTCTCAAAGTTGTGGTGGTTGAGAGGGATCCAACATATTCCCGTGCTTCTACAGTCCTGTCAGCTGGAGGCATACGTCAGCAGTTCTCCCGACCAGAAAATATCCAGATGTCACTCTTCTCAGCCCAGTTTCTGCGTAATATCAATGAACATCTTGGTGTTGTCAATGAGGATCGTATAGACATTCAGTTTAACCCATCAGGATACCTGTTCTTGGCCAGTGAGGAAGGTGCTACTATAATGGAAGAAAACTACAATGTGCAGAGAGAGTGTGGCGCCCAAGTGACCCTTATGTTGCCAGATCAGCTGAAGAAGAAGTTTCCATGGATCAATACAAATGGTGTGGCACTTGCTTCCTATGGTTTAGAGAACGAAGGTTGGTTTGACCCTTGGACACTGCTCAACGCTTTCCGGAGAAAAGCTCTCTCCATGGGGGTTTATCAATGCCATGGGGAGGTAACAGATTTTAGCACCGCCAAGCGTGAAATGATAACTGCAGATGGTGACCCAGTGACATTCAGTAGAATCGGACATGTTACTGTACAAATGCCAAATAGCCTGGAGAGCCAAAGTGTTGAGTGTTCCCTTGTGATCAATGCAGCGGGTGCCTGGTCCTCAAAGGTGGCTGAACTGGCTGGGATTGGTACTGGGCCGTCAAATTCATTAGAAGGGATAAAGCTTCCCGTGGAGCCCAAAAAGAGATATGTGTATGTCGTTCACTGCCCTAATGGGCCTGGATTGGACTGTCCATTGCTGATTGACAATTCAGGTGCCTATTTTAGAAGAGAAGGATTAGGAGGCAACTATATAGCAGGAAAATCTCCAGCAGAGGAGGAAGAGCCAGATATAAGTAACATGGAAGTAGATCATGATTTCTTCCAGGAGAAGGTATGGCCATTATTAGCACATCGGGTGCCAGCTTTTGAATCCTTAAAGGTTAAAACAGCGTGGGCTGGATATTATGACTACAATACCTATGATCAGAATGGTGTGGTGGGAATGCATCCACTAGTAAACAACTTATTCTTTGCCACTGGGTTCAGTGGCCATGGTTTGCAGCACTCACCGGCTGTGGGCCGTGCTGTGGCTGAGCTTATCGTGGATGGAGGCTTCAAGACGCTCAATCTCAGCTCTTTTTCCTTCAGGCGTTTCTGGTCACAGGAGCCTCTCTTGGAGCGGAACATAGTGTGAGCTACTGGAGAGTAGACCCGAAAAAACTGGTGAAAGAAGGCGCTAGCAGAGAGCCACAGCCAATTTCCTGCATTGAGAGACCAGTCTGCATCCAGGAGCAACTCCTCAAAGACCTAAAGCAAAGAAACAACACTTTGCTTTAATGAATTCTTCTATTTTCTGATAATAACCACATGGAGGGCGTTAGTGGTACCAAAATTGCCTTACTCTTTCTGGACTACAAAACAATTGACTTATCATTTATATTATTGCAGGCAGAGCAAATCACTGCCTGTATGTATATATTTTACATTGCATCTTGCATACACAGCACTGTACAATTTACAATTAATTTCAAAGTTTTAGGACAAAGAACATACAATAAATATACAAGATTACACCTCAAAAACAGAATAGGACAAGATGCACTGCTCTGATAAAGCACTTATGATTGATTTAGTACATCTGCTTTCTTAGGTTTATGATTTGCTGAAGTACCCTTACTTCCATCACCCACCAACAAGAATTTAAGAAGACTGGACAGCTGCAGACTGGATCAAAAAATCAAATCTGTCTTAAGATCGGTCTCATTGGAAGTCTTTACCAGATACAAATCAACCATTATTGTAGGAATGAGCAAAAAAAAAAAAATAACAGAAGAAAATATGCTCCATTTGAAATTATTAATCGAGTAGCTCTATCAAAAAATAAATCTAATGTAAGCTACCCATCACTGTTCAAAACTGGATTCTACCCGAAGAGAAGTTTTTGAGATATGATTTTTTTTTGTTAAAGGGATAGAAAGTGATTAGAATACTTTTGTTTTGATTTTGTTTTAATTTTCCATCAGACCTTTTGTCCTTCTTCCCATGATATCCAGAGGTCCCCTCGGGTGAGGAAACAAGCGACAGCATGTCGAGCTAAGTGGTGGATCCATCCTTCTGTTCTTAGCTGTGTCATAATTGCATCTATGAATGGGTAACCTGTTCTACCCTGAGAAGAAATGTTAAACAAAAGCAGGGTCAGATAAACTGCATAATATATTGTTCTTATGGCTGGGGTTTCTAACCTATGGGTTGAGACATTATTAAATTGGATCACAGTGACCAATTTACATTTTAATAACAGAAAACATATTGAACTAGCTACTTACAAAGTATTCCCTGCAGTTCCCCCCCC
>TRINITY_DN5299_c0_g1_i4 len=725 path=[0:0-415 2:416-425 3:426-483 4:484-485 5:486-487 7:488-582 8:583-614 10:615-648 12:649-661 14:662-673 17:674-681 18:682-685 19:686-688 22:689-724] apov1.1.S apovitellenin 1 gene 1 S homeolog
ACCCTGAGCCACCATTGATAAAACTTCTCTTATTACAAGGATGGATGGGAATGTTCACTTACTGGTGGTGATAGTGAAGAGGAGTCTTATATACGGTTGTGCCATTGTTACTTAAAGTTTCTCTGATCTCAAGCTGTAAACTGTTGCTCTGATGAGTGTGGCAGCTCCAACCAAATGAAGATAAATACACTTAACTGGACATCCCTCATCCTCTTTTTACTGATGATTAAGTAAAATGAAGTCAACCAACATGAAAAAATGTAACTTATTTCTCTTTTATCTTCTGCCCTTATGAAATGAACCAGAACAAGGGTCACATAATATGTTTTATTGGGCCTTTAATACCAGGTTCTGAACATTCAACTAATGACCATTACCCAGCCCTGATCTTCAAGGATCCCATATCTTGCTGATATGCTGATATATAGCTCTTCGACTCATGAAATGCGTATTCCACGGACAACAGTGACATTAACAGGAGGCGCCTGGCTTCTTGGAACCATTCCTTGTCACCCAAGTCCCTGAGCATCTCTCCAGCTTTGGGTGAAGCTGCATTCACTGCATTATACACATATGTCCCAATGTCATCGACTACAGGGAACCATCCGCGCATCGCGTGTCTCTTTGTCACCCCCTTTCTGTCACCCTCTGTTCTTGTGTTCTTCTGAGTTGGGAGCAATGGAGCCCCCGTACTCTGCGTTGATACCAAGCAGTGGTATCAACGC
>TRINITY_DN2367_c0_g1_i1 len=1851 path=[0:0-1850] krt15.2; keratin 15, type I gene 2
TTTTTTGCACATCCCATGAAGACCTTTATTTGAATTAAAGAGACAGAACATTCTCCATCAAGGGACAAGGTTTTGCAGATCACCATTGTGGTTGCTGTAGATATAAACCATTTTTTGTTCTTTATTTTAATAAAATCTCCTCTGCTTGATGGGCTTACACACAAAAACAGAGCTCTTGAGTCTTTATATTGCAACCTTGTTCCTCTATTTCAGTTAAAACCTCCTGCATGCCAATAGTCTGGGAATGGATAGACAACACAAAACCTTTCTTTATGTTTTTACATCAAGCACTTAGATATGAAAGAATGTTTAGTAATAGCTTTTCCAGTGTAATAACTATTGTGAACATTTGCACTACCTAATTATAGCCCAGACTAATTTAGAACACTGTTATTGTCTGCCCCTAAATCAGATGGTGGAGATGTCTCTGGTTGGCAGGAATATCTCAGAATGTCTTTAGAAGTGTGGAGGAGACAACTTTTCCACCAGATTCTTCTTTTATTATGGCATGTACACGGGTTGTTGTGCTTACACTGCTGGATCCTGATTTTGCATCTGCTCCAGCTGCAATTAGGCCTCCCTCTCCCTCCAGGAGCTGATGGTATTTGAGAATTTCCTGTTCAAGGCGACTCTTGATGTCTAATAAAGCCTTATACTCCTGGTTTTGTCGCTCAAGGTCTAAACGGAGGTCTGAAAACTGAGCCTCCATGCCACTAATGATGTTCTGAATCTGGAAGAGCTGTGCAGCATAACGACCTTCTGTATCTGCCAATGATGCCTCTAAACCTGCTTTAGTGCTTAGCAGAGTTTGAAGCTCCAACTCCAGACTCTGCAGTCCACGTCTCAATTCCGTGCTTTCACTCTTGGTGCTCTGTACTTGCTGTGTGCTGGTAACCACTTCACTCTGCAAACTTTCAACCTGACTGAGGAACCAAGCTTCAGCTTCCCTCCTGTTTTTATCAGCAATATGTTCATATTGTTCACGGAGATCATTCAGTGTGTTTAGAAGGCCAATTCCAGGAGCAGCATCCAACTCCACACTCACTGTTCCTGCTGCATGCTGCTTTTTCTCACTGACTTCCTCCTCATGATTTGCCTTCAAGTAAACCAGCTCTTCCTTCAGACTTTCAATCTGGAGCTCCAGGTCAGATCTGGTCATGGTCAAGTCATCCAAGATGTGTCTCAGGCCACTGATGTCATTTTCAACAGTCTGGCGCAGGGAAAACTCATTCTCATACTTGAGTTTGAAGTCATCAGCAGCCAATCTCGCATTGTCAATCTCCAGGACAAGTTTGGCATTGTTAATCGTGGCAGCTAAAATCTTACCACGTAGTTCATCAATAGTGGTGTATAAAGCTCCATAGCTTTGTTCACGGATAGCAACTGACCCTTGCTTAGCATACCAATCACGTATCTTGGCCTCCAGGTCAGTGTTGGCTTCCTCCAAAGACTTCACTTTAGCCAAGTAACTGGCAAGGCGGCCATTAAGGTTCTGCATGGTGTACTTCTCATTGCCAGAGAGCAGCCCTTCTCCACCACCTGAACCAAAGCTAGCACCAAACGAGGCCCCACCACCAAATCCAGAGCCAAAACCCGAGCCAAAGCCCATTCCTGACCCAGCTGCTGATGTGAACTTGGCACTGGACATGGAAATGCCAGATCCTCCATAACCACCATGCATACTAGCAGCTCTGTAGATTCCATCAGAGAGCCTGGAGGAGATCCCACTGCTCAGTCTAGATGAACTGCTCTGAGAGAAGCTATAACTCATTGTGATGGTGGGAGATGGTTAAAACAGTGACACGATTCCTACAGTGTTTTAGCCAGCAGTCCCGTACTCTGCGTTGATAC
>TRINITY_DN209_c5_g1_i3 len=2516 path=[1:0-1383 2:1384-1391 4:1392-1470 6:1471-1476 8:1477-1676 12:1677-1690 13:1691-1691 14:1692-2515] dnajc9.L; DnaJ heat shock protein family (Hsp40) member C9 
ACAGGTATCAACGCAGAGTACGGGAGATCGAAATCCGCCCAGTGTTGTCTGTCAGAGCGGAGGCTTGTGATTGGCGAGCAGGTGAAGAAGCTGATAAGTTCTCCTCATACAAGGGAATCAGGATGCCGGGCCTTCTGGAGAGCTGCGAGCGACACTTTGGTACGTCAGACCTTTACAAAGTGCTCGGGGTACGCAAGGAAGCGGGAGAAGGAGAGATACGTCGCGGGTATCACCGGGTGTCATTGAAAGTGCATCCCGACCGGGTTCAAGATGGAGAAAAAGAAACTGCCACTGCTGAGTTCCAGATACTTGGAAAGGTTTATGCAGTTCTCAGTGATAAGGAGCAGCGAGCTTTGTATGATGAGCAGGGAATTGTTGATGAAGAAACTGATACACTGAGTCAAGATCGCAACTGGGAAGAATATTGGAGGTTGCTGTTTAAAAAGATAACTGTGGAAGATATAAAGGCATATGAAGAAAAGTACAAAGGCTCAGAAGAAGAAAAGAATGATATTATATCCGCTTATATGGACTTCGAGGGAGACCTGGATGGCATCATGGAGTCTGTGCCATGTGCAGAGTTTGAAGATGAACCAAGAATTAGACAGATCATTCAGAAAGCCATTAAATCCAAAGAAATTCCTTCCTATGACACTTTTGTTAAGGAGACTAAAAAGAAACGTGAACAACGAAAGAAACGGGCTCATGAGGAAGCAAAAGAAGCTGAAGAAATGAAAAAGGAAATGGGCCTTGGTGATGATAATGATGACCTTAAGGCACTGATACAAAAAAGGCAGAGTGATCGAAAAAAAGAAATGGACGGTTTCTTTGATCAGTTGGAAGCAAAATACTGCAGTAATTCGAAAAAAGCACCAAGCAAATCCAAAGCACCAAAGCGAGGGAAGAAGTGATTTTTGTCCGATTGTGTTCCTCCGATTGTGTCGTCATTTTAAGGTCCATTAAACTTCATTAGAGAAGCCAAAAAGAAAAGTGAAAGTTTTACAAAAAGGGATTTCTATTTGTCCTTAGATTGTATCTTTATAGTAGTTACCCAAATGGGTTCATAATGTTGTGGATCTACTGAAGAGGGCACCTGTCTCGGAAGATTTAACAGCACTTACAAATGATGTAGATGTCCTTGGAGTTCAATTTTAATATTAGCTTGAACCAATTTAAAGGGTAAACAAAACACTTTATTTAGCATAGCATGTGTGGCAAAGCTAAGACTCCCTGCATGGCCAAAGCCATTTAGTTTATTGTACAGACATGATAAAAACAATAGGCAACAAAGATTTTGCTTTGTGCCGTCTAAGTTTGAGTTCAAAACCCTATTGAAAATGAATAACATTTTTGAACTCTGTGTACCTAGGGTATGAGTTTTTTTTTTTTTTTAAAGCAGGAAAAAACAGCTATGACCCATAAACTTATTGTAATGTTTCCTAAGCCTACAATTTTGCTTGTTGCTCATCCAATTCCAGCATTTCTCAGCTTTGTAGATCCAATAGGCACCTGTACTTCTTTTATAACATTCAGCTCCATGAGCTTTACTGGTCTCTCCTTAATATGAATCGGCAAAATGTATATTTGTAAATGTTGTGCTTTTTGTTGTAGATGCACTTACCTTTTAAATGCAAAATCCCAACCAATAAAATATTTTGTTAACTTAAAAAAAAAAAAAAAAAAAAAAAAAAACTCATACCCTAGGTACACAGAGTTCAAAAATGTTATTCATTTTCAATAGGGTTTTGAACTCAAACTTAGACGGCACAAAGCAAAATCTTTGTTGCCTATTGTTTTTATCATGTCTGTACAATAAACTAAATGGCTTTGGCCATGCAGGGAGTCTTAGCTTTGCCACACATGCTATGCTAAATAAAGTGTTTTGTTTACCCTTTAAATTGGTTCAAGCTAATATTAAAATTGAACTCCAAGGACATCTACATCATTTGTAAGTGCTGTTAAATCTTCCGAGACAGGTGCCCTCTTCAGTAGATCCACAACATTATGAACCCATTTGGGTAACTACTATAAAGATACAATCTAAGGACAAATAGAAATCCCTTTTTGTAAAACTTTCACTTTTCTTTTTGGCTTCTCTAATGAAGTTTAATGGACCTTAAAATGACGACACAATCGGAGGAACACAATCGGACAAAAATCACTTCTTCCCTCGCTTTGGTGCTTTGGATTTGCTTGGTGCTTTTTTCGAATTACTGCAGTATTTTGCTTCCAACTGATCAAAGAAACCGTCCATTTCTTTTTTTCGATCACTCTGCCTTTTTTGTATCAGTGCCTTAAGGTCATCATTATCATCACCAAGGCCCATTTCCTTTTTCATTTCTTCAGCTTCTTTTGCTTCCTCATGAGCCCGTTTCTTTCGTTGTTCACGTTTCTTTTTAGTCTCCTTAACAAAAGTGTCATAGGAAGGAATTTCTTTGGATTTAATGGCTTTCTGAATGATCTGTCTAATTCTTGGTTCATCTTCAAACTCTGCACATGGCACAGACTCCATGATG
>TRINITY_DN1063_c0_g1_i6 len=2714 path=[1:0-77 3:78-112 4:113-305 6:306-915 10:916-980 11:981-1013 13:1014-1076 14:1077-1105 16:1106-1754 17:1755-1806 19:1807-1848 20:1849-2607 21:2608-2646 23:2647-2673 24:2674-2713] ets2.S; v-ets avian erythroblastosis virus E26 oncogene homolog 2
CTATGAATCAAAATATGTTACAGAAATAGGTTTTCATAAAGGCACTTGTTAAAAAAATACAACTTCCTTTCAAAAGAAGAAAAAAAAAAACTGCAGTAAAAAAAGCTGCTACGCACACTTCTCCTGGCCTATGAAAATGAGGCAAATACATTGAAAAAAATATATATAATGAATTCCTATTATGAAATCAAAAGGCAATTAAAAACAAAAAGTACAACAGTGAACAAATTATGTACAGCAAAGGAGAACGGGACATTCTATCGGTTGTGCAATATTGTGAACAATCGGCATCAAATCACCCACATGGGAAACACAAGGCTTCATTTAAAACAAAAAAAAAATAGCGGTAATCAAGTCATGTACAAACAGTCTTTCCTAGCGCTCCTGTCTTAATGCTCTCGGCAAAGGAAACGTTACGATATTTACAACGAGAATCATATTTAAAATACAAATAATTAAAATATACATATACCTGGCAAAATAAAATGAGTTCTGTAAATATAAAACCTCATATATTATTTACAAGCGGTCTTGAAATGATGAGAATGCCATGATCACATTTTCAGTTCATTAAAAACAAAAATGAGGGTGTATATGGACAGCCTGGTGACTCGGTCCTGGTCAATATCTACATGACTTAGAGATGATGCCAATGAAGACTTTTGGAGAATTCTCAAAACTCCCAGATCCTTGGAATGTGGTTCAGCAACAGGCACCAGGCCTTGATCCGATATCCAAGCACGACAGTTTTAAACAAACTAAAACGCCTCTCTTCTCGTTCATGGACTGAAACAAGAAGTTCCTTAAGGCGGGACGCATAAAAGTAGCATACACGACGAGTTACTTATGTGCAATCTTCTTTTTCAAGTGAGCAACTGTATAACAATAAGCGCCATCCATCAGATCTTTACTCCGACCCAAGAGAATGGAGTAAATCACTAATTCATGTCACGGATTTGGTGTTTTTCCAGTCTTTTCAACTCATGCACTGTTCGGCCATCTTGGACAGCCCACCAAATGATCCTTTTCGCAGTACAATAATGTGCTCTTCAAGCCAACGATTTAGTTCGCACTCCAGCAGCTCATTCGTCTGTGTCGGGCTGTACGCCAAGCATGGCGTGCAATTCATCCGGTGTGTAACCCAACAGGTTGTGAAGGTCGCAAACGAATCGGTACACGTATCTCTTCCCCGACGTCTTGTGGATAATGTTCTTGTCGTAATAATACCTTAGTCCTCGACTAAGCTTCTCGTAGTTCATTTTTGGTTTGTTTTTCCTTTTGCCCCACCGACGGGCCACCTCATCTGGGTCAGTCAGCTTGAACTCCCATCCGTCCCCAGTCCAGCTGATAAACGACTGACATGACTTATCCGTTAAAAGTTCCAACAGGAATTGCCACAGCTGGATTGGTCCACTTCCTGTGAATCCAGCAAGTATAGATGCCGGAATAACAGGCTTTCCTAGTTCTGCTGGCTCACATCTGTCTTGGATGTAATCTTTAAAGGACATAGGTGGCTTATTTAGGCATAACGCCTGGCTGCCATCTTCCTCAAAGCCATCGTAAGAAGGAACTCGTTGCATGTCGACAAGAGAAGACTGGCTCGTCCAAGACTGTAGCAGAGATTCGGTGCTCTCGAAGCTTTCCGTGCCACTGTCTCCCGAGTCATAATCTCTCAACTTGCCAGAATTTAGAGAATTGAGCAGTGCATTCATATTGCTCCTGGCAAAGTCCTGGCTGGCGGGTGAGTAGTTTACCGCTCTGGATTTAAGGCAAGAGCTGGGGTATTGTTGGAACTCTTGTTTGGGGTGGAGGAGAGTTTGACCAGTAGGCACTGAGCACATTTCATTAAACATTCCATTTTTAGGGTAATTGTGTACTTGTGCTCCACACTGCAAAGGATCAGTGGTGAAATTTAAGGAATCTGCATTCATCCAATGGTTTATGGAGTCCTGGTTAGAATGGTCAATATACGGCTCCTGTGCCTTTTCTTGATGTTCTTTCATCATTTCTTCCAGATGCTCCCAAAGGATGTCCCCAACAAAGTCAGGAGCCAGAGCCAAAAACCTTTCCTTCCCGAGGCTGCACAGCTCGTGTCCATTCATGAGAAACTTCTGAAAATTCACGTTCTGCAAGGAGAATTCTTTGGCGGCCCACCAAAGCCACTGAAACACATTATTCTCATCCCAAAGCCACGGATTGCTGAGGATGCCAAGTCTGAATCGTTCCTTTGCAAAACCATTAAACGTATTTTTCAGCGCTTGACTCATCACAGCTTTGCTGCACGGGGTAAGCAGGGGTAACTCACAGCTGCTGGAATCATGAGAGTAAGAGTCTAGGCCCGTTGGTACCGCCTGTTCTTCTTCATAAGCAGAAAACAAGCCAGAGTAGAGGGAGGTCGGCACGCTGACATTGTCGAACGCTAACTGGCGCTTGAGCATTGCTCTGTGGCCGTTATACACCGGAGCCACTTGGTCCATGTTCCGAATTCCAAACTCTGTCATTGGCCCTGAGCAGTTGGGGAGAAGCCCTTCCCGGGTATATCTATTTATATCCACAAATAATCCGAAGGGAGTCGTTAGCGACACTGAGATTTATTAGTATCCCGCGATTCCCAGGCTGCTGAGCTCCCGTGAGCGTTATCTCTATATTCGGCCGGAGCTCCGATACGAAGCGCTTTCAGGGCTGACTGTGCTCCCGTACTCTGCGTTGATACC
>TRINITY_DN3375_c0_g1_i6 len=2423 path=[2:0-455 3:456-524 5:525-569 8:570-610 9:611-621 11:622-636 13:637-735 15:736-768 16:769-876 18:877-919 19:920-957 21:958-1012 22:1013-1047 24:1048-1075 25:1076-1104 27:1105-1138 28:1139-1164 30:1165-1195 31:1196-1287 33:1288-1320 34:1321-1434 36:1435-1484 37:1485-1815 39:1816-1834 40:1835-1857 42:1858-1869 43:1870-1949 45:1950-2422] cops3.L;COP9 signalosome subunit 3
ATATAATGTACAACATTTCTAGCTATTTCTTTAGTTAGACTTTAGTTCTCCTTTAAGTACTAGAAGGGGTAAAGGCTATTTATTACACATGTAGACTGCTGCCCAATGATTTAAAAGGGAAACATGTTGTTTATGGGAGTTGTTTCTGTGCTGGGGGAGAACCAATGACCTTGTGCCTCAGTGACACGCCCCTCCCAACAAAGGGATCTACCAATAGACCAACCATTATCCATGTAACACAAAACTGATTCACTGTAAACAAAATAAAAGCTACTGGATCCTCAAATGAATTTATCATTCATAGGAAAAGTCAAAGAGACCACTCCACACAGCGTAACACCAAGATCTTCATATACACGAAGGATGAAAGTGAGAGAAATTCACAACAGGTTTTTTTTTTTTTTTATTGTTTCTTGAAGATTGTTTTTCCTTGACGTGTATCTGTGTAACGCATAGAGAGCAGAGAGAAAACGTTATAAATATTTGCTTTCGTTCCCATGTAAATGTCCAGGCAGAAACGATCATCTCCCTCTCTGCTCTTTCTTATTTCCGAAGACGTATTTTTATTGTGGGGATTCTGAAGGTTTTTGCATTGCTGATTATCCACCTTTGTTGGGAATACATAGCGGAATCTTGTGTTTCAAGAGTAGCTGGATGGTTTGCTTCCTGAATCGTCCTCCTGAGACCCCATGCTCTTTTGAACAAACTGAGGGTTCACAGTGATCTCTTGATCCATCGCTTTCAACCTATCGTCAAGGTCTATACATCGCAGCATCTCTTGATCAATATTGTGAAGCATCGCAGGATTGTTATATTTCTCTGGGTTATCATGAAAACACACCATGCCGTCCTTTTGATTAATACTTGCAAATATTTCCCCATCTTCTATCATGTGCAGCACATATTTCTCTGCCTCCTGGGCACCAGATAACTGAACGCGACTGGCCATGTCTTGTAATGACAGCGTTAAAAATGTTTTGGTGAGTCTCTGAATGTTCTTCTTGTAGAGGGAAGACAGACATTGCTTCACCAGCCCCATGTTGTTGTCTCTGGTAAAGGTTTCGTTGTGTTTGCTCACAAGGTTCCGCAGCTCAGCTGGGTTATTAGTTGAATAGACTTGTGCCAGCTCATGGTAAGCGTTGCTAAGAGGCTTGATAAATCTACCAACAATCTGAGAGGTATATTTAGGCAGCTGTTGAACCTTTCCATGGAGTATCAAGGAGACCAGGATATACTTTTTATAGGCTTCCAGCATTATGTGACTAACAGCCATGGCTGGAGTTGTGATAGCCTGCTCGTAAAAGTACAGCGCTCTCTCGAGGTTCTTCAAGCCTGTGTAAATCATCCCTCCGTAATAGTAGTAACACAGGAAGGGCTTGGCATCATAAGCTCCATTCTCTTTGCAAATGTCCATCATATCTACATCTAGGTACGCCAGAGCGGGCTTGAAACATCTTGCTAATAGACAGAGCTGGCAGAGATCTGCGTGGATTGAGGTCAGCTGGTTTGCATTCATTTGCATCTTGTCTATGGCTTGTCTGATAACACAGATTCCACGCAAGGGCTGCTTCCTTTCCACAAGTGCATTTGTTAACTGATGGCAGAGACCAGCAAAAGTATCTGTCGCGTATCTAATGTGCTCTCCGTTGCAAGTGCTGATAAAGAGTTGAACTTGAGAAAATAAAGTCTCAAAATCGGGAATACTGGGCATGGAAAACTTCACAAACAAAACAGCTAAGACCCCGAGAGAATGCTCCTGAACATCCAGCGCCCCCAGAACAGTGTCCAGATGAGAGAGGTTCTTTGCTAGTAGCTCCCCACTCTTGTTAATTAGTTCACATAGCTGTGTCATTTGCCCTTGGGAAGACAGCTGTCTCACGCTGTTCACGAACTGCTCCAGTGCGGAAGCCATGCCGCCTGGATCCCGCCTCCCTCCGCGGGGAACTGCAGCCGGGAGCGGAGCGGGAGGAACAGAGTGGAGCACAGGCCAGGCTCACTCACACAGCGCTTAATCTCGCGAGAGCCCACACCGGCCCGCGCGCTGTCACTGGTGTAACATCTCGCGAGAATCTCCGGATCACTAGCTATCGCGAGATTGCTGTTGGCGCTGAAGGCTTGGACCTGGAGTGTTCTCGCGAGAACTTGCCGGGGGACAGCTGGAATGTGTGGGCGTAACTTGTTTCTGCTGAGTCGTCGAAGGCGGATGCGTAGATCTTAAAGGTGCCGCTTACTCGAAAAAGGGGGAGCCAAAGGGATGATGTCTACGTAATAAGGCGTTTGGTTTATAGGAAACCTGCAGACACCATGTAATATTGACCAATGTACACCCTTAACAACCAATTAGATACTTGTACACCCGTAGCAACCAATCAGATACTTGTACACCCGTAGCAACCAATCAGATACTTGTACACCCGTAGCAA
>TRINITY_DN492_c0_g2_i2 len=2345 path=[0:0-647 2:648-676 3:677-755 5:756-1512 7:1513-1698 8:1699-1729 10:1730-2344] unclear hit: eukaryotic translation initiation factor 4E binding protein 2 (eif4ebp2)
AATTTTCTGTCCAGCCCAATCGACTAACTGATATCGGTCGGCTCATTTCCCACCATACACACACCGAATATTGTACGAAAATTCATTTTATACAACATTATCTGTGCGTTTATAGCCACCTTTGGCCATAGACCATGCTGGTGGCCACAGACGTACAAAGTTGCATATCACCAAAGGGTCGCAAACATGGCAGCTCCCTATGAAGAGACATGATTTTAGCTGCATTTAGCCCAGGGCACTCGCATAAATCATAATGCCAACAAGTGTTTGCTGTCTTTTGTTTGCCCTTTATTCACAACGGTCTTTGATAACAAGCACATTTGTGCAGCTATGTTTCCTTGTTACACAGCCTTATTATTCTTGAGAATTAGCCATGCCACATTATAAAAACATATGCCAAGTATTTGATTATCTAAATATTACAATGCTGTCAATATTACAAAATGCCTTGAGGAAAGAACACTTTATTTCCCTTTTGGTTGCATTTCCCAATATTATGATGATTAATTGCTGACACCTACAGGTATTTCTTCTGCAACTGCAGAAAATCCCACATCAGAAAACGATACTTGATATCTGCACAAAGATGTTTATTATTGCATTGCAATGGAAATGTCCTATTTCAATTGGTTTGCCAATTAATAATAAAAAAAAATAGTTCTTCTTTCAGACCAAATTAAGCCTAGCGCAAAACCATATGCATTCAGAAAGAATTTAGATTTTTGGCTGCTAGGCAAATCTGAAAATATTCAAAGTAAAAAAAAAAAAATGAACCAATGCCCACATATACAGGTCATTCAAATATGATTGGTGACAGGAGAGGGCATAGGATAAGTTCTCCATGCAGTGTATTAAGACTCAGGGCTTGTAAAGGCAAACATCTGGAGATATATTACATCTCCAAAATGAAGGGACAGCCAAACAACAAAAGTATGAGTGGACAATGGTTGTCTGCTTTTAGCCTAAAGGCACAGATAAATATGCTTTTTCCATATTTCAATCTTCTCCATACCCAGTGAGCATAGTCTACATGCCCTAATTCTTTTATAGAGTAATGACACTAGCAGGCAAGAAATTCTAGTTTCAACTGGGATGGATTCCTTGAGCAAAGTAATACTCTCATTAAAGGAAAAGTACATTTGTGCAGAGGTAGAAGCAAGGCCAGATGAAGTTTAACAGAATTATGACAAAAAAAAGGCAGCTTTTTAATTTTAATGTAACATCTAAAAATCAAACTGCTGATGACTGAGTGGTAGCCAACATGAAAGGAATGGTATAGGACATGCTGGGAGAGCAGCCCTTCCCAGACTAAAGCAGTAAATACTCCCAATCATGCCCACACACATTTTCTATGTGCCGAGGAACATACTCCTAGTAGTAGATTATAAAAGTATGCCCTTAAATTGCTTCCACAGTGCCTAGTCACATATTCTCAAAATTCCTTTAAAAAAAAGTTATTTCTTAAAGCACGCAAACTATATTTACTAGCCAATTATAGCATGATTGTACATTTCTGCTTCTTTAAAAGGACGTAAGGGGTAAAGTGGAGGCCTGATTTGTATTTACCCCCCATTCATAGGTAAAGGATCCATGTATTTTGAGTCCCTTCATCAAAATATTACTGTTCAAAAGAAGAAGTTGAAGAAAGTGAATTGATAAGGGCAGATATTAGAAAGCACTAGGTGCCCTTAGTACACTGTCTGGTGCACTTGTACATTAGGGTTTATATAAAAAAAAAAATATCCAAAAGCAAAAATCAAGGTTGGTGACCTAAATCGGGAGCCCACAGTTTGCTACACAATTTCCTCAAACTTAATTAAGAAAAAAAGCAAAACGTAAAGGACTTTATGAAAGTTCAAGTAAAAGTAGATGATCCGCTTTTCCTTTTCAGCAATGAAATGAGAGCGAGAAGTTAGTTCTTATAACAAACAGCGCGCTTCCATGCTGCAAACTTGATGAGCACCATATGATCTTTTCGGTCCACAAAACCAACTAAATTTGATCAAGAATGGCTATTAATATTTGCTCTAAACATTTCGTTAGTTGTACGTTTTAGACGACTTGACTCAAAGTGGTTTGACGTCAGAGCATTTTCTTCAGAAATGTTAAATCTGTATATAGACCAGACCAAATTAATTGTGCAGTGCCATTCAGCTGGTACAGCCACAATCAAATGGGTCCGTCACGCAACCGTGGGATATGACAACGTCTAGATGCAGCTTATTTAAGGATTAGTAACCTGTGAAATATATGGCCCAGATTTTAAGTCATCAGGAACTAACTTAAAAGTGGCCCATATTTAAAGCTGAATATCACTGGGGAAACATTTGCCTGCATTAGAGCTC
>TRINITY_DN3927_c0_g1_i6 len=1810 path=[0:0-31 2:32-94 4:95-134 6:135-179 7:180-263 9:264-367 10:368-400 11:401-443 13:444-748 15:749-815 16:816-859 18:860-886 19:887-1183 20:1184-1258 22:1259-1580 24:1581-1685 25:1686-1714 27:1715-1809] vim.L vimentin
TGTATAAGAGACAGGCTAAAGATTTATTGAAACACTGGCTGCACATTGCCTGGAAATAAAAAAAGCTTTTTTTTTCTTGGGAATGTAAAAACAGTTTTGCGGCACAAGATCATTGCTATTACAGTCTGTAAACTAACAGAAGGAAGCTCGCTGTCTGTCCTTCTCAAGGGCTGTCAATCTAACTCCTAACTGTCAGCTAGCTTATCCCTATCTTGCGCTCTCCAATTGACTGCAAAAGGCACTTGAAAGCTGGTTTTCTTCAAAGTTTAAAACACACCAATATTATGGTGCTGGCAAAGGTTCTCTCTTATGTTCCAGATATTGCAATTTCACTCAAAGTCATCGTGGTGCTGAGAGCTTTCATTAACAACCTGTCCATCTCTTGTCTCCACAGTCTTGATTAGCACAGTGCGCTTGGAATGAGTTTCCGCTGGGTGAGAATCAAGGTTGGTTTCTCTCAGGCTCATTGTAGAAAAGGAATGTACAGGAAGGGAAATTCTGCTCTCCTCTCCCTCCAGGAGTTTCCTGTAAGTAGCAATCTCAATATCAAGAGCCATCTTAACATTGAGCAGATCCTGGTACTCTCGCAAGTGCCGAGCCATCTCCTCCTTCATGTTCTGGATCTCCTCTTGCAGGCGTTGAATAGTGTCCTGATAATTAGCAGCTTCAATGGCAAAGTTCTCCTCCATTTCTCGCATTTGGCGTTCGTAAGACTCGTTAGATCCTTTCATTGCATCAATCTCGCAGGTGAGAGTCTGGATCTGTCTGCGGAAGTCACTGGTCTCCTGCTTAGCTTGACGCAGGGCATCATTGTTGCGATTAGCAGCTTCCGACAGGTCAGCAAACTTAGACTTGTACCATTCTTCGGCATCTGAAAGATTTTTAGCAGCCACATTCTCATACTGCTGGCGAACATCACGTAGGGCAGCAGTGAGATCTGGTTTCGACACATCCATGTCAACCTGGATATGTGACTCCTGAATTTGAAGCTGAAGTTCTCGGATTTCCTCATCATGGAGTTTCTTTAGGAAAGCAATCTCCTCTTGTAAGGACTCAACCTTCCTCTCCAAGTCAATGCGTGCCAGAGAGGCATTATCCACATCCTGTCTGAATGACTGCAGGTTGCCTTCTGCTTCTTCTTTTTGGATCATTTCATCTTGAAGCTTTTCTCTCAGTCGCTGGAGATCATCGGCCAGATTGTCCCTGTCTACTTCTACTCGGGCTTTGTCGTTGGTCGCTTGGTCCAGTTGTCTGCGGAGCTCCCTCATCTCCTCCTCATAGAGATCCCCTATCCGTGATGTGCCTTTCCCTTTCAGCTGCTCCAGCTCCGCCACCAGGATCTTGTTCTGCTGCTCCAGGAATCGCACCTTGTCGATGAAGTTGGCGAACCTGTCGTTGAGCTCGATCATCTCCGCCTTCTCGTTGGTCCTGTTGGCTTTAAACTCGAGGTTGACGGCATCTGCCAGGGCGAAGTCTACGGAGTCTGCCATTCTGGCCGGGGGGAGGCTGCTCCTCAGTCTGACCGAGGACGACTTGAACACCGCGGGGGAGGAAGATGTGGAATAGACCATTCTGCTGCTGGTGCTCGGGCGCATGGCACTGCCCAGGGTGTACCTAGTGCTGGAGGTCGCATAGCGATTTCCAGAGGAGCTCGACCTAGGGTTTCCCCCGAAAATCCTTCTGTAGGATGACTTGGTTGTTGCCATGTTGGTGCTTCCTTTGTAATTCCTTGAGTTTATCTGCGAGCGGCTTCCCTTGACCTCACCCAGCCAGTACTTTACCAGGACTCCCCGTACTCTGCGTTGATACC
>TRINITY_DN5058_c0_g1_i1 len=1492 path=[2:0-1356 3:1357-1374 4:1375-1375 6:1376-1380 8:1381-1384 11:1385-1491] rhog.S ras homolog family member G S homeolog
TTTTTTTAAGCATATAAAAGATTTTATTGGCCGCCCGAGGGCACAGTTCTGCGGGTATAAAACATTTCTCACAGGCAGTTACATCACAAAAATCTGTTACTTTCACAACACACAATAAAAAAGCCAGAGCAGTGGAGCTGCCATACAGCAGGAGCGAGCACACGGGTCCATCTCGCTGCTGTACCCCAGCAGTTACAGGGTTAATGGGACACTTGATCATTAGTAGGTGAAAGGCAATTACCACACGGATCCCAACCACACAGACTTGTTTATTTGGAACCTCGGCCAATAAGATTCAACCATATTCAGAGCTATCTGCTCCCTGGGGGCGGCCCCTGCACTGCGACCACATCTCCTCGGGAATCATGCACAGCACCAAACTGGCCAATGGAGTGCAGAGACTTGGCCCCGGCATGAATATGTTTGACCCCCACTAGGAATTCTGGGATCAAACATGGAGCTGCACAAGCTGAACTCAACTCCAACAAATCACCTGACAGCTGGAGGGCCCCTCAGGAACAAACAGGAAATCTGCTGCTTATGTGTGACAATAAGAGCTGGGCGTTTCCTCTGATTGGCCGGTCATTGGCAAGAAGAAAATTTAGAAATTTCCTATAAATAATCTTCTGACTAATCAGCCACGTTTGCATATTTCCACCCTGGCTCATAACGAGGTGCCAGGCAAGTGGAAAACACTGTTCCTCCAAGATCCTCACGTTCATTGCTCAATGGAGAGCTCTTGTGATCACCCACAAAAGGAACTGACACCCTATTTGTTAACCCCTCCCCCTACAATTCCTCGTTCTAATGTGCTGCTACAGTATATGTTTTCTTATTTCTTCCGGATGTTTGTCAGTGTTATTGGCTCCTCCTGGGCGTGGCCAAGATGCGATGATACAAGTGAGTGTTAAAATAGCAGTGTAGTTGGGATGAAAAACCTGTAAATAAGAATCAGTGACTTGCTAAGATATCATATACATTGTAACAGTTATTAGAGGACCTGAGTATAAGGACAGCCCCTGAGGAACCCCCCCGAGTGACGTAATCCTGAGTGGTTGCTCTAACTGGAAGCTTTAACAATAAAACCATTAGATATCCCAGGACCCAACTTCCAGGGGAAACTTACAGCCAGAATTAGGCATTTTTAAATAAAATGTAATATTCGGGAAATTGGCCCCTGGGGGTTACATCGATTTTTAGAAATAAAAGTGGCAATATACACCTGTGCTCAACATGTGGGAAATAGCATTTTTGTTTGTTATGCATTATTTTGTAGTCGTTAGAGCTAAATAGTAATTAAAGGGATCCTGTCATCAAAAAACATGTTTTTTTCAAAACACATCAGTTAATAGTGCTACTCCAGCAGAATTCTGCACTGAAATCCACCTCTCAAAAGAGCAAACAGACCTCTATATATTCAGTTTTGAAATCTGCCTGCTCCCTTGAGAAATGGACCCCAGTGAAAAACTCCGCTGGAGCAGCACTACTAACC
>TRINITY_DN1524_c0_g1_i5 len=2354 path=[0:0-465 2:466-466 3:467-1139 5:1140-1192 6:1193-1265 7:1266-1297 9:1298-1361 11:1362-1402 12:1403-1586 14:1587-1684 16:1685-1690 17:1691-2353] glt8d1.L; glycosyltransferase 8 domain containing 1
CATTATACCTTATAATACATGAGTGATACAGATTTTTATCTCAGAACAATGTCCAACTCAGAAACATTATCATTTAATAGATTTAATATATAAAACACACAATAAGAAAACTAGGGGGCGCTGCTCCATGTTGCAAGTAACAGATTTCCTCGGGGAGAGACAAATTACCAATTGTGATTGGCCCAGTCTCTAAACATAAATTTATTTAGGCAGCAACACTATGGCAAGAAGGAGGCGGGGTTCCTGCGGTTATGACTAGTGGGAGGAGCAATAGTTGTGTCTCAGTCCCGCCTCTAGTGGTGACCAAGTCAGCTAGACATTCCCTGCAATGGCAGGAGGGCCCCAACCCTGTATTGTGATTGGCTCTACTGTGAATTCTGACTTTGAGGAGTAAGTTGGGGGAAACGTGAGCAGTTTATAGGGGTGCGAGGCCTTGGCCCCCACAGGAACAATGAGACCCCCCCCCCAATGTCACAATAGGAACAAGTGATCGTCCCCCGTCTGTCTCTTCCCTGAGTTCCACCGACGTCCGGGAGTCGCAGACTCATCACTCTCCGTGTCTCCTGATTAGTGAGAACTGCCCCGAGGGGTCTGGCAGGAACCACTTCTCCCAAATCTCAGGGAACGAAGATGTTCTGCCCCAAGGTTTAAAATGTCCGTTCCAATGGAGAAGTTTAGCGGCTTTCACAAACTGCGGCGAGTAACGTTTTCCTGTACTTGAACCCAGGTGTCGGACGTGCCACAGAGGGTTAATGTTGGAGTAAAGTCTATAGAAGACGATGAGGAGGGGGGGTGCAGCGATGTTGCCAGATAACGACTTGCTGTAAAGTTCCTCAGTCACATCCAGCTCCATCCACTTCTCTAACTGGCGGGTGATGTTCTGCCGCCGCCACTCGGTCAGGTTGGCCACAAAGACCCCCGGGTTGAAGGAGCAGGTGTTGGGCTTTATTCCCAGGCTCCGGACCCGCTCCTTCTTGTAGTCCAGGAAACCAACGTAATTGTACTGATTAGCGCCCCCTCTGACGGGGAATTTGGAGGTGACGGAGTCACAATCCTCAGAGAAGGCAGCCGCATGGCCGGGGCTGATTGGGGTGTTGTACAGTTGCACAATGTCATCTTGTACAATAACATCATCATCCAAGTAAATGACTTTTTTGGCTCCAGGCAGCAAACTGGGCAGGTAGAACCGGGCAAATGTCATTGGTTTCACTGGCTCCGCCCCGGCATCAACCCGCACTTTCCCATCAAGAACACGGGCATCAAATGGCAGCAGCTTATAGGCCACTCGTTTCAGGTCTGTGCCATCCAGCCAAGAACTGATGTGTTTCTTGGTGTCATTGGTAGTGATGATATAGAAGACAATGTTAGACTTTGTGTTGCTGGAGATACTGTTAATTGTGGCGATGAGGCCTCCAAGTCGTTCCTCCACCCCTGGAATCACCACAGCAATCTCCTCCCCGTGTCCTTGTTCTGGGCCAATGTCAGGGGCATCAGGAAGGGATTCGAGTTGGTGGAACACTAAGGGGCCAGTGTCTGAGTTTTGTCTCTTCAGGATGTCGCTGAGGCCCAGGATATTGTGGTGCAGAATCAGCAGGAAGATTACGGCAGACAATAGGATCACGGCTACATGAACTTTGCGCAACGTCATCTTCTCATCTGCTCCTCAACTCAACGACACAAAGATCCTGTAATCTTCTCTCTGCCCCCGAAACCCGCAGCTACTCCGCCGCCAACTCAGCTCTAGCCGCGCGCGTCACGTGACCGGGAGACTGACCGTCGCCGCCATTTTGAATTCGGGCAGAAATTCCCACTTGTCTCTTGTTTCAGCGACACCCGCGGGGGAGTCCCCAGTGAGAGCATTTGTGTCCCTCCCCTTCTGCTCTACGTCACACGGATATTCTGGCTGCTGCTCCGCCCCCGCTATTACCTCCCCGGCAAGACCCACAATTTATTTGTTTGACTACAAAGTGCGACTTCCGGTCACGCATGCGCGGACACTGCAGGGAAAGAACCTGAGTCTGTAACTTCCGCTCATCCTGAGACAAGAGTCTGACACGCTGCGCCTCTCCCTGGTGTGTCACTTTGTACAGGCACAGCCCGGGCAGCCTCTGGCTTCTCTTCCGCCCCACGTGATTCTGTTTTATTCACACGTCACATTAATCCCAGGGCCCGGCGAGTAAGAGCCCAGAATCCATCAGGCACAGACCGGACCGCTATTGGGTCCTTCTCGCTGTCAGCTGACGTCGCTGTGAGTCGGCGGCACTTCCGGTGGAGCAGCTGTTACTGGCGCGATCATGTTGGGAATATTCAGCTCCATCCCGACCCAGATGGACTTCAAGGGTCAGAAACTCGCC
```

# Checking specific transcripts
I'd like to figure out which transcript(s) correspond to dmw, dmrt1S, and dmrt1L.  I can use blast to query the assembly with these seqs to get the transcript ID, and then I can use R to print out the counts in each individual.

```
module load nixpkgs/16.09  gcc/7.3.0 blast+/2.9.0
makeblastdb -in dmw_trinity_assembly_all_batches.Trinity.fasta -dbtype nucl -out dmw_trinity_assembly_all_batches.Trinity.fasta_blastable
blastn -query dmw_mRNA_NM_001114842.1_ex4_only.fasta -db dmw_trinity_assembly_all_batches.Trinity.fasta_blastable -out dmw_ex4only_to_dmw_assemb.out
```