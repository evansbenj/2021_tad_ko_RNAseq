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
