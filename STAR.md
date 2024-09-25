# Path
```
/home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq
```
for Jade the directory would be something like this:


```
cd /home/ben/projects/rrg-ben/froglady/
mkdir XL_v10_refgenome
cd XL_v10_refgenome
gunzip XENLA_10.1_genome.fa.gz
```

# Download the XL genome seq from xenbase and the gtf file
```
wget https://download.xenbase.org/xenbase/Genomics/JGI/Xenla10.1/XENLA_10.1_genome.fa.gz
wget https://download.xenbase.org/xenbase/Genomics/JGI/Xenla10.1/latest/XENLA_10.1_Xenbase.gtf.gz
```

Now make a directory for your scripts
```
cd ..
mkdir jade_scripts
cd jade_scripts
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
--genomeDir /home/ben/projects/rrg-ben/ben/2021_XL_v10_refgenome/ \
--genomeFastaFiles /home/ben/projects/rrg-ben/ben/2021_XL_v10_refgenome/XENLA_10.1_genome.fa \
--sjdbGTFfile /home/ben/projects/rrg-ben/ben/2021_XL_v10_refgenome/XENLA_10.1_GCF_XBmodels.gtf \
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

STAR --genomeDir /home/ben/projects/rrg-ben/ben/2021_XL_v10_refgenome/ \
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
featureCounts -T 4 -s 0 -p --countReadPairs -C \
  -a /home/ben/projects/rrg-ben/ben/2021_XL_v10_refgenome/XENLA_10.1_GCF_XBmodels.gtf \
  -o ${2} \
  ${1}

```


