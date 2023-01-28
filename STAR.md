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
