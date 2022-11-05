# Path:
```
/home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/XL_v10_transcriptome
```

# Creating a decoy file

There are some dependencies for the pipeline to do this.  One is mashmap and this is here:
```
/home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/bin/MashMap/mashmap-Linux64-v2.0/mashmap
```

Another is bedtools:
```
module load StdEnv/2020 bedtools/2.30.0
```

Salmon:
```
module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 salmon/1.7.0
```
# Make decoy list 
I'm following directions from here: https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
```
grep "^>" ../../2021_XL_v10_refgenome/XENLA_10.1_genome.fa | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
```


Concat transcriptome and genome for making decoy file:
```
cat XENLA_10.1_GCF_XBmodels.transcripts.fa.gz ../../2021_XL_v10_refgenome/XENLA_10.1_genome.fa.gz > XL_v10_transcriptome_and_genome.fa.gz
```

# make the decoy genome
In this directory:
```
/home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/XL_v10_transcriptome
```
I ran this script first to generate the decoyed ref transcriptome:
```
#!/bin/sh
#SBATCH --job-name=salmon
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=64gb
#SBATCH --output=salmon.%J.out
#SBATCH --error=salmon.%J.err
#SBATCH --account=def-ben

module load StdEnv/2020 bedtools/2.30.0
module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 salmon/1.7.0
../bin/SalmonTools/scripts/generateDecoyTranscriptome.sh -m /home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAse
q/bin/MashMap/mashmap-Linux64-v2.0/mashmap -a /home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/XL_v10_trans
criptome/XENLA_10.1_GCF_XBmodels.gff3 -g /home/ben/projects/rrg-ben/ben/2021_XL_v10_refgenome/XENLA_10.1_genome
.fa -t /home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/XL_v10_transcriptome/XENLA_10.1_GCF_XBmodels.transc
ripts.fa -o /home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/XL_v10_transcriptome/salmonindex_withannotatio
n
```

# Index the decoyed transcriptome
```
#!/bin/sh
#SBATCH --job-name=salmon
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=64gb
#SBATCH --output=salmon.%J.out
#SBATCH --error=salmon.%J.err
#SBATCH --account=def-ben

module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 salmon/1.7.0
salmon index -t /home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/XL_v10_transcriptome/XL_v10_transcriptome_
and_genome.fa.gz -i /home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/XL_v10_transcriptome/salmonindex_witha
nnotation -d /home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/XL_v10_transcriptome/decoys.txt -p 2 -k 31
```

# Quantify

```
#!/bin/sh
#SBATCH --job-name=salmon
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=64gb
#SBATCH --output=salmon.%J.out
#SBATCH --error=salmon.%J.err
#SBATCH --account=def-ben

module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 salmon/1.7.0
salmon index -t /home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/XL_v10_transcriptome/XL_v10_transcriptome_
and_genome.fa.gz -i /home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/XL_v10_transcriptome/salmonindex_witha
nnotation -d /home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/XL_v10_transcriptome/decoys.txt -p 2 -k 31
```
