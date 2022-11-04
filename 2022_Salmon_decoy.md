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
# Make decoy list (following directions from here: https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/)
```
grep "^>" ../../2021_XL_v10_refgenome/XENLA_10.1_genome.fa | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
```


Concat transcriptome and genome for making decoy file:
```
cat XENLA_10.1_GCF_XBmodels.transcripts.fa.gz ../../2021_XL_v10_refgenome/XENLA_10.1_genome.fa.gz > XL_v10_transcriptome_and_genome.fa.gz
```
