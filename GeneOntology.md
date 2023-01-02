# Gene ontology analysis

The main problem with doing an analysis of differential expression in Xenopus is that the names of most transcripts is silly and does not permit analysis if function.  For example, what does "XBXL10_1g43753" do?  

To begin to deal with this, I am going to blast each differentially expressed Xenopus transcript against a human transcriptome and identify the top hit.  I will use the acronym of these genes to perform a gene ontology analysis in humanns using this tool: http://geneontology.org/

# Human transcriptome
In this directory:
```
/home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/human_transcriptome
```

I downloaded a human transcriptome here: https://www.gencodegenes.org/human/

This transcriptome has 252416 transcripts of which 115846 are protein coding.

# Make a blast db out of the human transcriptome
```
makeblastdb -in gencode.v42.transcripts.fa -dbtype nucl -out gencode.v42.transcripts.fa_blastable
```

# Query seqs
Get the exact header for each file. First at a bar `|` to the end of each seq to prevent extra matches:
`emacs -nw id.file`

```
XBXL10_1g43753|
XBXL10_1g26238|
XBXL10_1g569|
```
Now use this to grep the headers:
```
for i in `cat ./id.file `; do grep -i $i XENLA_10.1_GCF_XBmodels.transcripts.fa >> ID_in_sequence_file;done
```
Now remove the greater than sign:
```
sed -i 's/>//g' ID_in_sequence_file
```
Now use seqtk to extract the fasta entries
```
seqtk subseq XENLA_10.1_GCF_XBmodels.transcripts.fa ID_in_sequence_file > output.fasta
```

# Get best alignment of blast results (based on bit score)
I am using dc-megablast, which is for somewhat similar sequences (here I am comparing human and frog protein coding).  I should also consider using `-task blastn` which allows for even more divergent matches.  The default is `-task megablast` which is for highly similar seqs.
```
blastn -task dc-megablast -query scanw_kallisto_edgeR_de.fasta -db ../human_transcriptome/gencode.v42.transcripts.fa_blastable -outfmt 6 | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > scanw_kallisto_edgeR_de_to_human_best_single_hits.blastn
```
# Extract acronyms from best blast hit
```
cat scanw_kallisto_edgeR_de_to_human_best_single_hits.blastn | cut -f2 | cut -f6 -d\|
```
