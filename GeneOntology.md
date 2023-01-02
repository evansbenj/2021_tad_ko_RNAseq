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
Make a header for each file. 

First add a bar `|` to the end of each line of a text file containing unique identifers to prevent extra matches:
`sed -i 's/$/\|/' ccdc_kallisto_edgeR.file`

```
XBmRNA21528|
XBmRNA44681|
XBmRNA48999|
XBmRNA4722|
XBmRNA50613|
XBmRNA51371|
XBmRNA54222|
XBmRNA54223|
XBmRNA56258|
XBmRNA66942|
XBmRNA68936|
XBmRNA68976|
XBmRNA69753|
XBmRNA74882|
XBmRNA77377|
XBmRNA82371|
XBmRNA82372|
XBmRNA82995|
XBmRNA9690|
XBmRNA596|
```
Now use this to grep the headers:
```
for i in `cat ./ccdc_kallisto_edgeR.file `; do grep -i $i XENLA_10.1_GCF_XBmodels.transcripts.fa >> ccdc_kallisto_edgeR_sequence_file;done
```
Now remove the greater than sign:
```
sed -i 's/>//g' ccdc_kallisto_edgeR_sequence_file
```
Now use seqtk to extract the fasta entries
```
module load StdEnv/2020 seqtk/1.3
seqtk subseq XENLA_10.1_GCF_XBmodels.transcripts.fa ccdc_kallisto_edgeR_sequence_file > ccdc_kallisto_edgeR_sequence_file_output.fasta
```

# Get best alignment of blast results (based on bit score)
I am using dc-megablast, which is for somewhat similar sequences (here I am comparing human and frog protein coding).  I should also consider using `-task blastn` which allows for even more divergent matches.  The default is `-task megablast` which is for highly similar seqs.
```
blastn -task dc-megablast -query ccdc_kallisto_edgeR_sequence_file_output.fasta -db ../human_transcriptome/gencode.v42.transcripts.fa_blastable -outfmt 6 | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > ccdc_kallisto_edgeR_de_to_human_best_single_hits.blastn
```
# Extract acronyms from best blast hit
```
cat ccdc_kallisto_edgeR_de_to_human_best_single_hits.blastn | cut -f2 | cut -f6 -d\|
```
use this list for the GO analysis here:  http://geneontology.org/
