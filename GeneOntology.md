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
Get the exact header for each file. First at a bar to the end of each seq to prevent extra matches:
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





According to this link: https://edwards.sdsu.edu/research/perl-one-liner-to-extract-sequences-by-their-identifer-from-a-fasta-file/ Get fasta entries like this:

```
perl -ne 'if(/^>(\S+)/){$c=grep{/^$1$/}qw(42732422)}print if $c' AO245_newtrim_scaffolds.fa
```
or with a file:
```
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' ids.file fasta.file
```
awk -v seq="42732422 10004 289064 42567259-,...,42597465-" -v RS='>' '$1 == seq {print RS $0}' AO245_newtrim_scaffolds.fa
