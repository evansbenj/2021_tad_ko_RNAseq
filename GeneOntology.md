# Gene ontology analysis

The main problem with doing an analysis of differential expression in Xenopus is that the names of most transcripts is silly and does not permit analysis if function.  For example, what does "XBXL10_1g43753" do?  

To begin to deal with this, I am going to blast each differentially expressed Xenopus transcript against a human transcriptome and identify the top hit.  I will use the acronym of these genes to perform a gene ontology analysis in humanns using this tool: http://geneontology.org/

# DE genes
I am writing all the DE genes to csv files.

Cat them and skip first line
```
awk 'FNR>1' MF_STAR_ccdc_DE_DeSeq2.csv MF_STAR_dmrt1L_DE_DeSeq2.csv MF_STAR_dmrt1S_DE_DeSeq2.csv wtko_STAR_ccdc_DE_DeSeq2.csv wtko_STAR_dmw_DE_DeSeq2.csv wtko_STAR_scan_DE_DeSeq2.csv > All_STAR_DE_DeSeq2.csv
```
get names only
```
cut -f1 -d "," All_STAR_DE_edgeR.csv | cut -f2 -d "\"" > All_STAR_DE_edgeR.names
```
or
```
cut -f1 -d "," All_Kallisto_EdgeR.csv | cut -f2 -d "\"" | cut -f2 -d '|' > All_Kallisto_EdgeR.names
```

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

First add a bar `|` to the end of each line of a text file containing unique identifers to prevent extra matches (must be done on graham because does not work locally):
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
for i in `cat ./ccdc_kallisto_edgeR.file `; do grep -i $i ../XL_v10_transcriptome/XENLA_10.1_GCF_XBmodels.transcripts.fa >> ccdc_kallisto_edgeR_sequence_file;done
```
```
for i in `cat ./All_STAR_DE_DeSeq2.names `; do grep -i $i ../XL_v10_transcriptome/XENLA_10.1_GCF_XBmodels.transcripts.fa >> All_STAR_DE_DeSeq2.full_names;done
```
Now remove the greater than sign:
```
sed -i 's/>//g' ccdc_kallisto_edgeR_sequence_file
```
Now use seqtk to extract the fasta entries
```
module load StdEnv/2020 seqtk/1.3
seqtk subseq ../XL_v10_transcriptome/XENLA_10.1_GCF_XBmodels.transcripts.fa ccdc_kallisto_edgeR_sequence_file > ccdc_kallisto_edgeR_sequence_file_output.fasta
```

# Get best alignment of blast results (based on bit score)
I am using dc-megablast, which is for somewhat similar sequences (here I am comparing human and frog protein coding).  I should also consider using `-task blastn` which allows for even more divergent matches.  The default is `-task megablast` which is for highly similar seqs.
```
module load nixpkgs/16.09 gcc/7.3.0 blast+/2.10.1 
blastn -task dc-megablast -query ccdc_kallisto_edgeR_sequence_file_output.fasta -db ../human_transcriptome/gencode.v42.transcripts.fa_blastable -outfmt 6 | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > ccdc_kallisto_edgeR_de_to_human_best_single_hits.blastn
```
# extract acronyms of successful queries
```
cat dmw_kallisto_deseq2_de_to_human_best_single_hits.blastn| cut -f1 | cut -f3 -d\|
```
# Extract acronyms from best blast hit
```
cat ccdc_kallisto_edgeR_de_to_human_best_single_hits.blastn | cut -f2 | cut -f6 -d\|
```

# Check if how many orthologs of the queries have orthologs
```
wc -l ccdc_kallisto_edgeR.file
cat ccdc_kallisto_edgeR_de_to_human_best_single_hits.blastn | cut -f2 | cut -f6 -d\| | wc -l
```

use this list for the GO analysis here:  http://geneontology.org/

# Colate results:
```perl
#!/usr/bin/perl
use warnings;
use strict;
use List::MoreUtils qw(uniq);

# This program reads in three files and colates the results
# the first one is a list of DE genes for a particular analsis
# the order of the first list will be the same as the order of the output
# the second file has the fasta headers for each of these genes;
# we need to harvest the xenbase gene acronym from this one
# the third one has blastn output; we need to harvest human 
# gene acronyms from this one when available

# Colate_for_gene_ontology.pl All_STAR_DE_DeSeq2.full_names All_STAR_DE_DeSeq2_to_human_best_single_hits.blastn All_STAR_DE_DeSeq2.out

my $input1 = $ARGV[0];
my $input2 = $ARGV[1];
my $outputfile = $ARGV[2];

unless (open DATAINPUT1, $input1) {
	print "Can not find the input file.\n";
	exit;
}
unless (open DATAINPUT2, $input2) {
	print "Can not find the input file.\n";
	exit;
}
unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile\n";
	exit;
}

my @temp;
my @temp1;
my %gene_hash;
my $counter=0;
# first load up the gene IDs and xenbase acronyms
while ( my $line = <DATAINPUT1>) {
	@temp = split(" ",$line);
	@temp1 = split(/\|/,$temp[0]);
	$gene_hash{$counter}[0] = $temp1[1];
	$gene_hash{$counter}[1] = $temp[1];
	$gene_hash{$counter}[2] = "";
	$counter+=1;
}	
close DATAINPUT1;
# now get the acronym of the human ortholog if present
my $y;
my @temp2;
while ( my $line = <DATAINPUT2>) {
	@temp = split(" ",$line);
	@temp1 = split(/\|/,$temp[0]);	
	@temp2 = split(/\|/,$temp[1]);	
 	for ($y = 0 ; $y < $counter ; $y++ ) {
		if($temp1[1] eq $gene_hash{$y}[0]){
	 		$gene_hash{$y}[2] = $temp2[5];
	 	}
 	}
}	

for ($y = 0 ; $y < $counter ; $y++ ) {
	print  $gene_hash{$y}[0],"\t",$gene_hash{$y}[1],"\t",$gene_hash{$y}[2],"\n";
	print OUTFILE $gene_hash{$y}[0],"\t",$gene_hash{$y}[1],"\t",$gene_hash{$y}[2],"\n";
}		
close DATAINPUT2;
close OUTFILE;
```
