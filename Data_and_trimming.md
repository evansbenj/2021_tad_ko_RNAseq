# RNAseq

Protocol for RNAseq: Clontech/Takara SMARTer v4 cDNA conversion kit followed by Illumina Nextera XT library prep. This protocol does not generate direction poly(A) mRNA-seq data and hence you cannot infer the direction of the transcript from the data

# Path

```
/home/ben/projects/rrg-ben/ben/2021_XL_ko_tad_RNAseq/
```

# trimming

Apparently this kit uses the same seq on the F and R read, so the only seq to trim is:
```
CTGTCTCTTATACACATCT
```

# trimmomatic
```
#!/bin/sh
#SBATCH --job-name=trimmomatic
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=36:00:00
#SBATCH --mem=8gb
#SBATCH --output=trimmomatic.%J.out
#SBATCH --error=trimmomatic.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this
# sbatch ./2020_trimmomatic.sh ../raw_data/dmrt1


module load StdEnv/2020
module load trimmomatic/0.39

#v=1
#  Always use for-loop, prefix glob, check if exists file.
for file in $1/*R1_001.fastq.gz ; do         # Use ./* ... NEVER bare *
  if [ -e "$file" ] ; then   # Check whether file exists.
  	#if [[ $v -eq 1 ]]
	#then # if/then branch
	java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE ${file::-15}R1_001.fastq.gz ${file::-15}R2_001.fastq.gz 
${file::-15}R1_trim_001.fastq.gz ${file::-15}R1_trim_single_001.fastq.gz ${file::-15}R2_trim_001.fastq.gz ${file::-15
}R2_trim_single_001.fastq.gz ILLUMINACLIP:Nextera_XT_adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 M
INLEN:36
	#	  v=0
	#else # else branch
  	#	v=1
	#fi
  fi
done 
```
