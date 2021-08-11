# Assembly with trinity (de novo)

I may not need this but was having issues on cedar so before I executed the script below I activated a python environment with numpy:
```
module load nixpkgs/16.09 python/3.7.4
virtualenv /home/$USER/my_venv
source /home/$USER/my_venv/bin/activate
pip install numpy
module load StdEnv/2020 scipy-stack/2021a
avail_wheels numpy --version "1.15*"
pip install numpy --no-index
```

here's the trinity assemnbly script:
```
#!/bin/sh
#SBATCH --job-name=trinity
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --exclusive
#SBATCH --time=240:00:00
#SBATCH --mem=0
#SBATCH --output=trinity.%J.out
#SBATCH --error=trinity.%J.err
#SBATCH --account=def-ben


module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3 trinity/2.12.0

Trinity --seqType fq --samples_file ccdc_fastq_list.txt --CPU "${SLURM_CPUS_PER_TASK}" \
   --full_cleanup --max_memory 110G --min_kmer_cov 2 --bflyCalculateCPU \
   --include_supertranscripts --output ccdc_trinity_assembly_all_batches
```
