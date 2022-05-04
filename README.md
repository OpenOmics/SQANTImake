# SQANTImake
Pipeline to generate custom SQANTI GTF &amp; modify it for use in RNA-seek

```
git clone https://github.com/tomh1lll/SQANTImake.git
cd SQANTImake

module load singularity
singularity pull docker://quay.io/biocontainers/agat:0.8.0--pl5262hdfd78af_0
```
Dry run of the pipeline
To perform a dry run a step of the pipeline (changing the directory to the directory containing the snakemake file & dependencies):

```
sh pipeline_submit.sh npr /data/NCBR/projects/NIAMS-11/SQANTImake
```

This should also install any missing dependencies within the SQANTImake directory

Actual run of the pipeline
Once everything seems to work, to perform a full run of a step of the pipeline, submit:

```
sbatch --partition=norm --gres=lscratch:500 --time=10-00:00:00 --mail-type=BEGIN,END,FAIL pipeline_submit.sh process /data/NCBR/projects/NIAMS-11/SQANTImake
```
