#! /bin/bash

module load python

fastq=$1
sam=$2
prefix=$3
cupcake=$4

source /data/$USER/conda/etc/profile.d/conda.sh
conda activate SQANTI3.env
export PYTHONPATH=$PYTHONPATH:${cupcake}
collapse_isoforms_by_sam.py --input ${fastq}  --fq -s ${sam} --dun-merge-5-shorter -o ${prefix} -c 0.9 -i 0.95
