#! /bin/bash

script_dir=$1
GTF=$2
ref_gtf=$3
ref_fa=$4
cupcake=$5

module load python
source /data/$USER/conda/etc/profile.d/conda.sh
conda activate SQANTI3.env
export PYTHONPATH=$PYTHONPATH:${cupcake}

python ${script_dir}/sqanti3_qc.py ${GTF} ${ref_gtf} ${ref_fa}
