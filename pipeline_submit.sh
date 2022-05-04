#! /bin/bash
set -e

###################
#
# Launching shell script for SQANTI pipeline job submit
#
###################
module load python/3.7
module load snakemake/5.13.0

cd $SLURM_SUBMIT_DIR

R=$2
echo $R

mkdir -p $R/snakejobs
mkdir -p $R/reports

if [ -d "${R}/gffcompare" ] 
then
    echo "gffcompare already installed." 
else
    echo "Installing gffcompare."
    git clone https://github.com/gpertea/gffcompare.git
    cd ${R}/gffcompare
    make release
    cd ..
fi

if [ -d "${R}/gffread" ]
then
    echo "gffread already installed."
else
    echo "Installing gffread."
    git clone https://github.com/gpertea/gffread.git
    cd ${R}/gffread
    make release
    cd ..
fi

if [ -d "${R}/cDNA_Cupcake" ]
then
    echo "cDNA_Cupcake already installed."
else
    echo "Installing cDNA_Cupcake."
    git clone https://github.com/Magdoll/cDNA_Cupcake.git
fi

if [ -d "${R}/SQANTI3" ]
then
    echo "SQANTI3 already installed."
else
    echo "Installing SQANTI3."
    git clone https://github.com/ConesaLab/SQANTI3.git
fi

##
## Test commandline arguments
##
if [ $# -ne 2 ]; then
    echo " "
    echo "Requires a single commandline argument: npr or process"
    echo "Also requires path to working directory"
    echo " "
    exit
fi

if [ $1 != "npr" ] && [ $1 != "process" ] ; then
    echo " "
    echo "Invalid commandline option: $1"
    echo "Valid commandline options include: npr or process"
    echo " "
    exit
fi

SCRIPT=`realpath $0`
SCRIPTPATH=`dirname $SCRIPT`

echo $SCRIPT
echo $SCRIPTPATH

##
## Run snakemake
##
echo "Run snakemake"

CLUSTER_OPTS="sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name={params.rname} -e snakejobs/slurm-%j_{params.rname}.out -o snakejobs/slurm-%j_{params.rname}.out"

if [ $1 == "npr" ]
then
    snakemake -npr --snakefile $R/Snakefile --configfile $R/config.yaml -j 1
fi

if [ $1 == "process" ]
then
    snakemake --latency-wait 120  -s $R/Snakefile -d $R --printshellcmds  --configfile $R/config.yaml --cluster-config $R/cluster.json --keep-going --restart-times 1 --cluster "$CLUSTER_OPTS" -j 500 --rerun-incomplete --stats $R/reports/snakemake.stats | tee -a $R/reports/snakemake.log
fi
