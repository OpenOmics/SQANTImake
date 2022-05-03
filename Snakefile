#####################################################################################################
# SQANTI analysis workflow
#
# This Snakefile is designed to perform the generation of a custom annotation 
# with SQANTI to a custom genome (though hg38 human genome by default)
# uses minimap to align, followed by SQANTI, then GFFcompare and other tools
#
# Created: May 3, 2022
# Contact details: Tom Hill (tom.hill@nih.gov)
#
# Last Modified: May 3, 2022
#
#####################################################################################################

from os.path import join
from snakemake.io import expand, glob_wildcards
from snakemake.utils import R
from os import listdir
import pandas as pd

##
## Locations of working directories and reference genomes for analysis
##
sample_file= config["samples"]
rawdata_dir= config["rawdir"]
working_dir= config["workdir"]
ref_fa= config["ref_fa"]
ref_gtf= config["ref_gtf"]
sqanti= config["sqanti"]
gffread= config["gffread"]
gffcomp= config["gffcomp"]
cupcake= config["cupcake"]


df = pd.read_csv(sample_file, header=0, sep='\t')
SAMPLES=list(set(df['sample'].tolist()))

rule All:
    input:
      expand(join(working_dir, "pacbio/{samples}.sort.bam"), samples=SAMPLES),
      expand(join(working_dir, "SQANTI/{samples}.collapsed.gff"), samples=SAMPLES),
      expand(join(working_dir, "SQANTI/{samples}.collapsed_corrected.gtf.cds.gff"), samples=SAMPLES),
      expand(join(working_dir, "SQANTI/{samples}.collapsed_corrected.fasta"),samples=SAMPLES),
      expand(join(working_dir, "SQANTI/{samples}.collapsed_classification.filtered_lite.gtf"),samples=SAMPLES),
      join(working_dir,"SQANTI/total.combined.gtf"),
      join(working_dir,"SQANTI/total_filter.combined.gtf"),

rule minimap:
  input:
    FQ=join(rawdata_dir, "{samples}.ccs.fastq"),
  output:
    SAM=temp(join(working_dir, "pacbio/{samples}.ccs.sam")),
  params:
    rname="minimap",
    dir=directory(join(working_dir, "pacbio")),
    genome=ref_fa
  shell:
    """
    mkdir -p {params.dir}
    module load minimap2
    minimap2 -ax splice -t 8 -uf --secondary=no -C5 {params.genome} {input.FQ} > {output.SAM}
    """

rule sort:
  input:
    SAM=join(working_dir, "pacbio/{samples}.ccs.sam"),
  output:
    SAM=temp(join(working_dir, "pacbio/{samples}.sort.sam")),
    BAM=join(working_dir, "pacbio/{samples}.sort.bam"),
  params:
    rname="sort",
  shell:
    """
    module load samtools
    sort -k 3,3 -k 4,4n {input.SAM} > {output.SAM}
    samtools view -hb -@ 4 {output.SAM} > {output.BAM}
    """

rule collapse_isoforms:
  input:
    FQ=join(rawdata_dir, "{samples}.ccs.fastq"),
    SAM=join(working_dir, "pacbio/{samples}.sort.sam"),
  output:
    GTF=join(working_dir, "SQANTI/{samples}.collapsed.gff"),
  params:
    rname="collapse_isoforms",
    prefix=join(working_dir,"SQANTI/{samples}"),
    dir=directory(join(working_dir, "SQANTI")),
    cupcake=cupcake,
  shell:
    """
    module load python
    mkdir {params.dir}
    source /data/$USER/conda/etc/profile.d/conda.sh
    conda activate SQANTI3.env
    export PYTHONPATH=$PYTHONPATH:{params.cupcake}
    collapse_isoforms_by_sam.py --input {input.FQ}  --fq -s {input.SAM} --dun-merge-5-shorter -o {params.prefix} -c 0.9 -i 0.95
    """

rule SQANTI:
  input:
    GTF=join(working_dir, "SQANTI/{samples}.collapsed.gff"),
  output:
    GTF=join(working_dir, "SQANTI/{samples}.collapsed_corrected.gtf"),
    GFF=join(working_dir, "SQANTI/{samples}.collapsed_corrected.gtf.cds.gff"),
    FA=join(working_dir, "SQANTI/{samples}.collapsed_corrected.fasta"),
    TXT=join(working_dir, "SQANTI/{samples}.collapsed_classification.txt"),
  params:
    rname="SQANTI",
    script_dir=directory(join(working_dir, "SQANTI3")),
    ref_gtf=ref_gtf,
    ref_fa=ref_fa,
  shell:
    """
    module load python
    source /data/$USER/conda/etc/profile.d/conda.sh
    conda activate SQANTI3.env
    python {params.script_dir}/sqanti3_qc.py {input.GTF} {params.ref_gtf} {params.ref_fa}
    """

rule SQANTIfilter:
  input:
    SAM=join(working_dir, "pacbio/{samples}.sort.sam"),
    FA=join(working_dir, "SQANTI/{samples}.collapsed_corrected.fasta"),
    GFF=join(working_dir, "SQANTI/{samples}.collapsed_corrected.gtf.cds.gff"),
    TXT=join(working_dir, "SQANTI/{samples}.collapsed_classification.txt"),
  output:
    GFF=join(working_dir, "SQANTI/{samples}.collapsed_classification.filtered_lite.gtf"),
  params:
    rname="SQANTIfilter",
    script_dir=directory(join(working_dir, "SQANTI3")),
    dir=directory(join(working_dir, "SQANTI")),
  shell:
    """
    module load python
    source /data/$USER/conda/etc/profile.d/conda.sh
    conda activate SQANTI3.env
    cd {params.dir}
    python {params.script_dir}/sqanti3_RulesFilter.py --sam {input.SAM} {input.TXT} {input.FA} {input.GFF}
    """

rule mergeGTF:
  input:
    GFF=expand(join(working_dir, "SQANTI/{samples}.collapsed_corrected.gtf.cds.gff"),samples=SAMPLES),
    GTF=expand(join(working_dir, "SQANTI/{samples}.collapsed_classification.filtered_lite.gtf"),samples=SAMPLES),
  output:
    GTF=join(working_dir,"SQANTI/total.combined.gtf"),
    GTFF=join(working_dir,"SQANTI/total_filter.combined.gtf"),
  params:
    rname="mergeGTF",
    gffcompare_dir=gffcomp,
    ref_gtf=ref_gtf,
    prefix_filter=join(working_dir,"SQANTI/total_filter"),
  shell:
    """
    {params.gffcompare_dir}/gffcompare -r {params.ref_gtf} -o {params.prefix} {input.GFF}
    {params.gffcompare_dir}/gffcompare -r {params.ref_gtf} -o {params.prefix_filter} {input.GTF}
    """

rule AGAT:
  input:
    GTF=join(working_dir,"SQANTI/total.combined.gtf"),
    GTFF=join(working_dir,"SQANTI/total_filter.combined.gtf"),
  output:
    GTF=join(working_dir,"SQANTI/total.AGAT.gtf"),
    GTFF=join(working_dir,"SQANTI/total_filter.AGAT.gtf"),
  params:
    rname="AGAT",
    dir=working_dir,
  shell:
    """
    module load singularity python
    singularity exec -B {params.dir} agat_0.8.0--pl5262hdfd78af_0.sif agat_convert_sp_gff2gtf.pl --gff {input.GTF} -o {output.GTF}
    singularity exec -B {params.dir} agat_0.8.0--pl5262hdfd78af_0.sif agat_convert_sp_gff2gtf.pl --gff {input.GTFF} -o {output.GTFF}
    """

rule cleanUp:
  input:
    GTF=join(working_dir,"SQANTI/total.AGAT.gtf"),
    GTFF=join(working_dir,"SQANTI/total_filter.AGAT.gtf"),
  output:
    GTF=join(working_dir,"SQANTI/total.AGAT.clean.gtf"),
    GTFF=join(working_dir,"SQANTI/total_filter.AGAT.clean.gtf"),
  params:
    rname="cleanUp",
    dir=working_dir,
  shell:
    """
    module load python
    {params.dir}/clean_gtf.py {input.GTF} > {output.GTF}
    {params.dir}/clean_gtf.py {input.GTFF} > {output.GTFF}
    """
