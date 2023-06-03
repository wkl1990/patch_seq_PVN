#!/bin/bash

#PBS -N fastqc
#PBS -l nodes=1:ppn=4,walltime=1:00:00
#PBS -V
#PBS -m ae

#usage: ls data/rawdata | while read id ; do file=`basename $id`; \
#    qsub -v sample_id=$id -o ${id}.fastqc.log -e ${id}.fastqc.err fastqc.sh; done

# * export
export PATH=~/miniconda3/bin/:$PATH
source activate fastqc

# * get infor.
path=~
sample_id=$sample_id 

# * qc
if [ ! -d ${path}/qc/${sample_id} ]; then mkdir -p ${path}/qc/${sample_id}; fi
ls ${path}/data/rawdata/${sample_id}/*fq.gz | xargs fastqc -t 4 -o ${path}/qc/${sample_id}

