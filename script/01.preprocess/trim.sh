#!/bin/bash

#PBS -N trim
#PBS -l nodes=1:ppn=4,walltime=1:00:00
#PBS -V
#PBS -m ae

#usage: ls data/rwadata/RNA_SubLib_*/*L1_1.fq.gz | while read id ; do file=`basename $id`; \
#	 path=`echo ${id%/*}`; sample_id=`echo ${path##*/}`; file_id=`echo ${file%_*}`; \ 
#    qsub -v sample_id=$sample_id,file_id=$file_id -o ${file_id}.trimming.log -e ${file_id}.trimming.err trim.sh; done

# * export
export PATH=~/miniconda3/bin/:$PATH
source activate trim

# * get infor.
path=~
sample_id=$sample_id 
file_id=$file_id 

# * trim
if [ ! -d ${path}/trim/${sample_id}/${file_id} ]; then mkdir -p ${path}/trim/${sample_id}/${file_id}; fi
trimmomatic PE -threads 4 ${path}/data/rawdata/${sample_id}/${file_id}_1.fq.gz ${path}/data/rawdata/${sample_id}/${file_id}_2.fq.gz -baseout ${path}/trim/${sample_id}/${file_id}/${file_id}.trimmed.fq.gz \
	ILLUMINACLIP:~/adapters/TruSeq3-PE-2.fa:2:30:10:2:True SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:40 HEADCROP:2 CROP:55 -trimlog ${path}/trim/${sample_id}/${file_id}/${file_id}.trimming.log

