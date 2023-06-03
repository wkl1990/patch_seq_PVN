#!/bin/bash

#PBS -N dedup
#PBS -l nodes=1:ppn=1,walltime=2:00:00
#PBS -V
#PBS -m ae

#usage: ls bam/2_pass/final | while read id ; do file=`basename $id`; \
#    qsub -v cell=${file} -o ${file}.duplicate.log -e ${file}.duplicate.err star_duplicate.sh; done

# * export
export PATH=~/miniconda3/bin/:$PATH
source activate dedup

# * get infor.
path=~
cell=${cell} 

# * bamRemoveDuplicates
STAR --runThreadN 1 \
     --limitBAMsortRAM 10000000000 \
     --runMode inputAlignmentsFromBAM \
     --bamRemoveDuplicatesType UniqueIdenticalNotMulti \
     --inputBAMfile ${path}/bam/2_pass/final/${cell}/${cell}_Aligned.sortedByCoord.out.bam \
     --outFileNamePrefix ${path}/bam/2_pass/final/${cell}/${cell}_Aligned.markDup. 

samtools view -b -F 0x400 ${path}/bam/2_pass/final/${cell}/${cell}_Aligned.markDup.Processed.out.bam > ${path}/bam/2_pass/final/${cell}/${cell}_Aligned.markDup.removeDupl.bam

