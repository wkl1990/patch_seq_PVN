#!/bin/bash

#PBS -N mapping2
#PBS -l nodes=1:ppn=1,walltime=2:00:00
#PBS -V
#PBS -m ae

#usage: ls manifest/*txt | while read id ; do file=`basename $id`; \
#    qsub -v manifest=${file} -o ${file}.mapping2.log -e ${file}.mapping2.err mapping_2pass.sh; done

# * export
export PATH=~/miniconda3/bin/:$PATH
source activate mapping

# * get infor.
path=~
manifest=$manifest
cell=${manifest%.*} 
starindex=~/data/annotation/star_v23

# * align
## 1. include multiple files in one run
##    by using multiple fastq files in one time
##    - we use --readFilesManifest param for more samples, and add RG tag using --outSAMattributes
## 2.  in two pass Mode
## 3. PCR duplicates were masked and removed using STAR option bamRemoveDuplicates

cat ${path}/bam/1_pass/polyA/RNA_SubLib_*/RNA_SubLib_*_SJ.out.tab | awk '($5 > 0 && $7 > 2 && $6==0)' | cut -f1-6 | sort | uniq > ${path}/bam/2_pass/SJ_out_filtered.tab

# smartseq and 0.5 score mapping, using vM23
STAR --genomeDir $starindex \
     --readFilesManifest ${path}/manifest/${manifest} \
     --runThreadN 1 \
     --sjdbFileChrStartEnd ${path}/bam/2_pass/SJ_out_filtered.tab \
     --outFileNamePrefix ${path}/bam/2_pass/final/${cell}/${cell}_ \
     --outSAMtype BAM SortedByCoordinate \
     --limitBAMsortRAM 10000000000   \
     --genomeLoad NoSharedMemory \
     --outSAMattributes NH HI AS nM RG \
     --readFilesCommand zcat \
     --outSAMmultNmax -1 \
     --quantMode GeneCounts \
     --outWigType bedGraph \
     --outFilterScoreMinOverLread 0.5 \
     --outFilterMatchNminOverLread 0.5 \
     --soloType SmartSeq \
     --soloUMIdedup Exact NoDedup \
     --soloStrand Unstranded \
     --soloFeatures GeneFull \
     --clip3pAdapterSeq AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA --clip3pAdapterMMp 0.1 0.1

