#!/bin/bash

#PBS -N mapping.statistic
#PBS -l nodes=1:ppn=1,walltime=1:00:00
#PBS -V
#PBS -m ae

# * generate mapping statistic
path=~

echo -e "CellID\tInputReads\tInputReadLength\tUniquelyMappedReads\tUniquelyMappedRate\tMappedLength\tMultipleMappedReads\tMultipleMappedRate" > ${path}/statistics/2_pass.mapping.statistics.txt

ls ${path}/bam/2_pass/final/ | while read id 
do
     lib=`basename $id`
     if [[ -f ${path}/bam/2_pass/final/${lib}/${lib}_Log.final.out ]]; then
          InputReads=`cat ${path}/bam/2_pass/final/${lib}/${lib}_Log.final.out | grep "Number of input reads" | cut -f 2`
          InputReadLength=`cat ${path}/bam/2_pass/final/${lib}/${lib}_Log.final.out | grep "Average input read length" | cut -f 2`
          UniquelyMappedReads=`cat ${path}/bam/2_pass/final/${lib}/${lib}_Log.final.out | grep "Uniquely mapped reads number" | cut -f 2`
          UniquelyMappedRate=`cat ${path}/bam/2_pass/final/${lib}/${lib}_Log.final.out | grep "Uniquely mapped reads %" | cut -f 2`
          MappedLength=`cat ${path}/bam/2_pass/final/${lib}/${lib}_Log.final.out | grep "Average mapped length" | cut -f 2`
          MultipleMappedReads=`cat ${path}/bam/2_pass/final/${lib}/${lib}_Log.final.out | grep "Number of reads mapped to multiple loci" | cut -f 2`
          MultipleMappedRate=`cat ${path}/bam/2_pass/final/${lib}/${lib}_Log.final.out | grep "% of reads mapped to multiple loci" | cut -f 2`
          echo -e "${lib}\t${InputReads}\t${InputReadLength}\t${UniquelyMappedReads}\t${UniquelyMappedRate}\t${MappedLength}\t${MultipleMappedReads}\t${MultipleMappedRate}" >> ${path}/statistics/2_pass.mapping.statistics.txt
     else
          echo -e "${lib}\tError!\terror!" >> ${path}/statistics/2_pass.mapping.statistics.txt
     fi
done

