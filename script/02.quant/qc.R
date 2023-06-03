#!/usr/bin/R

suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("viridis"))
suppressPackageStartupMessages(library("reshape2"))


map_stat <- read.table("statistics/2_pass.mapping.statistics.txt", header=TRUE)
map_stat %>% mutate(depth=MappedLength*UniquelyMappedReads/1000/1000, UniquelyMappedPercent=as.numeric(str_remove(UniquelyMappedRate, "%")), 
    MultipleMappedPercent=as.numeric(str_remove(MultipleMappedRate, "%")), Sample.Number=str_remove(CellID, "RNA_SubLib_"), MappingRate=UniquelyMappedPercent+MultipleMappedPercent) -> map_stat
patch_condition <- read.csv("data/meta/condition_patch_seq.csv")
map_stat$condition <- patch_condition$condition[match(map_stat$Sample.Number, patch_condition$Sample.Number)]
map_stat %>% mutate(QC=case_when(MappingRate<50 | UniquelyMappedReads<500000 ~ "Fail", .default="Pass")) -> sample_information
write.csv(sample_information, file="data/meta/sample_information.csv", row.names=FALSE, quote=FALSE)

