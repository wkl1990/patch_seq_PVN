#!/usr/bin/R

suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("viridis"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("GenomicAlignments"))
suppressPackageStartupMessages(library("edgeR"))
suppressPackageStartupMessages(library("rtracklayer"))
suppressPackageStartupMessages(library("GenomicRanges"))

# * gtf annotation
gtfFile <- file.path("data/resource//modified_gencode.vM23.annotation.gtf")
gtf0 <- import(gtfFile)
idx_exon <- mcols(gtf0)$type == "exon"
gtf_exon <- gtf0[idx_exon]
idx_intron <- mcols(gtf0)$type == "intron"
gtf_intron <- gtf0[idx_intron]
exon_genes <- split(gtf_exon, mcols(gtf_exon)$gene_name)
intron_genes <- split(gtf_intron, mcols(gtf_intron)$gene_name)

# * bam files
fls <- list.files("bam/2_pass/final", recursive=TRUE, pattern="*markDup.removeDupl.bam$", full=TRUE)
names(fls) <- str_split(fls,"/",simplify=TRUE)[,11]
flag <- scanBamFlag(isUnmappedQuery=FALSE, isDuplicate=FALSE, isNotPassingQualityControls=FALSE, isPaired=NA)
param <- ScanBamParam(flag=flag)
bamlst <- BamFileList(fls)

# * get counts
exon_genehits <- summarizeOverlaps(exon_genes, bamlst, mode="IntersectionNotEmpty", param=param)
intron_genehits <- summarizeOverlaps(intron_genes, bamlst, mode="IntersectionNotEmpty", param=param)
exon_count <- assay(exon_genehits) %>% as.data.frame %>% rownames_to_column("gene_name")
intron_count <- assay(intron_genehits) %>% as.data.frame %>% rownames_to_column("gene_name")

gene_count <- bind_rows(exon_count, intron_count) %>% group_by(gene_name) %>% summarise_all(., sum, na.rm=TRUE)
gene_count_mtx <- gene_count %>% column_to_rownames("gene_name") %>% as.matrix
gene_cpm <- cpm(gene_count_mtx)

# * detected genes
gene_dection <- colSums(gene_cpm>0)
sample_information <- readRDS(file="data/meta/sample_information.rds")
sample_information$gene_detection <- gene_dection[match(sample_information$CellID, names(gene_dection))]

# * save data
saveRDS(gene_count_mtx, file="data/meta/gene_count_mtx.rds")
saveRDS(gene_cpm, file="data/meta/gene_cpm.rds")






