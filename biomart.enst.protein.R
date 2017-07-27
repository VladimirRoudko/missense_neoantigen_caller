#!/usr/bin/env Rscript

library("biomaRt")

args = commandArgs(trailingOnly=TRUE)

mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
data <- read.table(args[1], sep=" ")
ENST_ID <- t(data$V4)
seq = getSequence(id=ENST_ID, type="ensembl_transcript_id", seqType="peptide", mart = mart, verbose=FALSE)
exportFASTA(seq, args[2])
