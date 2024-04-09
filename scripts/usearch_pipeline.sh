#!/bin/bash
##### USEARCH PIPELINE TO MERGE READS, REMOVE DUPLICATES, FIND zOTUs #####

# set project ID for files
PRJ=$1

if [[ -z $PRJ ]];
then echo "Missing input to run command. Please provide an ID to put at front of output files"; exit;
else echo "All variables provided. Proceeding to USEARCH analysis."; fi

# merge read pairs
echo "merging read pairs"
date
usearch -fastq_mergepairs *R1*.fq -relabel @ -fastq_maxdiffs 10 -fastq_pctid 80 -fastqout ${PRJ}_merged.fastq

# filter
echo "filtering fastq"
date
usearch -fastq_filter ${PRJ}_merged.fastq -fastq_maxee 1.0 -fastaout ${PRJ}_filtered.fa

# dereplicate
echo "dereplicate filtered sequences"
date
usearch -fastx_uniques ${PRJ}_filtered.fa -fastaout ${PRJ}_uniques.fasta -sizeout -relabel ${PRJ}-uniq

# denoise/find zero-radius OTUs
echo "denoise/find ASVs"
date
usearch -unoise3 ${PRJ}_uniques.fasta -zotus ${PRJ}_zotus.fa -minsize 2

# map merged reads to zOTUs to get count table
echo "generating zOTU count table"
date
usearch -otutab ${PRJ}_merged.fastq -zotus ${PRJ}_zotus.fa -otutabout ${PRJ}_zotutab.txt -mapout ${PRJ}_zmap.txt -sample_delim .
echo "usearch analysis complete"
date

