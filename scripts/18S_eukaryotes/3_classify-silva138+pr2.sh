#!/bin/bash

# classify 18S
module load mothur
mothur "#classify.seqs(fasta=./SBNMS21-18S_zotus.fa, reference=/data/app/databases/silva-138.1/silva.nr_v138_1.align, taxonomy=/data/app/databases/silva-138.1/silva.nr_v138_1.tax)"
mv SBNMS21-18S_zotus.0_SSU_mothur.wang.taxonomy SBNMS21-18S_silva138_mothur.wang.taxonomy
mv SBNMS21-18S_zotus.0_SSU_mothur.wang.tax.summary SBNMS21-18S_silva138_mothur.wang.tax.summary

mothur "#classify.seqs(fasta=../2_usearch/SBNMS21-18S_zotus.fa, reference=/data/app/databases/18s/pr2_version_5.0.0_SSU_mothur.fasta, taxonomy=/data/app/databases/18s/pr2_version_5.0.0_SSU_mothur.tax)"
mv SBNMS21-18S_zotus.0_SSU_mothur.wang.taxonomy SBNMS21-18S_PR2v5_mothur.wang.taxonomy
mv SBNMS21-18S_zotus.0_SSU_mothur.wang.tax.summary SBNMS21-18S_PR2v5_mothur.wang.tax.summary
