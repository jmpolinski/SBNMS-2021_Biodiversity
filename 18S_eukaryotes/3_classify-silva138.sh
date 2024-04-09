#!/bin/bash

# classify 18S
module load mothur
mothur "#classify.seqs(fasta=./SBNMS21-18S_zotus.fa, reference=/data/app/databases/silva-138.1/silva.nr_v138_1.align, taxonomy=/data/app/databases/silva-138.1/silva.nr_v138_1.tax)"


