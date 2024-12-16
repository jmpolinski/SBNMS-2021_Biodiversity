#!/bin/bash

#  working in directory /data/prj/ecosystem-diversity/Stellwagen/sbnms-2021/18S/Analyses/1_trimming/
# raw reads (combined from 3 sequencing runs) are in /data/prj/ecosystem-diversity/Stellwagen/sbnms-2021/18S/Data/combined-runs/

# create text file "samples" with samples IDs to use in a loop
ls ../../Data/combined-runs > samples
sed -i -e 's/_R.*._allRaw.fastq//g' samples
uniq samples > tmp && mv tmp samples

# Trim beginning of sequence to remove primer and end to remove low-quality bases (improves read merging)
for i in `cat samples`;  
do java -jar /data/app/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 5 -phred33 -summary sbnms-16s_trimmomatic.txt ../../Data/combined-runs/${i}_R1_allRaw.fastq ../../Data/combined-runs/${i}_R2_allRaw.fastq ${i}_R1_qc.fq ${i}_R1.orphan.fq ${i}_R2_qc.fq ${i}_R2.orphan.fq  HEADCROP:15 CROP:250 TRAILING:15 MINLEN:220; 
done

# remove unpaired reads, as they won't be used
rm *orphan.fq
