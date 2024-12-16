# Data analysis to find and classify ASVs in 18S SSU data

1. Trimmomatic (v0.39) to remove the primer from the beginning of the read and to trim lower quality bases at the end of the read.  

2. USEARCH (v11) for ASV identification.  

3. Classification based on the SILVA 138.1 and PR2 v5 databases using mothur classify.seqs.  
The two sets of classifications were manually curated for downstream analyses.  

*Analyses done by J. M. Polinski*
