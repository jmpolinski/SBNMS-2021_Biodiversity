# Stellwagen Bank National Marine Sanctuary 2021

12S amplicon sequencing data from SBNMS 2021 ecosystem team trip. See https://github.com/emmastrand/GMGI_Notebook/blob/main/posts/2024-02-27%20eDNA%20workflow%2012S.md for more details on the programs used.

### Workflow steps 

0. FastQC (`00-fastqc.sh`)  
1. Nf-core Ampliseq: Implements CutAdapt and DADA2 pipelines (`01-ampliseq.sh`).     
2. Taxonomic identification with blastn and taxonkit (`02-tax_ID.sh`).  
3. Taxonomic assignment and data table preparation (`03-datatable_prep_report_generation.Rmd`)

## 00. FastQC

I previously set up conda environments with multiqc installed. I should reorganize these environments, but for now use `haddock_methylation` by:  
- `source ~/../../work/gmgi/miniconda3/bin/activate`  
- `conda activate haddock_methylation`

To run a slurm array, create a rawdata list `ls -d /work/gmgi/ecosystem-diversity/Stellwagen/sbnms-2021/12S/raw_data/*.gz > rawdata`. There are 356 samples within this raw data folder. 

`00-fastqc.sh`:  

```
#!/bin/bash
#SBATCH --error=fastqc_output/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=fastqc_output/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --time=23:00:00
#SBATCH --job-name=SB_fastqc
#SBATCH --mem=5GB
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=2

module load OpenJDK/19.0.1 ## dependency on NU Discovery cluster 
module load fastqc/0.11.9

raw_path="/work/gmgi/ecosystem-diversity/Stellwagen/sbnms-2021/12S/raw_data"
dir="/work/gmgi/ecosystem-diversity/Stellwagen/sbnms-2021/12S/QC/raw_fastqc"

## File name based on rawdata list
mapfile -t FILENAMES < ${raw_path}/rawdata
i=${FILENAMES[$SLURM_ARRAY_TASK_ID]}

## FastQC program
fastqc ${i} --outdir ${dir}
```

To run slurm array = `sbatch --array=0-356 00-fastqc.sh`.

Once complete, `cat *output.* > ../fastqc_output.txt` to create one file with all the output. The length of this file should be 356. 

### Multiqc 

`00-mutliqc.sh`: 

```
#!/bin/bash
#SBATCH --error="%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output="%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --time=23:00:00
#SBATCH --job-name=multiqc
#SBATCH --mem=5GB
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=2

dir="/work/gmgi/ecosystem-diversity/Stellwagen/sbnms-2021/12S/QC/raw_fastqc"
multiqc_dir="/work/gmgi/ecosystem-diversity/Stellwagen/sbnms-2021/12S/QC/multiqc"

multiqc --interactive ${dir} -o ${multiqc_dir} --filename multiqc_raw.html
```

Check this output before proceeding.

## 01. Ampliseq 

**Create samplesheet**:

This was created with a custom R script on NU cluster. 

**Primer information**:

MiFish 12S amplicon F: ACTGGGATTAGATACCCC    
MiFish 12S amplicon R: TAGAACAGGCTCCTCTAG  

**Run ampliseq script**:

`01-ampliseq.sh`:

```
#!/bin/bash
#SBATCH --error="%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output="%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --time=20:00:00
#SBATCH --job-name=ampliseq_SB
#SBATCH --mem=50GB
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=2

# singularity module version 3.10.3
module load singularity/3.10.3

# nextflow module loaded on NU cluster is v23.10.1
module load nextflow/23.10.1

#set paths 
metadata="/work/gmgi/ecosystem-diversity/Stellwagen/sbnms-2021/12S/metadata" 

cd /work/gmgi/ecosystem-diversity/Stellwagen/sbnms-2021/12S

nextflow run nf-core/ampliseq -resume \
   -profile singularity \
   --input ${metadata}/samplesheet.csv \
   --FW_primer "ACTGGGATTAGATACCCC" \
   --RV_primer "TAGAACAGGCTCCTCTAG" \
   --outdir ./results \
   --trunclenf 100 \
   --trunclenr 100 \
   --trunc_qmin 25 \
   --max_len 200 \
   --max_ee 2 \
   --min_len_asv 100 \
   --max_len_asv 115 \
   --sample_inference pseudo \
   --skip_taxonomy \
   --ignore_failed_trimming
```

With 50 GB, 356 samples took ~ minutes.

# 2. Taxonomic Identification

We use GMGI, Mitofish, and NCBI databases. 

The following script takes the `ASV_seqs.len.fasta` file from DADA2 output (`asv_length_filter/`) and uses `blastn` to compare sequences to three databases: Blast nt, Mitofish, and our in-house GMGI database.

Make a directory called `BLASToutput` within the project directory.

`02-taxonomicID.sh`: 

```
#!/bin/bash
#SBATCH --error=output_messages/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=output_messages/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --time=20:00:00
#SBATCH --job-name=tax_ID
#SBATCH --mem=30GB
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=2

## load modules needed
module load ncbi-blast+/2.13.0

## set internet connection
export http_proxy="http://10.99.0.130:3128"
export https_proxy="http://10.99.0.130:3128"
export ftp_proxy="http://10.99.0.130:3128"

# set path to ASV_seqs.fasta
## needs to be changed for every project
fasta="/work/gmgi/ecosystem-diversity/Stellwagen/sbnms-2021/12S/results/asv_length_filter"
out="/work/gmgi/ecosystem-diversity/Stellwagen/sbnms-2021/12S/BLASToutput"
db="/work/gmgi/Fisheries/databases/12S/reference_fasta"
taxonkit="/work/gmgi/Fisheries/databases/taxonkit"

#### DATABASE QUERY ####
### NCBI database 
blastn -remote -db nt \
   -query ${fasta}/ASV_seqs.len.fasta \
   -out ${out}/BLASTResults_NCBI.txt \
   -max_target_seqs 10 -perc_identity 100 -qcov_hsp_perc 95 \
   -outfmt '6  qseqid   sseqid   sscinames   staxid pident   length   mismatch gapopen  qstart   qend  sstart   send  evalue   bitscore'

## Mitofish database 
blastn -db ${db}/Mitofish.fasta \
   -query ${fasta}/ASV_seqs.len.fasta \
   -out ${out}/BLASTResults_Mito.txt \
   -max_target_seqs 10 -perc_identity 100 -qcov_hsp_perc 95 \
   -outfmt '6  qseqid   sseqid  pident   length   mismatch gapopen  qstart   qend  sstart   send  evalue   bitscore'

## GMGI database 
blastn -db ${db}/GMGIVertRef.fasta \
   -query ${fasta}/ASV_seqs.len.fasta \
   -out ${out}/BLASTResults_GMGI.txt \
   -max_target_seqs 10 -perc_identity 98 -qcov_hsp_perc 95 \
   -outfmt '6  qseqid   sseqid   pident   length   mismatch gapopen  qstart   qend  sstart   send  evalue   bitscore'

############################

#### TAXONOMIC CLASSIFICATION #### 
## creating list of staxids from all three files 
awk -F $'\t' '{ print $4}' ${out}/BLASTResults_NCBI.txt | sort -u > ${out}/NCBI_sp.txt

## annotating taxid with full taxonomic classification
cat ${out}/NCBI_sp.txt | ${taxonkit}/taxonkit reformat -I 1 -r "Unassigned" > NCBI_taxassigned.txt
```

### Output 

- `BLASTResults_GMGI.txt`: GMGI's curated database results    
- `BLASTResults_Mito.txt`: MitoFish database results    
- `BLASTResults_NCBI.txt`: NCBI results    
- `NCBI_taxassigned.txt`: full taxonomic classification for all staxids found in NCBI output. 

Files to use within `/work/gmgi/Fisheries/databases/12S`: `taxonomic_classification_fishbase.csv` and `taxonomic_classification_mitofish.csv`. 

Switch to `Datatable preparation` Rmarkdown. 