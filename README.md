
# MEMDS output sorting, SLURM version

The pipeline is designed to run on distributed computing systems (clusters) managed by SLURM
(https://slurm.schedmd.com/overview.html)

## Prerequisite programs:
The pipeline requires Conda package managment program to run

#### Install miniconda: (https://docs.conda.io/en/latest/miniconda.html)
#### Install bioconda: (https://bioconda.github.io)

#### Use Conda to create environment with the following programs:
conda create -n modules3  python=2.7  bwa cutadapt fastqc pear perl picard pysam biopython samtools seqtk trimmomatic

#### Check that all the programs are present in the newly created environment:
bwa  0.7.17
samtools 1.9
cutadapt 1.18
fastqc 0.11.8
pear 0.9.6
python 2.7.15
pysam 0.15.2
biopython 1.73
trimmomatic 0.39
perl 5.26.2
seqtk 1.3

## Pipeline usage:

#### a) Navigate to the folder containing pipeline scripts: 
cd path/to/project/scripts
#### b) Prepare the parameter files needed to run the pipeline, as outlined in the guide **"Parameter_file_preparation.pdf"**
-------
#### 1) Joining partial \*.fastq files (if data is separated across several lanes):
a) **Prepare:** more_scripts/samples_table_0.txt\
b) **Run:** *bash concatenate_partfiles.sh*\
c) **Check output file:** "more_scripts/samples_table_0.sh.concat.sh"; make sure it doesn't include comments indicating errors\
d) **Run**: *srun bash more_scripts/samples_table_0.sh.concat.sh*\
e) **Check the concatenated output files**\

#### 2) Organizing parameter files:
a) **Prepare in the "scripts/design" folder**: params_1.sh; samples_table.txt; factors_table.txt\
b) In the **sequences directory** (as listed in **params_1.sh**), prepare reference fasta files. Each reference file should contain single reference sequence. All reference files should have the extension **.fa**.\
c) **Run**: \
   *bash setting_1-PE.sh* (Paired-end data) **or**\
   *bash setting_1-SE.sh* (Single-end data)\
d) **Check**: output file **design/samples_table.sh**. Make sure it contains all parameter data and doesn't include comments indicating errors.\

#### 3) Quality control and merging (for paired-end data):\
a) **Run**: \
    *bash filter-PE4.sh $i* ($i in 1...5 (5 is optional)) for paired-end data **or**\
    *bash filter-SE4.sh $i* ($i in 1...4 (4 is optional)) for single-end data.\
   **Note**: All steps are running on Slurm. The steps should be run in a sequential manner; **don't start** a new step before completion of the previous one!\
b) In each step check that the jobs are completed on the grid using **"sacct"** command.\
c) **Check** the output at each step:\
   **step1:** Check Fastqc \*.html output for quality analysis of the raw data.\
   **step2: For paired-end data only.** Check Pear results (merged fastq files). Check **\*.assembled.info** file for the percentage of reads that were merged.\
   **step3:** Check quality-filtered fastq files - Cutadapt and Trimmomatic output and log files (Each filtered fastq can be viewed with the tool "less" in linux).\
   **step4:** Check Fastqc \*.html output for quality analysis of the filtered data.\
   **step5: Optional.** Check subsampled reads' data for "sanity checks" in **tests** directory.

#### 4) Separating barcode sequences and gene sequence:
a) **Run:** *bash trim7.sh  1*\
b) **Check** that the jobs are completed on the grid using **"sacct"** command.\
c) **Check** in the **filtered** directory, trimming information files: **\*.trimmed.barcodes** files (correctly trimmed), **\*.wrongId.barcodes** files (unexpected identifiers), and **\*.trimmed.log** files (stating percentage of reads with correct identifiers).\ 
   Make sure that correct barcode sequences are recognized by the script and that most reads in the data contain correct barcodes!\ 

#### 5) Sorting paralogous genes based on unique sequence signature
a) **Run:** *bash sort2.sh 1*\
b) **Check** that the jobs are completed on the grid using **"sacct"** command.\
c) **Check** sorted fastq files in the **"sorted"** directory.\ 
   The **\*.withIndels.log** files contain all the correctly matched reads with a shift in the frame.\ 
   The **\*.others.fastq** files contain reads tha couldn't be matched to either reference sequence.

#### 6) Mapping filtered, trimmed and sorted reads to reference:
a) **Run:** *bash bwa9.sh 1*\
b) **Check** that dictionary files are generated for each **\*.fa** file in the **reference sequence directory**\
c) **Run:** *bash bwa9.sh 2*\
d) **Check** that the jobs are completed on the grid using **"sacct"** command.\
e) **Check** that read-reference alignment files (**\*.sam** and **\*.bam**) are created in the **mapping** directory for each sorted fastq file.

#### 7) Adjust data format so the mapping results can be viewed with IGV viewer.
bash create_dummy_genome5.sh 1

#### after this stage, all files for IGV are ready (visual inspection of the alignments) in the 'mapping' directory: *.bam, *.bam.bai, *.barcode.fasta

#### 8) detection of potential mutations
bash sam_to_mutation-table_5.0.sh 1
#### outputs are directories tables.BC3cutoff1 and tables.BC3cutoff2, with  mutationFrequncyPerBarcode.txt 
#### tables.BC3cutoff1 uses all tail barcodes for 3'BC the criterion and tables.BC3cutoff2 uses only tail barcodes that appear more than once in a family for the 3'BC criterion)

#### 9) detection of consensus mutations
#### $i in 1...2
bash consensus_15.sh $i

#### then check the final mutation tables

bash sam_to_mutation-list-3.sh 1
