
## BC sorting, SLURM version
(https://slurm.schedmd.com/overview.html)

## instalation instuctions:

#### install miniconda (https://docs.conda.io/en/latest/miniconda.html)
#### install bioconda (https://bioconda.github.io)
#### then:
conda create -n modules3  python=2.7  bwa cutadapt fastqc pear perl picard pysam biopython samtools seqtk multiqc trimmomatic

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

## Run instructions 

cd path/to/project/scripts

#### 1) join fastq partfiles, if needed
#### prepare:
#### more_scripts/samples_table_0.txt

bash concatenate_partfiles.sh

#### then check output file: "more_scripts/samples_table_0.sh.concat.sh", make sure it doesn't include comments indicating errors

srun bash more_scripts/samples_table_0.sh.concat.sh

#### then check the concatenated output files

#### 2) 
#### prepare:
#### design/params_1.sh
#### design/samples_table.txt
#### design/factors_table.txt
#### in the sequences directory, prepare reference fasta files, with the extension *.fa
bash setting_1-PE.sh
#### or:
bash setting_1-SE.sh

#### then check output file design/samples_table.sh, make sure it doesn't include comments indicating errors

#### 3) merge and quality filter

#### for paired-end data:
#### $i in 1...5 (5 is optional). All steps are running on Slurm.
bash filter-PE4.sh $i

#### or for single end data:
#### $i in 1...4 (4 is optional). All steps are running on Slurm.
bash filter-SE4.sh $i

#### in this and all following steps: check the jobs are completed on the grid using "sacct"
#### then check the quality and read-count results of the raw and filtered/assembled *.fastq files, for example in the Fastqc output HTML files, in "filtered" directory.
#### step1 : check  also fastqc *.html graphs for raw fastq
#### step2 : check Pear results (merged fastq files, look at assembled.info)
#### step3 : check quality-filtered fastq files (cutadapt + trimmomatic output and log files. Each filtered fastq can be viewed with the tool "less" in linux)
#### step4: check  also fastqc *.html graphs for filtered fastq
#### step5: check subsampled reads data for "sanity checks" in tests directory

#### 4) Select the region of interest
bash trim7.sh  1

#### then check, in filtered directory, trimming quality information files *.trimmed.barcodes files (correctly trimmed), *.wrongId.barcodes files (unexpected identifiers), and *.trimmed.log files.
#### if needed, further test the result fastq files.

#### 5) Sorting of paralogous genes based on unique sequence signature
bash sort2.sh 1
#### check "sorted" directory. The indels files contain all the correctly matched reads with a shift in the frame.

#### 6) Reads mapping
#### Creates sorted.bam files with the read-reference alighnment (to be used by the pipline and IGV), and sorted.bam.bai files (to be used by IGV) 
#### $i in 1...2
bash bwa9.sh  $i

#### 7) Adjust data format so the mapping results can be viewed with IGV viewer.
bash create_dummy_genome5.sh 1

#### after this stage, all files for IGV are ready (visual inspection of the alignments) in the 'mapping' directory: *.bam, *.bam.bai, *.barcode.fasta

#### 8)
bash sam_to_mutation-table_5.0.sh 1
#### outputs are directories tables.BC3cutoff1 and tables.BC3cutoff2, with  mutationFrequncyPerBarcode.txt 
#### tables.BC3cutoff1 uses all tail barcodes for 3'BC the criterion and tables.BC3cutoff2 uses only tail barcodes that appear more than once in a family for the 3'BC criterion)

#### 9)
#### $i in 1...2
bash consensus_15.sh $i

#### then check the final mutation tables

bash sam_to_mutation-list-3.sh 1
