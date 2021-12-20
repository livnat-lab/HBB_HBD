
### Activate Conda module ############
. /path/to/miniconda2/etc/profile.d/conda.sh
conda activate modules3
#############################

# Set work directories
params_adapters_1="$PWD/wildcard_adapters_1.fa"
params_dir_out_1='/data/home/My_experiment/out1'
params_dir_reference='/data/home/My_experiment/seqs'

# Set run parameters
params_ninimum_fastq_size_1=80
is_SE=0
TSS=1
