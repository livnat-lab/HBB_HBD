# A wrapper script for paired-end raw ".fastq" file quality analysis, merging and trimming

params_1="$PWD/design/params_1.sh"
params_2="$PWD/design/samples_table.sh"

###########################################
# Gather pipeline parameter data; check that relevant parameters are non-empty strings and create output folders

function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }

. "$params_1"
assert_
. "$params_2"
assert_

if [ ! -d "$params_dir_out_1" ] || [ ! -f "$params_adapters_1" ] || \
   [ ! -n "${params_ninimum_fastq_size_1+1}" ]; then
    exit "Error: not all parameters/files exit"
fi


echo "Out dir: $params_dir_out_1"
if [ ! -d "$params_dir_out_1/filtered" ]; then mkdir "$params_dir_out_1/filtered"; assert_; fi
if [ ! -d "$params_dir_out_1/fastqc" ]; then mkdir "$params_dir_out_1/fastqc"; assert_; fi
if [ ! -d "$params_dir_out_1/tests" ]; then mkdir "$params_dir_out_1/tests"; assert_; fi

#########
# Execute QC jobs on each analyzed sample
for i in ${!title[@]}; do
    echo $i
    
    title1="${title[i]}"_"idx$i"
    f1="${for1[i]}"
    f2="${rev1[i]}"
    
    echo "forward reads = $f1"
    echo "reverse reads = $f2"
    echo "title = $title1"
    
    out_prefix1="$params_dir_out_1/filtered/$title1"
    echo "out = $out_prefix1"
    
    # Validate existence of the input files
    if [ -f "$f1" ] && [ -f "$f2" ]; then
        
        # Run FastQC to assess data quality
        if [ $1 -eq 1 ]; then
            sbatch -N1 -n1 --ntasks-per-node=1 --time=120 -o "$params_dir_out_1/fastqc/$title1.f.out" -e "$params_dir_out_1/fastqc/$title1.f.err" \
                -p hive1d,hive7d,hiveunlim,queen,ckptqueen,ckpthp1d,ckpthp7d,ckpthp31d,ckptdell1d,ckptdell7d,ckptdell31d,preempt1d,preempt7d,preempt31d \
                --wrap "fastqc -o \"$params_dir_out_1/fastqc\" \"$f1\""
            sbatch -N1 -n1 --ntasks-per-node=1 --time=120 -o "$params_dir_out_1/fastqc/$title1.r.out" -e "$params_dir_out_1/fastqc/$title1.r.err" \
                -p hive1d,hive7d,hiveunlim,queen,ckptqueen,ckpthp1d,ckpthp7d,ckpthp31d,ckptdell1d,ckptdell7d,ckptdell31d,preempt1d,preempt7d,preempt31d \
                --wrap "fastqc -o \"$params_dir_out_1/fastqc\" \"$f2\""
        
        # Merge forward and reverse reads
        elif [ $1 -eq 2 ]; then
            sbatch --output="$out_prefix1.assembled.out" --error="$out_prefix1.assembled.err" \
                   -N1 -n1 --ntasks-per-node=1 --time=180 -p hive1d,hive7d,hiveunlim,queen,ckptqueen,ckpthp1d,ckpthp7d,ckpthp31d,ckptdell1d,ckptdell7d,ckptdell31d,preempt1d,preempt7d,preempt31d \
                    --wrap "pear -f \"$f1\" -r \"$f2\" -o \"$out_prefix1\" -j 1 | tee \"$out_prefix1.assembled.info\""
        
        # Trim merged reads to remove contaminants and low quality data
        elif [ $1 -eq 3 ]; then
            if [ -f "$out_prefix1.assembled.fastq" ]; then
                sbatch --output="$out_prefix1.assemb.filtered.out" --error="$out_prefix1.assemb.filtered.err" \
					-N1 -n1 --ntasks-per-node=1 --time=720 -p hive1d,hive7d,hiveunlim,queen,ckptqueen,ckpthp1d,ckpthp7d,ckpthp31d,ckptdell1d,ckptdell7d,ckptdell31d,preempt1d,preempt7d,preempt31d \
                    "$PWD/fastq-filter_job_3.sh" "SE" "$params_adapters_1" 30 3 \
                    "$params_ninimum_fastq_size_1" "$out_prefix1.assembled.fastq" \
                    "$out_prefix1.assembled.filtered.fastq"
            fi
        
        # Run FastQC to assess data quality after merging and trimming
        elif [ $1 -eq 4 ]; then
            if [ -f "$out_prefix1.assembled.filtered.fastq" ]; then
                sbatch -N1 -n1 --ntasks-per-node=1 --time=120 -o "$params_dir_out_1/fastqc/$title1.assembled.fastqc.out" \
                    -p hive1d,hive7d,hiveunlim,queen,ckptqueen,ckpthp1d,ckpthp7d,ckpthp31d,ckptdell1d,ckptdell7d,ckptdell31d,preempt1d,preempt7d,preempt31d \
                    --wrap "fastqc -o \"$params_dir_out_1/fastqc\" \"$out_prefix1.assembled.filtered.fastq\""
            fi
        
        # Create small subsamples of analyzed ".fastq" files for manual analyses, if needed:
        elif [ $1 -eq 5 ]; then
            if [ -f "$out_prefix1.assembled.fastq" ]; then
                sbatch --output="$params_dir_out_1/tests/$title1.assembled.subsample1000.out" --error="$params_dir_out_1/tests/$title1.assembled.subsample1000.err" \
					-N1 -n1 --ntasks-per-node=1 --time=60 -p hive1d,hive7d,hiveunlim,queen,ckptqueen,ckpthp1d,ckpthp7d,ckpthp31d,ckptdell1d,ckptdell7d,ckptdell31d,preempt1d,preempt7d,preempt31d \
                    --wrap "seqtk sample \"$out_prefix1.assembled.fastq\" 10000 > \"$params_dir_out_1/tests/$title1.assembled.subsample1000.fastq\""
            fi
            if [ -f "$f1" ] && [ -f "$f2" ]; then
                sbatch --output="$params_dir_out_1/tests/$title1.f.subsample1000.out" --error="$params_dir_out_1/tests/$title1.f.subsample1000.err" \
					-N1 -n1 --ntasks-per-node=1 --time=60 -p hive1d,hive7d,hiveunlim,queen,ckptqueen,ckpthp1d,ckpthp7d,ckpthp31d,ckptdell1d,ckptdell7d,ckptdell31d,preempt1d,preempt7d,preempt31d \
                    --wrap "seqtk sample -s100 \"$f1\" 10000 > \"$params_dir_out_1/tests/$title1.f.subsample1000.fastq\""
                sbatch --output="$params_dir_out_1/tests/$title1.r.subsample1000.out" --error="$params_dir_out_1/tests/$title1.r.subsample1000.err" \
					-N1 -n1 --ntasks-per-node=1 --time=60 -p hive1d,hive7d,hiveunlim,queen,ckptqueen,ckpthp1d,ckpthp7d,ckpthp31d,ckptdell1d,ckptdell7d,ckptdell31d,preempt1d,preempt7d,preempt31d \
                    --wrap "seqtk sample -s100 \"$f2\" 10000 > \"$params_dir_out_1/tests/$title1.r.subsample1000.fastq\""
            fi
        fi
    else
        echo 'At least one of the input files does not exist'
    fi
    
    echo '----------------'
done
