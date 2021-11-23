
params_1="$PWD/design/params_1.sh"
params_2="$PWD/design/samples_table.sh"

###########################################

function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }

. "$params_1"
assert_
. "$params_2"
assert_

if [ ! -d "$params_dir_out_1" ] || [ ! -f "$params_adapters_1" ] || \
   [ ! -n "${params_ninimum_fastq_size_1+1}" ]; then
    exit "error: not all parameters/files exit"
fi
#[ ! -f "$params_fastqc_1" ]  ||

echo "out dir: $params_dir_out_1"
if [ ! -d "$params_dir_out_1/filtered" ]; then mkdir "$params_dir_out_1/filtered"; assert_; fi
if [ ! -d "$params_dir_out_1/fastqc" ]; then mkdir "$params_dir_out_1/fastqc"; assert_; fi
if [ ! -d "$params_dir_out_1/tests" ]; then mkdir "$params_dir_out_1/tests"; assert_; fi

for i in ${!title[@]}; do
	echo $i
    title1="${title[i]}"_"idx$i"
    f1="${for1[i]}"
    echo "f1 = $f1"
	echo "title = $title1"
    out_prefix1="$params_dir_out_1/filtered/$title1"
    echo "out = $out_prefix1"
    if [ -f "$f1" ]; then
        if [ $1 -eq 1 ]; then
            sbatch -N1 -n1 --ntasks-per-node=1 --time=120 -o "$params_dir_out_1/fastqc/$title1.f.out" \
                -p hive1d,hive7d,hiveunlim,queen,ckptqueen,ckpthp1d,ckpthp7d,ckpthp31d,ckptdell1d,ckptdell7d,ckptdell31d,preempt1d,preempt7d,preempt31d \
                --wrap "fastqc -o \"$params_dir_out_1/fastqc\" \"$f1\""
        elif [ $1 -eq 2 ]; then
            sbatch --output="$out_prefix1.filtered.out" --error="$out_prefix1.filtered.err" \
					-N1 -n1 --ntasks-per-node=1 --time=720 -p hive1d,hive7d,hiveunlim,queen,ckptqueen,ckpthp1d,ckpthp7d,ckpthp31d,ckptdell1d,ckptdell7d,ckptdell31d,preempt1d,preempt7d,preempt31d \
                    "$PWD/fastq-filter_job_3.sh" "SE" "$params_adapters_1" 30 3 \
                    "$params_ninimum_fastq_size_1" "$f1" "$out_prefix1.filtered.fastq"
        elif [ $1 -eq 3 ]; then
            if [ -f "$out_prefix1.filtered.fastq" ] && [ -d "$params_dir_out_1/fastqc" ]; then
                sbatch -N1 -n1 --ntasks-per-node=1 --time=120 --output="$params_dir_out_1/fastqc/$title.fastqc.out" --error="$params_dir_out_1/fastqc/$title.fastqc.err" \
                    -p hive1d,hive7d,hiveunlim,queen,ckptqueen,ckpthp1d,ckpthp7d,ckpthp31d,ckptdell1d,ckptdell7d,ckptdell31d,preempt1d,preempt7d,preempt31d \
                    --wrap "fastqc -o \"$params_dir_out_1/fastqc\" \"$out_prefix1.filtered.fastq\""
            fi
        elif [ $1 -eq 4 ]; then
            # to test fastq subsamples manually, if needed:
            if [ -f "$out_prefix1.filtered.fastq" ]; then
                sbatch --output="$params_dir_out_1/tests/$title1.subsample1000.out" \
					-N1 -n1 --ntasks-per-node=1 --time=60 -p hive1d,hive7d,hiveunlim,queen,ckptqueen,ckpthp1d,ckpthp7d,ckpthp31d,ckptdell1d,ckptdell7d,ckptdell31d,preempt1d,preempt7d,preempt31d \
                    --wrap "head -n400 \"$out_prefix1.filtered.fastq\" > \"$params_dir_out_1/tests/$title1.filtered.head100.fastq\""
                    
                    sbatch --output="$params_dir_out_1/tests/$title1.subsample1000.out" \
					-N1 -n1 --ntasks-per-node=1 --time=60 -p hive1d,hive7d,hiveunlim,queen,ckptqueen,ckpthp1d,ckpthp7d,ckpthp31d,ckptdell1d,ckptdell7d,ckptdell31d,preempt1d,preempt7d,preempt31d \
                    --wrap "seqtk sample -s 100 \"$out_prefix1.filtered.fastq\" 10000 > \"$params_dir_out_1/tests/$title1.filtered.random-subsample1000.fastq\""
            fi
            if [ -f "$f1" ]; then
                sbatch --output="$params_dir_out_1/tests/$title1.f.subsample1000.out" \
					-N1 -n1 --ntasks-per-node=1 --time=60 -p hive1d,hive7d,hiveunlim,queen,ckptqueen,ckpthp1d,ckpthp7d,ckpthp31d,ckptdell1d,ckptdell7d,ckptdell31d,preempt1d,preempt7d,preempt31d \
                    --wrap "head -n400 \"$f1\" > \"$params_dir_out_1/tests/$title1.head100.fastq\""
                    
                sbatch --output="$params_dir_out_1/tests/$title1.f.subsample1000.out" \
					-N1 -n1 --ntasks-per-node=1 --time=60 -p hive1d,hive7d,hiveunlim,queen,ckptqueen,ckpthp1d,ckpthp7d,ckpthp31d,ckptdell1d,ckptdell7d,ckptdell31d,preempt1d,preempt7d,preempt31d \
                    --wrap "seqtk sample -s 100 \"$f1\" 10000 > \"$params_dir_out_1/tests/$title1.random-subsample1000.fastq\""
            fi
        fi
    else
        echo 'at least one of the files do not exist'
    fi
	echo '----------------'
done
