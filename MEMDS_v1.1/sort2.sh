# A wrapper script to run jobs sorting reads by their origin gene

params_1="$PWD/design/params_1.sh"
params_2="$PWD/design/samples_table.sh"

########################################################
# Gather pipeline parameter data and check that relevant parameters are non-empty strings and create output folders

function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }

. "$params_1"
assert_
. "$params_2"
assert_

if [ ! -d "$params_dir_out_1" ] || [ ! -d "$params_dir_out_1/filtered" ] || [ ! -n "${is_SE+1}" ]; then
    exit "Error: not all parameters/files exit"
fi
if [ ! -d  "$params_dir_out_1/sorted" ]; then
    mkdir "$params_dir_out_1/sorted"
    assert_
fi

#########
# Execute sorting jobs on each analyzed sample

for i in ${!title[@]}; do

    # Read in parameters need to sort the reads by their gene of origin
    title1="${title[i]}"_"idx$i"
    assert_
    prefix1="$params_dir_out_1/filtered/$title1"
    if [ $is_SE == 0 ]; then
        f1=$prefix1.assembled.filtered.fastq.trimmed.fastq
    else
        f1=$prefix1.filtered.fastq.trimmed.fastq
    fi
    
    pos="${sort_pos[$i]}"
    nucl="${sort_nucl[$i]}"
    refs="${sort_refs[$i]}"
    ref="${sort_ref[$i]}"
    match_per="${sort_match[$i]}"
    
    echo "fastq = $f1"
    echo "refs = $refs"
    echo "ref = $ref"
    echo "pos = $pos"
    echo "nucl = $nucl"
    echo "match_percent = $match_per"
    
    # Check that input file exists
    if [ -f "$f1" ]; then
        echo 'Input files are found'
       
        # Run the jobs, if '1' is specified; otherwise run the wrapper w/o job execution, to test that input parameters are correct
        if [ $1 -eq 1 ]; then
            echo "Sorting $title1"
            assert_
            sbatch --output="$params_dir_out_1/sorted/$title1.out" --error="$params_dir_out_1/sorted/$title1.err" -N1 -n1 --ntasks-per-node=1 \
                -p hive1d,hive7d,hiveunlim,queen,ckptqueen,ckpthp1d,ckpthp7d,ckpthp31d,ckptdell1d,ckptdell7d,ckptdell31d,preempt1d,preempt7d,preempt31d \
                --wrap "python sort2.py --fastq \"$f1\" --pos \"$pos\" --nucl \"$nucl\" --refs \"$refs\" --ref \"$ref\" --dirOut \"$params_dir_out_1/sorted\" --title \"$title1\" --match_percentage \"$match_per\""
        fi
    else
        echo 'Not all files exist'
    fi
    
    echo '-------------'
done

