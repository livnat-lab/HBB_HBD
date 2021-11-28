# A wrapper script to run BWA alignment jobs

params_1="$PWD/design/params_1.sh"
params_2="$PWD/design/samples_table.sh"

#####################################
# Gather pipeline parameter data and check that relevant parameters are non-empty strings and create output folders

function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }

. $params_1
assert_
. $params_2
assert_

if [ ! -d "$params_dir_out_1" ] || [ ! -d "$params_dir_out_1/filtered" ] || [ ! -n "${is_SE+1}" ] || \
    [ ! -n "${params_dir_reference+1}" ] || [ ! -n "${sort_refs+1}" ] || [ ! -n "${sort_ref+1}" ] || [ ! -n "${reference_size+1}" ]; then
    exit "Error: not all parameters/files exit"
fi


if [ ! -d "$params_dir_out_1/mapping" ]; then mkdir "$params_dir_out_1/mapping"; assert_; fi

#########
# Create dictionary files from reference sequences and run sequence alignment jobs
if [ $1 -eq 1 ]; then
    # for ref1 in  $(echo "${reference[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '); do # unique (non-redundant) names only
    for ref1 in "$params_dir_reference/"*.fa; do
        echo $ref1
        if [ ! -f $ref1.indices.OK ] && [ $1 -eq 1 ] && [ -f "$ref1" ]; then
            bwa index $ref1
            assert_
            samtools faidx $ref1
            assert_
            #java -jar $params_picardDir/picard.jar CreateSequenceDictionary REFERENCE=$ref1 OUTPUT=$ref1.dict 
			picard CreateSequenceDictionary REFERENCE=$ref1 OUTPUT=$ref1.dict 
            assert_
            touch $ref1.indices.OK
            assert_
        fi
    done
fi

# Run read aligning jobs on the sorted reads
if [ $1 -eq 2 ]; then
    for i in ${!title[@]}; do
        title1="${title[i]}"_"idx$i"; assert_
        prefix1="$params_dir_out_1/sorted/$title1";
        refs1=$(echo "${sort_ref[i]}"";""${sort_refs[i]}" | tr ";" "\n")
        
        echo "title = $title1"
        echo "prefix = $prefix1"
        
        # Iterate over sorted read files to map the reads against their reference
        k=0
        for r1 in $refs1; do
            let k=$k+1
            echo $k
            r2="$r1"
            if [ $k -eq 1 ]; then r2="$r2.others"; fi # first item in loop refers to "others"
            
            f1="$prefix1.$r2.fastq"
            fout1="$params_dir_out_1/mapping/$title1.$r2.bwa"
            ref1="$params_dir_reference"/"$r1".fa
            ref1_length="${reference_size[$i]}"
            
            echo 'len = '$ref1_length
            echo 'ref = '$ref1
            echo 'f1 = '$f1
            echo 'fout = '$fout1
            
            # Check that all required parameters are OK before executing the jobs
            if [ -f "$f1" ] && [ -d "$params_dir_out_1/mapping" ] && [ -f "$ref1" ] && [ -f "$ref1.indices.OK" ] && [ $ref1_length != "" ]; then
                echo 'OK'
                if [ $1 -eq 2 ]; then
                    sbatch  -N1 -n20 --ntasks-per-node=20 \
                            -p hive1d,hive7d,hiveunlim,queen,ckptqueen,ckpthp1d,ckpthp7d,ckpthp31d,ckptdell1d,ckptdell7d,ckptdell31d,preempt1d,preempt7d,preempt31d \
                            -o "$fout1".out -e "$fout1".err \
                            --export=ref1="$ref1",f1="$f1",fout1="$fout1",threads=20,ref1_length="$ref1_length",params_1="$params_1" \
                            bwa_job9.sh
                fi
            else
                echo 'Error: not all input files or parameters exist'
            fi
        done
        echo '--------------'
    done
fi

