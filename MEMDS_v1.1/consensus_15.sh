# A wrapper script to run jobs summarizing mutation information per read family

params_1="$PWD/design/params_1.sh"
params_2="$PWD/design/samples_table.sh"

###########################################
# Gather pipeline parameter data

function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }

. $params_1
assert_
. $params_2
assert_

#########
# Run consensus tables generating jobs
for BC3_min_cutoff in 1; do
	
    # Check input directory existence and create output directory
    if [ ! -d "$params_dir_out_1/tables.BC3cutoff$BC3_min_cutoff" ]; then
    	exit "Error: not all parameters/files exit"
    fi

    if [ ! -d "$params_dir_out_1/tables_consensus.BC3cutoff$BC3_min_cutoff" ]; then 
	mkdir "$params_dir_out_1/tables_consensus.BC3cutoff$BC3_min_cutoff"; assert_; 
    fi
	
    # Execute the jobs on each analyzed sample
    for i in ${!title[@]}; do
        # Gather input parameters
        size_r1="${size_r[i]}"
        title1="${title[i]}"_"idx$i"
        refs1=$(echo "${sort_ref[i]}"";""${sort_refs[i]}" | tr ";" "\n")
        limit_start="${limit_starts[$i]}"
        limit_end="${limit_ends[$i]}"
        seq_pos="${seq_pos[i]}"
        seq_mut="${seq_mut[i]}"
        seq_action="${seq_action[i]}"
		
        echo "limits = $limit_start - $limit_end"
        echo "BC3_min_cutoff = $BC3_min_cutoff"
        
        # Iterate over sorted read mutation data and run consensus jobs	
        k=0
        for r1 in $refs1; do
            let k=$k+1
            echo "k = "$k
            r2="$r1"
            if [ $k -eq 1 ]; then r2="$r2.others"; fi # first item in loop refers to "others"
            
            bam_name="$title1.$r2.bwa.sorted.bam"
            ref1="$params_dir_reference"/"$r1".fa
            echo "bam_name = $bam_name"
            echo "ref1 = $ref1"
            
            for sam_BX_tag in '' 'BC3_missing'; do # To include only reads having specific BC3s
                echo "sam_BX_tag = \"$sam_BX_tag\""
                tag1=""
                if [ "$sam_BX_tag" != "" ]; then tag1="tag-$sam_BX_tag."; fi
                
                # Define I/O variables
                table_in1="$params_dir_out_1/tables.BC3cutoff$BC3_min_cutoff/$bam_name.""$tag1""mutationFrequncyPerBarcode.txt"
                table_out2_dir="$params_dir_out_1/tables_consensus.BC3cutoff$BC3_min_cutoff/$bam_name.""$tag1"".cutoffs-update"
                table_out3_dir="$params_dir_out_1/tables_consensus.BC3cutoff$BC3_min_cutoff/$bam_name.""$tag1"
                
                table_out1="$params_dir_out_1/tables_consensus.BC3cutoff$BC3_min_cutoff/$bam_name.""$tag1""consensus.txt"
                table_out2="$table_out2_dir/$bam_name"
                table_out3="$table_out3_dir/$bam_name"
                
                echo "table_in1 = $table_in1"
                echo "table_out1 = $table_out1"
                echo "table_out2 = $table_out2"
                echo "size_r1 = $size_r1"
                
                # Check that all required parameters are OK before executing the jobs
                if [ -f "$table_in1" ] && [ -d "$params_dir_out_1/tables_consensus.BC3cutoff$BC3_min_cutoff" ]; then
                    if [ ! -z "$table_in1" ] && [ ! -z "$table_out1" ] && [ ! -z "$table_out2" ] && [ ! -z "$table_out3" ] && [ ! -z "$size_r1" ] && [ ! -z "$title1" ] && [ ! -z "$BC3_min_cutoff" ]; then
                        echo 'All files and parameters are found'
                        
                        c=$(cat $table_in1 | wc -l)
                        echo "c="$c
                        if [ $c -gt 1 ]; then
                            echo 'Input table is non-empty'
                            
                            # Run jobs to create consensus mutation profile per read family
                            if [ $1 -eq 1 ]; then
								sleep 2
                                sbatch -o "$table_out1.out" -e "$table_out1.1.err" \
									-p hive1d,hive7d,hiveunlim,queen,ckptqueen,ckpthp1d,ckpthp7d,ckpthp31d,ckptdell1d,ckptdell7d,ckptdell31d,preempt1d,preempt7d,preempt31d \
                                    --wrap "python consensus_15.py \"$table_in1\" \"$table_out1\" \"$TSS\""
                            
                            # Run jobs to compare consensus mutation profile against cutoff criteria sets and determine "true" mutations and ambigous positions
                            elif [ $1 -eq 2 ]; then
                                if [ -e "$table_out1" ]; then
                                    echo 'Input consensus table is found'
									
									if [ ! -d "$table_out2_dir" ]; then mkdir "$table_out2_dir"; assert_; fi
									sleep 2
                                    sbatch -o "$table_out1.cutoffs-update.2.out" -e "$table_out1.cutoffs-update.2.err" --mem=128000 \
										-p hive1d,hive7d,hiveunlim,queen \
										--wrap "python consensus_cutoffs_15.py \"$table_out1\" \"$table_out2\" \"$size_r1\" 1"
									
									if [ ! -d "$table_out3_dir" ]; then mkdir "$table_out3_dir"; assert_; fi
									sleep 2
                                    sbatch -o "$table_out1.2.out" -e "$table_out1.2.err" --mem=128000 \
										-p hive1d,hive7d,hiveunlim,queen \
										--wrap "python consensus_cutoffs_15.py \"$table_out1\" \"$table_out3\" \"$size_r1\" 0"
                                fi
                            fi
                        fi
                    else
                        echo "Error: some of the required parameters are empty"
                        exit
                    fi
                else
                    echo "Error: not all I/O data exists"
                    exit
                fi
            done
            echo '----------------'
        done
    done
done

conda deactivate
