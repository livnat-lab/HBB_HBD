
params_1="$PWD/design/params_1.sh"
params_2="$PWD/design/samples_table.sh"

run_on_slurm=1

#######################################

function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }

. $params_1
assert_
. $params_2
assert_

if [ ! -d "$params_dir_out_1/mapping" ]; then
    exit "error: not all parameters/files exit"
fi

for BC3_min_cutoff in 1; do

	outdir="$params_dir_out_1/tables.BC3cutoff$BC3_min_cutoff"
	if [ ! -d "$outdir" ]; then mkdir "$outdir"; assert_; fi
    
	for i in ${!title[@]}; do
		title1="${title[i]}"_"idx$i"
		refs1=$(echo "${sort_ref[i]}"";""${sort_refs[i]}" | tr ";" "\n")
		limit_start="${limit_starts[$i]}"
		limit_end="${limit_ends[$i]}"
		seq_pos="${seq_pos[i]}"
		seq_mut="${seq_mut[i]}"
		seq_action="${seq_action[i]}"
        echo "limits = $limit_start - $limit_end"
        echo "BC3_min_cutoff = $BC3_min_cutoff"
        k=0
        for r1 in $refs1; do
            echo "r1 = $r1"
            let k=$k+1
            echo "k = "$k
            r2="$r1"
            if [ $k -eq 1 ]; then r2="$r2.others"; fi # first item in loop refers to "others"
            bam="$params_dir_out_1/mapping/$title1.$r2.bwa.sorted.bam"
            ref1="$params_dir_reference"/"$r1".fa
            echo "bam = $bam"
            echo "ref1 = $ref1"
            if [ -f "$bam" ] && [ -f $ref1 ]; then
                echo 'files found'
                if [ ! -z "$limit_start" ] && [ ! -z "$limit_end" ] && [ ! -z "$seq_pos" ] && [ ! -z "$seq_mut" ] && [ ! -z "$seq_action" ] && [ ! -z "$title1" ] && [ ! -z "$bam" ] && [ ! -z "$ref1" ] && [ ! -z "$BC3_min_cutoff" ]; then
                    echo 'parameters found'
                for sam_BX_tag in '' 'BC3_missing'; do # to include only specific sam BX tags
                    echo "sam_BX_tag = \"$sam_BX_tag\""
					echo "title = $title1"
					echo "outlog = ""$outdir"/"$title1.$r2.bwa.sorted.bam.$sam_BX_tag.out"
                    #if [ "$bam" == "/data/home/livnat/dmelamed/daniel/run11_1and2/out1/mapping/contLL035_idx0.HBB.bwa.sorted.bam" ] && [ "$sam_BX_tag" == "" ]; then
                    #    echo xxxxxxxxx
                    if [ $1 -eq 1 ]; then
						sleep 2
                        
						sbatch --mem=20000  -N 1 -n 1 --ntasks-per-node=1 \
                                -p hive1d,hive7d,hiveunlim,queen,ckptqueen,ckpthp1d,ckpthp7d,ckpthp31d,ckptdell1d,ckptdell7d,ckptdell31d,preempt1d,preempt7d,preempt31d \
                                -o "$outdir"/"$title1.$r2.bwa.sorted.bam.$sam_BX_tag.out" -e "$outdir"/"$title1.$r2.bwa.sorted.bam.$sam_BX_tag.err" \
                                --wrap "python sam_to_mutation-table_5.0.py \
                                            --sam_files_path=\"$bam\" \
                                            --out_dir=\"$outdir\" \
                                            --ref_fasta=\"$ref1\" \
                                            --limit_start=\"$limit_start\" \
                                            --limit_end=\"$limit_end\" \
                                            --queryAlnMustStartAt0=1 \
                                            --include_BX_tag=\"$sam_BX_tag\" \
                                            --BC3_min_cutoff=\"$BC3_min_cutoff\" \
                                            --checkMut_name=\"$seq_mut\" \
                                            --checkMut_pos=\"$seq_pos\" \
                                            --checkMut_action=\"$seq_action\""
                    fi
                    #fi
                done
                else
                    echo 'error: not all parameters OK'
                    exit
                fi
            else
                echo 'error: not all input files exist'
                exit
            fi
        done
        echo '---------'
	done
done

