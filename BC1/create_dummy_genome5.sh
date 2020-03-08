
params_1="$PWD/design/params_1.sh"
params_2="$PWD/design/samples_table.sh"

#####################################

function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }

. $params_1
assert_
. $params_2
assert_

if [ ! -d "$params_dir_out_1/mapping" ]; then
    exit "error: not all parameters/files exit"
fi

for i in ${!title[@]}; do
	
    title1="${title[i]}"_"idx$i"; assert_
    refs1=$(echo "${sort_ref[i]}"";""${sort_refs[i]}" | tr ";" "\n")
    
    k=0
    for r1 in $refs1; do
		let k=$k+1
		echo $k
		r2="$r1"
        if [ $k -eq 1 ]; then r2="$r2.others"; fi # first item in loop refers to "others"
        bam1="$params_dir_out_1/mapping/$title1.$r2.bwa.sorted.bam"
        ref1="$params_dir_reference"/"$r1".fa
        header="$params_dir_out_1/mapping/$title1.$r2.bwa.sam.header"
        echo "ref1 = $ref1"
        echo "header = $header"
        echo "bam = $bam1"
        if [ -f  "$header" ] && [ -f "$ref1" ] && [ -f "$bam1" ] && [ -f "$ref1.indices.OK" ]; then
            echo 'OK'
            if [ $1 -eq 1 ]; then
                if [ ! -f $header.barcode ]; then
                    cat $header | perl -ne 'chomp; if($_ =~ m/SN:(\w+)/){printf("%s\n",$1)};' > $header.barcode
                    assert_
                    echo $header.barcode
                    python create_dummy_genome5.py "$ref1" "$header.barcode"
                    assert_
                else
                    echo 'error: overiding exisiting $header.barcode'
                    exit
                fi
            fi
        fi
        echo '-----------------'
    done
done

