#!/bin/sh

function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }

. $params_1
assert_

bwa mem -M -t $threads "$ref1" "$f1" > "$fout1.0.sam"
assert_

# detecting the barcode at the start of the query name, making it the referece name, removing it from the query name, filtering barcodes with Ns, filtering low MAPQ
cat "$fout1.0.sam" \
  | perl -ne 'chomp; @x=split /\t/,$_; if($x[0] =~ m/^(\w+)_(.*)/){$x[2]=$1; $x[0]=$2; if((!($x[2] =~ m/.*N.*/))and($x[4]>0)){$j=join "\t",@x; printf("%s\n",$j)}}' \
  | perl -ne 'chomp; @x=split /\t/,$_; if($x[2] =~ m/^(.*?)x(.*)$/){$bc5=$1; $bc3=$2; if(length($bc5)<1){$bc5="BC5_missing"}; if(length($bc3)<1){$bc3="BC3_missing"}; $x[2]=$bc5; push @x,"XB:Z:$bc3"; $j=join "\t",@x; printf("%s\n",$j)}' \
  > "$fout1.sam"
assert_

# prepare sam header
cut -f3 "$fout1.sam" | sort -u | perl -ne 'chomp; printf("\@SQ\tSN:%s\tLN:%s\n",$_,'$ref1_length')' > "$fout1.sam.header"
assert_

cat "$fout1.sam.header" "$fout1.sam" > "$fout1.1.sam"
assert_

samtools view -Sb  "$fout1.1.sam" > "$fout1.bam"
assert_

samtools sort -T "$fout1.sorted.temp" -o "$fout1.sorted.bam" -O bam "$fout1.bam"
assert_

samtools index "$fout1.sorted.bam"
assert_

rm "$fout1.1.sam"
assert_
rm "$fout1.bam"
assert_


