in1=design/samples_table.txt
in2=design/factors_table.txt
out1=design/samples_table.sh

python mapper_write-SE-or-PE-list_11.py --input_table $in1 --output_script $out1 --factors_table $in2 --SE
