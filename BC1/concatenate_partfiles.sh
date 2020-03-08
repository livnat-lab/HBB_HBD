in1=more_scripts/samples_table_0.txt
out1=more_scripts/samples_table_0.sh

python mapper_write-SE-or-PE-list_11.py --input_table "$in1" --output_script "$out1" --concatenate_partfiles --output_dir "../raw_concatenated"

