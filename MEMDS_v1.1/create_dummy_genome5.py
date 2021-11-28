from Bio import SeqIO
import glob
import sys

# This script creates a reference file needed for manual inspection of read alignments against the reference sequence by IGV
# Reference file contains an entry for each read family, with entry name being the unique BC5 identifying the family and entry sequence being the reference sequence used in the analysis. 

f_ref=str(sys.argv[1])
path_barcodes=str(sys.argv[2])

barcode_files = glob.glob(path_barcodes)
ref1 = open(f_ref,'rU') # Open reference sequence file and extract reference entry from it
for record1 in SeqIO.parse(ref1, "fasta"): pass
print(record1.id)
for i in range(0,len(barcode_files)):
    print(barcode_files[i])
    bc = open(barcode_files[i],'r')
    out1 = open("%s.fasta"%(barcode_files[i]),'w') # Open an empty "fasta" file to write reference sequences with updated names
    for j,r in enumerate(bc):
        if j>0:
            barcode = r.rstrip().split("\t")[0]
            record1.id = barcode # Update name of the reference sequence as the unique BC5 sequence
            SeqIO.write([record1], out1, "fasta") # Output reference sequence with the updated name
    out1.close()
    bc.close()
ref1.close()
