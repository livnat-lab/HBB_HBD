
from Bio import SeqIO
import glob
import sys

f_ref=str(sys.argv[1])
path_barcodes=str(sys.argv[2])

barcode_files = glob.glob(path_barcodes)
ref1 = open(f_ref,'rU')
for record1 in SeqIO.parse(ref1, "fasta"): pass
print(record1.id)
for i in range(0,len(barcode_files)):
    print(barcode_files[i])
    bc = open(barcode_files[i],'r')
    out1 = open("%s.fasta"%(barcode_files[i]),'w')
    for j,r in enumerate(bc):
        if j>0:
            barcode = r.rstrip().split("\t")[0]
            record1.id = barcode
            SeqIO.write([record1], out1, "fasta")
    out1.close()
    bc.close()
ref1.close()
