import sys
import re
from Bio import SeqIO

# A script to trim barcodes from the analyzed reads and filter reads having wrong barcode sequences

'''
# Fastq record example:

@M00654:137:000000000-AURH5:1:1101:15842:1877 1:N:0:1
TTACGTTATAGCTATGTGTTCACTAGCAACCTCAAACAGACACCATGGTGCATCTGACTCCAGGAGAAGTCTGCCGTTACTGCCCTGTGGTGCAAGGTGAACGTG
+
>A?>3BFBAB5DGBF5GGGGGGHHHHHHHHHHHCF4FHHHHHHGFGBBGBEGHHHGFH5FHHHHHA1F3F5FHGHGGGGFH5EFFCCGEGGEGGHHEDGHHGFHF
 
'''

#### Functions ####

# Check sequencing quality of the barcodes, return 'N' (unknown base) at low quality positions
def barcodeQuality_phred33(barcode1,barcode_qual,quality_cutoff=20):
    # Illumina 1.8+ Phred+33
    if len(barcode1) != len(barcode_qual): sys.exit('Mismatch between length of the barcode and length of the quality score array') 
    phred33 = [(ord(x)-33) for x in barcode_qual]
    for i in range(len(phred33)):
       if phred33[i] < quality_cutoff:
        z = list(barcode1)
        z[i] = 'N'
        barcode1 = "".join(z)
    return barcode1

# Barcode filtering and trimming function
def trimming(fastq_record,out1,out2,out1b,out2b,size_f,size_r,log1,total1,read_seq,read_pos,read_action):
    
    seq1 = fastq_record[1] # Read sequence
    qual1 = fastq_record[3] # Read sequencing quality string
    header1 = fastq_record[0] # Read fastq header
    
    excluded=[]
    for iii in range(len(read_seq)):
        read_seq1=read_seq[iii].split("|") # "|" is used to separate multiple variants of "correct" barcode signatures at same position 
        read_pos1=map(int,read_pos[iii].split("|"))
        include1=False
        
        # Check that barcode sequence at identified positions in the read matches expected signature
        for s in range(len(read_seq1)):
            if read_pos1[s] >= 0:  x = seq1[(read_pos1[s]-1):(read_pos1[s]-1+len(read_seq1[s]))] # Check 5' barcode
            else:                  x = seq1[(len(seq1)+read_pos1[s]-len(read_seq1[s])+1):(len(seq1)+read_pos1[s]+1)] # Check 3' barcode
            if (x == read_seq1[s]) and (read_action[iii] == 'include'):
                include1=True
        if not(include1):
            excluded.append("at %i expected '%s' observed '%s'"%(read_pos1[s],read_seq1[s],x))
    if len(excluded)>0: log1['wrong_identifier']+=1
    
    # Separate barcode and gene sequences
    seq3=seq1[:]
    qual3=qual1[:]
    
    # Extract barcode sequences and their associated quality
    barcodeF0 =       seq3[:size_f]
    barcodeF0_qual = qual3[:size_f]
    barcodeR0 =       seq3[(len(seq3)-size_r):]
    barcodeR0_qual = qual3[(len(seq3)-size_r):]
    
    # Substitute low quality positions in the 5' barcode by 'N'; separate its sequence from gene sequence in the read    
    if size_f>0:
        barcodeF = barcodeQuality_phred33(barcodeF0,barcodeF0_qual)
        seq3 =   seq3[size_f:]
        qual3 = qual3[size_f:]
    else:
        barcodeF = ""
    
    # Substitute low quality positions in the 3' barcode by 'N'; separate its sequence from gene sequence in the read 
    if size_r>0:
        barcodeR = barcodeQuality_phred33(barcodeR0,barcodeR0_qual)
        len3 = len(seq3)
        seq3 =   seq3[:(len3-size_r)]
        qual3 = qual3[:(len3-size_r)]
    else:
        barcodeR = ""
    
    # Create new fastq header = "@<BC5>x<BC3>_<original header>"
    header2 = "%s%sx%s_%s"%(header1[:1],barcodeF,barcodeR,header1[1:])
    
    total1+=1
    if total1==1: newline1=""
    else: newline1="\n"
    
    # Write output, correct barcode signature; out2 - list of extracted barcodes
    if len(excluded) == 0:
        out1.write("%s%s"%(newline1,"\n".join([header2,seq3,'+',qual3])))
        out2.write("%s\t%s\n"%(barcodeF0,barcodeR0))
    
    # Write output, wrong barcode signature; out2 - list of extracted barcodes
    else:
        out1b.write("%s%s"%(newline1,"\n".join([header2,seq3,'+',qual3])))
        out2b.write("%s\t%s\n"%(barcodeF0,barcodeR0))
    
    return(log1,total1)

#### Main ####
# Gather trimming parameters
f1=sys.argv[1]
outIndex1=sys.argv[2]
size_f = int(sys.argv[3])
size_r = int(sys.argv[4])
read_seq = sys.argv[5].split(",") # argv[5], argv[6], argv[7] - "," separates 5' and 3' barcode info
read_pos = sys.argv[6].split(",")
read_action = sys.argv[7].split(",")
if (len(read_seq) != len(read_pos)): sys.exit("err read_pos size")
if (len(read_seq) != len(read_action)): sys.exit("err read_action size")
print(read_seq)
print(read_pos)
print(read_action)

# Read input data and create output files 
h1=open(f1,'r')

out1=open("%s.fastq"%(outIndex1),'w')
out2=open("%s.barcodes"%(outIndex1),'w')
out2.write("seq5\tseq3\tmatch5\tmatch3\n")
out1b=open("%s.wrongId.fastq"%(outIndex1),'w')
out2b=open("%s.wrongId.barcodes"%(outIndex1),'w')
out2b.write("seq5\tseq3\tmatch5\tmatch3\n")

out3=open("%s.log"%(outIndex1),'w')
temp1=[]
print1=False
log1={'wrong_identifier':0}
log1_description={'wrong_identifier':"wrong 5' or 3' identifier"}
total1=0

# Iterate over input data and submit read entries for trimming
for i,r in enumerate(h1):
    r=r.rstrip()
    j=(i+1)%4
    if j==1: # fastq header
        if len(temp1)!=0: sys.exit('Data array should be empty when reading in new read')
        if r[0:1]!='@': sys.exit('Wrong fastq header format, expected \"@\" at the start')
        temp1.append(r)
    elif j==2: # fastq DNA sequence
        if len(temp1)!=1: sys.exit('Data array should contain only header before reading read sequence')
        temp1.append(r)
    elif j==3: # fastq "+" tag
        if len(temp1)!=2: sys.exit('Data array should contain only header and read sequence before fastq \"+\" tag')
        if r[0:1]!='+': sys.exit('Unexpected tag; \"+\" tag is expected as a third line of fastq read entry')
        temp1.append(r)
    elif j==0: # fastq quality squence
        if len(temp1)!=3: sys.exit('Data array should contain only header, read sequence and \"+\" tag before reading read sequencing quality')
        temp1.append(r)
        log1,total1 = trimming(temp1,out1,out2,out1b,out2b,size_f,size_r,log1,total1,read_seq,read_pos,read_action)
        temp1=[]
    if print1:
        if (i+1)%4000000==0: print str((i+1)/4)

# Output log data on percentage of removed reads
for k in sorted(log1.keys()):
    out3.write("%s = %i / %i (%f) (%s)\n"%(k,log1[k],total1,100*(float(log1[k])/float(total1)),log1_description[k]))
out3.write("5' info expected size = %i\n"%(size_f))
out3.write("3' info expected size = %i\n"%(size_r))

h1.close()
out1.close()
out2.close()
out3.close()
