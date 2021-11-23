import sys
import re
from Bio import SeqIO

f1=sys.argv[1]
outIndex1=sys.argv[2]
#seqRef0=sys.argv[3]
size_f = int(sys.argv[3])
size_r = int(sys.argv[4])
read_seq = sys.argv[5].split(",")
read_pos = sys.argv[6].split(",")
read_action = sys.argv[7].split(",")
if (len(read_seq) != len(read_pos)): sys.exit("err read_pos size")
if (len(read_seq) != len(read_action)): sys.exit("err read_action size")
print(read_seq)
print(read_pos)
print(read_action)

'''
 # fasq record:

 @M00654:137:000000000-AURH5:1:1101:15842:1877 1:N:0:1
 TTACGTTATAGCTATGTGTTCACTAGCAACCTCAAACAGACACCATGGTGCATCTGACTCCAGGAGAAGTCTGCCGTTACTGCCCTGTGGTGCAAGGTGAACGTG
 +
 >A?>3BFBAB5DGBF5GGGGGGHHHHHHHHHHHCF4FHHHHHHGFGBBGBEGHHHGFH5FHHHHHA1F3F5FHGHGGGGFH5EFFCCGEGGEGGHHEDGHHGFHF
 
'''

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

def barcodeQuality_phred33(barcode1,barcode_qual,quality_cutoff=20):
    # Illumina 1.8+ Phred+33
    if len(barcode1) != len(barcode_qual): sys.exit('err1') 
    phred33 = [(ord(x)-33) for x in barcode_qual]
    for i in range(len(phred33)):
       if phred33[i] < quality_cutoff:
        z = list(barcode1)
        z[i] = 'N'
        barcode1 = "".join(z)
    return barcode1

# seqRef1, fivePrimeWin=7,threePrimeWin=7,fivePrimeWinMismatch=1,threePrimeWinMismatch=1
def trimming(fastq_record,out1,out2,out1b,out2b,size_f,size_r,log1,total1,read_seq,read_pos,read_action):
    
    seq1 = fastq_record[1] # fastq DNA sequnce
    qual1 = fastq_record[3] # fastq quality string
    header1 = fastq_record[0] # fastq header
    
    excluded=[]
    for iii in range(len(read_seq)):
        read_seq1=read_seq[iii].split("|")
        read_pos1=map(int,read_pos[iii].split("|"))
        include1=False
        for s in range(len(read_seq1)):
            if read_pos1[s] >= 0:  x = seq1[(read_pos1[s]-1):(read_pos1[s]-1+len(read_seq1[s]))]
            else:                 x = seq1[(len(seq1)+read_pos1[s]-len(read_seq1[s])+1):(len(seq1)+read_pos1[s]+1)]
            if (x == read_seq1[s]) and (read_action[iii] == 'include'):
                include1=True
        if not(include1):
            excluded.append("at %i expected '%s' observed '%s'"%(read_pos1[s],read_seq1[s],x))
    if len(excluded)>0: log1['wrong_identifier']+=1
    
    # extract BC/identifier sequences
    
    seq3=seq1[:]
    qual3=qual1[:]
    barcodeF0 =       seq3[:size_f]
    barcodeF0_qual = qual3[:size_f]
    barcodeR0 =       seq3[(len(seq3)-size_r):]
    barcodeR0_qual = qual3[(len(seq3)-size_r):]
    
    if size_f>0:
        barcodeF = barcodeQuality_phred33(barcodeF0,barcodeF0_qual)
        seq3 =   seq3[size_f:]
        qual3 = qual3[size_f:]
    else:
        barcodeF = ""
    
    if size_r>0:
        barcodeR = barcodeQuality_phred33(barcodeR0,barcodeR0_qual)
        len3 = len(seq3)
        seq3 =   seq3[:(len3-size_r)]
        qual3 = qual3[:(len3-size_r)]
    else:
        barcodeR = ""
    
    header2 = "%s%sx%s_%s"%(header1[:1],barcodeF,barcodeR,header1[1:]) # new fastq header = "@<5' barcode>x<3' barcode>_<original header>"
    total1+=1
    if total1==1: newline1=""
    else: newline1="\n"
    
    if len(excluded) == 0:
        out1.write("%s%s"%(newline1,"\n".join([header2,seq3,'+',qual3])))
        out2.write("%s\t%s\n"%(barcodeF0,barcodeR0))
    else:
        out1b.write("%s%s"%(newline1,"\n".join([header2,seq3,'+',qual3])))
        out2b.write("%s\t%s\n"%(barcodeF0,barcodeR0))
    
    return(log1,total1)

for i,r in enumerate(h1):
    r=r.rstrip()
    j=(i+1)%4
    if j==1: # fastq header
        if len(temp1)!=0: sys.exit('err11')
        if r[0:1]!='@': sys.exit('err1')
        temp1.append(r)
    elif j==2: # fastq DNA sequence
        if len(temp1)!=1: sys.exit('err22')
        temp1.append(r)
    elif j==3: # fastq "+" tag
        if len(temp1)!=2: sys.exit('err33')
        if r[0:1]!='+': sys.exit('err3')
        temp1.append(r)
    elif j==0: # fastq quality squence
        if len(temp1)!=3: sys.exit('err44')
        temp1.append(r)
        log1,total1 = trimming(temp1,out1,out2,out1b,out2b,size_f,size_r,log1,total1,read_seq,read_pos,read_action)
        temp1=[]
    if print1:
        if (i+1)%4000000==0: print str((i+1)/4)

for k in sorted(log1.keys()):
    out3.write("%s = %i / %i (%f) (%s)\n"%(k,log1[k],total1,100*(float(log1[k])/float(total1)),log1_description[k]))
out3.write("5' info expected size = %i\n"%(size_f))
out3.write("3' info expected size = %i\n"%(size_r))

h1.close()
out1.close()
out2.close()
out3.close()
