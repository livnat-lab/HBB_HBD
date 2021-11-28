import csv
import sys
import re
import os
import argparse
from Bio import SeqIO

# Main function for sorting reads based on their gene of origin
def main1(print1=True):
    # Return float(x) if x in defined range
    def float1(x):
        x=float(x)
        if (x < 0.0) or (x > 1.0): raise argparse.ArgumentTypeError("%f not in range [0.0, 1.0]"%(x))
        return x
    
    # Gather input arguments    
    parser = argparse.ArgumentParser(description='ort',formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--fastq', action="store", dest="fastq",required=True)
    parser.add_argument('--pos', action="store", dest="pos",required=True)
    parser.add_argument('--nucl', action="store", dest="nucl",required=True)
    parser.add_argument('--refs', action="store", dest="refs",required=True)
    parser.add_argument('--ref', action="store", dest="ref",required=True)
    parser.add_argument('--dirOut', action="store", dest="dirOut",required=True)
    parser.add_argument('--title', action="store", dest="title",required=True)
    parser.add_argument('--match_percentage',action="store", dest="match_percentage",required=True,type=float1)
    p=parser.parse_args()
    print p
    
    # Read the input
    fq0=open(p.fastq,'r')
    fq1=filter(lambda row: row!="\n",fq0)
    temp1=[]
    fastq_handles={}
    
    
    for i,r in enumerate(fq1):
        r=r.rstrip()
        j=(i+1)%4
        if j==1: # Fastq header
            if len(temp1)!=0: sys.exit('Data array should be empty when reading in new read')
            if r[0:1]!='@': sys.exit('Wrong fastq header format, expected \"@\" at the start')
            temp1.append(r)
        elif j==2: # Fastq DNA sequence
            if len(temp1)!=1: sys.exit('Data array should contain only header before reading read sequence')
            temp1.append(r)
        elif j==3: # Fastq "+" tag
            if len(temp1)!=2: sys.exit('Data array should contain only header and read sequence before fastq \"+\" tag')
            if r[0:1]!='+': sys.exit('Unexpected tag; \"+\" tag is expected as a third line of fastq read entry')
            temp1.append(r)
        elif j==0: # fastq quality squence
            if len(temp1)!=3: sys.exit('Data array should contain only header, read sequence and \"+\" tag before reading read sequencing quality')
            temp1.append(r)
            sort1(p,temp1,fastq_handles)
            temp1=[]
        if print1:
            if (i+1)%4000000==0: print str((i+1)/4)
        #if i>1000: break # Uncomment this line to perform test run on the first 1000 entries of the ".fastq" file
        
    fq0.close()
    for k in fastq_handles.keys():
        print "closing \'%s\'"%(k)
        fastq_handles[k].close()

# Sort input reads based on their gene of origin
def sort1(p,fastq_record,fastq_handles):

    seq1 = fastq_record[1] # Read sequence
    qual1 = fastq_record[3] # Read sequencing quality string
    header1 = fastq_record[0] # Read fastq header
    
    # Gather sorting parameters - position and identity of sorting nucleotides
    nucl0 = p.nucl.split(";") # ";" - gene separator, here and in all parameters below
    pos0  = p.pos.split(";")
    refs = p.refs.split(";")
    ref = p.ref.split(";")
    match_percentage = p.match_percentage
    
    if ((len(nucl0)!=len(pos0))or(len(nucl0)!=len(refs))or(len(ref)!=1)) or (len(pos0)==0): sys.exit("Some of the sorting parameters are missing")
    
    sort_to_ref=ref[0]+".others" # Reads would be sorted to "others" file if gene of origin is not identified below
    put_in_logfile=False
    test_sort_to_ref=[]
    test_put_in_logfile=[]
    
    # Perform read sorting
    # Iterate over different sorting genes to check if given read matches to any of them
    for i in range(len(pos0)): 
        nucl1 = nucl0[i].split("&") # Separate sorting nucleotides of given gene by haplotype separator
        pos1  = pos0[i].split("&")  # Separate sorting positions of each haplotype by haplotype separator
        if (len(nucl1)!=len(pos1)) or (len(pos1)==0): sys.exit("Mismatch in number of haplotypes between sorting nucleotides and sorting positions")
        z=0
        haplotypes_with_offsets=0
        
        # Visit a haplotype (here haplotype means one or more alleles with know bp distance between each other)
        for j in range(len(pos1)): 
            nucl2 = nucl1[j].split(",") # Separate individual sorting nucleotides of given haplotype
            pos2  = map(int,pos1[j].split(",")) # Separate individual sorting positions of given haplotype
            if (len(nucl2)!=len(pos2)) or (len(pos2)==0): sys.exit("Mismatch between number of sorting nucleotides and number of sorting positions")
            
            # Check presence of sorting nucleotides at sorting positions to identify gene origin
            # Offset is used to account for potential indels before the sorting region
            # Read is considered as originating from gene X if amount of matched sorting nucleotides in it is equal or larger than user defined match percent
            for offset in [0,-1,1,-2,2,-3,3]:
                zz=0
                for k in range(len(pos2)):
                    pos3=pos2[k]
                    nucl3=nucl2[k]
                    s = seq1[(pos3-1+offset):(pos3+offset)]
                    if s==nucl3: zz+=1
                
                if (float(zz)/float(len(pos2))) >= float(match_percentage):
                    z+=1
                    if offset!=0: haplotypes_with_offsets+=1
                    break
        
        # Indicate potential origin genes of the read            
        if z==len(pos1):
            test_sort_to_ref.append(refs[i])
            if haplotypes_with_offsets>0: test_put_in_logfile.append(True)
            else: test_put_in_logfile.append(False)
            
    if len(test_sort_to_ref)==1: # Only if read sequence matches a single reference gene we know where to sort it
        sort_to_ref = test_sort_to_ref[0]
        put_in_logfile=test_put_in_logfile[0]
        
    title2="%s.%s"%(p.title,sort_to_ref)
    
    ## Write sorted reads to the output ##
    # Open new output file if it doesn't exist
    if fastq_handles.get(title2,"")=="":
        fastq_handles[title2] = open("%s/%s.fastq"%(p.dirOut,title2),'w')
        fastq_handles[title2+".withIndels.log"] = open("%s/%s.withIndels.log"%(p.dirOut,title2),'w')
        newline1=""
    else: newline1="\n"
    
    # Write sorted read to the relevant output file
    fastq_handles[title2].write("%s%s"%(newline1,"\n".join([header1,seq1,'+',qual1])))
    
    # Create a separate log of reads where match position was shifted due to indels
    if put_in_logfile:
        fastq_handles[title2+".withIndels.log"].write("%s%s"%(newline1,"\n".join([header1,seq1,'+',qual1])))

#### Main ####
main1()
