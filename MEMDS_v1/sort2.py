import csv
import sys
import re
import os
import argparse
from Bio import SeqIO

def main1(print1=True):
    def float1(x):
        x=float(x)
        if (x < 0.0) or (x > 1.0): raise argparse.ArgumentTypeError("%f not in range [0.0, 1.0]"%(x))
        return x
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
    fq0=open(p.fastq,'r')
    fq1=filter(lambda row: row!="\n",fq0)
    temp1=[]
    fastq_handles={}
    for i,r in enumerate(fq1):
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
            sort1(p,temp1,fastq_handles)
            temp1=[]
        if print1:
            if (i+1)%4000000==0: print str((i+1)/4)
        #if i>1000: break
    fq0.close()
    for k in fastq_handles.keys():
        print "closing \'%s\'"%(k)
        fastq_handles[k].close()

def sort1(p,fastq_record,fastq_handles):
    seq1 = fastq_record[1] # fastq DNA sequnce
    qual1 = fastq_record[3] # fastq quality string
    header1 = fastq_record[0] # fastq header
    
    nucl0 = p.nucl.split(";") # gene separator
    pos0  = p.pos.split(";")
    refs = p.refs.split(";")
    ref = p.ref.split(";")
    match_percentage = p.match_percentage
    if ((len(nucl0)!=len(pos0))or(len(nucl0)!=len(refs))or(len(ref)!=1)) or (len(pos0)==0): sys.exit("err1")
    sort_to_ref=ref[0]+".others"
    put_in_logfile=False
    test_sort_to_ref=[] #
    test_put_in_logfile=[] #
    for i in range(len(pos0)): # visit different genes, where for each gene is represented by one or more haplotypes
        nucl1 = nucl0[i].split("&") # haplotypes separator
        pos1  = pos0[i].split("&")
        if (len(nucl1)!=len(pos1)) or (len(pos1)==0): sys.exit("err2")
        z=0
        haplotypes_with_offsets=0
        for j in range(len(pos1)): # visit a haplotype (here haplotype means one or more alleles with know bp distance between each other)
            nucl2 = nucl1[j].split(",")
            pos2  = map(int,pos1[j].split(","))
            if (len(nucl2)!=len(pos2)) or (len(pos2)==0): sys.exit("err3")
            for offset in [0,-1,1,-2,2,-3,3]:
                zz=0
                for k in range(len(pos2)):
                    pos3=pos2[k]
                    nucl3=nucl2[k]
                    s = seq1[(pos3-1+offset):(pos3+offset)]
                    if s==nucl3: zz+=1
                #if zz==len(pos2):
                if (float(zz)/float(len(pos2))) >= float(match_percentage):
                    z+=1
                    if offset!=0: haplotypes_with_offsets+=1
                    break
        if z==len(pos1):
            #sort_to_ref = refs[i]
            #if haplotypes_with_offsets>0: put_in_logfile=True
            #break
            test_sort_to_ref.append(refs[i])
            if haplotypes_with_offsets>0: test_put_in_logfile.append(True)
            else: test_put_in_logfile.append(False)
    if len(test_sort_to_ref)==1: # only if the sequence matchs a single gene then we know how to sort
        sort_to_ref = test_sort_to_ref[0]
        put_in_logfile=test_put_in_logfile[0]
    title2="%s.%s"%(p.title,sort_to_ref)
    if fastq_handles.get(title2,"")=="":
        fastq_handles[title2] = open("%s/%s.fastq"%(p.dirOut,title2),'w')
        fastq_handles[title2+".withIndels.log"] = open("%s/%s.withIndels.log"%(p.dirOut,title2),'w')
        newline1=""
    else: newline1="\n"
    fastq_handles[title2].write("%s%s"%(newline1,"\n".join([header1,seq1,'+',qual1])))
    if put_in_logfile:
        fastq_handles[title2+".withIndels.log"].write("%s%s"%(newline1,"\n".join([header1,seq1,'+',qual1])))

main1()