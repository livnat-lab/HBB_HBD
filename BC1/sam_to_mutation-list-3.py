import sys
import re
import os
import glob
import argparse
import pysam
import collections
from collections import Counter
import itertools
from itertools import groupby
from Bio import SeqIO
import copy
# classes whose files should be stored in the same directory of this script:
from tb2 import tb
from tb2 import tables

parser = argparse.ArgumentParser(description='parser')
parser.add_argument('--sam_files_path', action="store", dest="path1")
parser.add_argument('--out_dir', action="store", dest="outDir1")
parser.add_argument('--ref_fasta', action="store", dest="seq_target_file")
parser.add_argument('--test1', action='store_true',dest="test", default=False)
parser.add_argument('--limit_start', action="store", dest="limit_start", default='-1')
parser.add_argument('--limit_end',   action="store", dest="limit_end"  , default="%i"%(pow(10,10)))
parser.add_argument('--queryAlnMustStartAt0', action='store',dest="queryAlnMustStartAt0",type=int)
parser.add_argument('--include_BX_tag', action='store',dest='include_BX_tag',default="")
parser.add_argument('--BC3_min_cutoff', action='store',dest="BC3_min_cutoff",type=int,default=2)
parser.add_argument('--checkMut_name', action='store',dest='checkMut_name',default="")
parser.add_argument('--checkMut_pos', action='store',dest='checkMut_pos',default="")
parser.add_argument('--checkMut_action', action='store',dest='checkMut_action',default="")
p=parser.parse_args()
p.limit_start = map(int,p.limit_start.split(","))
p.limit_end =   map(int,p.limit_end.split(","))
if (len(p.limit_start)!=len(p.limit_end)) or len(p.limit_start)<1: sys.exit('error: wrong input limits')
if p.queryAlnMustStartAt0 <=0: p.queryAlnMustStartAt0=False
else: p.queryAlnMustStartAt0=True
if(len(p.checkMut_name)>0):
    p.checkMut_name = p.checkMut_name.split(",")
    p.checkMut_pos = p.checkMut_pos.split(",")
    p.checkMut_action = p.checkMut_action.split(",")
    if(len(p.checkMut_name)!=len(p.checkMut_pos)): sys.exit("len(p.checkMut_name)!=len(p.checkMut_pos)")
    if(len(p.checkMut_name)!=len(p.checkMut_action)): sys.exit("len(p.checkMut_name)!=len(p.checkMut_action)")
    
# read reference sequence in fasta, should fit the sam file reference:
seqRef0_ = open(p.seq_target_file, "rU")
seqRef1 = ''
for i,r1 in enumerate(SeqIO.parse(seqRef0_, "fasta")): seqRef1=str(r1.seq); break
seqRef0_.close()
if seqRef1 == '': sys.exit("reference sequence is missing in file \"%\""%(seqRef0))
p.seq_target = seqRef1

# # functions ###############

def cigar_to_aligedSequence(cigar_string,pos,seq,isNotReference):
    if isNotReference:
        mut='D'  # deletion in read (add gap to read)
        clip='H' # hard clipping (add gap to read, usually after/before the read end/start)
    else:
        mut='I'  # insertion in read (add gap to reference)
        clip='S' # soft clipping (add gap to reference, usually after/before the reference end/start)
    cigar = re.findall(r'\d+[MIDNSHPX=]',cigar_string)
    aln = seq
    for c in range(len(cigar)):
        m1 = re.match(r'^(\d+)([MIDNSHPX=])$',cigar[c])
        if m1:
            type1   = m1.group(2)
            length1 = int(m1.group(1))
            if type1 == clip:
                aln = aln[0:(pos)] + '-'*length1 + aln[(pos):len(aln)]  # 'x'
            elif type1 == mut :
                aln = aln[0:(pos)] + '-'*length1 + aln[(pos):len(aln)] 
            pos = pos + length1
        else: sys.exit("eee")
    return(aln)

def indentify_mutations(aln_target,aln_read,qual,limit_start,limit_end):
    start_p=-1
    start_s=-1
    t = list(aln_target)
    q = list(aln_read)
    q = q + list('-'*(len(t)-min(len(q),len(t))))  # 'x'
    qual=list(qual)
    acc_p = start_p
    t_index=[] # indices of bases in the target seq
    for x in range(len(t)):
        if((t[x]!='-')and(t[x]!='x')): acc_p=acc_p+1
        t_index.append(acc_p)
    acc_s = start_s
    q_index=[] # indices of bases in the query seq
    for x in range(len(q)):
        if((q[x]!='-')and(q[x]!='x')): acc_s=acc_s+1
        q_index.append(acc_s)
    
    # mutation calls:
    tq = []; tq_idx = []; t_idx = []
    for x in range(len(t)):
        if((t[x]!=q[x])and(t[x]!='x')and(q[x]!='x')):
            tq.append("%i%s%s"%(t_index[x]+1,t[x],q[x])) # mutation info
            tq_idx.append(q_index[x]) # mutation position in quey
            t_idx.append(t_index[x]) # mutation position in target
        else:
            tq.append(None)
            tq_idx.append(None)
            t_idx.append(None)
    sel1 = [ False for i in range(len(tq))]
    for k in range(len(limit_start)):
        sel1 = [((not(tq[i] is None) and ((limit_start[k]-1) <= t_idx[i]) and ((limit_end[k]-1) >= t_idx[i])) or (sel1[i])) for i in range(len(tq))]
    mut=[]
    mutidx=[] # index of mutation in the query sequece
    if(sum(sel1)>0):
        mut =    [ tq[i]     for i in range(len(tq)) if sel1[i] ]
        mutidx = [ tq_idx[i] for i in range(len(tq)) if sel1[i] ]
    else:
        mut    = ['WT']
        mutidx = [None]
    
    # collapse consecutive read's insertions positions
    if(len(mut)>1):
        for k in range(len(mut)):
            if(k>0):
                if(mutidx[k-1]+1 == mutidx[k]):
                    if(re.match(r'^\d+-\w+',mut[k])):
                        if(re.match(r'^\d+-\w+',mut[k-1])):
                            mut[k] =   mut[k-1] + re.sub(r'^\d+-','',mut[k])
                            mut[k-1] = None
        mut =    [mut[i]    for i in range(len(mut)) if not(mut[i] is None)]
        mutidx = [mutidx[i] for i in range(len(mut)) if not(mut[i] is None)]
    
    qaul_mut=[]
    for z in range(len(mut)):
        #if((mut[z] != 'WT') and not(re.match(r'.*-|x.*',mut[z]))):
        if((mut[z] != 'WT') and not(re.match(r'x.*',mut[z]))):  #***#
            if re.match(r'^\d+-\w+$',mut[z]): # insertion  #***#
                asciiChar = qual[mutidx[z]]
                phred  = ord(asciiChar)-33
                qaul_mut.append(phred)
            elif re.match(r'^\d+\w+-$',mut[z]): # delition  #***#
                qaul_mut.append('NA')
            elif re.match(r'^\d+\w\w$',mut[z]): # mismatch
                asciiChar = qual[mutidx[z]]
                phred  = ord(asciiChar)-33
                qaul_mut.append(phred)
            else:
                sys.exit('err: unexpected mutation name "%s\"'%(mut[z]))
        else:
            qaul_mut.append('NA')
    
    qaul_mut1=[]
    all_mut1=[]
    for z in range(len(mut)):
        if not(re.match(r'.*x.*',mut[z])):
            all_mut1.append(mut[z])
            qaul_mut1.append(qaul_mut[z])
    
    return(all_mut1,qaul_mut1)

def summarize_mutations(sam1,seq_target,limit_start,limit_end):
    c={'RNAME':0,'POS':1,'CIGAR':2,'SEQ':3,'QUAL':4,'BC3':5,'READ_NAME':6}
    alns = {}
    mutations=[]
    for s in range(len(sam1)):
        bc3_check=sam1[s][c['BC3']] # @@@@@@@@@@@@
        seq_read = sam1[s][c['SEQ']]
        qual     = sam1[s][c['QUAL']]
        bc5 = sam1[s][c['RNAME']]
        bc3 = sam1[s][c['BC3']]
        read_name= sam1[s][c['READ_NAME']]
        
        if alns.get(seq_read,"")=="":
            cigar_str = sam1[s][c['CIGAR']]
            pos       = sam1[s][c['POS']]
            aln_target = cigar_to_aligedSequence(cigar_str,pos,seq_target,False) # reference sequence
            aln_read   = cigar_to_aligedSequence(cigar_str,0,seq_read,True)      # read sequence
            alns[seq_read] = [aln_target,aln_read]
        else:
            aln_target = alns[seq_read][0]
            aln_read   = alns[seq_read][1]
        
        all_mut0,qaul_mut0 = indentify_mutations(aln_target,aln_read,qual,limit_start,limit_end)
        
        if len(all_mut0)!=len(qaul_mut0): sys.exit("err")
        for i in range(len(all_mut0)):
            m1=re.match('^(\d+)(-|\w+?)(\w+|-)$',all_mut0[i])
            from1='NA'
            to1='NA'
            pos1='NA'
            if m1:
                from1=m1.group(2)
                to1=m1.group(3)
                pos1=m1.group(1)
            h=collections.OrderedDict([('read_name',read_name),('bc5',bc5),('bc5_count',len(sam1)),('bc3',bc3),('quality',qaul_mut0[i]),('pos',pos1),('from',from1),('to',to1),('allele',all_mut0[i])])
            
            mutations.append(h)
    
    return(mutations)

### main ########################################################################

files1 = glob.glob(p.path1)
# sam columns: 'QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','MRNM','MPOS','ISIZE','SEQ','QUAL','TAGs']

for i in range(len(files1)):
    m1=re.match('.*\\/(.*)',files1[i])
    if m1:
        title1=m1.group(1)
        if p.include_BX_tag != "": title1=title1+".tag-"+p.include_BX_tag
        samfile = pysam.AlignmentFile(files1[i], "rb")
        rname_prev=''
        rnames={}
        visited_ref=0
        visited_read=0
        aln=[]
        outfile1 = "%s/%s.mutations.txt"%(p.outDir1,title1)
        outfile2 = "%s/%s.mutations.log.txt"%(p.outDir1,title1)
        out1=open(outfile1,'w')
        out1_visited=False
        out2=open(outfile2,'w')
        #out1.write(("\t".join(out1_cols))+"\n")
        log={'unmapped':0,'ref_start_ne0':0,'query_start_ne0':0,'start_ne0':0,'total':0}
        for read in samfile.fetch(until_eof=True):
            rname=samfile.getrname(read.reference_id)
            summarize_mutations1=[]
            if not(read.is_unmapped) and (read.reference_start==0) and ((read.query_alignment_start==0) or (p.queryAlnMustStartAt0==False)):
                if read.has_tag('XB'):
                    tag1=read.get_tag('XB')
                else: tag1='NA'
                if (p.include_BX_tag == "") or ((tag1 != 'NA') and (p.include_BX_tag != "") and (p.include_BX_tag == tag1)):
                    if rname != rname_prev and rname_prev != '':
                        if rnames.get(rname,"")!="": sys.exit("expected the target name to be sorted\n")
                        visited_ref+=1
                        if visited_ref % 1000 == 0:
                            sys.stdout.write('\r' + "visited ref=%i read=%i alnSize=%i title=%s      "%(visited_ref,visited_read,len(aln),title1))
                            sys.stdout.flush()
                        summarize_mutations1 = summarize_mutations(aln,p.seq_target,p.limit_start,p.limit_end)
                        for j in range(len(summarize_mutations1)):
                            row1=[summarize_mutations1[j][k] for k in summarize_mutations1[j].keys()]
                            if not(out1_visited):
                                out1_visited=True
                                out1.write("%s\n"%("\t".join(summarize_mutations1[j].keys())))
                            # rname_prev
                            out1.write("%s\n"%("\t".join(map(str,row1))))
                        aln=[]
                    
                    aln.append([rname,read.reference_start,read.cigarstring,read.query_sequence,read.qual,tag1,read.query_name]) # reference_start = 0-based leftmost coordinate; sam file fields: 'RNAME','POS','CIGAR','SEQ','QUAL'
                    rnames[rname]=True
                    rname_prev = rname
                    visited_read +=1
            
            if read.is_unmapped: log['unmapped'] +=1
            else:
                if read.reference_start !=0: log['ref_start_ne0'] +=1 # reference leftmost coordinate
                if read.query_alignment_start !=0: log['query_start_ne0'] +=1  # start index of the aligned query portion
                if (read.reference_start !=0) or (read.query_alignment_start !=0): log['start_ne0'] +=1
            log['total'] += 1
            
            if len(summarize_mutations1)>0:
                pass
            
            #if visited_read>30: break
        
        if len(aln)>0:
            summarize_mutations1 = summarize_mutations(aln,p.seq_target,p.limit_start,p.limit_end)
            for j in range(len(summarize_mutations1)):
                row1=[summarize_mutations1[j][k] for k in summarize_mutations1[j].keys()]
                if not(out1_visited):
                    out1_visited=True
                    out1.write("%s\n"%("\t".join(summarize_mutations1[j].keys())))
                out1.write("%s\n"%("\t".join(map(str,row1))))
            aln=[]
        
        for k in log.keys():
            out2.write("reads %s\t%i\n"%(k,log[k]))
        out2.write("reads used\t%i\n"%(visited_read))
        out2.write("queryAlnMustStartAt0\t%s\n"%(p.queryAlnMustStartAt0))
        out2.write("bam\t%s\n"%(files1[i]))
        out2.write("fasta\t%s\n"%(p.seq_target_file))
        
        samfile.close()
        out1.close()
        out2.close()
    print "\nok"


