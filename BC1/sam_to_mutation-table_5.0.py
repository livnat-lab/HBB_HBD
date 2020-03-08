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

def aggregate_(keys_,values_,i=0,j=1):
    groups1=[]
    if len(keys_)!=len(values_): sys.exit('aggregate_() err1')
    arr_of_arr = sorted(zip(keys_,values_), key=lambda x: x[i])
    for key1, group1 in groupby(arr_of_arr, lambda x: x[i]):
        groups1.append([key1,[item1[j] for item1 in group1]])
    return(groups1)

def indentify_mutations(aln_target,aln_read,qual,min_phred33,limit_start,limit_end):
    start_p=-1
    start_s=-1
    t = list(aln_target)
    q = list(aln_read)
    q = q + list('-'*(len(t)-min(len(q),len(t))))  # 'x'
    qual=list(qual)
    #print 'aln:'
    #print "".join(q)
    #print "".join(qual)
    #print "".join(t)
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
    
    # exclude mutations matching low quality read regions
    high_qaul_mut=[]
    for z in range(len(mut)):
        # if((mut[z] != 'WT') and not(re.match(r'.*-|x.*',mut[z]))): #***#
        if((mut[z] != 'WT') and not(re.match(r'x.*',mut[z]))):  #***#
            if re.match(r'^\d+-\w+$',mut[z]): # insertion  #***#
                asciiChar = qual[mutidx[z]]
                phred  = ord(asciiChar)-33
                if(phred >= min_phred33): high_qaul_mut.append(mut[z])
            elif re.match(r'^\d+\w+-$',mut[z]): # delition  #***#
                high_qaul_mut.append(mut[z])
            elif re.match(r'^\d+\w\w$',mut[z]): # mismatch
                asciiChar = qual[mutidx[z]]
                phred  = ord(asciiChar)-33
                if(phred >= min_phred33): high_qaul_mut.append(mut[z])
            else:
                sys.exit('err: unexpected mutation name "%s\"'%(mut[z]))
        else:
            high_qaul_mut.append(mut[z])
    
    high_qaul_mut = [high_qaul_mut[i] for i in range(len(high_qaul_mut)) if not(re.match(r'.*x.*',high_qaul_mut[i]))]
    all_mut = [mut[i] for i in range(len(mut)) if not(re.match(r'.*x.*',mut[i]))]
    
    # query sequence qualities in all positions of the target sequence
    isHighQual = []
    for x in range(len(t)):
        qual1 = qual[q_index[x]]
        if(t[x]!='x') and (t[x]!='-'): # skipping insertion quality
            if((q[x]!='-')and(q[x]!='x')):
                isHighQual.append((ord(qual1)-33)>min_phred33)
            elif(q[x]=='-'):
                isHighQual.append(True)
            else: isHighQual.append(False)
        else: isHighQual.append(None)
    isHighQual = [isHighQual[i] for i in range(len(isHighQual)) if not(isHighQual[i] is None)]
    
    return(high_qaul_mut,all_mut,isHighQual)

def summarize_mutations(sam1,seq_target,limit_start,limit_end,
                        checkMut_name,checkMut_pos,checkMut_action,
                        c={'RNAME':0,'POS':1,'CIGAR':2,'SEQ':3,'QUAL':4,'BC3':5},
                        min_phred33=28,BC3_min_cutoff=2):
    high_qaul_mut = [] # HQ mutation names
    high_qaul_mut_bc3 = []
    all_mut = [] # mutation names
    all_mut_bc3 = []
    isHighQual = [] # rows: reads, columns: quality in all aligned target positions (True = high , False = low)
    tables1=tables()
    barcode3 = []
    
    barcode3_check = {}
    for s in range(len(sam1)):
        bc3_check=sam1[s][c['BC3']]
        if barcode3_check.get(bc3_check,'')=='': barcode3_check[bc3_check]=1
        elif bc3_check != 'X': barcode3_check[bc3_check] +=1
    
    alns = {}
    for s in range(len(sam1)):
        bc3_check=sam1[s][c['BC3']]
        if barcode3_check[bc3_check] >= BC3_min_cutoff:
            seq_read = sam1[s][c['SEQ']]
            qual     = sam1[s][c['QUAL']]
            
            if alns.get(seq_read,"")=="":
                cigar_str = sam1[s][c['CIGAR']]
                pos       = sam1[s][c['POS']]
                aln_target = cigar_to_aligedSequence(cigar_str,pos,seq_target,False) # reference sequence
                aln_read   = cigar_to_aligedSequence(cigar_str,0,seq_read,True)      # read sequence
                alns[seq_read] = [aln_target,aln_read]
            else:
                aln_target = alns[seq_read][0]
                aln_read   = alns[seq_read][1]
            
            high_qaul_mut0,all_mut0,isHighQual0 = indentify_mutations(aln_target,aln_read,qual,min_phred33,limit_start,limit_end)
            
            # exclude reads with specific mutations
            all_mut0_filtered = all_mut0[:]
            for zz in range(len(checkMut_name)):
                if checkMut_action[zz] == "exclude":
                    all_mut0_filtered   = [all_mut0_filtered[z]  for z in range(len(all_mut0_filtered)) if all_mut0_filtered[z] != "%s%s"%(checkMut_pos[zz],checkMut_name[zz])]
            
            # concatenate all mutations from surviving reads
            if len(all_mut0) == len(all_mut0_filtered):
                barcode3.append(sam1[s][c['BC3']]) #@@@
                all_mut += all_mut0
                high_qaul_mut += high_qaul_mut0
                isHighQual.append(isHighQual0)
                all_mut_bc3 =       all_mut_bc3       + [barcode3[len(barcode3)-1] for z in range(len(all_mut0))]
                high_qaul_mut_bc3 = high_qaul_mut_bc3 + [barcode3[len(barcode3)-1] for z in range(len(high_qaul_mut0))]
     
    if(len(all_mut)==0):
        mut3=[]
    else:
        c1=Counter(all_mut)
        mut1b= [ [k,c1[k]] for k in c1.keys()]
        if(len(high_qaul_mut)>0): # HQ mutations were found
            c2 = Counter(high_qaul_mut)
            mut1a = [ [k,c2[k]] for k in c2.keys()]
            mut2 = tables1.merge(mut1b,mut1a,ncolx0=2,ncoly0=2) # mut2 : [[mutation_name,mutation_count,mutation_count_HQ], ...]
        else:
            mut1a = []
            mut2 = [[mut1b[i][0],mut1b[i][1],0] for i in range(len(mut1b))]
        
        isHighQual_ = map(list, zip(*isHighQual)) # transpose quality arrays
        isHighQual_sum = [sum(isHighQual_[i]) for i in range(len(isHighQual_))] # count of high quality reads along the target sequence positions
        
        for i in range(len(mut2)):
            totalCount = len(isHighQual_[0])
            mut2[i].append(totalCount) # total_count
            #if re.match(r'^\d+.*',mut2[i][0]): #**#
            if re.match(r'^\d+\w\w$',mut2[i][0]): #**#
                mut2[i].append( isHighQual_sum[int(re.match(r'^(\d+).*',mut2[i][0]).group(1))-1] ) # total_counts_high_qual
                if float(mut2[i][4]) != 0.0: mut2[i].append(float(mut2[i][2])/float(mut2[i][4])) # HQ mutation frequency
                else: mut2[i].append(0.0)
            elif re.match(r'^\d+\w+-$',mut2[i][0]) or re.match(r'^\d+-\w+$',mut2[i][0]): #**#
                mut2[i].append(float(totalCount)) # total_counts_high_qual #**#
                mut2[i].append(float(mut2[i][2])/float(totalCount)) # HQ mutation frequency #**#
            elif mut2[i][0] == 'WT':
                mut2[i].append(len(isHighQual_[0]))
                mut2[i].append(float(mut2[i][2])/len(isHighQual_[0]))
            else:
                sys.exit('unexpected mutation type: \"%s\"'%(mut2[i][0]))
                #mut2[i].append(len(isHighQual_[0])) #**#
                #mut2[i].append(1.0) # HQ mutation frequency #**#
        
        mut2 = sorted(mut2, key=lambda x: x[0])
        
        # aading 3' barcoded information:
        bc3_per_HQmuts = aggregate_(keys_=high_qaul_mut,values_=high_qaul_mut_bc3) # [[mut1,[bcX,bcY,bcX,..]], [mut2,[bcX,bcX,bcZ,..]], ...]
        bc3_summary = []
        for w in range(len(bc3_per_HQmuts)):
            bc3_mut = bc3_per_HQmuts[w][0]
            bc3_counter = Counter(bc3_per_HQmuts[w][1])
            bc3_names =  ";".join(bc3_counter.keys())
            bc3_counts =  ";".join(map(str,bc3_counter.values()))
            bc3_names_unique = map(str,set(bc3_counter.keys()))
            bc3_count_nonMissing = len(bc3_names_unique)-1 if 'BC3_missing' in bc3_names_unique else len(bc3_names_unique)
            bc3_count_missing    = 1 if 'BC3_missing' in bc3_names_unique else 0
            bc3_summary.append([bc3_mut,bc3_names,bc3_counts,bc3_count_nonMissing,bc3_count_missing])
        mut3 = tables1.merge(mut2,bc3_summary,ncolx0=6,ncoly0=5)
        for i in range(len(mut3)):
            hqTotal=float(mut3[i][4])
            mutHQcount = int(mut3[i][2])
            if mutHQcount > 0:
                mutHQCountsBc3 = map(float,mut3[i][7].split(";"))
                if hqTotal>0: mutHQFreqBc3 = map(lambda z: z/hqTotal,mutHQCountsBc3)
                else: mutHQFreqBc3 = [0.0 for z in range(len(mutHQCountsBc3))]
                mut3[i].append(";".join(map(str,mutHQFreqBc3)))
            else:
                mut3[i].append(0)
        
        # in mutation positions, find the BC3s of non-mutated reads:
        
        cols1={'Mutation':0,'Mutation_count':1,'HQ_mutation_count':2,'Total_count':3,'HQ_total_count':4,'Mutation_freq':5,
                   'barcode3_IDs':6,'barcode3_HQ_mut_counts':7,'barcode3_existing':8,'barcode3_missing':9,'barcode3_HQ_mut_freq':10}
        bc3_total_counts = Counter(barcode3)
        bc3_total_counts1 = dict(bc3_total_counts)
        bc3_total_counts_str = ";".join(map(str,bc3_total_counts.values()))
        bc3_ids_str = ";".join(bc3_total_counts.keys())
        table1=tb(cols1)
        mut3 = table1.ins(mut3,'barcode3_all',[ bc3_ids_str for w in range(len(mut3))])
        mut3 = table1.ins(mut3,'barcode3_all_counts',[ bc3_total_counts_str for w in range(len(mut3))])
        mutations_pos1=[]
        for x in high_qaul_mut:
            mtch=re.match(r"^(\d+)([ACTG-]+)",x)
            if mtch: mutations_pos1.append(int(mtch.group(1)))
            else: mutations_pos1.append(-1)
        bc3s_perPos_inMuts = aggregate_(keys_=map(str,mutations_pos1),values_=high_qaul_mut_bc3)
        bc3Counts_perPos_inMuts={}
        for pos in bc3s_perPos_inMuts:
            bc3Counts_perPos_inMuts[pos[0]] = dict(Counter(pos[1]))
        bcs3WTcounts = []
        bcs3WTtypes = []
        bcs3WTtypesMissing = []
        for x in range(len(mut3)):
            mutation_name1=mut3[x][cols1['Mutation']]
            mtch=re.match(r"^(\d+)([ACTG-]+)",mutation_name1)
            if mtch:
                pos=mtch.group(1)
                if bc3Counts_perPos_inMuts.get(pos,"") != "":
                    bcs3_inMuts1= bc3Counts_perPos_inMuts[pos]
                    bcs3 = mut3[x][table1.cols['barcode3_all']].split(';')
                    bcs3counts = mut3[x][table1.cols['barcode3_all_counts']].split(';')
                    bcs3WTcounts0 = []
                    for y in range(len(bcs3)):
                        if(bcs3_inMuts1.get(bcs3[y],"")!=""):
                            d = int(bc3_total_counts1[bcs3[y]]) - int(bcs3_inMuts1[bcs3[y]])
                            if d<0: sys.exit("errrr3")
                            bcs3WTcounts0.append(str(d))
                        else:
                            bcs3WTcounts0.append(str(bc3_total_counts1[bcs3[y]]))
                    bcs3WTcounts.append(";".join(bcs3WTcounts0))
                    bcs3WTtypes.append(       sum([int(bcs3WTcounts0[y])>0 and (bcs3[y] != 'BC3_missing') for y in range(len(bcs3WTcounts0))]))
                    bcs3WTtypesMissing.append(sum([int(bcs3WTcounts0[y])>0 and (bcs3[y] == 'BC3_missing') for y in range(len(bcs3WTcounts0))]))
                else:
                    bcs3WTcounts.append('NA')
                    bcs3WTtypes.append('NA')
                    bcs3WTtypesMissing.append('NA')
            elif mutation_name1 == 'WT':
                bcs3WTcounts.append("NA")
                bcs3WTtypes.append("NA")
                bcs3WTtypesMissing.append("NA")
            else:
                sys.exit("errrr2")
        
        mut3 = table1.ins(mut3,'barcode3_all_WTcounts',bcs3WTcounts)
        mut3 = table1.ins(mut3,'barcode3_all_WTcount',bcs3WTtypes)
        mut3 = table1.ins(mut3,'barcode3_missing_all_WTcount',bcs3WTtypesMissing)
       
    # mut3 columns with iformation about barcode3 WT positions:
    # barcode3_all            = barcodes 3' names of all reads in a barcode 5' family
    # barcode3_all_counts     = the corresponding  barcodes 3' counts
    # barcode3_all_WTcounts   = the corresponding barcodes 3' counts for WT alleles at the tested position of mutation
    # barcode3_all_WTcount    = the counts of barcodes 3' types, not including 'BC3_missing', assoicated with WT at the tested position, 
    
    return(mut3) # mut3 : [[mutation_names,mutation_count,mutation_count_HQ,total_count,total_counts_high_qual,HQ_mut_freq,...], ...]


### main ########################################################################

files1 = glob.glob(p.path1)
# sam columns: 'QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','MRNM','MPOS','ISIZE','SEQ','QUAL','TAGs']

out1_cols=['ID','Mutation','Mutation_count','HQ_mutation_count',
    'Total_count','HQ_total_count','Mutation_freq',
    'barcode3_IDs','barcode3_HQ_mut_counts','barcode3_existing','barcode3_missing','barcode3_HQ_mut_freq',
    'barcode3_all','barcode3_all_counts','barcode3_all_WTcounts','barcode3_all_WTcount','barcode3_missing_all_WTcount']

for i in range(len(files1)):
    m1=re.match('.*\\/(.*)',files1[i])
    if m1:
        title1=m1.group(1)
        print "file1 = "+str(files1[i])
        print "title1 = "+str(title1)
        if p.include_BX_tag != "": title1=title1+".tag-"+p.include_BX_tag
        samfile = pysam.AlignmentFile(files1[i], "rb")
        rname_prev=''
        rnames={}
        visited_ref=0
        visited_read=0
        aln=[]
        outfile1 = "%s/%s.mutationFrequncyPerBarcode.txt"%(p.outDir1,title1)
        outfile2 = "%s/%s.mutationFrequncyPerBarcode.log.txt"%(p.outDir1,title1)
        out1=open(outfile1,'w')
        out2=open(outfile2,'w')
        out1.write(("\t".join(out1_cols))+"\n")
        log={'unmapped':0,'ref_start_ne0':0,'query_start_ne0':0,'start_ne0':0,'total':0}
        for read in samfile.fetch(until_eof=True):
            rname=samfile.getrname(read.reference_id)
            
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
                        #if (rname_prev == "AAGGACGTCGTGCACCGC") or (rname_prev == "AATCTAAATGATCGCCGC") or (rname_prev == "ACACATCTTGGCAGCCGC") or (rname_prev == "ACTGCACCGTAAGACCGC"):
                        if True:
                            summarize_mutations1 = summarize_mutations(aln,p.seq_target,p.limit_start,p.limit_end,p.checkMut_name,p.checkMut_pos,p.checkMut_action,BC3_min_cutoff=p.BC3_min_cutoff)
                            for j in range(len(summarize_mutations1)):
                                out1.write("%s\t%s\n"%(rname_prev,"\t".join(map(str,summarize_mutations1[j]))))
                        aln=[]
                    aln.append([rname,read.reference_start,read.cigarstring,read.query_sequence,read.qual,tag1]) # reference_start = 0-based leftmost coordinate; sam file fields: 'RNAME','POS','CIGAR','SEQ','QUAL'
                    rnames[rname]=True
                    rname_prev = rname
                    visited_read +=1
            
            if read.is_unmapped: log['unmapped'] +=1
            else:
                if read.reference_start !=0: log['ref_start_ne0'] +=1 # reference leftmost coordinate
                if read.query_alignment_start !=0: log['query_start_ne0'] +=1  # start index of the aligned query portion
                if (read.reference_start !=0) or (read.query_alignment_start !=0): log['start_ne0'] +=1
            log['total'] += 1
        
        if len(aln)>0:
            summarize_mutations1 = summarize_mutations(aln,p.seq_target,p.limit_start,p.limit_end,p.checkMut_name,p.checkMut_pos,p.checkMut_action,BC3_min_cutoff=p.BC3_min_cutoff)
            for j in range(len(summarize_mutations1)):
              out1.write("%s\t%s\n"%(rname_prev,"\t".join(map(str,summarize_mutations1[j]))))
        
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


