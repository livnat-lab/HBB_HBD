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

# Classes whose files should be stored in the same directory as this script:
from tb2 import tb
from tb2 import tables

# A script for identifying mutations in the aligned reads

# Gather input arguments
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

if (len(p.limit_start)!=len(p.limit_end)) or len(p.limit_start)<1: sys.exit('Error: wrong region of interest limits')
if p.queryAlnMustStartAt0 <=0: p.queryAlnMustStartAt0=False
else: p.queryAlnMustStartAt0=True

if(len(p.checkMut_name)>0):
    p.checkMut_name = p.checkMut_name.split(",")
    p.checkMut_pos = p.checkMut_pos.split(",")
    p.checkMut_action = p.checkMut_action.split(",")
    if(len(p.checkMut_name)!=len(p.checkMut_pos)): sys.exit("Mismatch between number of control mutation positions and control mutation names")
    if(len(p.checkMut_name)!=len(p.checkMut_action)): sys.exit("For each control mutation \"seq_action\" should be indicated")
    
# Read reference sequence fasta file, should fit to reference used for read alignment
seqRef0_ = open(p.seq_target_file, "rU")
seqRef1 = ''
for i,r1 in enumerate(SeqIO.parse(seqRef0_, "fasta")): seqRef1=str(r1.seq); break
seqRef0_.close()
if seqRef1 == '': sys.exit("reference sequence is missing in file \"%\""%(seqRef0))
p.seq_target = seqRef1

#### Functions ####

# A function to add indel positions to the reference and query sequences, to ensure that same coordinate refers to same position in both the reference and the query
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

# A function to aggregate values from input arrays by common term in the first column
def aggregate_(keys_,values_,i=0,j=1):
    groups1=[]
    if len(keys_)!=len(values_): sys.exit('Aggregate function err: mismatch in array sizes')
    arr_of_arr = sorted(zip(keys_,values_), key=lambda x: x[i])
    for key1, group1 in groupby(arr_of_arr, lambda x: x[i]):
        groups1.append([key1,[item1[j] for item1 in group1]])
    return(groups1)

# A function to identify mutated positions in the query reads
def indentify_mutations(aln_target,aln_read,qual,min_phred33,limit_start,limit_end):
    start_p=-1
    start_s=-1
    t = list(aln_target)
    q = list(aln_read)
    q = q + list('-'*(len(t)-min(len(q),len(t))))  # 'x'
    qual=list(qual)
    acc_p = start_p
    
    # Indices of bases in the reference sequence, indels get the same postion as the last non-indel base before them
    t_index=[]
    for x in range(len(t)):
        if((t[x]!='-')and(t[x]!='x')): acc_p=acc_p+1
        t_index.append(acc_p)
    acc_s = start_s
    
    # Indices of bases in the query sequence, indels get the same postion as the last non-indel base before them
    q_index=[] # indices of bases in the query seq
    for x in range(len(q)):
        if((q[x]!='-')and(q[x]!='x')): acc_s=acc_s+1
        q_index.append(acc_s)
    
    # Mutation calls
    tq = []; tq_idx = []; t_idx = []
    for x in range(len(t)):
        if((t[x]!=q[x])and(t[x]!='x')and(q[x]!='x')):
            tq.append("%i%s%s"%(t_index[x]+1,t[x],q[x])) # Mutation info
            tq_idx.append(q_index[x]) # Mutation position in the query
            t_idx.append(t_index[x]) # Mutation position in the reference
        else:
            tq.append(None)
            tq_idx.append(None)
            t_idx.append(None)
            
    # Check that any of the identified mutations fall in the ROI (region of interest)        
    sel1 = [ False for i in range(len(tq))]
    for k in range(len(limit_start)):
        sel1 = [((not(tq[i] is None) and ((limit_start[k]-1) <= t_idx[i]) and ((limit_end[k]-1) >= t_idx[i])) or (sel1[i])) for i in range(len(tq))]
   
    mut=[]
    mutidx=[] # Index of mutation in the query sequece
   
    if(sum(sel1)>0):
        mut =    [ tq[i]     for i in range(len(tq)) if sel1[i] ]
        mutidx = [ tq_idx[i] for i in range(len(tq)) if sel1[i] ]
    else:
        mut    = ['WT']
        mutidx = [None]
    
    # Collapse representation of the consecutive insertion positions (e.g. "40--/CC" -> "40-/CC")
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
    
    # Check and exclude mutations matching positions with low sequencing quality. Deletions are considered high quality by default
    high_qaul_mut=[]
    for z in range(len(mut)):
        if((mut[z] != 'WT') and not(re.match(r'x.*',mut[z]))):  
            if re.match(r'^\d+-\w+$',mut[z]): # Insertion
                asciiChar = qual[mutidx[z]]
                phred  = ord(asciiChar)-33
                if(phred >= min_phred33): high_qaul_mut.append(mut[z])
            elif re.match(r'^\d+\w+-$',mut[z]): # Deletion
                high_qaul_mut.append(mut[z])
            elif re.match(r'^\d+\w\w$',mut[z]): # Mismatch
                asciiChar = qual[mutidx[z]]
                phred  = ord(asciiChar)-33
                if(phred >= min_phred33): high_qaul_mut.append(mut[z])
            else:
                sys.exit('Err: unexpected mutation name "%s\"'%(mut[z]))
        else:
            high_qaul_mut.append(mut[z])
    
    high_qaul_mut = [high_qaul_mut[i] for i in range(len(high_qaul_mut)) if not(re.match(r'.*x.*',high_qaul_mut[i]))]
    all_mut = [mut[i] for i in range(len(mut)) if not(re.match(r'.*x.*',mut[i]))]
    
    # Check sequencing quality along all positions of the reference, excluding insertions
    # t[x] - reference position at coordinate x; q[x] - read position at coordinate x.
    isHighQual = []
    for x in range(len(t)):
        qual1 = qual[q_index[x]]
        if(t[x]!='x') and (t[x]!='-'): # Skipping insertion quality
            if((q[x]!='-')and(q[x]!='x')):
                isHighQual.append((ord(qual1)-33)>min_phred33)
            elif(q[x]=='-'):
                isHighQual.append(True)
            else: isHighQual.append(False)
        else: isHighQual.append(None)
    isHighQual = [isHighQual[i] for i in range(len(isHighQual)) if not(isHighQual[i] is None)]
    
    return(high_qaul_mut,all_mut,isHighQual)

# A function to collect information on all mutations in agiven read family
def summarize_mutations(sam1,seq_target,limit_start,limit_end,
                        checkMut_name,checkMut_pos,checkMut_action,
                        c={'RNAME':0,'POS':1,'CIGAR':2,'SEQ':3,'QUAL':4,'BC3':5},
                        min_phred33=28,BC3_min_cutoff=2):
                        
    high_qaul_mut = [] # HQ (high-quality) mutation names
    high_qaul_mut_bc3 = []
    all_mut = [] # All mutation names
    all_mut_bc3 = []
    isHighQual = [] # Position quality table; rows: reads, columns: quality in aligned read positions (True = high , False = low)
    tables1=tables()
    barcode3 = []
    
    # Check number of reads associated with different BC3 in the family
    barcode3_check = {}
    for s in range(len(sam1)):
        bc3_check=sam1[s][c['BC3']]
        if barcode3_check.get(bc3_check,'')=='': barcode3_check[bc3_check]=1
        elif bc3_check != 'X': barcode3_check[bc3_check] +=1
    
    # Gather and concatenate information on mutations across reads in the BC5 read family
    alns = {}
    for s in range(len(sam1)):
        bc3_check=sam1[s][c['BC3']]
        if barcode3_check[bc3_check] >= BC3_min_cutoff: # Only BC3 groups within the read family associated with number of reads >= cutoff are reported in the mutation summary table
            seq_read = sam1[s][c['SEQ']]
            qual     = sam1[s][c['QUAL']]
            
            # Align reference and read sequences to ensure that same coordinate points to same position in both
            if alns.get(seq_read,"")=="":
                cigar_str = sam1[s][c['CIGAR']]
                pos       = sam1[s][c['POS']]
                aln_target = cigar_to_aligedSequence(cigar_str,pos,seq_target,False) # Add indel positions to the reference sequence
                aln_read   = cigar_to_aligedSequence(cigar_str,0,seq_read,True)      # Add indel positions to the read
                alns[seq_read] = [aln_target,aln_read]
            else:
                aln_target = alns[seq_read][0]
                aln_read   = alns[seq_read][1]
            
            # Collect mutation info
            high_qaul_mut0,all_mut0,isHighQual0 = indentify_mutations(aln_target,aln_read,qual,min_phred33,limit_start,limit_end)
            
            # Exclude reads with specific mutations
            all_mut0_filtered = all_mut0[:]
            for zz in range(len(checkMut_name)):
                if checkMut_action[zz] == "exclude":
                    all_mut0_filtered   = [all_mut0_filtered[z]  for z in range(len(all_mut0_filtered)) if all_mut0_filtered[z] != "%s%s"%(checkMut_pos[zz],checkMut_name[zz])]
            
            # Concatenate mutations from surviving reads
            if len(all_mut0) == len(all_mut0_filtered):
                barcode3.append(sam1[s][c['BC3']])
                all_mut += all_mut0
                high_qaul_mut += high_qaul_mut0
                isHighQual.append(isHighQual0)
                all_mut_bc3 =       all_mut_bc3       + [barcode3[len(barcode3)-1] for z in range(len(all_mut0))]
                high_qaul_mut_bc3 = high_qaul_mut_bc3 + [barcode3[len(barcode3)-1] for z in range(len(high_qaul_mut0))]
     
    if(len(all_mut)==0):
        mut3=[]
    else:
        c1=Counter(all_mut)
        mut1b= [ [k,c1[k]] for k in c1.keys()] # Return counts of all mutations
        
        if(len(high_qaul_mut)>0): # HQ mutations were found
            c2 = Counter(high_qaul_mut)
            mut1a = [ [k,c2[k]] for k in c2.keys()] # Return counts of high-quality (HQ) mutations
            mut2 = tables1.merge(mut1b,mut1a,ncolx0=2,ncoly0=2) # mut2 : [[mutation_name,mutation_count,mutation_count_HQ], ...]; Merge all and HQ mutation counts into a single table
        else:
            mut1a = []
            mut2 = [[mut1b[i][0],mut1b[i][1],0] for i in range(len(mut1b))]
        
        # Count number of high-quality reads associated with each reference position
        isHighQual_ = map(list, zip(*isHighQual)) # Transpose quality arrays
        isHighQual_sum = [sum(isHighQual_[i]) for i in range(len(isHighQual_))] # Count number of HQ reads for each position
        
        # Add read counts and mutation frequencies to mutation summary table
        # Existing columns: {'Mutation':0,'Mutation_count':1,'HQ_mutation_count':2}
        for i in range(len(mut2)):
            totalCount = len(isHighQual_[0]) # Add total read count in the family
            mut2[i].append(totalCount) 
            
            # Case - substitution
            if re.match(r'^\d+\w\w$',mut2[i][0]): 
                mut2[i].append( isHighQual_sum[int(re.match(r'^(\d+).*',mut2[i][0]).group(1))-1] ) # Add counts of HQ reads at mutated positions
                if float(mut2[i][4]) != 0.0: mut2[i].append(float(mut2[i][2])/float(mut2[i][4])) # Calculate HQ mutation frequency (HQ_mut_count/HQ_read_count)
                else: mut2[i].append(0.0)
            
            # Case - indel (insertion/deletion)    
            elif re.match(r'^\d+\w+-$',mut2[i][0]) or re.match(r'^\d+-\w+$',mut2[i][0]):
                mut2[i].append(float(totalCount)) # For indels all reads are considered HQ - adding total read count as HQ read count
                mut2[i].append(float(mut2[i][2])/float(totalCount)) # Calculate indel mutation frequency (HQ_mut_count/total_read_count)
            
            # Case - WT    
            elif mut2[i][0] == 'WT':
                mut2[i].append(len(isHighQual_[0])) # As for indels - total read count is added
                mut2[i].append(float(mut2[i][2])/len(isHighQual_[0])) # WT freq - WT_read_count/total_read_count
                
            else:
                sys.exit('Unexpected mutation type found: \"%s\"'%(mut2[i][0]))
        
        mut2 = sorted(mut2, key=lambda x: x[0])
        
        # Add 3' barcode information
        # Existing columns: {'Mutation':0,'Mutation_count':1,'HQ_mutation_count':2,'Total_count':3,'HQ_total_count':4,'Mutation_freq':5}
        bc3_per_HQmuts = aggregate_(keys_=high_qaul_mut,values_=high_qaul_mut_bc3) # Associate HQ mutations with BC3 sequences of reads in which they were found: [[mut1,[bcX,bcY,bcX,..]], [mut2,[bcX,bcX,bcZ,..]], ...]
        bc3_summary = []
        
        for w in range(len(bc3_per_HQmuts)):
            bc3_mut = bc3_per_HQmuts[w][0]
            
            # Calculate counts of reads with particular mutation in different BC3 groups
            bc3_counter = Counter(bc3_per_HQmuts[w][1]) 
            bc3_names =  ";".join(bc3_counter.keys())
            bc3_counts =  ";".join(map(str,bc3_counter.values()))
            
            # Check if any of reads lacks 3' barcode
            bc3_names_unique = map(str,set(bc3_counter.keys()))
            bc3_count_nonMissing = len(bc3_names_unique)-1 if 'BC3_missing' in bc3_names_unique else len(bc3_names_unique)
            bc3_count_missing    = 1 if 'BC3_missing' in bc3_names_unique else 0
            
            bc3_summary.append([bc3_mut,bc3_names,bc3_counts,bc3_count_nonMissing,bc3_count_missing]) # Aggregate information on 3' barcodes associated with each mutation
            
        mut3 = tables1.merge(mut2,bc3_summary,ncolx0=6,ncoly0=5) # Merge BC3 and mutation count information into a single table
        
        # Calculate and add HQ mutation frequency per each 3BC for each mutation
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
        
        # Add information on 3' barcodes of all reads and 3' barcodes of non-mutated reads
        # Existing columns: {'Mutation':0,'Mutation_count':1,'HQ_mutation_count':2,'Total_count':3,'HQ_total_count':4,'Mutation_freq':5,
        #                    'barcode3_IDs':6,'barcode3_HQ_mut_counts':7,'barcode3_existing':8,'barcode3_missing':9,'barcode3_HQ_mut_freq':10}
        
        cols1={'Mutation':0,'Mutation_count':1,'HQ_mutation_count':2,'Total_count':3,'HQ_total_count':4,'Mutation_freq':5,
                   'barcode3_IDs':6,'barcode3_HQ_mut_counts':7,'barcode3_existing':8,'barcode3_missing':9,'barcode3_HQ_mut_freq':10}
       
        # Identify and count all 3' barcodes found in read family
        bc3_total_counts = Counter(barcode3)
        bc3_total_counts1 = dict(bc3_total_counts)
        bc3_total_counts_str = ";".join(map(str,bc3_total_counts.values()))
        bc3_ids_str = ";".join(bc3_total_counts.keys())
        
        # Add all 3' barcode info
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
                    # Gather total and variant-associated read count per BC3
                    bcs3_inMuts1= bc3Counts_perPos_inMuts[pos]
                    bcs3 = mut3[x][table1.cols['barcode3_all']].split(';')
                    bcs3counts = mut3[x][table1.cols['barcode3_all_counts']].split(';')
                    bcs3WTcounts0 = []
                    
                    # Calculate number of WT reads associated with each BC3
                    for y in range(len(bcs3)):
                        if(bcs3_inMuts1.get(bcs3[y],"")!=""):
                            d = int(bc3_total_counts1[bcs3[y]]) - int(bcs3_inMuts1[bcs3[y]])
                            if d<0: sys.exit("Variant-associated read count can not be larger than total read count per BC3")
                            bcs3WTcounts0.append(str(d))
                        else:
                            bcs3WTcounts0.append(str(bc3_total_counts1[bcs3[y]]))
                    
                    # Calculate number of different 3' barcode types and number of 3' barcode missing reads for each variant-associated BC3 group        
                    bcs3WTcounts.append(";".join(bcs3WTcounts0))
                    bcs3WTtypes.append(       sum([int(bcs3WTcounts0[y])>0 and (bcs3[y] != 'BC3_missing') for y in range(len(bcs3WTcounts0))]))
                    bcs3WTtypesMissing.append(sum([int(bcs3WTcounts0[y])>0 and (bcs3[y] == 'BC3_missing') for y in range(len(bcs3WTcounts0))]))
                
                # Pass on variants below sequencing quality cutoff  
                else:
                    bcs3WTcounts.append('NA')
                    bcs3WTtypes.append('NA')
                    bcs3WTtypesMissing.append('NA')
            
            # For WT families WT info fields are "NA", since variant-associated BC3 fields already reflect this information 
            elif mutation_name1 == 'WT':
                bcs3WTcounts.append("NA")
                bcs3WTtypes.append("NA")
                bcs3WTtypesMissing.append("NA")
            else:
                sys.exit("Unexpected variant name was found in \"Mutation\" field")
        
        # Add BC3 information om WT reads in variant-associated BC3 groups to the data frame 
        mut3 = table1.ins(mut3,'barcode3_all_WTcounts',bcs3WTcounts)
        mut3 = table1.ins(mut3,'barcode3_all_WTcount',bcs3WTtypes)
        mut3 = table1.ins(mut3,'barcode3_missing_all_WTcount',bcs3WTtypesMissing)
       
    # Columns with information about barcode3 WT positions:
    # barcode3_all            = list of all 3' barcode sequences associated with a 5' barcode family
    # barcode3_all_counts     = the corresponding  3' barcode counts
    # barcode3_all_WTcounts   = the corresponding barcodes 3' counts for WT alleles in BC3 groups containing variants
    # barcode3_all_WTcount    = the count of 3' barcode types, not including 'BC3_missing', assoicated with WT at the tested position, 
    
    return(mut3) # mut3 : [[mutation_names,mutation_count,mutation_count_HQ,total_count,total_counts_high_qual,HQ_mut_freq,...], ...]


#### Main ####
# Read and parse read alignment files (in BAM format)
# Alignment file columns: 'QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','MRNM','MPOS','ISIZE','SEQ','QUAL','TAGs'

files1 = glob.glob(p.path1)

out1_cols=['ID','Mutation','Mutation_count','HQ_mutation_count',
    'Total_count','HQ_total_count','Mutation_freq',
    'barcode3_IDs','barcode3_HQ_mut_counts','barcode3_existing','barcode3_missing','barcode3_HQ_mut_freq',
    'barcode3_all','barcode3_all_counts','barcode3_all_WTcounts','barcode3_all_WTcount','barcode3_missing_all_WTcount']

for i in range(len(files1)):
    m1=re.match('.*\\/(.*)',files1[i])
    if m1:
        # Read in input and create output files
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
        
        # Go over reads in the input file and parse their alignment information
        for read in samfile.fetch(until_eof=True):
            rname=samfile.getrname(read.reference_id)
            
            # Check that read is mapped to the reference and alignment starts at the first position of both the read and the reference
            if not(read.is_unmapped) and (read.reference_start==0) and ((read.query_alignment_start==0) or (p.queryAlnMustStartAt0==False)):
                if read.has_tag('XB'):
                    tag1=read.get_tag('XB')
                else: tag1='NA'
                
                # Process reads having user specified BC3 sequences (or all reads with 3BC, if there is no specification)
                if (p.include_BX_tag == "") or ((tag1 != 'NA') and (p.include_BX_tag != "") and (p.include_BX_tag == tag1)):
                    if rname != rname_prev and rname_prev != '':
                        if rnames.get(rname,"")!="": sys.exit("Expected reads to be sorted by their BC5 identifier\n")
                        
                        visited_ref+=1
                        if visited_ref % 1000 == 0:
                            sys.stdout.write('\r' + "visited ref=%i read=%i alnSize=%i title=%s      "%(visited_ref,visited_read,len(aln),title1))
                            sys.stdout.flush()
                        
                        # Output mutation info
                        if True:
                            summarize_mutations1 = summarize_mutations(aln,p.seq_target,p.limit_start,p.limit_end,p.checkMut_name,p.checkMut_pos,p.checkMut_action,BC3_min_cutoff=p.BC3_min_cutoff)
                            for j in range(len(summarize_mutations1)):
                                out1.write("%s\t%s\n"%(rname_prev,"\t".join(map(str,summarize_mutations1[j]))))
                        aln=[]
                        
                    aln.append([rname,read.reference_start,read.cigarstring,read.query_sequence,read.qual,tag1]) # Reference_start = 0-based leftmost coordinate; sam file fields: 'RNAME','POS','CIGAR','SEQ','QUAL'
                    rnames[rname]=True
                    rname_prev = rname
                    visited_read +=1
            
            if read.is_unmapped: log['unmapped'] +=1 # Log cases unmapped reads
            else:
                if read.reference_start !=0: log['ref_start_ne0'] +=1 # Log cases whereby alignment doesn't start at the first position of the reference
                if read.query_alignment_start !=0: log['query_start_ne0'] +=1  # Log cases whereby aligned portion of the query doesn't start from the first query position.
                if (read.reference_start !=0) or (read.query_alignment_start !=0): log['start_ne0'] +=1 # Log cases whereby either reference or query alignment doesn't start from the first position
            log['total'] += 1 # Log total number of processed reads
        
        # Output mutation info for the last family in the alignment file
        if len(aln)>0:
            summarize_mutations1 = summarize_mutations(aln,p.seq_target,p.limit_start,p.limit_end,p.checkMut_name,p.checkMut_pos,p.checkMut_action,BC3_min_cutoff=p.BC3_min_cutoff)
            for j in range(len(summarize_mutations1)):
              out1.write("%s\t%s\n"%(rname_prev,"\t".join(map(str,summarize_mutations1[j]))))
        
        # Output log data
        for k in log.keys():
            out2.write("reads %s\t%i\n"%(k,log[k]))
        out2.write("reads used\t%i\n"%(visited_read))
        out2.write("queryAlnMustStartAt0\t%s\n"%(p.queryAlnMustStartAt0))
        out2.write("bam\t%s\n"%(files1[i]))
        out2.write("fasta\t%s\n"%(p.seq_target_file))
        
        samfile.close()
        out1.close()
        out2.close()
    print "\nOK"


