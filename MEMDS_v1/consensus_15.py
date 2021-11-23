import re
import re
import sys
from tb2 import tb

''' input table format:

    Each BC5 family can be represented by one or more lines, where each line stores statistics of a single mutation, or about WT reads (reads with no mutation at all).
    If specific BC5 group includes both reads with and without mutations, then the WT statistics are stored in one of the lines, and the other line/s store mutation/s statistics.
    
    For lines that store mutation statistics:
   
   ID                               = BC5 string
   * Mutation                       = mutation string
   * Mutation_count                 = reads counts with mutations
   * HQ_mutation_count              = hiqh-quality reads counts with mutations
   * Total_count                    = total reads counts
   * HQ_total_count                 = total hiqh-quality reads counts
   * Mutation_freq                  = Mutation frequency calcluated as HQ_mutation_count/Total_HQ_reads,
                                      where Total_HQ_reads denotes the total count of non-mutated reads in specific position.
   
   **  barcode3_IDs                 = BC3 strings types, linked to mutations, semicolon separated
   **  barcode3_HQ_mut_counts       = high quality reads counts (semicolon separated, corresponding to barcode3_IDs)
   **  barcode3_existing            = number of different types of BC3 strings, linked to mutations
   **  barcode3_missing             = Deprecated
   **  barcode3_HQ_mut_freq         = freqeuncy of mutations in each BC3 string, (semicolon separated, corresponding to barcode3_IDs)
   
   *** barcode3_all                 = all BC3 strings types
   *** barcode3_all_counts          = reads counts (semicolon separated, corresponding to barcode3_all)
   *** barcode3_all_WTcounts        = reads counts linked to WT positions, or reads without mutations, (semicolon separated, corresponding to barcode3_all)
   *** barcode3_all_WTcount         = number of different types of BC3 strings linked to WT positions, or reads without mutations.
   *** barcode3_missing_all_WTcount = Deprecated
   
   For lines that store WT statistics (reads with no mutations at all):
  
   last 3 columns are NA, and their equivalents are as foolows: 
   * barcode3_all_WTcounts         ->  barcode3_HQ_mut_counts
   * barcode3_all_WTcount          ->  barcode3_existing
   * barcode3_missing_all_WTcount  =  Deprecated
   
'''

def float1(x):
    if x=='NA': x=0
    return(float(x))

def wt_freq(positions,freqs):
    pos_freq={}
    if len(positions)==len(freqs):
        for i in range(len(positions)):
            if pos_freq.get(positions[i],"")=="": pos_freq[positions[i]]= float1(freqs[i])
            else: pos_freq[positions[i]] = pos_freq[positions[i]] + float1(freqs[i])
        wt_freq1=[ 1.0 - pos_freq[positions[i]] for i in range(len(positions))]
    else: sys.exit('expected equal arrays')
    return wt_freq1

def casesAboveMinimunCount(countsStr,minCount):
    def int1(x): return(int(x) if str(x)!="NA" else 0) # expecting "NA" when no HQ mutation is found
    counts1 = map(int1,countsStr.split(";"))
    return sum([c > minCount for c in counts1])    

def consensus(m_,translationStartSite):
    t=tb(cols)
    m_1 = t.ne(m_,'Mutation',"WT")
    m_2 = t.eq(m_,'Mutation',"WT")
    HQmutations_found=False
    if(len(m_1) > 0):
        if sum(map(int,t.sel(m_1,'HQ_mutation_count'))) > 0:
            HQmutations_found=True
    if HQmutations_found:
        positions_to_sort=[int(re.match(r'.*?(\d+).*?',x).group(1)) for x in t.sel(m_1,'Mutation')]
        m_1 = t.ins(m_1,'pos',positions_to_sort)
        m_1 = t.sort(m_1,'pos')
        mut = [(re.match(r'.*?([A-Za-z\-]+).*?',x).group(1)) for x in t.sel(m_1,'Mutation')]
        pos=t.sel(m_1,'pos')
        checkConflicts=[True if (int(pos[c])==int(pos[c-1])) else False for c in range(1,len(pos))] # count of cases where different types of mutations occur at the same position
        if sum(checkConflicts) != 0:
             print "two types of mutations at the same postion in \"%s\""%(t.sel(m_1,'ID')[0])
        procede1=True
        if procede1:
            h = ";".join(t.sel(m_1,'Mutation'))
            total_count = t.sel(m_1,'Total_count')[0]
            counts_hq = ";".join(map(str,t.sel(m_1,'HQ_mutation_count')))
            total_hq =  ";".join(map(str,t.sel(m_1,'HQ_total_count')))
            
            #freqs_hq0 = [ c/float(total_count) for c in map(float,t.sel(m_1,'HQ_mutation_count'))] #* *#
            #freqs_hq = ";".join(map(str,freqs_hq0)) #* *#
            
            freqs_hq0 = t.sel(m_1,'Mutation_freq')
            freqs_hq =  ";".join(map(str,freqs_hq0))
            wt_freqs1 = ";".join(map(str,wt_freq(t.sel(m_1,'pos'),freqs_hq0)))
            
            if sum([True for zzz in freqs_hq0 if ((float(zzz) < 0.0) or (float(zzz) > 1.0))])>0:
                print "out of limits freq values for ID=\"%s\""%(t.sel(m_1,'ID')[0]) 
            
            barcode3_existing =   ";".join(map(str,t.sel(m_1,'barcode3_existing'))) # total count of different types of barcodes in mutation position (or for all reads if the reads are all WT)
            barcode3_missing =    ";".join(map(str,t.sel(m_1,'barcode3_missing'))) 
            barcode3_all_WTcount = ";".join(map(str,t.sel(m_1,'barcode3_all_WTcount')))
            barcode3mut_above1reads = ";".join([str(casesAboveMinimunCount(c,1)) for c in t.sel(m_1,'barcode3_HQ_mut_counts')]) # barcode3mut_above1reads variable equals to the total count of different types of barcodes represented by above 1 read (at least 2), for reads that include a mutation in each tested position
            barcode3WT_above1reads  = ";".join([str(casesAboveMinimunCount(c,1)) for c in t.sel(m_1,'barcode3_all_WTcounts')])  # barcode3WT_above1reads  variable equals to the total count of different types of barcodes represented by above 1 read (at least 2), for reads that do not include a mutation in each tested position
            barcode3_above1reads    = ";".join([str(casesAboveMinimunCount(c,1)) for c in t.sel(m_1,'barcode3_all_counts')])    # barcode3WT_above1reads  variable equals to the total count of different types of barcodes represented by above 1 read (at least 2), for all reads
            
            # add pos relative to TSS (translation start site), and i,d labels (insertions/deletions):
            pos_tss = [int(p) - translationStartSite + 1 for p in t.sel(m_1,'pos')]
            for pos_tss_0 in pos_tss:
                if pos_tss_0 <= 0: sys.exit("error: expected positive TSS, found: \'%i\'"%(pos_tss_0))
            t.ins(m_1,'pos_TSS',pos_tss)
            haplotype_TSS=[]
            for x in range(len(mut)):
                pos0 = pos_tss[x]
                mut0 = mut[x]
                if re.match(r'[A-Za-z]+-',mut0):   mut0=re.sub(r'-','d',mut0)
                elif re.match(r'-[A-Za-z]+',mut0): mut0=re.sub(r'-','i',mut0)
                haplotype_TSS.append("%s%s"%(pos0,mut0))
            h_tss = ";".join(haplotype_TSS)
            
            # df columns: 'Barcode', 'Consensus', 'Consensus_TSS', 'Mutation_freqs', 'Mutation_counts', 'total_count','Positions',
            #             'WT_freqs_in_mutations','barcode3_existing','barcode3mut_above1reads','barcode3WT_above1reads','barcode3_above1reads','barcode3_all_WTcount'
            # df columns to be added in the next section here: 'unmutatedReads_counts','unmutatedReads_freq'
            df=[t.sel(m_1,'ID')[0],h,h_tss,freqs_hq,counts_hq,total_count,len(m_1),
                wt_freqs1,barcode3_existing,barcode3mut_above1reads,barcode3WT_above1reads,barcode3_above1reads,barcode3_all_WTcount]
        else:
            df=[]
            print "skipped ID=\"%s\" - ambiguous mutations"%(t.sel(m_1,'ID')[0])
    elif(len(m_2) > 0):
        # if the data refers to a WT group, then in the input tables: mutations statistics (frequency,count) should include zero values. The total count column should include the total count of reads.
        barcode3_existing =        ";".join(map(str,t.sel(m_2,'barcode3_existing')))
        barcode3_missing =         ";".join(map(str,t.sel(m_2,'barcode3_missing')))
        barcode3mut_above1reads =  ";".join([str(casesAboveMinimunCount(c,1)) for c in t.sel(m_2,'barcode3_HQ_mut_counts')]) # barcode3_above1reads variable equals to the total count of different types of barcodes represented by above 1 read (at least 2)
        barcode3_above1reads =     ";".join([str(casesAboveMinimunCount(c,1)) for c in t.sel(m_2,'barcode3_all_counts')])
        
        total_count = t.sel(m_2,'Total_count')[0] #* *#
        #freqs_hq0 = [ c/float(total_count) for c in map(float,t.sel(m_2,'HQ_mutation_count'))] #* *#
        #freqs_hq = ";".join(map(str,freqs_hq0)) #* *#
        
        freqs_hq =  ";".join(map(str,t.sel(m_2,'Mutation_freq')))
        
        df=[t.sel(m_2,'ID')[0],'WT','WT',0,0,total_count,0,
            freqs_hq,barcode3_existing,'NA',barcode3_above1reads,barcode3_above1reads,barcode3_existing]
    else:
        df=[]
        print "skipped ID=\"%s\""%(t.sel(m_1,'ID')[0] if(len(m_1) > 0) else t.sel(m_2,'ID')[0])
        
    if len(df) > 0:
        if(len(m_2) > 0):
            if len(m_2)!=1: sys.exit('err3 in \'%s\''%(df[0]))
            freq_hq2 =  t.sel(m_2,'Mutation_freq')[0]
            count_hq2 = t.sel(m_2,'HQ_mutation_count')[0]
            total_count2 = t.sel(m_2,'Total_count')[0]
            if HQmutations_found:
                if (int(total_count2) != int(df[5])):
                    sys.exit('err4 in \'%s\''%(df[0]))
            df = df + [count_hq2,freq_hq2]
        else:
            df = df + [0,0]
    
    return(df)

if __name__ == '__main__':
    
    cols={'ID':0,'Mutation':1,'Mutation_count':2,'HQ_mutation_count':3,
     'Total_count':4,'HQ_total_count':5,'Mutation_freq':6,
     'barcode3_IDs':7, 'barcode3_HQ_mut_counts':8,'barcode3_existing':9,'barcode3_missing':10,'barcode3_HQ_mut_freq':11,
     'barcode3_all':12,'barcode3_all_counts':13,'barcode3_all_WTcounts':14,'barcode3_all_WTcount':15,'barcode3_missing_all_WTcount':16}
    out_cols=['Barcode', 'Consensus', 'Consensus_TSS', 'Mutation_freqs', 'Mutation_counts', 'total_count','Positions',
              'WT_freqs_in_mutations','barcode3_existing','barcode3mut_above1reads','barcode3WT_above1reads','barcode3_above1reads','barcode3_all_WTcount',
              'unmutatedReads_counts','unmutatedReads_freq']
    
    in1 = sys.argv[1]
    out1 = sys.argv[2]
    translationStartSite = int(sys.argv[3])
    
    t1=open(in1,'r')
    o1=open(out1,'w')
    prev_id=''
    barcode_group=[]
    visited={}
    o1.write("\t".join(out_cols)+"\n")
    for i,r in enumerate(t1):
        if i==0:
            x=r.rstrip().split("\t")
            for y in range(len(x)):
                if cols.get(x[y],"")=="": sys.exit('wrong header 1')
                elif str(cols[x[y]]) != str(y): sys.exit('wrong header 2')
        else:
            x=r.rstrip().split("\t")
            if len(x)!=len(cols): sys.exit('err1')
            if (prev_id != x[cols['ID']]) and (len(barcode_group)>0):
                c = consensus(barcode_group,translationStartSite)
                if len(c) > 0:
                    o1.write("\t".join(map(str,c))+"\n")
                else:
                    pass # print "ID=\"%s\" was skipped"%(prev_id)
                barcode_group=[]
                if visited.get(x[cols['ID']],"")!="": sys.exit('err2')
                visited[x[cols['ID']]]=True
                #print '-------'
            barcode_group.append(x)
            prev_id = x[cols['ID']]
        #if i > 100000: break
    c = consensus(barcode_group,translationStartSite)
    if len(c) > 0:
        o1.write("\t".join(map(str,c))+"\n")
        barcode_group=[]
        t1.close()
        o1.close()
    else:
        pass #print "ID=\"%s\" was skipped"%(prev_id)
