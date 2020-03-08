
import sys
import re
import math

f1=      sys.argv[1]
outIdx=  sys.argv[2]
second_barcode_size = int(sys.argv[3])
BC3above1_or_BC3atLeast3 = True if int(sys.argv[4]) == 1 else False

# col      = {'Barcode':0,'Consensus':1,'Consensus_TSS':2,'Mutation_freqs':3,'Mutation_counts':4  ,'total_count':5   ,
#             'Positions':6,'WT_freqs_in_mutations':7,'barcode3_existing':8,'barcode3_above1reads':9,'barcode3_missing':10,
#               'barcode3_all_WTcount':11,'barcode3_missing_all_WTcount':12,'unmutatedReads_counts':13,'unmutatedReads_freq':14}

col = {'Barcode':0, 'Consensus':1, 'Consensus_TSS':2, 'Mutation_freqs':3, 'Mutation_counts':4,
       'total_count':5,'Positions':6,'WT_freqs_in_mutations':7,'barcode3_existing':8, 'barcode3mut_above1reads':9,
       'barcode3WT_above1reads':10,'barcode3_above1reads':11,'barcode3_all_WTcount':12,'unmutatedReads_counts':13,'unmutatedReads_freq':14}

# barcode3_above1reads  -> barcode3mut_above1reads
# barcode3_all_WTcount -> 
# barcode3_missing ->deleted
# barcode3_missing_all_WTcount ->deleted

''' input table format:

* = corresponding to mutation positions in BC3

Barcode                   = CB5
* Consensus               = semicolon separared mutations
* Consensus_TSS           = semicolon separared mutations        
* Mutation_freqs	      = semicolon separared HQ mutation frequencies
* Mutation_counts         = semicolon separared HQ mutation read counts
total_count               = total read counts
* Positions               = number of mutations
* WT_freqs_in_mutations   = semicolon separared frequency of WT in positions of mutations
* barcode3_existing       = semicolon separared number of different types of BC3 in mutation position (or among all reads, if Barcode is WT)
* barcode3mut_above1reads = semicolon separared number of different types of BC3 liked to mutations, in mutation position, represented by least 2 reads each
* barcode3WT_above1reads  = semicolon separared number of different types of BC3 linked to WT in mutation position, represented by least 2 reads each
* barcode3_above1reads    = semicolon separared number of different types of BC3 linked to any kind of read, represented by least 2 reads each
* barcode3_all_WTcount	  = semicolon separared number of different types of BC3 linked to WT, in mutation position
unmutatedReads_counts	  = count of reads with no mutations at all
unmutatedReads_freq       = frequency of reads with no mutations at all from all reads

'''


col_out = col

col_names     = [z[0] for z in sorted(col.items(),     key=lambda x: x[1])]
col_names_out = [z[0] for z in sorted(col_out.items(), key=lambda x: x[1])]

min_freq=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
min_count=[0,1,2,3,4,5,6,7,8,9,10,25,50] #[0,3,10,25,50]
if second_barcode_size > 0:
    if BC3above1_or_BC3atLeast3:
        minCountBarcode = [0,1,2,3,4,5]
        bc3groupCountOK=3
        bc3groupsWithReadsAbove1_OK=2
        print "bc3groupCountOK=%i bc3groupsWithReadsAbove1_OK=%i"%(bc3groupCountOK,bc3groupsWithReadsAbove1_OK)
    else:
        minCountBarcode = [0,1,2,3,4,5]
        bc3groupCountOK=0
        bc3groupsWithReadsAbove1_OK=0
        print "bc3groupCountOK=%i bc3groupsWithReadsAbove1_OK=%i"%(bc3groupCountOK,bc3groupsWithReadsAbove1_OK)
else:
    minCountBarcode = [0]
    bc3groupCountOK=0
    bc3groupsWithReadsAbove1_OK=0

out_files={}
out_files_names={}
out_files2={}
rejected_cases={}
for f in range(len(min_freq)):
    for c in range(len(min_count)):
        for bc in range(len(minCountBarcode)):
            name1="mutFreq%s_readCount%i_BC3WithMut%i_BC3above%i"%(str(round(min_freq[f],2)),min_count[c],minCountBarcode[bc],bc3groupsWithReadsAbove1_OK)
            out_files_names[name1] = "%s.%s.txt"%(outIdx,name1)
            out_files[name1]=open("%s.%s.txt"%(outIdx,name1),'w')
            out_files[name1].write(("\t".join(col_names_out))+"\n")
            out_files2[name1]=open("%s.%s.cons-count.txt"%(outIdx,name1),'w')
            out_files2[name1].write("consensus\tbarcodes_count\n")
            rejected_cases["%s.low_count"%(name1)]=0
            rejected_cases["%s.other"%(name1)]=0
            rejected_cases["%s.from_WThaplotype"%(name1)]=0

h1=open(f1,'r')

def float1(x):
    if x=='NA': x=0
    return(float(x))

def select1(mut_name_,x,toMut,toN,writeN=False,add=True,maxVal=False): # select mutations and their info from array x, when accepted or ambiguous
    mut_val=[]
    mut_name=[]
    pos1=positions(mut_name_)
    #print mut_name_
    #print pos1
    pos_accepted={} # tells which positions have accepted mutations, since in such positions ambiguous mutations should not be reported
    for i in range(len(mut_name_)):
        pos_accepted[pos1[i]]=False
    for i in range(len(mut_name_)):
        if toMut[i]: pos_accepted[pos1[i]]=True
    for i in range(len(x)):
        if toMut[i]: # an accepted mutation
            mut_val.append(x[i]) # accepted mutation information
            mut_name.append(mut_name_[i])
        elif toN[i] and not(pos_accepted[pos1[i]]): # if the mutation is not accepted, yet WT postion is not accepted either
            mut_name0 = (re.sub('[A-Za-z\-]+','N',mut_name_[i]))
            if writeN: x0=re.sub('[A-Za-z\-]+','N',x[i])
            else:      x0=x[i]
            
            if (len(mut_name)==0) or (mut_name[len(mut_name)-1] != mut_name0):
                mut_val.append(x0)
                mut_name.append(mut_name0)
            elif add:
                if maxVal:
                    pass
                else:
                    mut_val[len(mut_val)-1] += x0
    return mut_val

def positions(mut_name_):
    positions1=[]
    for i in range(len(mut_name_)):
        m1=re.match('.*?(\d+).*?',mut_name_[i])
        if m1: positions1.append(int(m1.group(1)))
        else: sys.exit('unexpected positions !')
    return(positions1)

def intIfNotNA(x):
    if x=='NA': return('NA')
    else: return(int(x))

for i,r in enumerate(h1):
    r1=r.rstrip().split('\t')
    if i==0:
        if col_names != r1: sys.exit('unexpected column names\n')
        for y in range(len(r1)):
            if col.get(r1[y],"")=="": sys.exit('wrong header 1')
            elif str(col[r1[y]]) != str(y): sys.exit('wrong header 2')
    else:
        freqs0=   r1[col['Mutation_freqs']] # HQ frequencies
        freqs1=   map(float1,freqs0.split(";"))
        
        counts0=  r1[col['Mutation_counts']] # HQ counts
        counts1 = map(int,counts0.split(";"))
        
        count =   int(r1[col['total_count']])
        
        wt_freq_in_mut0 = r1[col['WT_freqs_in_mutations']]
        wt_freq_in_mut1 = map(float1,wt_freq_in_mut0.split(";"))
        
        barcode3_existing0 = r1[col['barcode3_existing']]
        barcode3_existing1 = map(int,barcode3_existing0.split(";"))
        
        barcode3mut_above1reads0 = r1[col['barcode3mut_above1reads']]
        barcode3mut_above1reads1 = map(intIfNotNA,barcode3mut_above1reads0.split(";"))
        
        barcode3WT_above1reads0 = r1[col['barcode3WT_above1reads']]
        barcode3WT_above1reads1 = map(int,barcode3WT_above1reads0.split(";"))
        
        barcode3_above1reads0   = r1[col['barcode3_above1reads']]
        barcode3_above1reads1   = map(int,barcode3_above1reads0.split(";"))
        
        barcode3_all_WTcount0 = r1[col['barcode3_all_WTcount']]
        barcode3_all_WTcount1 = map(intIfNotNA,barcode3_all_WTcount0.split(";"))
        
        aboveMin=bc3groupsWithReadsAbove1_OK
        groupsMin=bc3groupCountOK
        
        if (len(freqs1)!=len(counts1)) or (len(freqs1)!=len(wt_freq_in_mut1)) or (len(freqs1)!=len(barcode3_existing1)) \
            or (len(freqs1)!=len(barcode3mut_above1reads1) or (len(freqs1)!=len(barcode3WT_above1reads1)) \
            or (len(freqs1)!=len(barcode3_above1reads1)) or (len(freqs1)!=len(barcode3_all_WTcount1))):
            sys.exit('err1')
        
        if(min(counts1)<0):
            print("excluded:")
            print(r1)
        else:
            allels1 = r1[col['Consensus']].split(';')
            allels2 = r1[col['Consensus_TSS']].split(';')
            if (len(allels1)==len(freqs1)):
                for f in range(len(min_freq)):
                    for c in range(len(min_count)):
                        for bc in range(len(minCountBarcode)):
                            #name1 ="mutFreq%s_readCount%i_BC3WithMut%i"           %(str(round(min_freq[f],2)),min_count[c],minCountBarcode[bc])
                            name1 = "mutFreq%s_readCount%i_BC3WithMut%i_BC3above%i"%(str(round(min_freq[f],2)),min_count[c],minCountBarcode[bc],bc3groupsWithReadsAbove1_OK)
                            if (count >= min_count[c]):
                                if r1[col['Consensus']] != 'WT':
                                    is_mut_accepted0    = [ frq >= min_freq[f] for frq in freqs1]
                                    is_WT_accepted0     = [ frq >=  min_freq[f] for frq in wt_freq_in_mut1]
                                    
                                    is_secondBC_accepted_mut =    [ (barcode3_existing1[u]    >= minCountBarcode[bc]) and
                                                                    ((barcode3mut_above1reads1[u] >=  aboveMin) or (barcode3_existing1[u] >= groupsMin))
                                                                    for u in range(len(barcode3_existing1)) ]
                                    
                                    is_secondBC_accepted_wt  =    [ (barcode3_all_WTcount1[u]    >= minCountBarcode[bc]) and
                                                                    ((barcode3WT_above1reads1[u] >=  aboveMin) or (barcode3_all_WTcount1[u] >= groupsMin))
                                                                    for u in range(len(barcode3_all_WTcount1)) ]
                                    
                                    is_mut_accepted = [is_mut_accepted0[zzz] and  is_secondBC_accepted_mut[zzz] for zzz in range(len(is_mut_accepted0))]
                                    is_WT_accepted  = [is_WT_accepted0[zzz]  and  is_secondBC_accepted_wt[zzz]  for zzz in range(len(is_WT_accepted0))]
                                    is_WT_not_accepted = [not(zzz) for zzz in is_WT_accepted]
                                    
                                    # conditions of accepting mutations:
                                    # 1) Accepted mutation names, in positions with mutation above the cutoffs.
                                    # 2) 'N's, in positions with neither an accepted mutation, nor 'WT' nucleotide above cutoffs.
                                    # 3) 'WT', if no mutations or Ns are reported along the entire gene, and 'WT' above cutoff.
                                    # 4) if multiple mutations are accepted at the same position, then all are shown.
                                    # 5) if multiple 'N's exist at the same position, then a single 'N' is shown.
                                    
                                    if (sum(is_mut_accepted) > 0) or (sum(is_WT_not_accepted) > 0): # some mutation positions are accepted, or their WT cannot be accepted
                                        r2 = r1[:]
                                        r2[col['Consensus']] =             ";".join(             select1(allels1,allels1,                 is_mut_accepted,is_WT_not_accepted,writeN=True,add=False))
                                        r2[col['Consensus_TSS']] =         ";".join(             select1(allels1,allels2,                 is_mut_accepted,is_WT_not_accepted,writeN=True,add=False))
                                        r2[col['Mutation_freqs']] =        ";".join(map(str,     select1(allels1,freqs1,                  is_mut_accepted,is_WT_not_accepted)))
                                        r2[col['Mutation_counts']] =       ";".join(map(str,     select1(allels1,counts1,                 is_mut_accepted,is_WT_not_accepted)))
                                        r2[col['WT_freqs_in_mutations']] = ";".join(map(str,     select1(allels1,wt_freq_in_mut1,         is_mut_accepted,is_WT_not_accepted)))
                                        r2[col['Positions']] =                          str( len(select1(allels1,allels1,                 is_mut_accepted,is_WT_not_accepted)))
                                        r2[col['barcode3_existing']] =      ";".join(map(str,    select1(allels1,barcode3_existing1,      is_mut_accepted,is_WT_not_accepted)))
                                        r2[col['barcode3mut_above1reads']] =";".join(map(str,    select1(allels1,barcode3mut_above1reads1,is_mut_accepted,is_WT_not_accepted)))
                                        r2[col['barcode3_all_WTcount']] =   ";".join(map(str,    select1(allels1,barcode3_all_WTcount1,   is_mut_accepted,is_WT_not_accepted)))
                                        
                                        r2[col['barcode3mut_above1reads']] =   ";".join(map(str,    select1(allels1,barcode3mut_above1reads1,is_mut_accepted,is_WT_not_accepted)))
                                        r2[col['barcode3WT_above1reads']] =    ";".join(map(str,    select1(allels1,barcode3WT_above1reads1, is_mut_accepted,is_WT_not_accepted)))
                                        r2[col['barcode3_all_WTcount']] =      ";".join(map(str,    select1(allels1,barcode3_all_WTcount1,   is_mut_accepted,is_WT_not_accepted)))
                                        r2[col['barcode3_above1reads']] =      ";".join(map(str,    select1(allels1,barcode3_above1reads1,   is_mut_accepted,is_WT_not_accepted)))
                                        
                                        out_files[name1].write(("\t".join(r2))+"\n")
                                    elif sum(is_WT_not_accepted) == 0: # no positions with mutations are accepted, yet all their WT reads are accepted
                                        r2 = r1[:]
                                        r2[col['Consensus']] =            "WT"
                                        r2[col['Consensus_TSS']] =        "WT"
                                        r2[col['Mutation_freqs']] =       r2[col['unmutatedReads_freq']]
                                        r2[col['Mutation_counts']] =      r2[col['unmutatedReads_counts']]
                                        r2[col['WT_freqs_in_mutations']] = 'NA'
                                        r2[col['Positions']] =            str(0)
                                        r2[col['barcode3_existing']] =  r2[col['barcode3_all_WTcount']]
                                        r2[col['barcode3_above1reads']] = 'NA'
                                        out_files[name1].write(("\t".join(r2))+"\n")
                                    else:
                                        rejected_cases["%s.other"%(name1)] +=1
                                else:
                                    if(len(barcode3_existing1)!=1): sys.exit('error: unexpected length in \'barcode3_existing1\' column !')
                                    if(len(barcode3WT_above1reads1)!=1): sys.exit('error: unexpected length in \'barcode3WT_above1reads1\' column !')
                                    
                                    if (barcode3_existing1[0] >= minCountBarcode[bc]) and ((barcode3WT_above1reads1[0] >=  aboveMin) or (barcode3_existing1[0] >= groupsMin)):
                                        r2 = r1[:]
                                        r2[col['Consensus']] =            "WT"
                                        r2[col['Consensus_TSS']] =        "WT"
                                        r2[col['Mutation_freqs']] =       r2[col['unmutatedReads_freq']]
                                        r2[col['Mutation_counts']] =      r2[col['unmutatedReads_counts']]
                                        r2[col['WT_freqs_in_mutations']] = 'NA'
                                        r2[col['Positions']] =            str(1)
                                        r2[col['barcode3_above1reads']] = 'NA'
                                        out_files[name1].write(("\t".join(r2))+"\n")
                                    else:
                                        #print r1
                                        #print barcode3_existing1
                                        rejected_cases["%s.from_WThaplotype"%(name1)] +=1
                            else:
                                rejected_cases["%s.low_count"%(name1)] +=1
            else: sys.exit('err2')
        #if i >3000: break

h1.close()

for k in out_files.keys(): out_files[k].close()

for k in out_files_names.keys():
    hap_bcCount={}
    h2=open(out_files_names[k],'r')
    for i,r in enumerate(h2):
        if i > 0:
            r1=r.rstrip().split('\t')
            hap=r1[col['Consensus']]
            if hap_bcCount.get(hap,"") == "": hap_bcCount[hap]=1
            else: hap_bcCount[hap] += 1
    h2.close()
    hap_bcCount1=sorted(hap_bcCount.items(), reverse=True,key=lambda k: k[1])
    out_files2[k].write("%s\t%i\n"%('rejected.low_count',rejected_cases["%s.low_count"%(k)]))
    out_files2[k].write("%s\t%i\n"%('rejected.from_WThaplotype',rejected_cases["%s.from_WThaplotype"%(k)]))
    out_files2[k].write("%s\t%i\n"%('rejected.other',rejected_cases["%s.other"%(k)]))
    for hap in hap_bcCount1:
        out_files2[k].write("%s\t%i\n"%(hap[0],hap[1]))

