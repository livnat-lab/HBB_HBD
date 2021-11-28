import sys
import re
import math

# This script compares summary table of mutations per read family against a set of cutoff criteria to determine which mutations should be considered as "true" 

f1=      sys.argv[1]
outIdx=  sys.argv[2]
second_barcode_size = int(sys.argv[3])
BC3above1_or_BC3atLeast3 = True if int(sys.argv[4]) == 1 else False

col = {'Barcode':0, 'Consensus':1, 'Consensus_TSS':2, 'Mutation_freqs':3, 'Mutation_counts':4,
       'total_count':5,'Positions':6,'WT_freqs_in_mutations':7,'barcode3_existing':8, 'barcode3mut_above1reads':9,
       'barcode3WT_above1reads':10,'barcode3_above1reads':11,'barcode3_all_WTcount':12,'unmutatedReads_counts':13,'unmutatedReads_freq':14}


''' 
Input table format:

* = data corresponding to mutation positions in different BC3 groups within given read family

Barcode                   = BC5, a unique identifier of each family
* Consensus               = semicolon separared mutations
* Consensus_TSS           = semicolon separared mutations, with position relative to the translation start site        
* Mutation_freqs	   = semicolon separared HQ mutation frequencies
* Mutation_counts         = semicolon separared HQ mutation read counts
total_count               = total read counts
* Positions               = number of mutations
* WT_freqs_in_mutations   = semicolon separared frequency of WT in positions of mutations
* barcode3_existing       = semicolon separared number of different types of BC3 linked to each mutation position (or among all reads, if Barcode is WT)
* barcode3mut_above1reads = semicolon separared number of different types of BC3 linked to each mutation position and represented by at least 2 reads
* barcode3WT_above1reads  = semicolon separared number of different types of BC3 linked to WT reads in each mutation position and represented by least 2 reads
* barcode3_above1reads    = semicolon separared number of different types of BC3 linked to any kind of read, represented by least 2 reads each
* barcode3_all_WTcount	   = semicolon separared number of different types of BC3 linked to WT, in mutation position
unmutatedReads_counts	   = count of reads with no mutations at all
unmutatedReads_freq       = frequency of reads with no mutations at all from all reads

'''

#### Declaration of cutoff criteria ####
col_out = col

col_names     = [z[0] for z in sorted(col.items(),     key=lambda x: x[1])]
col_names_out = [z[0] for z in sorted(col_out.items(), key=lambda x: x[1])]

min_freq=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
min_count=[0,1,2,3,4,5,6,7,8,9,10,25,50]
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

#### Consensus table analysis ####

# Create output files
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

# Substitute 'N/A' values by '0'
def float1(x):
    if x=='NA': x=0
    return(float(x))

# Select mutations passing the cutoff criteria. Substitute mutations at ambigous positions by 'N's.
def select1(mut_name_,x,toMut,toN,writeN=False,add=True,maxVal=False): 
    mut_val=[]
    mut_name=[]
    pos1=positions(mut_name_)
    pos_accepted={} # Tells which positions have accepted mutations, since in such positions ambiguous mutations should not be reported
    
    for i in range(len(mut_name_)):
        pos_accepted[pos1[i]]=False
    for i in range(len(mut_name_)):
        if toMut[i]: pos_accepted[pos1[i]]=True
    for i in range(len(x)):
        if toMut[i]: # An accepted mutation
            mut_val.append(x[i]) # Append accepted mutation information
            mut_name.append(mut_name_[i])
        elif toN[i] and not(pos_accepted[pos1[i]]): # Ambigous position: both mutation and WT are not accepted at given position
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

# Extract mutation position from variant name (e.g. 39CT -> 39)
def positions(mut_name_):
    positions1=[]
    for i in range(len(mut_name_)):
        m1=re.match('.*?(\d+).*?',mut_name_[i])
        if m1: positions1.append(int(m1.group(1)))
        else: sys.exit('Err: unclear mutation position!')
    return(positions1)

# Check if x is 'NA' and return its value as integer if not 'NA'
def intIfNotNA(x):
    if x=='NA': return('NA')
    else: return(int(x))


#### Main ####
for i,r in enumerate(h1):
    r1=r.rstrip().split('\t')
    if i==0:
        if col_names != r1: sys.exit('Err: unexpected column names\n')
        for y in range(len(r1)):
            if col.get(r1[y],"")=="": sys.exit('Err: wrong header in the input table')
            elif str(col[r1[y]]) != str(y): sys.exit('Err: wrong header position in the input table')
    else:
        # Extract mutation information
        freqs0=   r1[col['Mutation_freqs']] # High quality mutation frequencies
        freqs1=   map(float1,freqs0.split(";"))
        
        counts0=  r1[col['Mutation_counts']] # High quality mutation counts
        counts1 = map(int,counts0.split(";"))
        
        count =   int(r1[col['total_count']])
        
        wt_freq_in_mut0 = r1[col['WT_freqs_in_mutations']]
        wt_freq_in_mut1 = map(float1,wt_freq_in_mut0.split(";"))
        
        barcode3_existing0 = r1[col['barcode3_existing']] # Count of different BC3 associated with each mutation
        barcode3_existing1 = map(int,barcode3_existing0.split(";"))
        
        barcode3mut_above1reads0 = r1[col['barcode3mut_above1reads']] # Count of different BC3 associated with each mutation and represented by at least 2 reads
        barcode3mut_above1reads1 = map(intIfNotNA,barcode3mut_above1reads0.split(";"))
        
        barcode3WT_above1reads0 = r1[col['barcode3WT_above1reads']] # Count of different BC3 associated with WT reads at each mutated position and represented by at least 2 reads
        barcode3WT_above1reads1 = map(int,barcode3WT_above1reads0.split(";"))
        
        barcode3_above1reads0   = r1[col['barcode3_above1reads']] # Total count of different BC3 at each mutated position (associated aither with mutated or WT reads at the position), represented by at least 2 reads
        barcode3_above1reads1   = map(int,barcode3_above1reads0.split(";"))
        
        barcode3_all_WTcount0 = r1[col['barcode3_all_WTcount']] # Count of different BC3 associated with WT reads at each mutated position 
        barcode3_all_WTcount1 = map(intIfNotNA,barcode3_all_WTcount0.split(";"))
        
        aboveMin=bc3groupsWithReadsAbove1_OK
        groupsMin=bc3groupCountOK
        
        if (len(freqs1)!=len(counts1)) or (len(freqs1)!=len(wt_freq_in_mut1)) or (len(freqs1)!=len(barcode3_existing1)) \
            or (len(freqs1)!=len(barcode3mut_above1reads1) or (len(freqs1)!=len(barcode3WT_above1reads1)) \
            or (len(freqs1)!=len(barcode3_above1reads1)) or (len(freqs1)!=len(barcode3_all_WTcount1))):
            sys.exit('Err: missing data on some mutation positions')
        
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
                            name1 = "mutFreq%s_readCount%i_BC3WithMut%i_BC3above%i"%(str(round(min_freq[f],2)),min_count[c],minCountBarcode[bc],bc3groupsWithReadsAbove1_OK)
                            if (count >= min_count[c]): # Check that family size is above cutoff criteria
                                if r1[col['Consensus']] != 'WT':
                                
                                    # Check that mutation and/or WT frequencies are above cutoff criteria
                                    is_mut_accepted0    = [ frq >= min_freq[f] for frq in freqs1]
                                    is_WT_accepted0     = [ frq >=  min_freq[f] for frq in wt_freq_in_mut1]
                                    
                                    # Check that number of different BC3 associated witn mutation or WT reads at each mutation position are above the cutoffs
                                    is_secondBC_accepted_mut =    [ (barcode3_existing1[u]    >= minCountBarcode[bc]) and
                                                                    ((barcode3mut_above1reads1[u] >=  aboveMin) or (barcode3_existing1[u] >= groupsMin))
                                                                    for u in range(len(barcode3_existing1)) ]
                                    
                                    is_secondBC_accepted_wt  =    [ (barcode3_all_WTcount1[u]    >= minCountBarcode[bc]) and
                                                                    ((barcode3WT_above1reads1[u] >=  aboveMin) or (barcode3_all_WTcount1[u] >= groupsMin))
                                                                    for u in range(len(barcode3_all_WTcount1)) ]
                                    
                                    # Check for each mutated position if mutation- or WT-linked reads are pssing the combined cutoff criteria
                                    is_mut_accepted = [is_mut_accepted0[zzz] and  is_secondBC_accepted_mut[zzz] for zzz in range(len(is_mut_accepted0))]
                                    is_WT_accepted  = [is_WT_accepted0[zzz]  and  is_secondBC_accepted_wt[zzz]  for zzz in range(len(is_WT_accepted0))]
                                    is_WT_not_accepted = [not(zzz) for zzz in is_WT_accepted]
                                    
                                    '''
                                    Output of the cutoff check:
                                    1) Mutation name in positions with mutation above the cutoff criteria (e.g. 39CT).
                                    2) 'N's in positions with neither a mutation, nor 'WT' nucleotide are above cutoff criteria.
                                    3) 'WT', if no mutations or 'N's are reported, and 'WT' nucleotides are above cutoff criteria in all positions with suspected mutations.
                                    4) If multiple mutations are accepted at the same position, then all are shown. 
                                    5) If multiple 'N's are declared at the same position, they are collapsed together (e.g.: instead of '49NN' -> '49N' would be shown).
                                    '''
                                    
                                    # Update mutation summary data to return accepted mutations or 'N's for ambigous positions
                                    if (sum(is_mut_accepted) > 0) or (sum(is_WT_not_accepted) > 0): 
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
                                        
                                        out_files[name1].write(("\t".join(r2))+"\n") # Write updated mutation data to the output file
                                        
                                    # Update mutation summary data to present read families as 'WT' if all the mutations in the family are rejected, but their associated WT reads are all accepted   
                                    elif sum(is_WT_not_accepted) == 0: 
                                        r2 = r1[:]
                                        r2[col['Consensus']] =            "WT"
                                        r2[col['Consensus_TSS']] =        "WT"
                                        r2[col['Mutation_freqs']] =       r2[col['unmutatedReads_freq']]
                                        r2[col['Mutation_counts']] =      r2[col['unmutatedReads_counts']]
                                        r2[col['WT_freqs_in_mutations']] = 'NA'
                                        r2[col['Positions']] =            str(0)
                                        r2[col['barcode3_existing']] =  r2[col['barcode3_all_WTcount']]
                                        r2[col['barcode3_above1reads']] = 'NA'
                                        out_files[name1].write(("\t".join(r2))+"\n") # Write updated mutation data to the output file
                                    else:
                                        rejected_cases["%s.other"%(name1)] +=1
                                else:
                                    if(len(barcode3_existing1)!=1): sys.exit('Err: for WT family single entry in \'barcode3_existing1\' column is expected!')
                                    if(len(barcode3WT_above1reads1)!=1): sys.exit('Err: for WT family single entry in \'barcode3WT_above1reads1\' column is expected!')
                                    
                                    # Update mutation summary data on the WT families passing the cutoff criteria
                                    if (barcode3_existing1[0] >= minCountBarcode[bc]) and ((barcode3WT_above1reads1[0] >=  aboveMin) or (barcode3_existing1[0] >= groupsMin)):
                                        r2 = r1[:]
                                        r2[col['Consensus']] =            "WT"
                                        r2[col['Consensus_TSS']] =        "WT"
                                        r2[col['Mutation_freqs']] =       r2[col['unmutatedReads_freq']]
                                        r2[col['Mutation_counts']] =      r2[col['unmutatedReads_counts']]
                                        r2[col['WT_freqs_in_mutations']] = 'NA'
                                        r2[col['Positions']] =            str(1)
                                        r2[col['barcode3_above1reads']] = 'NA'
                                        out_files[name1].write(("\t".join(r2))+"\n") # Write updated mutation data to the output file
                                    else:
                                        rejected_cases["%s.from_WThaplotype"%(name1)] +=1 # Reject WT family as rejected if any of cutoff criteria is not cleared 
                            else:
                                rejected_cases["%s.low_count"%(name1)] +=1 # Reject read family (mutation or WT) due to low count of total reads in the family
            else: sys.exit('Err: Missing frequency count for some mutations in the read family')
        #if i>3000: break # Uncomment to perform test run of the script on the first 3000 lines of the input

h1.close()

for k in out_files.keys(): out_files[k].close()

# Write output for consensus count files summing for each mutation number of families it was found in
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

