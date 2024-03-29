import csv
import sys
import re
import os
import argparse

class mapper_call:
    # Create output files and check input information for errors 
    def __init__(self,f1,fout1,isSE,dirout1,concat,f_factors):
        
        # Gather input parameters and create otuput files
        if(isSE): self.pairing='SE'
        else: self.pairing='PE'
        self.f1=f1
        self.f_factors=f_factors
        self.dirout1=dirout1
        self.fout1=fout1
        self.fout2=fout1+".concat.sh"
        self.fout3=fout1+".concatenated.sh"
        self.concat=concat
        self.out1=open(self.fout1,'w')
        if(self.concat):
            self.out2=open(self.fout2,'w')
            self.out3=open(self.fout3,'w')
            
        h1=open(f1,'r')
        h1_=filter(lambda row: row[0]!='#' and row[0]!='' and row[0]!='\t' and row!="\n",h1)
        d1 = csv.DictReader(h1_,delimiter='\t')
        
        # Create a grouping array for files coming from the same sample
        
        # PE: self.g -> {sample1} -> {pair1} -> [partfile1_f,partfile1_r]
        # SE: self.g -> {sample1} -> {pair1} -> [partfile1_f]
        self.g={}
        notFile=[]
        files={}; files_duplicates=[]
        for i,r1 in enumerate(d1):
            for x in ['sample','pair','file']:
                if r1[x] is None:
                    print r1
                    sys.exit("Wrong table format")
            if self.g.get(r1['sample'],"")=="": self.g[r1['sample']]={}
            if self.g[r1['sample']].get(r1['pair'],"")=="": self.g[r1['sample']][r1['pair']]=[]
            self.g[r1['sample']][r1['pair']].append(r1['file'])
            if not(os.path.isfile(r1['file'])): notFile.append(r1['file'])
            if files.get(r1['file'],"")=="": files[r1['file']]=True
            else: files_duplicates.append(r1['file'])
        h1.close()
        
        # PE: self.grp -> {sample1} -> [[partfile1_f,partfile2_f,partfile3_f],[partfile1_r,partfile2_r,partfile3_r]]
        # SE: self.grp -> {sample1} -> [[partfile1_f,partfile2_f,partfile3_f],[]]
        self.grp={}
        for i,group1 in enumerate(self.g.keys()):
            self.grp[group1]=[[],[]]
            for j,pair1 in enumerate(self.g[group1].keys()):
                if    self.pairing == 'SE' and len(self.g[group1][pair1])!=1: sys.exit('Expect SE input! %i'%(len(self.g[group1][pair1])))
                elif  self.pairing == 'PE' and len(self.g[group1][pair1])!=2: sys.exit('Expect PE input! %i'%(len(self.g[group1][pair1])))
                for k in range(len(self.g[group1][pair1])):
                    self.grp[group1][k].append(self.g[group1][pair1][k])
        
        # Check and report errors in the input            
        self.out1.write( '# Comments:\n' )
        self.out1.write( '# Script template created by %s from table %s\n'%(os.path.basename(__file__),self.f1) )
        for nf in notFile:
            self.out1.write( "# Error: didn't find file: %s\n"%(nf) )
        for fd in files_duplicates:
            self.out1.write( "# Error: duplicated file: %s\n"%(fd) )
        
        self.f_factors_ok=False
        if self.f_factors != "":
            self.factors={}
            self.f_factors_ok=True
            
            # self.factors -> {sample1} -> {factor1} -> string
            h2=open(self.f_factors,'r')
            h2_=filter(lambda row: row[0]!='#' and row[0]!='' and row[0]!='\t' and row!="\n",h2)
            d2 = csv.DictReader(h2_,delimiter='\t')
            self.factor_names=d2.fieldnames
            if self.factor_names[0] != 'sample': sys.exit("Found \"%s\" while expected \"sample\" in \"%s\""%(self.factor_names[0],self.f_factors))
            if len(self.factor_names) <2: sys.exit("Expected at least two columns in \"%s\""%(self.f_factors))
            for i in range(len(self.factor_names)):
                if not(re.match(r'^[A-Za-z]+[A-Za-z0-9_]*$',self.factor_names[i])): sys.exit('Error: unexpected symbol in \"factors_table\" column name \"%s\"'%(self.factor_names[i]))
            
            for i,r2 in enumerate(d2):
                if self.factors.get(r2['sample'],"")=="": self.factors[r2['sample']]={}
                else: sys.exit("Sample \"%s\" appears more than once in \"%s\""%(r2['sample'],self.f_factors))
                if self.grp.get(r2['sample'],"")=="":
                    self.out1.write( "# Error: sample \"%s\" is found in the factors input table, but not in the main input table\n"%(r2['sample']) )
                    self.f_factors_ok=False
                for j in range(1,len(self.factor_names)):
                    self.factors[r2['sample']][self.factor_names[j]] = r2[self.factor_names[j]]
            h2.close()
           
            for k in self.grp.keys():
                if self.factors.get(k,"")=="":
                    self.out1.write( "# Error: sample \"%s\" is found in the main input table, but not in the factors input table\n"%(k) )
                    self.f_factors_ok=False
    
    # Close output files    
    def __del__(self):
        self.out1.close()
        if(self.concat):
            self.out2.close()
            self.out3.close()
    
    # Check that all input ".fastq" files have same extension    
    def addExtension(self,files1):
        ext=''
        for i,f in enumerate(files1):
            m1=[    re.match(r'.*(\.fq)$',f),
                    re.match(r'.*(\.fastq)$',f),
                    re.match(r'.*(\.fastq\.gz)$',f),
                    re.match(r'.*(\.fq\.gz)$',f)]
            for m in m1:
                if m:
                    z=m.group(1)
                    if i==0: ext=z
                    elif ext!=z: sys.exit("error: mix of diffferent file extensions: \'%s\' \'%s\'"%(z,ext))
        return ext
    
    # Write input information as a table of bash arrays. if --concat-files - also create bash scripts for concatenating part files
    def write_bash_script(self,delim1=','):
        def o(x):
            return "\"%s\""%(x)
        
        cmd=[]
        cmd_cat=[]
        cmd_concatenated=[]

        for i,g in enumerate(sorted(self.grp.keys())): # Groups
            if self.pairing=='PE':
                if(len(self.grp[g])==2) and (len(self.grp[g][0])>0) and (len(self.grp[g][1])>0) and (len(self.grp[g][0])==len(self.grp[g][1])):
                    
                    # Gathering sample names and their paths for bash array file
                    for1=delim1.join(self.grp[g][0])
                    rev1=delim1.join(self.grp[g][1])
                    cmd1= \
                    "title[{i}]={sample}\nfor1[{i}]={for1}\nrev1[{i}]={rev1}" \
                    .format(for1=o(for1),rev1=o(rev1),i=i,sample=o(g))
                    cmd.append(cmd1)
                    
                    # Creating commands for joining (using cat command) all partial files and lanes, for each sample and pair:
                    for1_cat=' \\\n'.join(map(o,self.grp[g][0]))
                    rev1_cat=' \\\n'.join(map(o,self.grp[g][1]))
                    if len(self.grp[g][0])>1:
                        cmd1_cat= \
                        "# {sample}\ncat \\\n{for1_cat} \\\n   > $dirout_cat/{sample}_f{ext}\nassert_\ncat \\\n{rev1_cat} \\\n   > $dirout_cat/{sample}_r{ext}\nassert_" \
                        .format(for1_cat=for1_cat,rev1_cat=rev1_cat,sample=o(g),ext=self.addExtension(self.grp[g][0]+self.grp[g][1]))
                    else:
                        cmd1_cat= \
                        "# {sample}\ncp  \\\n{for1_cat} \\\n   $dirout_cat/{sample}_f{ext}\nassert_\ncp  \\\n{rev1_cat} \\\n   $dirout_cat/{sample}_r{ext}\nassert_" \
                        .format(for1_cat=for1_cat,rev1_cat=rev1_cat,sample=o(g),ext=self.addExtension(self.grp[g][0]+self.grp[g][1]))
                    cmd_cat.append(cmd1_cat)
                    
                    # Gathering concatenated file output data:
                    cmd1_concatenated = \
                    "title[{i}]=\"{sample}\"\nassert_\nfor1[{i}]=\"{dirout_cat}/{sample}_f{ext}\"\nassert_\nrev1[{i}]=\"{dirout_cat}/{sample}_r{ext}\"\nassert_" \
                    .format(sample=g,ext=self.addExtension(self.grp[g][0]+self.grp[g][1]),i=i,dirout_cat=self.dirout1)
                    cmd_concatenated.append(cmd1_concatenated)
                else:
                    self.out1.write( self.grp[g] +"\n")
                    sys.exit("%s cannot be PE, has %i records"%(g,len(self.grp[g])))
                    
            elif self.pairing=='SE':
                if(len(self.grp[g])==2) and (len(self.grp[g][0])>0) and (len(self.grp[g][1])==0):
                    
                    # Gathering sample names and their paths for bash array file
                    for1=delim1.join(self.grp[g][0])
                    cmd1= \
                    "title[{i}]={sample}\nfor1[{i}]={for1}" \
                    .format(for1=o(for1),i=i,sample=o(g))
                    cmd.append(cmd1)
                    
                    # Creating commands for joining (using cat command) all partial files and lanes, for each sample:
                    for1_cat=" \\\n".join(map(o,self.grp[g][0]))
                    if len(self.grp[g][0])>1:
                        cmd1_cat= \
                        "# {sample}\ncat \\\n{for1_cat} \\\n   > $dirout_cat/{sample}_f{ext}\nassert_" \
                        .format(for1_cat=for1_cat,sample=o(g),ext=self.addExtension(self.grp[g][0]))
                    else:
                        cmd1_cat= \
                        "# {sample}\ncp  \\\n{for1_cat} \\\n   $dirout_cat/{sample}_f{ext}\nassert_" \
                        .format(for1_cat=for1_cat,sample=o(g),ext=self.addExtension(self.grp[g][0]))
                    cmd_cat.append(cmd1_cat)
                    
                    # Gathering concatenated file output data:
                    cmd1_concatenated = \
                    "title[{i}]=\"{sample}\"\nassert_\nfor1[{i}]=\"{dirout_cat}/{sample}_f{ext}\"\nassert_" \
                    .format(sample=g,ext=self.addExtension(self.grp[g][0]),i=i,dirout_cat=self.dirout1)
                    cmd_concatenated.append(cmd1_concatenated)
                    
                else:
                    self.out1.write( self.grp[g] +"\n")
                    sys.exit("%s cannot be SE, has %i records"%(g,len(self.grp[g])))
            
            # Gathering "factors_table" information for bash array file        
            if self.f_factors!="":
                if self.f_factors_ok:
                    if self.factors.get(g,"")!="":
                        for z in range(1,len(self.factor_names)):
                            if(self.factors[g].get(self.factor_names[z],"")!=""):
                                cmd1 = "{name}[{i}]={value}  # {sample}" \
                                .format(name=self.factor_names[z], i=i,value=o(self.factors[g][self.factor_names[z]]),sample=o(g))
                                cmd.append(cmd1)
                            else: sys.exit('Err: missing factors_table information for one of the samples\n')
                    else: sys.exit('Err: missing sample file information\n')
       
        '''
        Writing output files:
        1) out1 - bash array file containing input file and "factors_table" information (latter if --concatenate_partfiles false)
        2) out2 - bash file containing commands for concatenating partial ".fastq" input
        3) out3 - bash array file containing names and locations of the concatenated ".fastq" files
        '''
                   
        self.out1.write( '\n# files:' + "\n" )
        if(self.concat): self.out2.write("dirout_cat=%s\n\n"%(o(self.dirout1)))
        self.out1.write( "function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }" +"\n")
        if(self.concat): self.out2.write( "function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }" +"\n")
        if(self.concat): self.out3.write( "function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }" +"\n")
        for i,c in enumerate(cmd):
            self.out1.write( c +"\n")
        for i,c in enumerate(cmd_cat):   
            if(self.concat): self.out2.write( c +"\n")
        for i,c in enumerate(cmd_concatenated):   
            if(self.concat): self.out3.write( c +"\n")
        self.out1.write( "\n#run:" +"\n")
        if self.pairing == 'PE':
            self.out1.write( "#for i in ${!title[@]}; do\n#\techo $i\n#\techo ${title[i]}\n#\techo ${for1[i]}\n#\techo ${rev1[i]}\n#\techo '----------------'\n#done" +"\n")
            if(self.concat):
                self.out3.write( "#for i in ${!title[@]}; do\n#\techo $i\n#\techo ${title[i]}\n#\techo ${for1[i]}\n#\techo ${rev1[i]}\n#\techo '----------------'\n#done" +"\n")
        elif self.pairing == 'SE':
            self.out1.write( "#for i in ${!title[@]}; do\n#\techo $i\n#\techo ${title[i]}\n#\techo ${for1[i]}\n#\techo '----------------'\n#done" +"\n")
            if(self.concat):
                self.out3.write( "#for i in ${!title[@]}; do\n#\techo $i\n#\techo ${title[i]}\n#\techo ${for1[i]}\n#\techo '----------------'\n#done" +"\n")

# Contents of 'help' output for --input_table and --factors_table parameters
def help_input_table():
    comment1=(
    "REQUIRED, tab-delimited, comment/blank lines are ignored\n"
    "pair	sample	file\n"
    "1      1B4     /home/1B4_R1_001.fastq\n"
    "2      1B4     /home/1B4_R1_002.fastq\n"
    "3      1B4     /home/1B4_R1_003.fastq\n"
    "1      1B4     /home/1B4_R2_001.fastq\n"
    "2      1B4     /home/1B4_R2_002.fastq\n"
    "3      1B4     /home/1B4_R2_003.fastq\n"
    "\n"
    "1      1B5     /home/1B5_R1_001.fastq\n"
    "2      1B5     /home/1B5_R1_002.fastq\n"
    "1      1B5     /home/1B5_R2_001.fastq\n"
    "2      1B5     /home/1B5_R2_002.fastq\n"
    )
    return(comment1)

def help_factor_table():
    comment1=(
        "OPTIONAL, tab-delimited\n"
        "Please see an example of factors_table.txt inside \"design\" folder\n"
    )
    return(comment1)

#### Main ####
parser = argparse.ArgumentParser(description='Function: store a table with SE/PE fastq information, including partfiles if exist, as bash arrays',formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--input_table', action="store", dest="f1",help=help_input_table(),required=True)
parser.add_argument('--output_script', action="store", dest="fout1",required=True,help='REQUIRED: destination to write output')
parser.add_argument('--factors_table', action="store", dest="f_factors",default="",help=help_factor_table())
parser.add_argument('--SE', action='store_true',dest="isSE", default=False,help='OPTIONAL, default: PE, paired-end')
parser.add_argument('--concatenate_partfiles', action='store_true',dest="concat", default=False,help='OPTIONAL: create bash script for partfile concatenation')
parser.add_argument('--output_dir', action="store", dest="dirout1",default="",help="(required with --concatenate_partfiles, to store concatenated files)")
p=parser.parse_args()
if (p.dirout1 == "") and p.concat: parser.error("--output_dir argument is required when using --concatenate_partfiles")

x=mapper_call(p.f1,p.fout1,p.isSE,p.dirout1,p.concat,p.f_factors)
x.write_bash_script()
