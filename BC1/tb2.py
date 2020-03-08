import copy
import collections
import itertools
import sys

class tables:
    def __init__(self):
        pass
    
    def merge(self,x,y,ncolx0,ncoly0,NA=0,x_notEmpty=True):
        # assuming the first column in each table x and y are the keys
        
        '''
        example of using this function:
        
        from tb2 import tables
        t=tables()
        x=[['x',1,2],['y',3,4],['z',5,6]]
        y=[['z',11,21,61],['x',31,41,21],['w',51,51,11],['q',01,31,11]]
        t.merge(x,y,3,4)
        
        output:
        [[['y', 3, 4, None, None, None, None]], ['x', 1, 2, 31, 41, 21], ['z', 5, 6, 11, 21, 61], [['q', None, None, 1, 31, 11]], [['w', None, None, 51, 51, 11]]]
        '''
       
        if (len(x) > 0) and (len(y) > 0):
            x_=copy.deepcopy(x) # [[x1,x2,..xn],[x1,x2,..xn],...]
            y_=copy.deepcopy(y) # [[y1,y2,..yn],[y1,y2,..yn],...]
            for i in range(len(x_)):
                if len(x_[i]) != ncolx0: sys.exit("err1: expected constant column length in x\n")
            for i in range(len(y_)):
                if len(y_[i]) != ncoly0: sys.exit("err2: expected constant column length in y\n")
            if ncolx0==ncoly0:
                u=1
                for i in range(len(x_)): x_[i].append(0)
            else: u=0
            
            result=collections.defaultdict(list)
            for item in itertools.chain(x_,y_):
                result[item[0]].append(item)
            r=[list(itertools.chain.from_iterable(value)) for value in result.values()] # array of: [key1,val1a,val1b..,key2,val2a,val2b]
            r2=[]
            for i in range(len(r)):
                if len(r[i])==(ncolx0+ncoly0+u): # [key1,val1a,val1b,...,key2,val2a,val2b]
                    r2.append(r[i][:ncolx0]+r[i][(ncolx0+1+u):])
                elif len(r[i])==ncolx0+u:
                    r2.append(r[i][:ncolx0]+[NA for j in range(ncoly0-1)])
                elif len(r[i])==ncoly0:
                    r2.append([r[i][0]]+[NA for j in range(ncolx0-1)]+r[i][1:])
                else:
                    print r[i]
                    sys.exit('err3: merge\n')
                
            return(r2)
        elif (len(x) > 0):
            x_=copy.deepcopy(x) # [[x1,x2,..xn],[x1,x2,..xn],...]
            for i in range(len(x_)):
                if len(x_[i]) != ncolx0: sys.exit("err1: expected constant column length in x\n")
            for i in range(len(x_)):
                for j in range(ncoly0-1): x_[i].append(NA)
            return x_
        elif (len(y) > 0):
            if x_notEmpty: sys.exit("err4: x table must not be empty\n")
            y_=copy.deepcopy(y) # [[y1,y2,..yn],[y1,y2,..yn],...]
            for i in range(len(y_)):
                if len(y_[i]) != ncoly0: sys.exit("err2: expected constant column length in y\n")
            y__ = []
            for i in range(len(y_)):
                x1 = [y_[i][0]]
                x3 = y_[i][1:]
                x2 = [ NA for j in range(ncolx0-1)]
                y__.append(x1+x2+x3)
            return(y__)
        else:
            return []

class tb:
    def __init__(self,cols):
        # cols should be a hash table of names and indices, that looks like: {'column_name_1':0,'column_name_2':1,'column_name_3':2, ...}
        self.cols=copy.deepcopy(cols)
        
    def ne(self,x,col,val1):
        y=[]
        for i in range(len(x)):
            if x[i][self.cols[col]] != val1:
                y.append(x[i])
        return(y)
    
    def eq(self,x,col,val1):
        y=[]
        for i in range(len(x)):
            if x[i][self.cols[col]] == val1:
                y.append(x[i])
        return(y)
    
    def sel(self,x,col):
        return([x[z][self.cols[col]] for z in range(len(x))])
    
    def ins(self,x,col,vals):
        self.cols[col]=len(self.cols)
        for i in range(len(x)):
            x[i].append(vals[i])
        return(x)
    
    def sort(self,x,col):
        y=self.sel(x,col)
        idxs=sorted(range(len(y)),key=y.__getitem__)
        z=[]
        for i in idxs:
           z.append(x[i])
        return(z)
