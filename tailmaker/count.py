import re
dict1={}
dict2={}
list1=[]
with open('fiber_spike_depolymerase_uniq_1774_noreid.fa','r') as f1:
    for line in f1.readlines():
        line=line.strip()
        #print(line)
        if re.match('>',line) == None:
            seq=line
            #print(line)
            dict1[seq]=ID
        else:
            ID=line.split(' ')[0].split('>')[-1]
            list1.append(ID)
            dict2[ID]={}
#print(len(dict1.keys()))
#list1=[]
#for k in dict1.keys():
#    list1.append(dict1[k])
 #   print(dict1[k])
 #   print(k)
with open('gene.list','r') as f2:
    for L in f2.readlines():
        L=L.strip()
        f3=open(L,'r')
        num=0
        for n in f3.readlines():
            n=n.strip()
            if re.match('>',n) == None:
                if n in dict1:
                    ID2=dict1[n]
                    
                    dict2[ID2][L.split('/')[-1].split('_genome.genes.fa')[0]]=1
        f3.close()  
#list3=sorted(list1)
list2=sorted(dict2.keys())
print(len(list1),len(list2))
#for n in range(0,1774):
#    if n < 1769:
#        print(list3[n],list2[n])
#    else:
#       print(list3[n])
#print(list1)
#print(list2)
#ret=list(set(list1).difference(set(list2)))
#print(ret)
import pandas as pd
import numpy as np
from pandas import DataFrame

frame2 = DataFrame(dict2)
#pd.options.display.max_columns = None
#pd.options.display.max_rows = None

#np.set_printoptions(threshold=np.inf)
#显示所有列
#pd.set_option('display.max_columns', None)
#显示所有行
#pd.set_option('display.max_rows', None)
#设置value的显示长度为100，默认为50
#pd.set_option('max_colwidth',100)

#print(frame2)
#with pd.option_context('display.max_rows', None, 'display.max_columns', 100):
#    print(frame2)
print(frame2.to_string())
#for k in dict2.keys():                 
#    print(dict2[k])
