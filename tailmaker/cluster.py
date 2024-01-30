#!/usr/bin/env python3
import re
dict1={}
with open('all_tail_changeID.protein.fa_70.clstr','r') as f1:
    for line in f1.readlines():
        line=line.strip()
        if re.match('>',line):
            ID=line.split('>')[-1]
            dict1[ID]=[] 
        else:
            phage=line.split(', >')[-1].split('...')[0]
            dict1[ID].append(phage)

#for key in sorted(dict1.keys()):
#    print(key,dict1[key])
search=[]
dict2={}
with open('phage_id.txt','r') as f2:
    for line in f2.readlines():
        line=line.strip()
        search.append(line)

for k in search:
    for key in sorted(dict1.keys()):
        for n in dict1[key]:
            if re.match(k,n):
                if k in dict2:
                    dict2[k].append(key)
                else:
                    dict2[k]=[]
                    dict2[k].append(key)
for l in sorted(dict1.keys()):
    if l != sorted(dict1.keys())[-1]:
        print(l,end='\t')
    else:
        print(l)
for k in search:
    
    print(k,end='\t')
    for l in sorted(dict1.keys()):
        if l != sorted(dict1.keys())[-1]:
            if l in dict2[k]:
                print('1',end='\t')
            else:
                print('0',end='\t')
        else:
            if l in dict2[k]:
                print('1',end='\n')
            else:
                print('0',end='\n')
