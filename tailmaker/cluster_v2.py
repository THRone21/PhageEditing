#!/usr/bin/env python3
import re
dict1={}
search=[]
with open('all_tail_changeID_v2.protein.fa_70.clstr','r') as f1:
    for line in f1.readlines():
        line=line.strip()
        if re.match('>',line):
            ID=line.split('>')[-1]
            dict1[ID]=[]
        else:
            phage=line.split(', >')[-1].split(';')[0]
            dict1[ID].append(phage)
            search.append(phage)

for l in sorted(dict1.keys()):
    if l != sorted(dict1.keys())[-1]:
        print(l,end='\t')
    else:
        print(l)
for k in search:

    print(k,end='\t')
    for l in sorted(dict1.keys()):
        if l != sorted(dict1.keys())[-1]:
            if k in dict1[l]:
                print('1',end='\t')
            else:
                print('0',end='\t')
        else:
            if k in dict1[l]:
                print('1',end='\n')
            else:
                print('0',end='\n')
