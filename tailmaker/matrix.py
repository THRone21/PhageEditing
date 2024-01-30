#!/usr/bin/env python3
dict1={}
bac=[]
with open('full_matrix_129x96.txt','r') as f1:
    for line in f1.readlines():
        line=line.strip().split(',') 
        bac.append(line[2])
        if line[0] in dict1:
            dict1[line[0]][line[2]]=line[-1]
        else:
            dict1[line[0]]={}
            dict1[line[0]][line[2]]=line[-1]

#print(dict1['PHKP4049'])
uniq_bac=sorted(list(set(bac)))
print('          ',end='\t')
for n in range(0,len(uniq_bac)):
    if n != len(uniq_bac) -1 :
        print(uniq_bac[n],end='\t')
    else:
        print(uniq_bac[n],end='\n')
for key in sorted(dict1.keys()):
    print(key,end='\t')
    for n in range(0,len(uniq_bac)):
        if n != len(uniq_bac) -1 :
            print(dict1[key][uniq_bac[n]],end='\t')
        else:
            print(dict1[key][uniq_bac[n]],end='\n')
