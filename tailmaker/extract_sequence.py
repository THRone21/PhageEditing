#!/usr/bin/env python3
dict1={}
import sys
with open(sys.argv[1],'r') as gff:
    for line in gff.readlines():
        line=line.strip().split('\t')
        start_end='{};{}'.format(line[3],line[4])
        if line[-1] != '0' :
            pro_id=line[-1].split(';')[0].split('protein_id=')[-1]
            ptoduct=line[-1].split(';')[1]
        else:
            pro_id='hypothesis_protein_{}_{}'.format(line[3],line[4])
            ptoduct='hypothesis_protein'
        dict1[start_end]=[]
        dict1[start_end].append(pro_id)
        dict1[start_end].append(ptoduct)
#print(dict1)
import re
dict2={}
with open(sys.argv[2],'r') as gene:
    for line in gene.readlines(): 
        line=line.strip()
        if re.match('>',line) :
            ID='{};{}'.format(line.split('#')[1].strip(),line.split('#')[2].strip())
            #if ID in dict1:
            #    print('>{}'.format("\t".join(dict1[ID])))
            #else:
            #    continue
         #       print("ERRO!",line)
            dict2[ID]=[] 
        else:
            dict2[ID].append(line)

for key in dict2:
    if key in dict1:
        print('>{}'.format(";".join(dict1[key])))
        print("".join(dict2[key]))

