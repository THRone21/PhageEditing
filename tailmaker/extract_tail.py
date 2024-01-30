#!/usr/bin/env python3
import re
#host={}
#with open('phage_host.list','r') as f0:
#    for line in f0.readlines():
#        line=line.strip().split('\t')
#        host[line[0]]=line[1]
f2=open('tail_fiber.fasta','w')
f3=open('tail_fiber_protein.fasta','w')
with open('pp_3.list','r') as f1:
    for line in f1.readlines():
        line=line.strip().split('\t')
        gff,fna,faa=line[0],line[1],line[2]
        fna_f=open(fna,'r')
        dict_fna={}
        for LINE in fna_f.readlines():
            LINE=LINE.strip()
            if re.match('>',LINE):
                ID=LINE.split(' # ')[1].strip()+'#'+LINE.split(' # ')[2].strip()
                dict_fna[ID]=[]
            else:
                dict_fna[ID].append(LINE)
        fna_f.close()
        faa_f=open(faa,'r')
        dict_faa={}
        for LINE in faa_f.readlines():
            LINE=LINE.strip()
            if re.match('>',LINE):
                ID=LINE.split(' # ')[1].strip()+'#'+LINE.split(' # ')[2].strip()
                dict_faa[ID]=[]
            else:
                dict_faa[ID].append(LINE)
        faa_f.close()
        gff_f=open(gff,'r')
        for LINE in gff_f.readlines():
            LINE=LINE.strip()
            new_id=LINE.split('\t')[3].strip() +'#'+LINE.split('\t')[4].strip()
            if re.search('tail fiber',LINE.split('\t')[-1]):
#                new_id=LINE.split('\t')[3].strip() +'#'+LINE.split('\t')[4].strip()
                print('>{}_{}_{}'.format(gff.split('_')[0],new_id,LINE.split('\t')[-1]),file=f2)
                print('>{}_{}_{}'.format(gff.split('_')[0],new_id,LINE.split('\t')[-1]),file=f3)
                print("".join(dict_fna[new_id]),file=f2)
                print("".join(dict_faa[new_id]),file=f3)
            elif re.search('spike',LINE.split('\t')[-1]):
                print('>{}_{}_{}'.format(gff.split('_')[0],new_id,LINE.split('\t')[-1]),file=f2)
                print('>{}_{}_{}'.format(gff.split('_')[0],new_id,LINE.split('\t')[-1]),file=f3)
                print("".join(dict_fna[new_id]),file=f2)
                print("".join(dict_faa[new_id]),file=f3)
        gff_f.close() 
