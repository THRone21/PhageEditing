dict0={}
dict1={}
dict2={}
with open('20210709/48/03.2_Spades/scaffolds.fasta.humanvirus.outfmt6.txt','r') as f0:
    for line in f0.readlines():
        line=line.strip().split('\t')
        seq_id=line[0].split('_cov')[0]
        if seq_id not in dict0:
            dict0[seq_id]={}
            dict0[seq_id][line[1]]=[]
            dict0[seq_id][line[1]].append(int(line[3])) 
        else:
            if line[1] not in dict0[seq_id]:
                dict0[seq_id][line[1]]=[]
                dict0[seq_id][line[1]].append(int(line[3])) 
            else:
                dict0[seq_id][line[1]].append(int(line[3]))

with open('20210709/48/03.2_Spades/scaffolds.fasta.mammalvirus.outfmt6.txt','r') as f1:
    for line in f1.readlines():
        line=line.strip().split('\t')
        seq_id=line[0].split('_cov')[0]
        if seq_id not in dict1:
            dict1[seq_id]={}
            dict1[seq_id][line[1]]=[]
            dict1[seq_id][line[1]].append(int(line[3]))
        else:
            if line[1] not in dict1[seq_id]:
                dict1[seq_id][line[1]]=[]
                dict1[seq_id][line[1]].append(int(line[3]))
            else:
                dict1[seq_id][line[1]].append(int(line[3]))

with open('20210709/48/03.2_Spades/scaffolds.fasta.vetervirus.outfmt6.txt','r') as f2:
    for line in f2.readlines():
        line=line.strip().split('\t')
        seq_id=line[0].split('_cov')[0]
        if seq_id not in dict2:
            dict2[seq_id]={}
            dict2[seq_id][line[1]]=[]
            dict2[seq_id][line[1]].append(int(line[3]))
        else:
            if line[1] not in dict2[seq_id]:
                dict2[seq_id][line[1]]=[]
                dict2[seq_id][line[1]].append(int(line[3]))
            else:
                dict2[seq_id][line[1]].append(int(line[3]))

all_tag=[]
with open('202107all_id.tag','r') as f4:
    for line in f4.readlines():
        line=line.strip()
        all_tag.append(line)
import re
with open('20210709/48/03.2_Spades/48.id.test.txt','r') as f3:
    for line in f3.readlines():
        line=line.strip()
        print(line,line.split('_')[3],'total length',sep='\t')
        li=line.split('_cov')[0]
        blast_len=0
        if li in dict0:
            for s in dict0[li]:
                for d in range(0,len(dict0[li][s])):
                    print(line,dict0[li][s][d],sep='\t',end='\t')
                    blast_len+=int(dict0[li][s][d])
                    for f in range(0,len(all_tag)):
                        if re.search(s,all_tag[f]):
                            print('{}:humanvirus;{};{}'.format(s,all_tag[f].split(',')[1].strip('"'),all_tag[f].split(',')[0].strip('"')))
        if li in dict1:
            for g in dict1[li]:                    
                for h in range(0,len(dict1[li][g])):
                    print(line,dict1[li][g][h],sep='\t',end='\t')
                    blast_len+=int(dict1[li][g][h])
                    for j in range(0,len(all_tag)):
                        if re.search(g,all_tag[j]):
                            print('{}:mammalvirus;{};{}'.format(g,all_tag[j].split(',')[1].strip('"'),all_tag[j].split(',')[0].strip('"')))
        if li in dict2:
            for z in dict2[li]:
                for x in range(0,len(dict2[li][z])):
                    print(line,dict2[li][z][x],sep='\t',end='\t')
                    blast_len+=int(dict2[li][z][x])
                    for c in range(0,len(all_tag)):
                        if re.search(z,all_tag[c]):
                            print('{}:vetervirus;{};{}'.format(z,all_tag[c].split(',')[1].strip('"'),all_tag[c].split(',')[0].strip('"')))     
        print(line,int(line.split('_')[3])-blast_len,'new length',sep='\t')
