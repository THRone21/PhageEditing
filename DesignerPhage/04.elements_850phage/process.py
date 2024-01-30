
dict_length = {}
with open('genome_length.txt','r') as f1:
    for line in f1.readlines():
        line=line.strip().split('\t')
        dict_length[line[0]] = line[1]

with open('blastn_kp_850_filter.txt','r') as f2:
    for line in f2.readlines():
        line=line.strip().split('\t')
        print(line[0],line[1],line[2],dict_length[line[0]],dict_length[line[1]],int(line[2])/((int(dict_length[line[0]])+int(dict_length[line[1]]))/2)*100,sep='\t')
