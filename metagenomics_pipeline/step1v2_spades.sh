##48.id.test.txt代表的是比对上3个病毒库的全部contig id
cat /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/limin/test/20210709/48/03.2_Spades/*outfmt6.txt |awk '{print $1}' |sort -n |uniq

##202107all_id.tag这个则是将当时0707、0709的样本比对上的全部病毒的具体信息给放在一起方便处理,为什么要取出来，因为整个大表太大了，取出来程序会比较快
##比对上的所有病毒id
cat /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/limin/test/2021070*/*/03.2_Spades/*outfmt6.txt |awk '{print $2}' |sort -n |uniq >/ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/limin/test/202107.id
##从大表里面取出这些病毒id号的具体信息
grep -f /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/limin/test/202107.id /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/limin/test/blast/viruses.csv >/ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/limin/test/202107all_id.tag


