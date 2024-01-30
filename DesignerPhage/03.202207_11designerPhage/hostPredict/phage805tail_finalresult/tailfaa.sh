#!/usr/bin/bash
for file in /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/limin/phage_Assembly/guizheng/kp_phage_GPH/terminase/repeat_joining/NCBI_KP/KP_phage_annotation/final_result/*.faa
    do
    fileID=`echo $file |awk -F / '{print $NF}'|awk -F . '{print $1"."$2"."$3}'|awk -F "_genome" '{print $NR}'`
    seqkit grep -f ./${fileID}_result_tail.list /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/limin/phage_Assembly/guizheng/kp_phage_GPH/terminase/repeat_joining/NCBI_KP/KP_phage_annotation/final_result/${fileID}_genome.proteins.faa >> ./${fileID}_tail.faa
done