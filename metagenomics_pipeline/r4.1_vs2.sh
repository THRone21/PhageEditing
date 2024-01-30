path=/zfssz6/BIGDATA_WAREHOUSE/bigdata_ods/for_cngb/cnsa/upload/liuwei8/cryosphere/Permafrost_29
cat $1|while read line
do
        ARR=($line)
        name=${ARR[0]}
        cd $name
        if [ ! -d "/04.1_vs2" ];then
                mkdir 04.1_vs2
        fi
	cd $path/$name/04.1_vs2
        cat>"j${name}.r4.1_vs2.sh"<<lv
conda activate vs2
/ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/liuwei/software/miniconda3/envs/vs2/bin/virsorter run -w result -i $path/$name/03.2_metaSpades_good_out/scaffolds.fasta --include-groups "dsDNAphage,ssDNA,lavidaviridae,NCLDV" --use-conda-off -j 4 all --latency-wait 60
lv
#        nohup sh j${name}.r4.1_vs2.sh
	qsub -cwd -l vf=1g,num_proc=4 -P P17Z10200N0246 -binding linear:4 -q st.q -V j${name}.r4.1_vs2.sh
done
