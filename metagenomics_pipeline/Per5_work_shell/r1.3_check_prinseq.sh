path=/ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/liuwei/cryosphere/Permafrost5_20211014/work
cat $1|while read line
do
	ARR=($line)
	name=${ARR[0]}
	cd 01.3_check_prinseq
	cat>"r1_check_prinseq${name}.sh"<<lv
#seqkit_soapnuke
/share/app/seqkit/0.14.0-dev/seqkit stats ${path}/${name}/01.3_prinseq/*good_out* >${name}_seqkit_stat.txt
lv
	qsub -cwd -l vf=10g,num_proc=2 -P P17Z10200N0246 -binding linear:2 -q st.q -V r1_check_prinseq${name}.sh
#	qsub -cwd -l vf=200g,num_proc=8 -P st_supermem -binding linear:8 -q st_supermem.q r2_kraken_${name}.sh
#	cd ..
done
