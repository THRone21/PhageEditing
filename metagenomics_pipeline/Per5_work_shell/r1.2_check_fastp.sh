cat $1|while read line
do
	ARR=($line)
	name=${ARR[0]}
	cd 01.2_check_fastp
	cat>"r1_check_fastp${name}.sh"<<lv
#seqkit_fastp
/share/app/seqkit/0.14.0-dev/seqkit stats ../${name}/01.1_fastp/${name}_fastp*.fq.gz >${name}_seqkit_stat_fastp.txt
lv
	qsub -cwd -l vf=10g,num_proc=2 -P P17Z10200N0246 -binding linear:2 -q st.q -V r1_check_fastp${name}.sh
#	qsub -cwd -l vf=200g,num_proc=8 -P st_supermem -binding linear:8 -q st_supermem.q r2_kraken_${name}.sh
	cd ..
done
