cat $1|while read line
do
	ARR=($line)
	name=${ARR[0]}
	cat>"r1_check_fastp_result${name}.sh"<<lv
#differences_between_1.fastp.result&seqkit_stat.txt
cat ../${name}/01.1fastp.result ./${name}_seqkit_stat_fastp.txt >check_fastp_${name}.txt
lv
	qsub -cwd -l vf=10g,num_proc=2 -P P17Z10200N0246 -binding linear:2 -q st.q -V r1_check_fastp_result${name}.sh
#	qsub -cwd -l vf=200g,num_proc=8 -P st_supermem -binding linear:8 -q st_supermem.q r2_kraken_${name}.sh
done
