cat $1|while read line
do
	ARR=($line)
	name=${ARR[0]}
	cat>"r1_check_soapnuke_basic${name}.sh"<<lv
#seqkit_soapnuke
#/share/app/seqkit/0.14.0-dev/seqkit stats ../${name}/01.2_SOAPnuke/${name}_SOAP.*.fq.gz >${name}_seqkit_stat.txt
cat ../${name}/01.2_SOAPnuke/Basic_Statistics_of_Sequencing_Quality.txt ${name}_seqkit_stat.txt >check_soapnuke_${name}.txt
lv
	qsub -cwd -l vf=10g,num_proc=2 -P P17Z10200N0246 -binding linear:2 -q st.q -V r1_check_soapnuke_basic${name}.sh
#	qsub -cwd -l vf=200g,num_proc=8 -P st_supermem -binding linear:8 -q st_supermem.q r2_kraken_${name}.sh
done
