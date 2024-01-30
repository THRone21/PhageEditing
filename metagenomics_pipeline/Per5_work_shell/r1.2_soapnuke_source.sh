cat $1|while read line
do
	ARR=($line)
	name=${ARR[0]}
	path=${ARR[1]}
	fq1=`find $path -name '*1.fq.gz'|awk 'NR==1 {print}'`
	fq2=${fq1/1.fq.gz/2.fq.gz}
	if [ ! -d $name ];then 
		mkdir $name
	fi
	cd $name
	if [ ! -d "/01.2_SOAPnuke" ];then
		mkdir 01.2_SOAPnuke
	fi
	cat>"r1_soapnuke_${name}.sh"<<lv
#soapnuke
source /home/liangtianzhu/.bashrc
/ldfssz1/ST_INFECTION/F16ZQSB1SY3030_Salmonella/USER/lindechun/bin/SOAPnuke filter -T 8 -1 ./01.1_fastp/${name}_fastp1.fq.gz  -2 ./01.1_fastp/${name}_fastp2.fq.gz -d -C ${name}_SOAP.1.fq -D ${name}_SOAP.2.fq -o ./01.2_SOAPnuke -Q 2
lv
	qsub -cwd -l vf=2g,num_proc=8 -P P17Z10200N0246 -binding linear:8 -q st.q -V r1_soapnuke_${name}.sh
	cd ..
done
