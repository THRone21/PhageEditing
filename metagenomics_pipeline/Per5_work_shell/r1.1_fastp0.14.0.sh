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
	if [ ! -d "/01.1_fastp" ];then
		mkdir 01.1_fastp
	fi
	cat>"r1_fastp_${name}.sh"<<lv
#fastp
source /home/liangtianzhu/.bashrc
/ldfssz1/ST_INFECTION/F16ZQSB1SY3030_Salmonella/USER/lindechun/software/fastp/fastp -i $fq1 -o ./01.1_fastp/${name}_fastp1.fq.gz -I $fq2 -O ./01.1_fastp/${name}_fastp2.fq.gz -5 -3 -q 20 -c -j ./01.1_fastp/fastp.json -h ./01.1_fastp/fastp.html -R ./01.1_fastp/out.prefix -l 30 1> 01.1fastp.result
lv
	qsub -cwd -l vf=40g,num_proc=2 -P P17Z10200N0246 -binding linear:2 -q st.q -V r1_fastp_${name}.sh
	cd ..
done
