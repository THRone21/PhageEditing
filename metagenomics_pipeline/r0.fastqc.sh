#与sample.list配合使用，命令sh r0.1.sh sample.list
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
	if [ ! -d "00.1_fastqc" ];then
		mkdir 00.1_fastqc
	fi
	cat>"r0.1_fastqc_${name}.sh"<<lv
#fastqc
/ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/limin/software/FastQC/fastqc -o 00.1_fastqc -t 4 -f fastq $fq1 $fq2 
lv
	qsub -cwd -l vf=40g,num_proc=2 -P P17Z10200N0246 -binding linear:2 -q st.q -V r0.1_fastqc_${name}.sh
	cd ..
done
