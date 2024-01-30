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
	if [ ! -d "/01.3_prinseq" ];then
                mkdir 01.3_prinseq
        fi
	cat>"r1_prinseq_${name}.sh"<<lv
#prinseq
/hwfssz5-tmp/ST_INFECTION/GlobalDatabase/share/software/Miniconda3/envs/prinseq++/bin/prinseq++ -threads 8 -fastq ./01.2_SOAPnuke/${name}_SOAP.1.fq.gz -fastq2 ./01.2_SOAPnuke/${name}_SOAP.2.fq.gz -lc_entropy=0.5 -lc_dust=0.5 -out_name ./01.3_prinseq/${name}
lv
	qsub -cwd -l vf=40g,num_proc=2 -P P17Z10200N0246 -binding linear:2 -q st.q -V r1_prinseq_${name}.sh
	cd ..
done
