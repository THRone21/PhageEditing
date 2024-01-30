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
	if [ ! -d "/01.2_SOAPnuke" ];then
		mkdir 01.2_SOAPnuke
	fi
	if [ ! -d "/01.3_prinseq" ];then
		mkdir 01.3_prinseq
	fi
	if [ ! -d "/02.1_minikraken" ];then
		mkdir 02.1_minikraken
	fi
	if [ ! -d "02.2_maxkraken" ];then
		mkdir 02.2_maxkraken
	fi
	if [ ! -d "03.5_megahit_edit_kmer" ];then
		mkdir 03.5_megahit_edit_kmer
	fi
	cat>"r1_filter_${name}.sh"<<lv
#fastp
/ldfssz1/ST_INFECTION/F16ZQSB1SY3030_Salmonella/USER/lindechun/software/fastp/fastp -i $fq1 -o ./01.1_fastp/${name}_fastp1.fq.gz -I $fq2 -O ./01.1_fastp/${name}_fastp2.fq.gz -5 -3 -q 20 -c -j ./01.1_fastp/fastp.json -h ./01.1_fastp/fastp.html -R ./01.1_fastp/out.prefix -l 30 1> 01.1fastp.result 2>r1.1fastp.sh.err
#soapnuke
/ldfssz1/ST_INFECTION/F16ZQSB1SY3030_Salmonella/USER/lindechun/bin/SOAPnuke filter -T 8 -1 ./01.1_fastp/${name}_fastp1.fq.gz  -2 ./01.1_fastp/${name}_fastp2.fq.gz -d -C ${name}_SOAP.1.fq -D ${name}_SOAP.2.fq -o ./01.2_SOAPnuke -Q 2 2>r1.2SOAPnuke.sh.err
#prinseq
/hwfssz5-tmp/ST_INFECTION/GlobalDatabase/share/software/Miniconda3/envs/prinseq++/bin/prinseq++ -threads 8 -fastq ./01.2_SOAPnuke/${name}_SOAP.1.fq.gz -fastq2 ./01.2_SOAPnuke/${name}_SOAP.2.fq.gz -lc_entropy=0.5 -lc_dust=0.5 -out_name ./01.3_prinseq/${name} 2>r1.3prinseq.sh.err
lv
	qsub -cwd -l vf=40g,num_proc=2 -P P17Z10200N0246 -binding linear:2 -q st.q -V r1_filter_${name}.sh
	cd ..
done
