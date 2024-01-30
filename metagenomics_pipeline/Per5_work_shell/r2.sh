cat $1|while read line
do
	ARR=($line)
	name=${ARR[0]}
	cd $name
	cat>"r2_kraken_${name}.sh"<<lv
#r2.1
source /home/liangtianzhu/.bashrc
/hwfssz5/ST_INFECTION/AntibioticResistance/P18Z10200N0164_Resistance_PN/PMseq/User/liangtianzhu/sortware/karken2/kraken2-2.0.8-beta/k2dir/kraken2  --db /hwfssz5/ST_INFECTION/AntibioticResistance/P18Z10200N0164_Resistance_PN/PMseq/User/liangtianzhu/sortware/karken2/kraken2-2.0.8-beta/db/minikraken2/minikraken2_v1_8GB/ -threads 10 --paired ./01.3_prinseq/${name}_good_out_R1.fastq ./01.3_prinseq/${name}_good_out_R2.fastq --output ./02.1_minikraken/${name}.kraken --report ./02.1_minikraken/${name}.kraken.report 2>r2.1kraken.sh.err
/ldfssz1/ST_META/share/User/zhujie/.conda/envs/bioenv/bin/bracken -d  /hwfssz5/ST_INFECTION/AntibioticResistance/P18Z10200N0164_Resistance_PN/PMseq/User/liangtianzhu/sortware/karken2/kraken2-2.0.8-beta/db/minikraken2/minikraken2_v1_8GB/ -i ./02.1_minikraken/${name}.kraken.report  -o ./02.1_minikraken/${name}.txt -w ./02.1_minikraken/${name}.report -r 100 -l F 2>r2.2brakenF.sh.err
#r2.2
/hwfssz5/ST_INFECTION/AntibioticResistance/P18Z10200N0164_Resistance_PN/PMseq/User/liangtianzhu/sortware/karken2/kraken2-2.0.8-beta/k2dir/kraken2 --db /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/luoyunzhe/database/kraken_max/maxikraken2_1903_140GB/ -threads 10 --paired ./01.3_prinseq/${name}_good_out_R1.fastq ./01.3_prinseq/${name}_good_out_R2.fastq --output 02.2_maxkraken/${name}.kraken --report 02.2_maxkraken/${name}.kraken.report 2>r2.1kraken.sh.err
lv
#	qsub -cwd -l vf=40g,num_proc=4 -P P17Z10200N0246 -binding linear:4 -q st.q -V r2_kraken_${name}.sh
	qsub -cwd -l vf=200g,num_proc=8 -P st_supermem -binding linear:8 -q st_supermem.q r2_kraken_${name}.sh
	cd ..
done
