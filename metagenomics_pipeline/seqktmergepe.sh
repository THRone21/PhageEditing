list=(10 28)
path=/ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/liuwei/cryosphere/Permafrost5_20211014/
for i in ${list[@]}
do
	name=$i
	mkdir $path/$name/01.4_blastn_decontaminant_H2O
	cd $path/$name/01.4_blastn_decontaminant_H2O
	/share/app/seqtk/1.3/seqtk mergepe $path/$name/01.3_prinseq/${name}_good_out_R1.fastq $path/$name/01.3_prinseq/${name}_good_out_R2.fastq>${name}.fastq
	/share/app/seqkit/0.14.0-dev/seqkit fq2fa ${name}.fastq -o ${name}.fasta
	/hwfssz5/ST_INFECTION/Phage/P17Z10200N0246_Phage_XMF/limin/software/ncbi-blast-2.9.0+/bin/blastn -query ${name}.fasta -db /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/liuwei/cryosphere/Permafrost5_20211014/db_blastn_contaminant/H2O_reads -out ./good_out.fa_H2O.outfmt6.txt -outfmt 6 -evalue 1e-5 -num_threads 8 -max_target_seqs 1
	awk '{print $1}' ./good_out.fa_H2O.outfmt6.txt>list
	/share/app/seqkit/0.14.0-dev/seqkit grep -v -f list ${name}.fasta >${name}_good_out_decontaminant_H2O.fasta
done
