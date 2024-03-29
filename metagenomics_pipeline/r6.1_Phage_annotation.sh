for i in $(ls /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/liuwei/cryosphere/Permafrost3/202107_06_bc48_pooling/04.2_blastn_phage/spades_phage_blastn/phage.scaffolds.fasta)
	do
		A=`echo $(basename $i .fa)`
		##目的路径 防止任务跑断 故判断续跑检查用
		address="/ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/liuwei/cryosphere/Permafrost3/202107_06_bc48_pooling/06.1_annotation$A"
		if [ -d $address ];
		then
		    echo "$adress had finished" >>annotation.log
		else
			#mkdir $adress
			echo "mkdir $address
chgrp -R Phage $address
chmod g+s $address
cd $address
mkdir ./01.prodigal
/hwfssz5/ST_INFECTION/Phage/P17Z10200N0246_Phage_XMF/share/02.annotation/software/Prodigal/prodigal -i $i -o "./01.prodigal/$A\_genome.prodigal.gff" -a "./01.prodigal/$A\_genome.proteins.faa" -d "./01.prodigal/$A\_genome.genes.fasta" -f gff -c -p meta
mkdir ./02.blastp
/hwfssz5/ST_INFECTION/Phage/P17Z10200N0246_Phage_XMF/limin/software/ncbi-blast-2.9.0+/bin/blastp -query "./01.prodigal/$A\_genome.proteins.faa" -db /hwfssz5/ST_INFECTION/Phage/P17Z10200N0246_Phage_XMF/share/02.annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out "./02.blastp/$A\_genome.proteins.blastp.txt" -outfmt \"6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore\"
mkdir ./03.phmmer
/hwfssz5/ST_INFECTION/Phage/P17Z10200N0246_Phage_XMF/share/02.annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt "./01.prodigal/$A\_genome.proteins.faa" /hwfssz5/ST_INFECTION/Phage/P17Z10200N0246_Phage_XMF/share/02.annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_INFECTION/Phage/P17Z10200N0246_Phage_XMF/share/02.annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt "./01.prodigal/$A\_genome.proteins.faa"  /hwfssz5/ST_INFECTION/Phage/P17Z10200N0246_Phage_XMF/share/02.annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_INFECTION/Phage/P17Z10200N0246_Phage_XMF/share/02.annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_INFECTION/Phage/P17Z10200N0246_Phage_XMF/share/02.annotation/db/Pfam/Pfam-A.hmm "./01.prodigal/$A\_genome.proteins.faa"
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_INFECTION/Phage/P17Z10200N0246_Phage_XMF/share/02.annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_INFECTION/Phage/P17Z10200N0246_Phage_XMF/share/02.annotation/software/get_result.py "./02.blastp/$A\_genome.proteins.blastp.txt" ./05.result/hmm.hit.uniq "./05.result/$A\_result.txt"
mkdir ./06.gff2genebank
/hwfssz5/ST_INFECTION/Phage/P17Z10200N0246_Phage_XMF/share/02.annotation/software/result2gff.py "./05.result/$A\_result.txt" "./01.prodigal/$A\_genome.prodigal.gff" "./06.gff2genebank/$A\_genome.result.gff"
seqret -sequence $i -feature -fformat gff -fopenfile "./06.gff2genebank/$A\_genome.result.gff" -osformat genbank -outseq "$address/$A\_genome.gbk" ">>run_annotation.sh
         fi
   done
#cp "./06.gff2genebank/$A\_genome.result.gb" /hwfssz5/ST_INFECTION/Phage/P17Z10200N0246_Phage_XMF/songwenchen/Denovo/phage_xiangya/annotation/
#/hwfssz5/ST_INFECTION/Phage/P17Z10200N0246_Phage_XMF/share/02.annotation/software/bcbb/gff/Scripts/gff/gff_to_genbank.py "./06.gff2genebank/$A\_genome.result.gff" ./spades/scaffolds.fasta
