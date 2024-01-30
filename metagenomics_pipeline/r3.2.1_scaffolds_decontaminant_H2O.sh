source /home/liangtianzhu/.bashrc
/hwfssz5/ST_INFECTION/Phage/P17Z10200N0246_Phage_XMF/limin/software/ncbi-blast-2.9.0+/bin/blastn -query ../03.2_Spades/scaffolds.fasta -db /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/liuwei/cryosphere/Permafrost5_20211014/db_blastn_contaminant/H2O_scaffolds -out ./scaffolds.fa_H2O.outfmt6.txt -outfmt 6 -evalue 1e-5 -num_threads 8 -max_target_seqs 1
awk '{print$1}' ./scaffolds.fa_H2O.outfmt6.txt>list
/share/app/seqkit/0.14.0-dev/seqkit grep -v -f list ../03.2_Spades/scaffolds.fasta >12_scaffolds_decontaminant_H2O.fasta 2>r3.2.1_scaffolds_decontaminant_H2O.sh.err
