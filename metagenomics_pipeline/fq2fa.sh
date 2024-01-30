/share/app/seqtk/1.3/seqtk mergepe /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/liuwei/cryosphere/Permafrost5_20211014/H2O_1/01.3_prinseq/H2O_1_good_out_R1.fastq /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/liuwei/cryosphere/Permafrost5_20211014/H2O_1/01.3_prinseq/H2O_1_good_out_R2.fastq >./H2O_1.fq
/share/app/seqtk/1.3/seqtk mergepe /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/liuwei/cryosphere/Permafrost5_20211014/H2O_2/01.3_prinseq/H2O_2_good_out_R1.fastq /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/liuwei/cryosphere/Permafrost5_20211014/H2O_2/01.3_prinseq/H2O_2_good_out_R2.fastq >./H2O_2.fq
cat *.fq >H2O.fastq
/share/app/seqkit/0.14.0-dev/seqkit fq2fa H2O.fastq -o H2O.fasta
