#mkdir 03.1_megahit_kmer_min
/ldfssz1/ST_META/share/User/zhujie/.conda/envs/bioenv/bin/megahit -1 ../01.3_prinseq/43_good_out_R1.fastq -2 ../01.3_prinseq/43_good_out_R2.fastq -o 03.1_megahit_kmer_step --k-min 21 --k-max 141 --k-step 2
