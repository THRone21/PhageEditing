source /home/liangtianzhu/.bashrc
/share/app/bowtie2/2.4.1/bowtie2 --very-sensitive -x H2O_reads -1 ../01.3_prinseq/12_good_out_R1.fastq -2 ../01.3_prinseq/12_good_out_R2.fastq 2> log.bowtie2_log | /share/app/samtools/1.11/bin/samtools fastq -N -c 5 -f 12 -F 256 -1 output.rm_r1.fastq -2 output.rm_r2.fastq
