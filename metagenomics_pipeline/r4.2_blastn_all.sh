/hwfssz5/ST_INFECTION/Phage/P17Z10200N0246_Phage_XMF/limin/software/ncbi-blast-2.9.0+/bin/blastn -query ../../03.2_Spades/scaffolds.fasta -db /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/liuwei/cryosphere/db_test/db_virus_host_blastn/Virus-host-dephage-db/Virus_host_dephage_db -out ./spades.scaffolds.fasta.virus.outfmt6.txt -outfmt 6 -evalue 1e-5 -num_threads 8 -max_target_seqs 1 &&\
cat ./*outfmt6.txt |awk '{print $1}' |sort -n |uniq >allid.txt &&\
cat ./*outfmt6.txt |awk '{print $2}' |sort -n |uniq >48.id &&\
grep -f 48.id /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/liuwei/cryosphere/Permafrost3/202107_06_bc43_pooling/04.4_blastn_Virus-host-dephage-db/04.4.1_blastn_spades/virus_host_dephage_db.header2.csv>id.tag &&\
python3 r4.6spades_plot.py
