#!/bin/sh
for filename in $(ls *.fa)
do
    blastn -query $filename -strand both -db KP_247phage -out ${filename}_result.txt -outfmt 6 -num_threads 12 -dust no -evalue 1e-5 -max_target_seqs 5;
done
