awk '{print $1}' spades.scaffolds.fasta.virus.outfmt6.txt >query.id.txt
seqkit grep -n -f query.id.txt ../../03.2_Spades/scaffolds.fasta >scaffolds.allvirus.fasta
