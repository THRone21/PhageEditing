mkdir /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_11292
cd /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_11292
mkdir ./01.prodigal
prodigal -i /home/data3/ziggy/annotation_pipeline/experiment_seq/KP_Generate_v6_11292.fa -o ./01.prodigal/KP_Generate_v6_11292_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_11292_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_11292_genome.genes.fasta -f gff -c -p meta
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_11292_genome.proteins.faa -db /home/data3/ziggy/annotation_pipeline/db/nr/nr_phage.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_11292_genome.proteins.blastp.txt -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_11292_genome.proteins.faa /home/data3/ziggy/annotation_pipeline/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_11292_genome.proteins.faa  /home/data3/ziggy/annotation_pipeline/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /home/data3/ziggy/annotation_pipeline/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_11292_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/home/data3/ziggy/annotation_pipeline/script/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/home/data3/ziggy/annotation_pipeline/script/get_result.py ./02.blastp/KP_Generate_v6_11292_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_11292_result.txt
mkdir ./06.gff2genebank
/home/data3/ziggy/annotation_pipeline/script/result2gff.py ./05.result/KP_Generate_v6_11292_result.txt ./01.prodigal/KP_Generate_v6_11292_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_11292_genome.result.gff
seqret -sequence /home/data3/ziggy/annotation_pipeline/experiment_seq/KP_Generate_v6_11292.fa -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_11292_genome.result.gff -osformat genbank -outseq /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_11292/KP_Generate_v6_11292_genome.gbk 
cp  /home/data3/ziggy/annotation_pipeline/experiment_seq/KP_Generate_v6_11292.fa /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_11292/KP_Generate_v6_11292_genome.gbk ./06.gff2genebank/KP_Generate_v6_11292_genome.result.gff ./05.result/KP_Generate_v6_11292_result.txt ./01.prodigal/KP_Generate_v6_11292_genome.proteins.faa ./01.prodigal/KP_Generate_v6_11292_genome.genes.fasta /home/data3/ziggy/annotation_pipeline/experiment_output/final_result 
echo /home/data3/ziggy/annotation_pipeline/experiment_output/final_result/KP_Generate_v6_11292 >>/home/data3/ziggy/annotation_pipeline/experiment_output/fa.list  
mkdir /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_12421
cd /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_12421
mkdir ./01.prodigal
prodigal -i /home/data3/ziggy/annotation_pipeline/experiment_seq/KP_Generate_v6_12421.fa -o ./01.prodigal/KP_Generate_v6_12421_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_12421_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_12421_genome.genes.fasta -f gff -c -p meta
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_12421_genome.proteins.faa -db /home/data3/ziggy/annotation_pipeline/db/nr/nr_phage.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_12421_genome.proteins.blastp.txt -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_12421_genome.proteins.faa /home/data3/ziggy/annotation_pipeline/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_12421_genome.proteins.faa  /home/data3/ziggy/annotation_pipeline/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /home/data3/ziggy/annotation_pipeline/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_12421_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/home/data3/ziggy/annotation_pipeline/script/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/home/data3/ziggy/annotation_pipeline/script/get_result.py ./02.blastp/KP_Generate_v6_12421_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_12421_result.txt
mkdir ./06.gff2genebank
/home/data3/ziggy/annotation_pipeline/script/result2gff.py ./05.result/KP_Generate_v6_12421_result.txt ./01.prodigal/KP_Generate_v6_12421_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_12421_genome.result.gff
seqret -sequence /home/data3/ziggy/annotation_pipeline/experiment_seq/KP_Generate_v6_12421.fa -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_12421_genome.result.gff -osformat genbank -outseq /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_12421/KP_Generate_v6_12421_genome.gbk 
cp  /home/data3/ziggy/annotation_pipeline/experiment_seq/KP_Generate_v6_12421.fa /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_12421/KP_Generate_v6_12421_genome.gbk ./06.gff2genebank/KP_Generate_v6_12421_genome.result.gff ./05.result/KP_Generate_v6_12421_result.txt ./01.prodigal/KP_Generate_v6_12421_genome.proteins.faa ./01.prodigal/KP_Generate_v6_12421_genome.genes.fasta /home/data3/ziggy/annotation_pipeline/experiment_output/final_result 
echo /home/data3/ziggy/annotation_pipeline/experiment_output/final_result/KP_Generate_v6_12421 >>/home/data3/ziggy/annotation_pipeline/experiment_output/fa.list  
mkdir /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_13944
cd /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_13944
mkdir ./01.prodigal
prodigal -i /home/data3/ziggy/annotation_pipeline/experiment_seq/KP_Generate_v6_13944.fa -o ./01.prodigal/KP_Generate_v6_13944_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_13944_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_13944_genome.genes.fasta -f gff -c -p meta
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_13944_genome.proteins.faa -db /home/data3/ziggy/annotation_pipeline/db/nr/nr_phage.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_13944_genome.proteins.blastp.txt -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_13944_genome.proteins.faa /home/data3/ziggy/annotation_pipeline/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_13944_genome.proteins.faa  /home/data3/ziggy/annotation_pipeline/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /home/data3/ziggy/annotation_pipeline/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_13944_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/home/data3/ziggy/annotation_pipeline/script/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/home/data3/ziggy/annotation_pipeline/script/get_result.py ./02.blastp/KP_Generate_v6_13944_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_13944_result.txt
mkdir ./06.gff2genebank
/home/data3/ziggy/annotation_pipeline/script/result2gff.py ./05.result/KP_Generate_v6_13944_result.txt ./01.prodigal/KP_Generate_v6_13944_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_13944_genome.result.gff
seqret -sequence /home/data3/ziggy/annotation_pipeline/experiment_seq/KP_Generate_v6_13944.fa -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_13944_genome.result.gff -osformat genbank -outseq /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_13944/KP_Generate_v6_13944_genome.gbk 
cp  /home/data3/ziggy/annotation_pipeline/experiment_seq/KP_Generate_v6_13944.fa /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_13944/KP_Generate_v6_13944_genome.gbk ./06.gff2genebank/KP_Generate_v6_13944_genome.result.gff ./05.result/KP_Generate_v6_13944_result.txt ./01.prodigal/KP_Generate_v6_13944_genome.proteins.faa ./01.prodigal/KP_Generate_v6_13944_genome.genes.fasta /home/data3/ziggy/annotation_pipeline/experiment_output/final_result 
echo /home/data3/ziggy/annotation_pipeline/experiment_output/final_result/KP_Generate_v6_13944 >>/home/data3/ziggy/annotation_pipeline/experiment_output/fa.list  
mkdir /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_17710
cd /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_17710
mkdir ./01.prodigal
prodigal -i /home/data3/ziggy/annotation_pipeline/experiment_seq/KP_Generate_v6_17710.fa -o ./01.prodigal/KP_Generate_v6_17710_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_17710_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_17710_genome.genes.fasta -f gff -c -p meta
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_17710_genome.proteins.faa -db /home/data3/ziggy/annotation_pipeline/db/nr/nr_phage.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_17710_genome.proteins.blastp.txt -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_17710_genome.proteins.faa /home/data3/ziggy/annotation_pipeline/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_17710_genome.proteins.faa  /home/data3/ziggy/annotation_pipeline/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /home/data3/ziggy/annotation_pipeline/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_17710_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/home/data3/ziggy/annotation_pipeline/script/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/home/data3/ziggy/annotation_pipeline/script/get_result.py ./02.blastp/KP_Generate_v6_17710_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_17710_result.txt
mkdir ./06.gff2genebank
/home/data3/ziggy/annotation_pipeline/script/result2gff.py ./05.result/KP_Generate_v6_17710_result.txt ./01.prodigal/KP_Generate_v6_17710_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_17710_genome.result.gff
seqret -sequence /home/data3/ziggy/annotation_pipeline/experiment_seq/KP_Generate_v6_17710.fa -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_17710_genome.result.gff -osformat genbank -outseq /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_17710/KP_Generate_v6_17710_genome.gbk 
cp  /home/data3/ziggy/annotation_pipeline/experiment_seq/KP_Generate_v6_17710.fa /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_17710/KP_Generate_v6_17710_genome.gbk ./06.gff2genebank/KP_Generate_v6_17710_genome.result.gff ./05.result/KP_Generate_v6_17710_result.txt ./01.prodigal/KP_Generate_v6_17710_genome.proteins.faa ./01.prodigal/KP_Generate_v6_17710_genome.genes.fasta /home/data3/ziggy/annotation_pipeline/experiment_output/final_result 
echo /home/data3/ziggy/annotation_pipeline/experiment_output/final_result/KP_Generate_v6_17710 >>/home/data3/ziggy/annotation_pipeline/experiment_output/fa.list  
mkdir /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_17888
cd /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_17888
mkdir ./01.prodigal
prodigal -i /home/data3/ziggy/annotation_pipeline/experiment_seq/KP_Generate_v6_17888.fa -o ./01.prodigal/KP_Generate_v6_17888_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_17888_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_17888_genome.genes.fasta -f gff -c -p meta
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_17888_genome.proteins.faa -db /home/data3/ziggy/annotation_pipeline/db/nr/nr_phage.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_17888_genome.proteins.blastp.txt -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_17888_genome.proteins.faa /home/data3/ziggy/annotation_pipeline/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_17888_genome.proteins.faa  /home/data3/ziggy/annotation_pipeline/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /home/data3/ziggy/annotation_pipeline/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_17888_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/home/data3/ziggy/annotation_pipeline/script/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/home/data3/ziggy/annotation_pipeline/script/get_result.py ./02.blastp/KP_Generate_v6_17888_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_17888_result.txt
mkdir ./06.gff2genebank
/home/data3/ziggy/annotation_pipeline/script/result2gff.py ./05.result/KP_Generate_v6_17888_result.txt ./01.prodigal/KP_Generate_v6_17888_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_17888_genome.result.gff
seqret -sequence /home/data3/ziggy/annotation_pipeline/experiment_seq/KP_Generate_v6_17888.fa -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_17888_genome.result.gff -osformat genbank -outseq /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_17888/KP_Generate_v6_17888_genome.gbk 
cp  /home/data3/ziggy/annotation_pipeline/experiment_seq/KP_Generate_v6_17888.fa /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_17888/KP_Generate_v6_17888_genome.gbk ./06.gff2genebank/KP_Generate_v6_17888_genome.result.gff ./05.result/KP_Generate_v6_17888_result.txt ./01.prodigal/KP_Generate_v6_17888_genome.proteins.faa ./01.prodigal/KP_Generate_v6_17888_genome.genes.fasta /home/data3/ziggy/annotation_pipeline/experiment_output/final_result 
echo /home/data3/ziggy/annotation_pipeline/experiment_output/final_result/KP_Generate_v6_17888 >>/home/data3/ziggy/annotation_pipeline/experiment_output/fa.list  
mkdir /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_20860
cd /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_20860
mkdir ./01.prodigal
prodigal -i /home/data3/ziggy/annotation_pipeline/experiment_seq/KP_Generate_v6_20860.fa -o ./01.prodigal/KP_Generate_v6_20860_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_20860_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_20860_genome.genes.fasta -f gff -c -p meta
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_20860_genome.proteins.faa -db /home/data3/ziggy/annotation_pipeline/db/nr/nr_phage.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_20860_genome.proteins.blastp.txt -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_20860_genome.proteins.faa /home/data3/ziggy/annotation_pipeline/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_20860_genome.proteins.faa  /home/data3/ziggy/annotation_pipeline/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /home/data3/ziggy/annotation_pipeline/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_20860_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/home/data3/ziggy/annotation_pipeline/script/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/home/data3/ziggy/annotation_pipeline/script/get_result.py ./02.blastp/KP_Generate_v6_20860_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_20860_result.txt
mkdir ./06.gff2genebank
/home/data3/ziggy/annotation_pipeline/script/result2gff.py ./05.result/KP_Generate_v6_20860_result.txt ./01.prodigal/KP_Generate_v6_20860_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_20860_genome.result.gff
seqret -sequence /home/data3/ziggy/annotation_pipeline/experiment_seq/KP_Generate_v6_20860.fa -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_20860_genome.result.gff -osformat genbank -outseq /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_20860/KP_Generate_v6_20860_genome.gbk 
cp  /home/data3/ziggy/annotation_pipeline/experiment_seq/KP_Generate_v6_20860.fa /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_20860/KP_Generate_v6_20860_genome.gbk ./06.gff2genebank/KP_Generate_v6_20860_genome.result.gff ./05.result/KP_Generate_v6_20860_result.txt ./01.prodigal/KP_Generate_v6_20860_genome.proteins.faa ./01.prodigal/KP_Generate_v6_20860_genome.genes.fasta /home/data3/ziggy/annotation_pipeline/experiment_output/final_result 
echo /home/data3/ziggy/annotation_pipeline/experiment_output/final_result/KP_Generate_v6_20860 >>/home/data3/ziggy/annotation_pipeline/experiment_output/fa.list  
mkdir /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_24504
cd /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_24504
mkdir ./01.prodigal
prodigal -i /home/data3/ziggy/annotation_pipeline/experiment_seq/KP_Generate_v6_24504.fa -o ./01.prodigal/KP_Generate_v6_24504_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_24504_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_24504_genome.genes.fasta -f gff -c -p meta
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_24504_genome.proteins.faa -db /home/data3/ziggy/annotation_pipeline/db/nr/nr_phage.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_24504_genome.proteins.blastp.txt -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_24504_genome.proteins.faa /home/data3/ziggy/annotation_pipeline/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_24504_genome.proteins.faa  /home/data3/ziggy/annotation_pipeline/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /home/data3/ziggy/annotation_pipeline/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_24504_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/home/data3/ziggy/annotation_pipeline/script/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/home/data3/ziggy/annotation_pipeline/script/get_result.py ./02.blastp/KP_Generate_v6_24504_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_24504_result.txt
mkdir ./06.gff2genebank
/home/data3/ziggy/annotation_pipeline/script/result2gff.py ./05.result/KP_Generate_v6_24504_result.txt ./01.prodigal/KP_Generate_v6_24504_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_24504_genome.result.gff
seqret -sequence /home/data3/ziggy/annotation_pipeline/experiment_seq/KP_Generate_v6_24504.fa -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_24504_genome.result.gff -osformat genbank -outseq /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_24504/KP_Generate_v6_24504_genome.gbk 
cp  /home/data3/ziggy/annotation_pipeline/experiment_seq/KP_Generate_v6_24504.fa /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_24504/KP_Generate_v6_24504_genome.gbk ./06.gff2genebank/KP_Generate_v6_24504_genome.result.gff ./05.result/KP_Generate_v6_24504_result.txt ./01.prodigal/KP_Generate_v6_24504_genome.proteins.faa ./01.prodigal/KP_Generate_v6_24504_genome.genes.fasta /home/data3/ziggy/annotation_pipeline/experiment_output/final_result 
echo /home/data3/ziggy/annotation_pipeline/experiment_output/final_result/KP_Generate_v6_24504 >>/home/data3/ziggy/annotation_pipeline/experiment_output/fa.list  
mkdir /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_5638
cd /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_5638
mkdir ./01.prodigal
prodigal -i /home/data3/ziggy/annotation_pipeline/experiment_seq/KP_Generate_v6_5638.fa -o ./01.prodigal/KP_Generate_v6_5638_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_5638_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_5638_genome.genes.fasta -f gff -c -p meta
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_5638_genome.proteins.faa -db /home/data3/ziggy/annotation_pipeline/db/nr/nr_phage.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_5638_genome.proteins.blastp.txt -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_5638_genome.proteins.faa /home/data3/ziggy/annotation_pipeline/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_5638_genome.proteins.faa  /home/data3/ziggy/annotation_pipeline/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /home/data3/ziggy/annotation_pipeline/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_5638_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/home/data3/ziggy/annotation_pipeline/script/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/home/data3/ziggy/annotation_pipeline/script/get_result.py ./02.blastp/KP_Generate_v6_5638_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_5638_result.txt
mkdir ./06.gff2genebank
/home/data3/ziggy/annotation_pipeline/script/result2gff.py ./05.result/KP_Generate_v6_5638_result.txt ./01.prodigal/KP_Generate_v6_5638_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_5638_genome.result.gff
seqret -sequence /home/data3/ziggy/annotation_pipeline/experiment_seq/KP_Generate_v6_5638.fa -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_5638_genome.result.gff -osformat genbank -outseq /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_5638/KP_Generate_v6_5638_genome.gbk 
cp  /home/data3/ziggy/annotation_pipeline/experiment_seq/KP_Generate_v6_5638.fa /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_5638/KP_Generate_v6_5638_genome.gbk ./06.gff2genebank/KP_Generate_v6_5638_genome.result.gff ./05.result/KP_Generate_v6_5638_result.txt ./01.prodigal/KP_Generate_v6_5638_genome.proteins.faa ./01.prodigal/KP_Generate_v6_5638_genome.genes.fasta /home/data3/ziggy/annotation_pipeline/experiment_output/final_result 
echo /home/data3/ziggy/annotation_pipeline/experiment_output/final_result/KP_Generate_v6_5638 >>/home/data3/ziggy/annotation_pipeline/experiment_output/fa.list  
mkdir /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_6935
cd /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_6935
mkdir ./01.prodigal
prodigal -i /home/data3/ziggy/annotation_pipeline/experiment_seq/KP_Generate_v6_6935.fa -o ./01.prodigal/KP_Generate_v6_6935_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_6935_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_6935_genome.genes.fasta -f gff -c -p meta
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_6935_genome.proteins.faa -db /home/data3/ziggy/annotation_pipeline/db/nr/nr_phage.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_6935_genome.proteins.blastp.txt -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_6935_genome.proteins.faa /home/data3/ziggy/annotation_pipeline/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_6935_genome.proteins.faa  /home/data3/ziggy/annotation_pipeline/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /home/data3/ziggy/annotation_pipeline/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_6935_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/home/data3/ziggy/annotation_pipeline/script/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/home/data3/ziggy/annotation_pipeline/script/get_result.py ./02.blastp/KP_Generate_v6_6935_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_6935_result.txt
mkdir ./06.gff2genebank
/home/data3/ziggy/annotation_pipeline/script/result2gff.py ./05.result/KP_Generate_v6_6935_result.txt ./01.prodigal/KP_Generate_v6_6935_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_6935_genome.result.gff
seqret -sequence /home/data3/ziggy/annotation_pipeline/experiment_seq/KP_Generate_v6_6935.fa -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_6935_genome.result.gff -osformat genbank -outseq /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_6935/KP_Generate_v6_6935_genome.gbk 
cp  /home/data3/ziggy/annotation_pipeline/experiment_seq/KP_Generate_v6_6935.fa /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_6935/KP_Generate_v6_6935_genome.gbk ./06.gff2genebank/KP_Generate_v6_6935_genome.result.gff ./05.result/KP_Generate_v6_6935_result.txt ./01.prodigal/KP_Generate_v6_6935_genome.proteins.faa ./01.prodigal/KP_Generate_v6_6935_genome.genes.fasta /home/data3/ziggy/annotation_pipeline/experiment_output/final_result 
echo /home/data3/ziggy/annotation_pipeline/experiment_output/final_result/KP_Generate_v6_6935 >>/home/data3/ziggy/annotation_pipeline/experiment_output/fa.list  
mkdir /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_7294
cd /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_7294
mkdir ./01.prodigal
prodigal -i /home/data3/ziggy/annotation_pipeline/experiment_seq/KP_Generate_v6_7294.fa -o ./01.prodigal/KP_Generate_v6_7294_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_7294_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_7294_genome.genes.fasta -f gff -c -p meta
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_7294_genome.proteins.faa -db /home/data3/ziggy/annotation_pipeline/db/nr/nr_phage.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_7294_genome.proteins.blastp.txt -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_7294_genome.proteins.faa /home/data3/ziggy/annotation_pipeline/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_7294_genome.proteins.faa  /home/data3/ziggy/annotation_pipeline/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /home/data3/ziggy/annotation_pipeline/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_7294_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/home/data3/ziggy/annotation_pipeline/script/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/home/data3/ziggy/annotation_pipeline/script/get_result.py ./02.blastp/KP_Generate_v6_7294_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_7294_result.txt
mkdir ./06.gff2genebank
/home/data3/ziggy/annotation_pipeline/script/result2gff.py ./05.result/KP_Generate_v6_7294_result.txt ./01.prodigal/KP_Generate_v6_7294_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_7294_genome.result.gff
seqret -sequence /home/data3/ziggy/annotation_pipeline/experiment_seq/KP_Generate_v6_7294.fa -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_7294_genome.result.gff -osformat genbank -outseq /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_7294/KP_Generate_v6_7294_genome.gbk 
cp  /home/data3/ziggy/annotation_pipeline/experiment_seq/KP_Generate_v6_7294.fa /home/data3/ziggy/annotation_pipeline/experiment_output/KP_Generate_v6_7294/KP_Generate_v6_7294_genome.gbk ./06.gff2genebank/KP_Generate_v6_7294_genome.result.gff ./05.result/KP_Generate_v6_7294_result.txt ./01.prodigal/KP_Generate_v6_7294_genome.proteins.faa ./01.prodigal/KP_Generate_v6_7294_genome.genes.fasta /home/data3/ziggy/annotation_pipeline/experiment_output/final_result 
echo /home/data3/ziggy/annotation_pipeline/experiment_output/final_result/KP_Generate_v6_7294 >>/home/data3/ziggy/annotation_pipeline/experiment_output/fa.list  
python3 /home/data3/ziggy/annotation_pipeline/script/zhengli.py /home/data3/ziggy/annotation_pipeline/experiment_output/fa.list /home/data3/ziggy/annotation_pipeline/experiment_output/final_stat.txt 
