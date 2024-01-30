mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_10064
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_10064
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_10064.fasta -o ./01.prodigal/KP_Generate_v6_10064_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_10064_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_10064_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_10064_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_10064_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_10064_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_10064_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_10064_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_10064_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_10064_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_10064_result.txt ./01.prodigal/KP_Generate_v6_10064_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_10064_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_10064.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_10064_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_10064/KP_Generate_v6_10064_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_10963
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_10963
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_10963.fasta -o ./01.prodigal/KP_Generate_v6_10963_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_10963_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_10963_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_10963_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_10963_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_10963_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_10963_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_10963_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_10963_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_10963_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_10963_result.txt ./01.prodigal/KP_Generate_v6_10963_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_10963_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_10963.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_10963_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_10963/KP_Generate_v6_10963_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_11097
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_11097
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_11097.fasta -o ./01.prodigal/KP_Generate_v6_11097_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_11097_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_11097_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_11097_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_11097_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_11097_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_11097_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_11097_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_11097_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_11097_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_11097_result.txt ./01.prodigal/KP_Generate_v6_11097_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_11097_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_11097.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_11097_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_11097/KP_Generate_v6_11097_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_11292
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_11292
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_11292.fasta -o ./01.prodigal/KP_Generate_v6_11292_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_11292_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_11292_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_11292_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_11292_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_11292_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_11292_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_11292_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_11292_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_11292_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_11292_result.txt ./01.prodigal/KP_Generate_v6_11292_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_11292_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_11292.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_11292_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_11292/KP_Generate_v6_11292_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_11349
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_11349
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_11349.fasta -o ./01.prodigal/KP_Generate_v6_11349_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_11349_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_11349_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_11349_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_11349_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_11349_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_11349_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_11349_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_11349_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_11349_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_11349_result.txt ./01.prodigal/KP_Generate_v6_11349_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_11349_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_11349.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_11349_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_11349/KP_Generate_v6_11349_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_12316
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_12316
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_12316.fasta -o ./01.prodigal/KP_Generate_v6_12316_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_12316_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_12316_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_12316_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_12316_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_12316_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_12316_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_12316_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_12316_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_12316_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_12316_result.txt ./01.prodigal/KP_Generate_v6_12316_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_12316_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_12316.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_12316_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_12316/KP_Generate_v6_12316_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_12421
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_12421
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_12421.fasta -o ./01.prodigal/KP_Generate_v6_12421_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_12421_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_12421_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_12421_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_12421_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_12421_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_12421_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_12421_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_12421_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_12421_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_12421_result.txt ./01.prodigal/KP_Generate_v6_12421_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_12421_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_12421.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_12421_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_12421/KP_Generate_v6_12421_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_13320
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_13320
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_13320.fasta -o ./01.prodigal/KP_Generate_v6_13320_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_13320_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_13320_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_13320_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_13320_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_13320_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_13320_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_13320_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_13320_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_13320_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_13320_result.txt ./01.prodigal/KP_Generate_v6_13320_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_13320_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_13320.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_13320_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_13320/KP_Generate_v6_13320_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_13944
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_13944
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_13944.fasta -o ./01.prodigal/KP_Generate_v6_13944_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_13944_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_13944_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_13944_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_13944_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_13944_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_13944_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_13944_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_13944_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_13944_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_13944_result.txt ./01.prodigal/KP_Generate_v6_13944_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_13944_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_13944.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_13944_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_13944/KP_Generate_v6_13944_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_14426
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_14426
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_14426.fasta -o ./01.prodigal/KP_Generate_v6_14426_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_14426_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_14426_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_14426_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_14426_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_14426_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_14426_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_14426_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_14426_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_14426_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_14426_result.txt ./01.prodigal/KP_Generate_v6_14426_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_14426_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_14426.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_14426_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_14426/KP_Generate_v6_14426_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_147
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_147
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_147.fasta -o ./01.prodigal/KP_Generate_v6_147_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_147_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_147_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_147_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_147_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_147_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_147_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_147_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_147_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_147_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_147_result.txt ./01.prodigal/KP_Generate_v6_147_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_147_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_147.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_147_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_147/KP_Generate_v6_147_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_15638
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_15638
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_15638.fasta -o ./01.prodigal/KP_Generate_v6_15638_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_15638_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_15638_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_15638_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_15638_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_15638_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_15638_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_15638_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_15638_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_15638_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_15638_result.txt ./01.prodigal/KP_Generate_v6_15638_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_15638_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_15638.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_15638_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_15638/KP_Generate_v6_15638_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_16190
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_16190
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_16190.fasta -o ./01.prodigal/KP_Generate_v6_16190_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_16190_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_16190_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_16190_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_16190_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_16190_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_16190_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_16190_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_16190_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_16190_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_16190_result.txt ./01.prodigal/KP_Generate_v6_16190_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_16190_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_16190.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_16190_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_16190/KP_Generate_v6_16190_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_16347
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_16347
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_16347.fasta -o ./01.prodigal/KP_Generate_v6_16347_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_16347_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_16347_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_16347_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_16347_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_16347_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_16347_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_16347_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_16347_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_16347_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_16347_result.txt ./01.prodigal/KP_Generate_v6_16347_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_16347_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_16347.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_16347_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_16347/KP_Generate_v6_16347_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_16530
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_16530
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_16530.fasta -o ./01.prodigal/KP_Generate_v6_16530_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_16530_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_16530_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_16530_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_16530_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_16530_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_16530_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_16530_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_16530_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_16530_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_16530_result.txt ./01.prodigal/KP_Generate_v6_16530_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_16530_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_16530.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_16530_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_16530/KP_Generate_v6_16530_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_16676
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_16676
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_16676.fasta -o ./01.prodigal/KP_Generate_v6_16676_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_16676_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_16676_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_16676_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_16676_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_16676_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_16676_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_16676_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_16676_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_16676_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_16676_result.txt ./01.prodigal/KP_Generate_v6_16676_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_16676_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_16676.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_16676_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_16676/KP_Generate_v6_16676_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_16697
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_16697
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_16697.fasta -o ./01.prodigal/KP_Generate_v6_16697_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_16697_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_16697_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_16697_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_16697_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_16697_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_16697_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_16697_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_16697_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_16697_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_16697_result.txt ./01.prodigal/KP_Generate_v6_16697_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_16697_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_16697.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_16697_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_16697/KP_Generate_v6_16697_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_17189
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_17189
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_17189.fasta -o ./01.prodigal/KP_Generate_v6_17189_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_17189_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_17189_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_17189_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_17189_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_17189_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_17189_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_17189_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_17189_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_17189_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_17189_result.txt ./01.prodigal/KP_Generate_v6_17189_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_17189_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_17189.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_17189_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_17189/KP_Generate_v6_17189_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_17309
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_17309
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_17309.fasta -o ./01.prodigal/KP_Generate_v6_17309_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_17309_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_17309_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_17309_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_17309_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_17309_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_17309_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_17309_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_17309_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_17309_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_17309_result.txt ./01.prodigal/KP_Generate_v6_17309_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_17309_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_17309.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_17309_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_17309/KP_Generate_v6_17309_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_17485
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_17485
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_17485.fasta -o ./01.prodigal/KP_Generate_v6_17485_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_17485_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_17485_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_17485_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_17485_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_17485_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_17485_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_17485_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_17485_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_17485_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_17485_result.txt ./01.prodigal/KP_Generate_v6_17485_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_17485_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_17485.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_17485_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_17485/KP_Generate_v6_17485_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_17507
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_17507
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_17507.fasta -o ./01.prodigal/KP_Generate_v6_17507_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_17507_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_17507_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_17507_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_17507_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_17507_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_17507_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_17507_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_17507_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_17507_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_17507_result.txt ./01.prodigal/KP_Generate_v6_17507_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_17507_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_17507.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_17507_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_17507/KP_Generate_v6_17507_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_17710
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_17710
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_17710.fasta -o ./01.prodigal/KP_Generate_v6_17710_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_17710_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_17710_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_17710_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_17710_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_17710_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_17710_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_17710_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_17710_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_17710_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_17710_result.txt ./01.prodigal/KP_Generate_v6_17710_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_17710_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_17710.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_17710_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_17710/KP_Generate_v6_17710_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_17714
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_17714
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_17714.fasta -o ./01.prodigal/KP_Generate_v6_17714_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_17714_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_17714_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_17714_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_17714_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_17714_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_17714_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_17714_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_17714_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_17714_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_17714_result.txt ./01.prodigal/KP_Generate_v6_17714_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_17714_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_17714.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_17714_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_17714/KP_Generate_v6_17714_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_18153
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_18153
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_18153.fasta -o ./01.prodigal/KP_Generate_v6_18153_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_18153_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_18153_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_18153_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_18153_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_18153_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_18153_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_18153_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_18153_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_18153_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_18153_result.txt ./01.prodigal/KP_Generate_v6_18153_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_18153_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_18153.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_18153_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_18153/KP_Generate_v6_18153_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_18396
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_18396
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_18396.fasta -o ./01.prodigal/KP_Generate_v6_18396_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_18396_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_18396_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_18396_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_18396_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_18396_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_18396_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_18396_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_18396_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_18396_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_18396_result.txt ./01.prodigal/KP_Generate_v6_18396_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_18396_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_18396.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_18396_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_18396/KP_Generate_v6_18396_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_18558
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_18558
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_18558.fasta -o ./01.prodigal/KP_Generate_v6_18558_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_18558_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_18558_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_18558_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_18558_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_18558_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_18558_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_18558_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_18558_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_18558_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_18558_result.txt ./01.prodigal/KP_Generate_v6_18558_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_18558_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_18558.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_18558_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_18558/KP_Generate_v6_18558_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_19330
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_19330
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_19330.fasta -o ./01.prodigal/KP_Generate_v6_19330_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_19330_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_19330_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_19330_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_19330_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_19330_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_19330_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_19330_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_19330_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_19330_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_19330_result.txt ./01.prodigal/KP_Generate_v6_19330_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_19330_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_19330.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_19330_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_19330/KP_Generate_v6_19330_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_19393
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_19393
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_19393.fasta -o ./01.prodigal/KP_Generate_v6_19393_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_19393_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_19393_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_19393_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_19393_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_19393_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_19393_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_19393_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_19393_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_19393_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_19393_result.txt ./01.prodigal/KP_Generate_v6_19393_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_19393_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_19393.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_19393_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_19393/KP_Generate_v6_19393_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_19805
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_19805
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_19805.fasta -o ./01.prodigal/KP_Generate_v6_19805_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_19805_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_19805_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_19805_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_19805_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_19805_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_19805_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_19805_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_19805_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_19805_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_19805_result.txt ./01.prodigal/KP_Generate_v6_19805_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_19805_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_19805.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_19805_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_19805/KP_Generate_v6_19805_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_20325
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_20325
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_20325.fasta -o ./01.prodigal/KP_Generate_v6_20325_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_20325_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_20325_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_20325_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_20325_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_20325_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_20325_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_20325_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_20325_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_20325_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_20325_result.txt ./01.prodigal/KP_Generate_v6_20325_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_20325_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_20325.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_20325_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_20325/KP_Generate_v6_20325_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_20478
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_20478
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_20478.fasta -o ./01.prodigal/KP_Generate_v6_20478_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_20478_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_20478_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_20478_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_20478_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_20478_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_20478_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_20478_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_20478_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_20478_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_20478_result.txt ./01.prodigal/KP_Generate_v6_20478_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_20478_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_20478.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_20478_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_20478/KP_Generate_v6_20478_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_2064
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_2064
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_2064.fasta -o ./01.prodigal/KP_Generate_v6_2064_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_2064_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_2064_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_2064_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_2064_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_2064_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_2064_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_2064_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_2064_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_2064_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_2064_result.txt ./01.prodigal/KP_Generate_v6_2064_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_2064_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_2064.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_2064_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_2064/KP_Generate_v6_2064_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_20705
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_20705
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_20705.fasta -o ./01.prodigal/KP_Generate_v6_20705_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_20705_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_20705_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_20705_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_20705_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_20705_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_20705_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_20705_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_20705_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_20705_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_20705_result.txt ./01.prodigal/KP_Generate_v6_20705_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_20705_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_20705.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_20705_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_20705/KP_Generate_v6_20705_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_20860
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_20860
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_20860.fasta -o ./01.prodigal/KP_Generate_v6_20860_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_20860_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_20860_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_20860_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_20860_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_20860_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_20860_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_20860_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_20860_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_20860_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_20860_result.txt ./01.prodigal/KP_Generate_v6_20860_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_20860_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_20860.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_20860_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_20860/KP_Generate_v6_20860_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_20979
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_20979
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_20979.fasta -o ./01.prodigal/KP_Generate_v6_20979_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_20979_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_20979_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_20979_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_20979_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_20979_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_20979_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_20979_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_20979_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_20979_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_20979_result.txt ./01.prodigal/KP_Generate_v6_20979_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_20979_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_20979.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_20979_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_20979/KP_Generate_v6_20979_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_21092
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_21092
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_21092.fasta -o ./01.prodigal/KP_Generate_v6_21092_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_21092_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_21092_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_21092_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_21092_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_21092_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_21092_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_21092_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_21092_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_21092_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_21092_result.txt ./01.prodigal/KP_Generate_v6_21092_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_21092_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_21092.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_21092_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_21092/KP_Generate_v6_21092_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_21333
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_21333
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_21333.fasta -o ./01.prodigal/KP_Generate_v6_21333_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_21333_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_21333_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_21333_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_21333_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_21333_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_21333_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_21333_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_21333_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_21333_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_21333_result.txt ./01.prodigal/KP_Generate_v6_21333_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_21333_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_21333.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_21333_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_21333/KP_Generate_v6_21333_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_21940
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_21940
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_21940.fasta -o ./01.prodigal/KP_Generate_v6_21940_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_21940_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_21940_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_21940_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_21940_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_21940_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_21940_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_21940_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_21940_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_21940_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_21940_result.txt ./01.prodigal/KP_Generate_v6_21940_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_21940_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_21940.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_21940_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_21940/KP_Generate_v6_21940_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_22439
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_22439
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_22439.fasta -o ./01.prodigal/KP_Generate_v6_22439_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_22439_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_22439_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_22439_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_22439_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_22439_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_22439_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_22439_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_22439_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_22439_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_22439_result.txt ./01.prodigal/KP_Generate_v6_22439_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_22439_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_22439.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_22439_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_22439/KP_Generate_v6_22439_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_22730
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_22730
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_22730.fasta -o ./01.prodigal/KP_Generate_v6_22730_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_22730_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_22730_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_22730_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_22730_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_22730_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_22730_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_22730_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_22730_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_22730_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_22730_result.txt ./01.prodigal/KP_Generate_v6_22730_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_22730_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_22730.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_22730_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_22730/KP_Generate_v6_22730_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_22859
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_22859
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_22859.fasta -o ./01.prodigal/KP_Generate_v6_22859_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_22859_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_22859_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_22859_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_22859_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_22859_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_22859_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_22859_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_22859_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_22859_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_22859_result.txt ./01.prodigal/KP_Generate_v6_22859_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_22859_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_22859.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_22859_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_22859/KP_Generate_v6_22859_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_23442
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_23442
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_23442.fasta -o ./01.prodigal/KP_Generate_v6_23442_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_23442_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_23442_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_23442_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_23442_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_23442_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_23442_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_23442_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_23442_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_23442_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_23442_result.txt ./01.prodigal/KP_Generate_v6_23442_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_23442_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_23442.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_23442_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_23442/KP_Generate_v6_23442_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_23594
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_23594
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_23594.fasta -o ./01.prodigal/KP_Generate_v6_23594_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_23594_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_23594_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_23594_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_23594_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_23594_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_23594_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_23594_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_23594_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_23594_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_23594_result.txt ./01.prodigal/KP_Generate_v6_23594_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_23594_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_23594.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_23594_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_23594/KP_Generate_v6_23594_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_23907
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_23907
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_23907.fasta -o ./01.prodigal/KP_Generate_v6_23907_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_23907_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_23907_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_23907_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_23907_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_23907_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_23907_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_23907_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_23907_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_23907_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_23907_result.txt ./01.prodigal/KP_Generate_v6_23907_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_23907_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_23907.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_23907_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_23907/KP_Generate_v6_23907_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_24010
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_24010
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_24010.fasta -o ./01.prodigal/KP_Generate_v6_24010_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_24010_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_24010_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_24010_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_24010_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_24010_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_24010_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_24010_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_24010_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_24010_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_24010_result.txt ./01.prodigal/KP_Generate_v6_24010_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_24010_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_24010.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_24010_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_24010/KP_Generate_v6_24010_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_24132
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_24132
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_24132.fasta -o ./01.prodigal/KP_Generate_v6_24132_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_24132_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_24132_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_24132_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_24132_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_24132_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_24132_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_24132_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_24132_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_24132_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_24132_result.txt ./01.prodigal/KP_Generate_v6_24132_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_24132_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_24132.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_24132_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_24132/KP_Generate_v6_24132_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_24222
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_24222
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_24222.fasta -o ./01.prodigal/KP_Generate_v6_24222_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_24222_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_24222_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_24222_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_24222_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_24222_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_24222_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_24222_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_24222_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_24222_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_24222_result.txt ./01.prodigal/KP_Generate_v6_24222_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_24222_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_24222.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_24222_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_24222/KP_Generate_v6_24222_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_24504
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_24504
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_24504.fasta -o ./01.prodigal/KP_Generate_v6_24504_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_24504_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_24504_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_24504_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_24504_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_24504_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_24504_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_24504_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_24504_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_24504_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_24504_result.txt ./01.prodigal/KP_Generate_v6_24504_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_24504_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_24504.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_24504_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_24504/KP_Generate_v6_24504_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_2566
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_2566
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_2566.fasta -o ./01.prodigal/KP_Generate_v6_2566_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_2566_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_2566_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_2566_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_2566_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_2566_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_2566_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_2566_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_2566_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_2566_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_2566_result.txt ./01.prodigal/KP_Generate_v6_2566_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_2566_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_2566.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_2566_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_2566/KP_Generate_v6_2566_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_2849
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_2849
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_2849.fasta -o ./01.prodigal/KP_Generate_v6_2849_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_2849_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_2849_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_2849_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_2849_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_2849_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_2849_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_2849_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_2849_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_2849_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_2849_result.txt ./01.prodigal/KP_Generate_v6_2849_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_2849_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_2849.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_2849_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_2849/KP_Generate_v6_2849_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_2894
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_2894
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_2894.fasta -o ./01.prodigal/KP_Generate_v6_2894_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_2894_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_2894_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_2894_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_2894_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_2894_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_2894_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_2894_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_2894_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_2894_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_2894_result.txt ./01.prodigal/KP_Generate_v6_2894_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_2894_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_2894.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_2894_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_2894/KP_Generate_v6_2894_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_2901
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_2901
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_2901.fasta -o ./01.prodigal/KP_Generate_v6_2901_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_2901_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_2901_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_2901_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_2901_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_2901_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_2901_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_2901_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_2901_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_2901_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_2901_result.txt ./01.prodigal/KP_Generate_v6_2901_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_2901_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_2901.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_2901_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_2901/KP_Generate_v6_2901_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_3227
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_3227
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_3227.fasta -o ./01.prodigal/KP_Generate_v6_3227_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_3227_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_3227_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_3227_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_3227_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_3227_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_3227_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_3227_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_3227_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_3227_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_3227_result.txt ./01.prodigal/KP_Generate_v6_3227_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_3227_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_3227.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_3227_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_3227/KP_Generate_v6_3227_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_3429
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_3429
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_3429.fasta -o ./01.prodigal/KP_Generate_v6_3429_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_3429_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_3429_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_3429_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_3429_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_3429_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_3429_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_3429_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_3429_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_3429_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_3429_result.txt ./01.prodigal/KP_Generate_v6_3429_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_3429_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_3429.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_3429_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_3429/KP_Generate_v6_3429_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_4016
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_4016
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_4016.fasta -o ./01.prodigal/KP_Generate_v6_4016_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_4016_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_4016_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_4016_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_4016_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_4016_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_4016_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_4016_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_4016_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_4016_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_4016_result.txt ./01.prodigal/KP_Generate_v6_4016_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_4016_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_4016.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_4016_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_4016/KP_Generate_v6_4016_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_4084
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_4084
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_4084.fasta -o ./01.prodigal/KP_Generate_v6_4084_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_4084_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_4084_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_4084_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_4084_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_4084_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_4084_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_4084_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_4084_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_4084_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_4084_result.txt ./01.prodigal/KP_Generate_v6_4084_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_4084_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_4084.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_4084_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_4084/KP_Generate_v6_4084_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_4769
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_4769
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_4769.fasta -o ./01.prodigal/KP_Generate_v6_4769_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_4769_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_4769_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_4769_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_4769_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_4769_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_4769_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_4769_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_4769_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_4769_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_4769_result.txt ./01.prodigal/KP_Generate_v6_4769_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_4769_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_4769.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_4769_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_4769/KP_Generate_v6_4769_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_5019
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_5019
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_5019.fasta -o ./01.prodigal/KP_Generate_v6_5019_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_5019_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_5019_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_5019_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_5019_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_5019_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_5019_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_5019_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_5019_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_5019_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_5019_result.txt ./01.prodigal/KP_Generate_v6_5019_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_5019_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_5019.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_5019_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_5019/KP_Generate_v6_5019_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_5360
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_5360
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_5360.fasta -o ./01.prodigal/KP_Generate_v6_5360_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_5360_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_5360_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_5360_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_5360_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_5360_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_5360_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_5360_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_5360_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_5360_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_5360_result.txt ./01.prodigal/KP_Generate_v6_5360_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_5360_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_5360.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_5360_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_5360/KP_Generate_v6_5360_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_5430
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_5430
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_5430.fasta -o ./01.prodigal/KP_Generate_v6_5430_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_5430_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_5430_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_5430_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_5430_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_5430_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_5430_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_5430_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_5430_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_5430_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_5430_result.txt ./01.prodigal/KP_Generate_v6_5430_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_5430_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_5430.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_5430_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_5430/KP_Generate_v6_5430_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_5475
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_5475
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_5475.fasta -o ./01.prodigal/KP_Generate_v6_5475_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_5475_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_5475_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_5475_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_5475_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_5475_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_5475_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_5475_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_5475_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_5475_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_5475_result.txt ./01.prodigal/KP_Generate_v6_5475_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_5475_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_5475.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_5475_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_5475/KP_Generate_v6_5475_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_5638
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_5638
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_5638.fasta -o ./01.prodigal/KP_Generate_v6_5638_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_5638_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_5638_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_5638_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_5638_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_5638_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_5638_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_5638_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_5638_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_5638_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_5638_result.txt ./01.prodigal/KP_Generate_v6_5638_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_5638_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_5638.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_5638_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_5638/KP_Generate_v6_5638_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_6181
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_6181
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_6181.fasta -o ./01.prodigal/KP_Generate_v6_6181_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_6181_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_6181_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_6181_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_6181_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_6181_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_6181_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_6181_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_6181_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_6181_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_6181_result.txt ./01.prodigal/KP_Generate_v6_6181_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_6181_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_6181.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_6181_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_6181/KP_Generate_v6_6181_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_6331
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_6331
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_6331.fasta -o ./01.prodigal/KP_Generate_v6_6331_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_6331_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_6331_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_6331_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_6331_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_6331_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_6331_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_6331_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_6331_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_6331_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_6331_result.txt ./01.prodigal/KP_Generate_v6_6331_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_6331_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_6331.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_6331_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_6331/KP_Generate_v6_6331_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_6935
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_6935
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_6935.fasta -o ./01.prodigal/KP_Generate_v6_6935_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_6935_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_6935_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_6935_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_6935_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_6935_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_6935_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_6935_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_6935_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_6935_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_6935_result.txt ./01.prodigal/KP_Generate_v6_6935_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_6935_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_6935.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_6935_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_6935/KP_Generate_v6_6935_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_6953
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_6953
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_6953.fasta -o ./01.prodigal/KP_Generate_v6_6953_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_6953_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_6953_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_6953_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_6953_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_6953_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_6953_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_6953_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_6953_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_6953_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_6953_result.txt ./01.prodigal/KP_Generate_v6_6953_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_6953_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_6953.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_6953_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_6953/KP_Generate_v6_6953_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_6958
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_6958
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_6958.fasta -o ./01.prodigal/KP_Generate_v6_6958_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_6958_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_6958_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_6958_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_6958_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_6958_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_6958_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_6958_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_6958_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_6958_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_6958_result.txt ./01.prodigal/KP_Generate_v6_6958_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_6958_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_6958.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_6958_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_6958/KP_Generate_v6_6958_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_6972
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_6972
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_6972.fasta -o ./01.prodigal/KP_Generate_v6_6972_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_6972_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_6972_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_6972_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_6972_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_6972_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_6972_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_6972_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_6972_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_6972_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_6972_result.txt ./01.prodigal/KP_Generate_v6_6972_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_6972_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_6972.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_6972_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_6972/KP_Generate_v6_6972_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_7037
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_7037
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_7037.fasta -o ./01.prodigal/KP_Generate_v6_7037_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_7037_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_7037_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_7037_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_7037_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_7037_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_7037_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_7037_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_7037_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_7037_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_7037_result.txt ./01.prodigal/KP_Generate_v6_7037_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_7037_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_7037.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_7037_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_7037/KP_Generate_v6_7037_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_7294
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_7294
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_7294.fasta -o ./01.prodigal/KP_Generate_v6_7294_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_7294_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_7294_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_7294_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_7294_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_7294_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_7294_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_7294_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_7294_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_7294_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_7294_result.txt ./01.prodigal/KP_Generate_v6_7294_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_7294_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_7294.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_7294_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_7294/KP_Generate_v6_7294_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_742
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_742
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_742.fasta -o ./01.prodigal/KP_Generate_v6_742_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_742_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_742_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_742_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_742_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_742_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_742_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_742_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_742_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_742_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_742_result.txt ./01.prodigal/KP_Generate_v6_742_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_742_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_742.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_742_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_742/KP_Generate_v6_742_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_7486
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_7486
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_7486.fasta -o ./01.prodigal/KP_Generate_v6_7486_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_7486_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_7486_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_7486_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_7486_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_7486_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_7486_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_7486_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_7486_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_7486_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_7486_result.txt ./01.prodigal/KP_Generate_v6_7486_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_7486_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_7486.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_7486_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_7486/KP_Generate_v6_7486_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_7490
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_7490
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_7490.fasta -o ./01.prodigal/KP_Generate_v6_7490_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_7490_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_7490_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_7490_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_7490_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_7490_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_7490_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_7490_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_7490_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_7490_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_7490_result.txt ./01.prodigal/KP_Generate_v6_7490_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_7490_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_7490.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_7490_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_7490/KP_Generate_v6_7490_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_7784
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_7784
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_7784.fasta -o ./01.prodigal/KP_Generate_v6_7784_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_7784_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_7784_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_7784_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_7784_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_7784_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_7784_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_7784_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_7784_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_7784_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_7784_result.txt ./01.prodigal/KP_Generate_v6_7784_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_7784_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_7784.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_7784_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_7784/KP_Generate_v6_7784_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_7810
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_7810
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_7810.fasta -o ./01.prodigal/KP_Generate_v6_7810_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_7810_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_7810_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_7810_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_7810_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_7810_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_7810_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_7810_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_7810_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_7810_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_7810_result.txt ./01.prodigal/KP_Generate_v6_7810_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_7810_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_7810.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_7810_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_7810/KP_Generate_v6_7810_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_786
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_786
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_786.fasta -o ./01.prodigal/KP_Generate_v6_786_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_786_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_786_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_786_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_786_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_786_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_786_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_786_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_786_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_786_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_786_result.txt ./01.prodigal/KP_Generate_v6_786_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_786_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_786.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_786_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_786/KP_Generate_v6_786_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_8440
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_8440
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_8440.fasta -o ./01.prodigal/KP_Generate_v6_8440_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_8440_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_8440_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_8440_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_8440_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_8440_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_8440_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_8440_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_8440_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_8440_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_8440_result.txt ./01.prodigal/KP_Generate_v6_8440_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_8440_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_8440.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_8440_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_8440/KP_Generate_v6_8440_genome.gbk 
mkdir //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_9014
cd //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_9014
mkdir ./01.prodigal
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/Prodigal/prodigal -i /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_9014.fasta -o ./01.prodigal/KP_Generate_v6_9014_genome.prodigal.gff -a ./01.prodigal/KP_Generate_v6_9014_genome.proteins.faa -d ./01.prodigal/KP_Generate_v6_9014_genome.genes.fasta -f gff -c
mkdir ./02.blastp
blastp -query ./01.prodigal/KP_Generate_v6_9014_genome.proteins.faa -db /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/NCBI_phage/NCBI_protien_phage_remove_hypothetical.fasta -evalue 1e-5 -max_target_seqs 1 -out ./02.blastp/KP_Generate_v6_9014_genome.proteins.blastp.txt -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"
mkdir ./03.phmmer
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniref.output.txt --tblout ./03.phmmer/uniref.tblout.txt ./01.prodigal/KP_Generate_v6_9014_genome.proteins.faa /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniref/uniref_phage.fasta
grep '#' -v ./03.phmmer/uniref.tblout.txt>./03.phmmer/uniref.tblout.txt.hit 
mkdir ./04.hmmscan
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/phmmer -E 1e-5 -o ./03.phmmer/uniprotkb.output.txt --tblout ./03.phmmer/uniprotkb.tblout.txt ./01.prodigal/KP_Generate_v6_9014_genome.proteins.faa  /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/uniprotkb/uniprotkb_phage.fasta
grep '#' -v ./03.phmmer/uniprotkb.tblout.txt>./03.phmmer/uniprotkb.tblout.txt.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer-3.1b2/bin/hmmscan -E 1e-5 -o ./04.hmmscan/pfam.output.txt --tblout ./04.hmmscan/pfam.tblout.txt /hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/db/Pfam/Pfam-A.hmm ./01.prodigal/KP_Generate_v6_9014_genome.proteins.faa
grep '#' -v ./04.hmmscan/pfam.tblout.txt>./04.hmmscan/pfam.tblout.txt.hit
mkdir ./05.result
cat ./04.hmmscan/pfam.tblout.txt.hit ./03.phmmer/uniref.tblout.txt.hit ./03.phmmer/uniprotkb.tblout.txt.hit>./05.result/hmm.hit
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/hmmer_hit_uniq.py ./05.result/hmm.hit ./05.result/hmm.hit.uniq
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/get_result.py ./02.blastp/KP_Generate_v6_9014_genome.proteins.blastp.txt ./05.result/hmm.hit.uniq ./05.result/KP_Generate_v6_9014_result.txt
mkdir ./06.gff2genebank
/hwfssz5/ST_HEALTH/P17Z10200N0246/USER/share/02_annotation/software/result2gff.py ./05.result/KP_Generate_v6_9014_result.txt ./01.prodigal/KP_Generate_v6_9014_genome.prodigal.gff ./06.gff2genebank/KP_Generate_v6_9014_genome.result.gff
seqret -sequence /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/KP_Generate_v6_9014.fasta -feature -fformat gff -fopenfile ./06.gff2genebank/KP_Generate_v6_9014_genome.result.gff -osformat genbank -outseq //ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/DesignerPhage/PhageDesigner_v6_79z_fasta_20220513/annotation/KP_Generate_v6_9014/KP_Generate_v6_9014_genome.gbk 
