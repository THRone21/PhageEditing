/share/app/bowtie2-2.2.5/bowtie2-build RS-T294/RS-T294.contigs.fa contigs/RS-T294.contigs.fa.index && /share/app/bowtie2-2.2.5/bowtie2 -p 10 -x contigs/RS-T294.contigs.fa.index -1 /hwfssz1/ST_AGRIC/P18Z10200N0112/USER/jincanzhi/projects/ruili/01.cleandata/G9/RS-T294/RS-T294.pair.1.fq.gz -2 /hwfssz1/ST_AGRIC/P18Z10200N0112/USER/jincanzhi/projects/ruili/01.cleandata/G9/RS-T294/RS-T294.pair.2.fq.gz -U /hwfssz1/ST_AGRIC/P18Z10200N0112/USER/jincanzhi/projects/ruili/01.cleandata/G9/RS-T294/RS-T294.single.fq.gz -S RS-T294.sam