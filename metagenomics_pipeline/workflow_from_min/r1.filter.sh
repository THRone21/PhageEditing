#fastp
/ldfssz1/ST_INFECTION/F16ZQSB1SY3030_Salmonella/USER/lindechun/software/fastp/fastp -i /zfssz6/CNGB_DATA/BGISEQ01/DIPSEQ/DIPSEQT1/P17Z10200N0246_Temp/LYZLib20210518-2/210704_SEQ053_DP8400021210TR_L01_PERMAFROST0518-2001/DP8400021210TR_L01_43_1.fq.gz -o ./01.1_fastp/43_fastp1.fq.gz -I /zfssz6/CNGB_DATA/BGISEQ01/DIPSEQ/DIPSEQT1/P17Z10200N0246_Temp/LYZLib20210518-2/210704_SEQ053_DP8400021210TR_L01_PERMAFROST0518-2001/DP8400021210TR_L01_43_2.fq.gz -O ./01.1_fastp/43_fastp2.fq.gz -5 -3 -q 20 -c -j ./01.1_fastp/fastp.json -h ./01.1_fastp/fastp.html -R ./01.1_fastp/out.prefix -l 30 1> 01.1fastp.result 2>r1.1fastp.sh.err
#soapnuke
/ldfssz1/ST_INFECTION/F16ZQSB1SY3030_Salmonella/USER/lindechun/bin/SOAPnuke filter -T 8 -1 ./01.1_fastp/43_fastp1.fq.gz  -2 ./01.1_fastp/43_fastp2.fq.gz -d -C 43_SOAP.1.fq -D 43_SOAP.2.fq -o ./01.2_SOAPnuke -Q 2 2>r1.2SOAPnuke.sh.err
#prinseq
/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/Miniconda3/envs/prinseq++/bin/prinseq++ -threads 8 -fastq ./01.2_SOAPnuke/43_SOAP.1.fq.gz -fastq2 ./01.2_SOAPnuke/43_SOAP.2.fq.gz -lc_entropy=0.5 -lc_dust=0.5 -out_name ./01.3_prinseq/43 2>r1.3prinseq.sh.err
