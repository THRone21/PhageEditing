helpdoc(){
    cat <<EOF
Description:

    This is a help document, please use metawrap-env by input source activate metawrap-env
    - metagenome binning

Usage:

    $0 -n <project name> -p <thread number> -i <raw data filefold> -m <kneaddata mm> -M <megahit mm>

Option:

    -i    followed rawdata filefold; rawdata format: sample_R[12].fastq
    -p    followed thread/processor number
    -n    followed project name
    -m    followed kneaddata memory
    -M    followed megahit max memory percent
    -g    followed metadata file
    -G    followed Category include all group in metadata file
EOF
}

# 参数传递
while getopts ":i:p:n:m:M:g:G:" opt
do
    case $opt in
        i)
        rawdata=`echo $OPTARG`
        ;;
        p)
        thread=`echo $OPTARG`
        ;;
        n)
        project=`echo $OPTARG`
        ;;
        m)
        memory=`echo $OPTARG`
        ;;
        M)
        megahit_mm=`echo $OPTARG`
        ;;
        g)
        metadata=`echo $OPTARG`
        ;;
        G)
        all_group=`echo $OPTARG`
        ;;
        ?)
        echo "未知参数"
        exit 1;;
    esac
done

# 若无指定任何参数则输出帮助文档
if [ $# = 0 ]
then
    helpdoc
    exit 1
fi

########################
##### 原始数据质控 #####
########################

echo -e "\033[32m解压原始数据: \033[0m"
for i in ./$rawdata/*_R[12].fastq.gz; do
    gunzip $i
    echo -e "\033[32m$i Done...\033[0m"
done


echo -e "\033[32m准备文件夹: \033[0m"
mkdir Bin_all
mkdir tmp

echo -e "\033[32m质控原始数据: \033[0m"
mkdir Bin_all/kneaddata

for i in ./$rawdata/*_R1.fastq; do
    r1=$i
    r2=${i%_R1.fastq}_R2.fastq
    kneaddata -t $thread --max-memory $memory -v \
-i $r1 \
-i $r2 \
-o Bin_all/kneaddata \
--trimmomatic /home/cheng/pipelines/Genome/softwares/Trimmomatic-0.39/ \
--trimmomatic-options "ILLUMINACLIP:/home/cheng/pipelines/MetaGenome/data/adaptor_Illumina.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50" \
--fastqc /home/cheng/softwares/fastqc/FastQC/ \
--run-fastqc-start \
--run-fastqc-end \
--remove-intermediate-output
    echo -e "\033[32m$i Done...\033[0m"
done

mkdir Bin_all/Clean_data
echo -e "\033[32m格式化clean数据文件名: \033[0m"
for i in Bin_all/kneaddata/*_R1_kneaddata.trimmed.1.fastq; do
    base=${i%_R1_kneaddata.trimmed.1.fastq}
    name=${base##*/}
    mv $i Bin_all/Clean_data/${name}_1.fastq
    echo -e "\033[32m$i Done...\033[0m"
done

for i in Bin_all/kneaddata/*_R1_kneaddata.trimmed.2.fastq; do
    base=${i%_R1_kneaddata.trimmed.2.fastq}
    name=${base##*/}
    mv $i Bin_all/Clean_data/${name}_2.fastq
    echo -e "\033[32m$i Done...\033[0m"
done

echo -e "\033[32mraw clean数据统计: \033[0m"

python3 /home/cheng/huty/softwares/script_bin/summary_fastq.py rawdata Bin_all/summary_rawdata.txt

python3 /home/cheng/huty/softwares/script_bin/summary_fastq.py Bin_all/Clean_data Bin_all/summary_cleandata.txt

python3 /home/cheng/huty/softwares/script_bin/table_to_html_summary_data.py Bin_all/summary_rawdata.txt Bin_all/summary_rawdata_html.txt

python3 /home/cheng/huty/softwares/script_bin/table_to_html_summary_data.py Bin_all/summary_cleandata.txt Bin_all/summary_cleandata_html.txt

mkdir Bin_all/Bin_align
echo -e "\033[32m合并、组装: \033[0m"
cat Bin_all/Clean_data/*_1.fastq > Bin_all/Bin_align/R1.fastq
cat Bin_all/Clean_data/*_2.fastq > Bin_all/Bin_align/R2.fastq
echo -e "\033[32m合并 Done...\033[0m"

# 自动创建结果文件夹
time megahit -1 Bin_all/Bin_align/R1.fastq -2 Bin_all/Bin_align/R2.fastq -o Bin_all/Assembly_out --out-prefix final -t $thread -m $megahit_mm
echo -e "\033[32m组装 Done...\033[0m"

########################
##### 宏基因组分箱 #####
########################
mkdir Bin_all/Bin

echo -e "\033[32m建索引: \033[0m"
time bowtie2-build --large-index -f Bin_all/Assembly_out/final.contigs.fa Bin_all/Assembly_out/final.contigs --threads $thread --quiet

echo -e "\033[32m1 align: \033[0m"
time bowtie2 --mm -1 Bin_all/Bin_align/R1.fastq -2 Bin_all/Bin_align/R2.fastq -p $thread -x Bin_all/Assembly_out/final.contigs -S Bin_all/Bin_align/final.sam

echo -e "\033[32m2 sam to bam : \033[0m"
time samtools view -@ $thread -b -S Bin_all/Bin_align/final.sam -o Bin_all/Bin_align/final.bam
if [ ! -f "Bin_all/Bin_align/final.bam" ]; then
    echo -e "\033[31merror: 文件不存在\033[0m"
else
    echo -e "\033[32m文件存在，正在删除\033[0m"
    rm Bin_all/Bin_align/final.sam
fi

echo -e "\033[32m3 bam to sorted.bam: \033[0m"
time samtools sort -@ $thread -l 9 -O BAM Bin_all/Bin_align/final.bam -o Bin_all/Bin_align/final.sorted.bam
if [ ! -f "Bin_all/Bin_align/final.sorted.bam" ]; then
    echo -e "\033[31merror: 文件不存在\033[0m"
else
    echo -e "\033[32m文件存在，正在删除\033[0m"
    rm Bin_all/Bin_align/final.bam
fi

echo -e "\033[32m4 计算深度: \033[0m"
time jgi_summarize_bam_contig_depths --outputDepth Bin_all/Bin_align/final.depth.txt Bin_all/Bin_align/final.sorted.bam
echo -e "\033[32m4 计算深度 Done...\033[0m"

echo -e "\033[32m5 分箱: \033[0m"
time metabat2 -m 1500 -t $thread -i Bin_all/Assembly_out/final.contigs.fa -a Bin_all/Bin_align/final.depth.txt -o Bin_all/Bin/bin -v
echo -e "\033[32m5 分箱 Done...\033[0m"

#############################
##### bin质量评估、筛选 #####
#############################
mkdir Bin_all/Bin_quality
mkdir Bin_all/Bin_pick

echo -e "\033[32mcheckm评估Bin质量: \033[0m"
mkdir Bin_all/tmp
time checkm lineage_wf -f Bin_all/Bin_quality/checkm.txt -t $thread -x fa Bin_all/Bin/ Bin_all/Bin_quality/ --tmpdir Bin_all/tmp/
rm -r Bin_all/tmp

echo -e "\033[32mcheckm评估Bin质量 Done...\033[0m"

echo -e "\033[32m整理checkm结果: \033[0m"
grep 'bin' Bin_all/Bin_quality/checkm.txt | sed 's/^  //' | awk '{print $1,$13,$14}' | sed 's/\ /\t/g'| sed 's/\./\t/' | sort -n -k 2 | sed 's/\t/./' | sort -k 2 -n -r | sed '1i\BinID\tCompleteness\tContamination' > Bin_all/Bin_quality/checkm_raw.txt

echo -e "\033[32m筛选clean bin: \033[0m"
grep 'bin' Bin_all/Bin_quality/checkm.txt | sed 's/^  //' | awk '{print $1,$13,$14}' | sed 's/\ /\t/g'| sed 's/\./\t/' | sort -n -k 2 | sed 's/\t/./' | awk '{if($2>=70 && $3<=10) print $0}' | sort -k 2 -n -r | sed '1i\BinID\tCompleteness\tContamination' > Bin_all/Bin_quality/checkm_pick.txt

echo -e "\033[32m获取clean bin ID: \033[0m"
awk '{print $1}' Bin_all/Bin_quality/checkm_pick.txt | sed '1d' > Bin_all/Bin_quality/checkm_pick_id.txt

echo -e "\033[32m复制clean bin: \033[0m"
for i in `cat Bin_all/Bin_quality/checkm_pick_id.txt`; do
    mv Bin_all/Bin/$i.fa Bin_all/Bin_pick/
done

################
##### 统计 #####
################
mkdir Bin_all/Bin_summary

echo -e "\033[32m统计Size和GC含量: \033[0m"
stat_bin_gc.py Bin_all/Bin_pick/ Bin_all/Bin_summary/bin.gc.txt
sed '1d' Bin_all/Bin_summary/bin.gc.txt | sed '1 iBinID\tSize\tGC_percent' > Bin_all/Bin_summary/bin.gc_out.txt

echo -e "\033[32m统计N50和N90: \033[0m"
for i in Bin_all/Bin_pick/*.fa; do
    base=${i##*/}
    name=${base%.fa}
    echo $name >> Bin_all/Bin_summary/N50.N90.txt
    perl /home/cheng/huty/softwares/script_bin/N50.N90.pl $i >> Bin_all/Bin_summary/N50.N90.txt
    echo -e "\033[32m$i Done...\033[0m"
done

sed -i 's/N50: //g' Bin_all/Bin_summary/N50.N90.txt
sed -i 's/N90: //g' Bin_all/Bin_summary/N50.N90.txt

echo -e "\033[32m排序N50和N90: \033[0m"
R CMD BATCH --args /home/cheng/huty/softwares/script_bin/N50.N90.sort.R

echo -e "\033[32m合并checkm n50 n90 size gc%: \033[0m"
R CMD BATCH --args /home/cheng/huty/softwares/script_bin/N50.N90.GC.size.R

mv *.Rout tmp

##################
##### 可视化 #####
##################
mkdir Bin_all/Bin_plot
echo -e "\033[32m计算每个Bin中每条contig的GC含量：\033[0m"
for i in Bin_all/Bin_pick/bin.*.fa; do
    file=${i##*/}
    fold=${file%.fa}
    mkdir Bin_all/Bin_plot/$fold
    python3 /home/cheng/huty/softwares/script_bin/fasta_gc_percent.py $i Bin_all/Bin_plot/$fold/${fold}.gc.txt
    echo -e "\033[32m\t $i gc percent Done...\033[0m"
done

echo -e "\033[32m获取contig深度数据：\033[0m"
cat Bin_all/Bin_align/final.depth.txt | awk -F"\t" 'BEGIN{OFS="\t"}{print $1, $3}' > Bin_all/Bin_plot/all_contig_avgdepth.txt

echo -e "\033[32mmerge contig gc depth：\033[0m"
for i in Bin_all/Bin_plot/bin.*; do
    fold=${i##*/}
    Rscript /home/cheng/huty/softwares/script_bin/merge_gc_depth.R $i/${fold}.gc.txt ../all_contig_avgdepth.txt
    echo -e "\033[32mmerge $i gc depth Done...\033[0m"
done

echo -e "\033[32mmerge all bin gc depth data: \033[0m"
touch Bin_all/Bin_plot/all_bin_gc_depth.txt
for i in Bin_all/Bin_plot/bin*; do
    fold=${i##*/}
    cat $i/${fold}.gc.depth.txt | sed '1d' >> Bin_all/Bin_plot/all_bin_gc_depth.txt
    echo -e "\033[32m\tadd $i gc depth data: \033[0m"
done
sed -i '1 icontig\tGC_percent\tDepth\tBin' Bin_all/Bin_plot/all_bin_gc_depth.txt

echo -e "\033[32mbin gc depth 可视化: \033[0m"
Rscript /home/cheng/huty/softwares/script_bin/bin_scatter.R Bin_all/Bin_plot/all_bin_gc_depth.txt all_bin_gc_depth.pdf

convert -density 400 -quality 200 Bin_all/Bin_plot/all_bin_gc_depth.pdf Bin_all/Bin_plot/all_bin_gc_depth.png