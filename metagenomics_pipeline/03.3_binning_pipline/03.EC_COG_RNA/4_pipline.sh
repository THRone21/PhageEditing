####################################
##### prokka：type gene EC COG #####
####################################
echo -e "\033[32mprokka基因组注释: \033[0m"
for i in Bin_all/Bin_pick/*.fa; do
    file=${i##*/}
    base=${file%.fa}
    prokka $i --outdir Bin_all/Bin_prokka/$base --prefix $base --metagenome --cpus $thread --kingdom Bacteria
    echo -e "\033[32m$i prokka Done...\033[0m"
done

mkdir Bin_all/Bin_prokka/prokka_out_table
mkdir Bin_all/Bin_prokka/prokka_map
mkdir Bin_all/Bin_prokka/prokka_map_table

echo -e "\033[32m集合prokka注释tsv文件：\033[0m"
for i in Bin_all/Bin_prokka/bin*; do
    fold=${i##*/}
    mv $i/$fold.tsv Bin_all/Bin_prokka/prokka_out_table  # prokka注释结果
done

echo -e "\033[32m集合prokka比对gff文件：\033[0m"
for i in Bin_all/Bin_prokka/bin*; do
    fold=${i##*/}
    mv $i/$fold.gff Bin_all/Bin_prokka/prokka_map  # prokka比对文件
done

echo -e "\033[32m集合提取prokka gff比对信息：\033[0m"
for i in Bin_all/Bin_prokka/prokka_map/bin.*.gff; do
    base=${i##*/}
    grep '^k' $i > Bin_all/Bin_prokka/prokka_map_table/${base}.txt
    sed -i '1 iseqid\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tattributes' Bin_all/Bin_prokka/prokka_map_table/${base}.txt
done

echo -e "\033[32m统计prokka注释tsv文件：\033[0m"
R CMD BATCH --args /home/cheng/huty/softwares/script_bin/prokka_fun_count.R

mv Bin_all/*.Rout tmp
mv *.Rout tmp

########################
##### cog注释、绘图 #####
########################
mkdir Bin_all/cog
echo -e "\033[32mextract cog: \033[0m"
for i in Bin_all/Bin_prokka/prokka_out_table/bin.*.tsv; do
    file=${i##*/}
    fold=${file%.tsv}
    mkdir Bin_all/cog/$fold
    cat $i | awk -F"\t" '{if($6 != "") printf("%s\t%s\n", $1, $6)}' > Bin_all/cog/$fold/${fold}_cog.txt
    echo -e "\033[32m\t $i cog extract Done...\033[0m"
done

echo -e "\033[32mannotate cog: \033[0m"
for i in Bin_all/cog/bin.*; do
    fold=${i##*/}
    Rscript /home/cheng/huty/softwares/script_genome/prokka_cog_annotation.R $i/${fold}_cog.txt /home/cheng/huty/databases/cog/cog_anno.txt ${fold}_cog_annotation.txt
    echo -e "\033[32m\t $i annotate cog Done...\033[0m"
done

echo -e "\033[32msummary annotated cog: \033[0m"
for i in Bin_all/cog/bin.*; do
    fold=${i##*/}
    cat $i/${fold}_cog_annotation.txt | awk -F"\t" 'BEGIN{OFS="\t"}{if($6 != "NA") print $3}' | sed '1d' | sort | uniq -c | awk 'BEGIN{OFS="\t"}{print $2,$1}' | sed '1 ifunc\tCount' > $i/${fold}_cog_sum.txt
    echo -e "\033[32m\t $i summary annotated cog Done...\033[0m"
done

echo -e "\033[32m分类 annotated cog: \033[0m"
for i in Bin_all/cog/bin.*; do
    fold=${i##*/}
    Rscript /home/cheng/huty/softwares/script_genome/prokka_cog_count_annotation.R $i/${fold}_cog_sum.txt /home/cheng/huty/databases/cog/cog_category.txt ${fold}_cog_sum_annotation.txt
    echo -e "\033[32m\t $i 分类 annotated cog Done...\033[0m"
done

echo -e "\033[32m分类可视化 annotated cog: \033[0m"
for i in Bin_all/cog/bin.*; do
    fold=${i##*/}
    Rscript /home/cheng/huty/softwares/script_genome/prokka_cog_count_annotation_bar.R $i/${fold}_cog_sum_annotation.txt ${fold}_cog_sum_annotation.pdf
    convert -density 400 -quality 200 $i/${fold}_cog_sum_annotation.pdf $i/${fold}_cog_sum_annotation.png
    echo -e "\033[32m\t $i 分类可视化 annotated cog Done...\033[0m"
done
