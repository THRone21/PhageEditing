##########################
##### kegg注释、绘图 #####
##########################
mkdir Bin_all/kegg

echo -e "\033[32m\nkofamscan注释KEGG：\033[0m"
for i in Bin_all/Bin_prokka/bin.*; do
    fold=${i##*/}
    mkdir Bin_all/kegg/$fold
    exec_annotation -f mapper -o Bin_all/kegg/$fold/${fold}_kegg_raw.txt $i/${fold}.faa
    echo -e "\033[32m\t$i KEGG Done...\033[0m"
done

for i in Bin_all/kegg/bin*; do
    fold=${i##*/}
    cat $i/${fold}_kegg_raw.txt | awk -F"\t" '{if($2 != "") print $2,$1}' | sed 's/ /\t/' | sed '1 ik_id\tgene_id' > $i/${fold}_kegg.txt
    echo -e "\033[32m\t$i KEGG注释结果整理 Done...\033[0m"
done

for i in Bin_all/kegg/bin*; do
    fold=${i##*/}
    Rscript /home/cheng/huty/softwares/script_genome/kofamscan_annotation.R $i ${fold}_kegg.txt ${fold}_kegg_pathway.txt /home/cheng/huty/databases/kofamkoala/KEGG_orthology2pathway.txt
    echo -e "\033[32m\t$i KEGG to pathway Done...\033[0m"
done

for i in Bin_all/kegg/bin*; do
    fold=${i##*/}
    cat $i/${fold}_kegg_pathway.txt | sed '1d' | awk -F"\t" '{print $3}' | sort | uniq -c | awk '{printf("%s\t%s\n", $2, $1)}' | sort -t $'\t' -k 2nr | sed '1 ipathway_id\tpathway_num' > $i/${fold}_kegg_pathway_num.txt
    echo -e "\033[32m\t$i KEGG pathway 统计 Done...\033[0m"
    Rscript /home/cheng/huty/softwares/script_bin/kegg_level_annotation.R $i ${fold}_kegg_pathway_num.txt ${fold}_kegg_pathway_num_annotation.txt /home/cheng/huty/databases/KEGG/pathway_annotation.txt    
    echo -e "\033[32m\t$i KEGG pathway 分类注释 Done...\033[0m"
    Rscript /home/cheng/huty/softwares/script_bin/kegg_pathway_level_bar.R ${i}/${fold}_kegg_pathway_num_annotation.txt ${fold}_kegg_pathway_num_annotation
    echo -e "\033[32m\t$i KEGG pathway barplot绘图 Done...\033[0m"
done

for i in Bin_all/kegg/bin*; do
    fold=${i##*/}
    Rscript /home/cheng/huty/softwares/script_genome/kofamscan_annotation.R $i ${fold}_kegg.txt ${fold}_kegg_annotation.txt /home/cheng/huty/databases/kofamkoala/ko_mapping.txt
    echo -e "\033[32m\t$i KEGG pathway 详细注释 Done...\033[0m"
done

for i in Bin_all/kegg/bin*; do
    fold=${i##*/}
    mkdir $i/Figure
    map=`cat $i/${fold}_kegg_pathway.txt | sed '1d' | awk -F"\t" '{print $3}' | sort | uniq`
    for j in $map; do
        cp /home/cheng/Databases/map/${j}.png $i/Figure/
    done
    echo -e "\033[32m\t$i Figure Done...\033[0m"
done

########################
##### GO注释、绘图 #####
########################
mkdir Bin_all/emapper

echo -e "\033[32m\n emapper注释GO：\033[0m"

for i in Bin_all/Bin_prokka/bin.*; do
    fold=${i##*/}
    mkdir Bin_all/emapper/$fold
    emapper.py -i $i/${fold}.ffn \
--output_dir Bin_all/emapper/$fold \
-o $fold \
--seed_ortholog_evalue 0.00001 \
--data_dir /media/cheng/disk2/Databases/eggnog \
-m diamond --cpu $thread
    echo -e "\033[32m\t$i emapper Done...\033[0m"
done

for i in Bin_all/emapper/bin.*; do
    fold=${i##*/}
    #cat $i/${fold}.emapper.annotations | awk -F"\t" '{print $7}' | sed '/^$/d' | sed 's/,/\n/g' | sed '1d' | sort | uniq | sed '1 iKEGG_KOs' > $i/${fold}_kegg.txt
    cat $i/${fold}.emapper.annotations | sed '/^#/d' | awk -F"\t" '{print $1,$6}' | awk '{if($2 != "") print $0}' | sed '1d'| sed '1 igene_id\tGO_terms' > $i/${fold}_go_map.txt
    cat $i/${fold}_go_map.txt | sed '1d' | awk '{print $2}' | sed 's/,/\n/g' | sort | uniq -c | awk '{print $2,$1}' | sed 's/ /\t/' | sed '1 igo_id\tgo_gene_num' > $i/${fold}_go.txt
    Rscript /home/cheng/huty/softwares/script_genome/emapper_go_annotation.R $i ${fold}_go.txt ${fold}_go_annotation.txt /home/cheng/huty/databases/GO/ncbi/go.annotation.txt
    Rscript /home/cheng/huty/softwares/script_genome/emapper_go_bar.R $i/${fold}_go_annotation.txt ${fold}_go_annotation.pdf
    convert -density 400 -quality 200 $i/${fold}_go_annotation.pdf $i/${fold}_go_annotation.png
    echo -e "\033[32m\t$i GO Done...\033[0m"
done
