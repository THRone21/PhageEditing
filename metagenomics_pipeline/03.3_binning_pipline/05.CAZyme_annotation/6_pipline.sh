############################
##### cazyme数据库注释 #####
############################
mkdir Bin_all/Bin_cazy

echo -e "\033[32mcazyme数据库注释：\033[0m"
for i in Bin_all/Bin_prokka/bin.*; do
    base=${i##*/}
    time /home/cheng/softwares/miniconda2/bin/diamond blastp --db /media/cheng/disk2/Databases/protein/CAZy/CAZyDB.07312018.dmnd --query $i/$base.faa -e 0.00001 --outfmt 6 --more-sensitive --max-target-seqs 1 --threads $thread --quiet --out Bin_all/Bin_cazy/$base.raw
    echo -e "\033[32m$i Done...\033[0m" 
done
echo -e "\033[32mcazyme数据库注释 Done...\033[0m"

echo -e "\033[32m统计cazy level 1 (group): \033[0m"
for i in Bin_all/Bin_cazy/bin.*.raw; do
    base=${i##*/}
    name=${base%.raw}
    awk '{print $1,$2}' $i | uniq | awk '{print $2}' | sed 's/|/ /g' | awk '{print $2}' | sed 's/[|_]/ /g' | awk '{print $1}' | sed 's/[0123456789]//g' > Bin_all/Bin_cazy/$name.clean
    echo -e "\033[32m$i Done...\033[0m"
done

echo -e "\033[32mcazyme group 统计: \033[0m"
R CMD BATCH --args /home/cheng/huty/softwares/script_bin/cazy_group.R
echo -e "\033[32mcazyme group 热图、箱图: \033[0m"
R CMD BATCH --args /home/cheng/huty/softwares/script_bin/cazy_group_pheatmap_boxplot.R

echo -e "\033[32mcazyme group 饼图: \033[0m"
R CMD BATCH --args /home/cheng/huty/softwares/script_bin/cazy_group_pie.R
echo -e "\033[32mcazyme group 堆叠图: \033[0m"
R CMD BATCH --args /home/cheng/huty/softwares/script_bin/cazy_group_stack.R

echo -e "\033[32m统计cazy level 2 (all) 统计: \033[0m"
for i in Bin_all/Bin_cazy/bin.*.raw; do
    base=${i##*/}
    name=${base%.raw}
    count_cazy.py $i > Bin_all/Bin_cazy/$name.sum
    echo -e "\033[32m$i Done...\033[0m"
done

echo -e "\033[32m统计cazy level 2 (all) 统计: \033[0m"
R CMD BATCH --args /home/cheng/huty/softwares/script_bin/cazy_all.R
echo -e "\033[32m统计cazy level 2 (all) 饼图: \033[0m"
R CMD BATCH --args /home/cheng/huty/softwares/script_bin/cazy_all_pie.R

echo -e "\033[32mcazyme pdf to png: \033[0m"
convert -density 400 -quality 200 Bin_all/Bin_cazy/cazy_group_pheatmap.pdf Bin_all/Bin_cazy/cazy_group_pheatmap.png
convert -density 400 -quality 200 Bin_all/Bin_cazy/cazy_group_boxplot.pdf Bin_all/Bin_cazy/cazy_group_boxplot.png
convert -density 400 -quality 200 Bin_all/Bin_cazy/cazy_group_stack.pdf Bin_all/Bin_cazy/cazy_group_stack.png

for i in Bin_all/Bin_cazy/cazy_all_pie/*.pdf; do
    base=${i%.pdf}
    convert -density 400 -quality 200 $i ${base}.png
    echo -e "\033[32m\t $i all pie Done...\033[0m"
done

for i in Bin_all/Bin_cazy/cazy_group_pie/*.pdf; do
    base=${i%.pdf}
    convert -density 400 -quality 200 $i ${base}.png
    echo -e "\033[32m\t $i group pie Done...\033[0m"
done

mv *.Rout tmp