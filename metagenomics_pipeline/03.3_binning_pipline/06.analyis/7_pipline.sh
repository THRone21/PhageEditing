#######################
##### bin进化分析 #####
#######################
mkdir Bin_all/Bin_phylo

echo -e "\033[32mphylophlan进化分析: \033[0m"
mkdir /home/cheng/huty/softwares/phylophlan/input/$project.denovo
mkdir /home/cheng/huty/softwares/phylophlan/input/$project.insert

for i in Bin_all/Bin_prokka/bin.*; do
    fold=${i##*/}
    cp $i/$fold.faa /home/cheng/huty/softwares/phylophlan/input/$project.denovo
    cp $i/$fold.faa /home/cheng/huty/softwares/phylophlan/input/$project.insert
done

temp=`pwd`
cd /home/cheng/huty/softwares/phylophlan
echo -e "\033[32mde novo phylophlan进化分析: \033[0m"
python phylophlan.py -u $project.denovo --nproc $thread
echo -e "\033[32minsert phylophlan进化分析: \033[0m"
python phylophlan.py -i -t $project.insert --nproc $thread
cd $temp

temp=`pwd`
cd /home/cheng/huty/softwares/phylophlan/output/$project.insert
echo -e "\033[32m合并insert分析物种注释结果: \033[0m"
cat imputed_conf_* | sort -t $'\t' -k 2r | sed '1 iid\ttax' > bin_taxonomy.txt
cd ../
mv /home/cheng/huty/softwares/phylophlan/output/$project.insert $temp/Bin_all/Bin_phylo 
mv /home/cheng/huty/softwares/phylophlan/output/$project.denovo $temp/Bin_all/Bin_phylo
cd $temp

rm -r /home/cheng/huty/softwares/phylophlan/input/$project.denovo
rm -r /home/cheng/huty/softwares/phylophlan/input/$project.insert
rm -r /home/cheng/huty/softwares/phylophlan/data/$project.denovo
rm -r /home/cheng/huty/softwares/phylophlan/data/$project.insert

echo -e "\033[32m绘制de novo phylophlan进化树: \033[0m"
temp=`pwd`
cd Bin_all/Bin_phylo/$project.denovo
sed '1d' ../$project.insert/bin_taxonomy.txt | sed 's/.p__/\t/g' | sed 's/.c__/\t/g' | awk '{print $1,$3}' | sed 's/ /\t/g' | sed '1 ibin\tPhylum' > map.txt
R CMD BATCH --args /home/cheng/huty/softwares/script_bin/tree_bins_circ.R  # 画图
convert -density 400 -quality 200 tree_bins.pdf tree_bins.png
cd $temp

mv *.Rout tmp


##########################
##### 绘制circos圈图 #####
##########################
mkdir Bin_all/Bin_circos
mkdir Bin_all/Bin_circos/Figure

bash /home/cheng/huty/softwares/script_circos/circos.sh

echo -e "\033[32mpdf to png: \033[0m"
for i in Bin_all/Bin_circos/Figure/*.circos.pdf; do
    base=${i%.pdf}
    convert -density 400 -quality 200 $i ${base}.png
done 

mv *.Rout tmp
mv .RData tmp

