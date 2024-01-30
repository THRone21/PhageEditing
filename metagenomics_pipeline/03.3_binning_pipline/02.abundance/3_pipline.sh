#######################
##### 计算bin丰度 #####
#######################
echo -e "\033[32mquant_bins计算bin丰度: \033[0m"
metawrap quant_bins -t $thread -o Bin_all/Bin_quant/ -b Bin_all/Bin_pick/ -a Bin_all/Assembly_out/final.contigs.fa Bin_all/Clean_data/*.fastq
echo -e "\033[32mquant_bins计算bin丰度 Done...\033[0m"

mv Bin_all/Bin_quant/bin_abundance_table.tab Bin_all/Bin_quant/bin_abundance_table.txt

echo -e "\033[32m绘制bin丰度热图、柱形图: \033[0m"
R CMD BATCH --args /home/cheng/huty/softwares/script_bin/quant.pheatmap.R
R CMD BATCH --args /home/cheng/huty/softwares/script_bin/quant.bar.R

convert -density 400 -quality 200 Bin_all/Bin_quant/bin_abundance_pheatmap.pdf Bin_all/Bin_quant/bin_abundance_pheatmap.png
convert -density 400 -quality 200 Bin_all/Bin_quant/bin_abundance_bar.pdf Bin_all/Bin_quant/bin_abundance_bar.png

echo -e "\033[32bin丰度lefse分析: \033[0m"
python3 /home/cheng/huty/softwares/script_bin/quant_lefse.py $metadata $all_group /home/cheng/huty/databases/group_color.list

