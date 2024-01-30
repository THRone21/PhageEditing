#!/usr/bin/env bash
echo -e "\033[32mcircos绘图: \033[0m"
# 准备文件
echo -e "\033[32m准备文件: \033[0m"
for i in Bin_all/Bin_prokka/prokka_map/bin.*.gff; do
    name=${i##*/}
    base=${name%.gff}
    mkdir Bin_all/Bin_circos/$base
    grep '^k141' $i > Bin_all/Bin_circos/$base/bin.gene
    # 取gene mapping信息    
    grep '^##sequence' $i | awk '{print $2,$3,$4}' | sed 's/ /\t/g' | sed '1 icontig\tstart\tend' > Bin_all/Bin_circos/$base/bin.contig
    # 取contig长度信息
done
echo "circos: prepare done"

# 0. 文件：chr
echo -e "\033[32m0. 文件：chr: \033[0m"
for i in Bin_all/Bin_circos/bin.*; do
    temp=`pwd`
    cd $i
    R CMD BATCH --args /home/cheng/huty/softwares/script_circos/conf.chr.R
    cd $temp
done
echo "circos: chr done"

echo -e "\033[32m1. 文件：+CDS -CDS tRNA rRNA tmRNA\033[0m"
# 1. 文件：+CDS -CDS tRNA rRNA tmRNA
for i in Bin_all/Bin_circos/bin.*; do
    awk '{if($3=="CDS" && $7=="+") print $1,$4,$5,$7}' $i/bin.gene > $i/bin.positive.cds
    awk '{if($3=="CDS" && $7=="-") print $1,$4,$5,$7}' $i/bin.gene > $i/bin.negative.cds
    awk '{if($3=="tRNA") print $1,$4,$5,$3}' $i/bin.gene > $i/bin.tRNA
    awk '{if($3=="rRNA") print $1,$4,$5,$3}' $i/bin.gene > $i/bin.rRNA
    awk '{if($3=="tmRNA") print $1,$4,$5,$3}' $i/bin.gene > $i/bin.tmRNA
done
echo "circos: +CDS -CDS tRNA rRNA tmRNA done"

echo -e "\033[32m2. 文件：GC GH CE AA PL CBM\033[0m"
# 2. 文件：GC GH CE AA PL CBM
for i in Bin_all/Bin_cazy/bin.*.raw; do
    base=${i##*/}
    bold=${base%.raw}
    awk '{print $1,$2}' $i | uniq | sed 's/[_|]/ /g' | awk '{print $1,$2,$4}' | sed 's/ /_/' | sed 's/[0-9]*$//' | sed '1 i id cazy' > Bin_all/Bin_circos/$bold/bin.cazy
    # cazy注释结果
done

for i in Bin_all/Bin_circos/bin.*; do
    sed 's/ID=//g' $i/bin.gene | sed 's/;/\t/g' | awk '{print $1,$4,$5,$9}' | sed '1 i contig start end id' > $i/bin.func
    # gff文件
done

for i in Bin_all/Bin_circos/bin.*; do
    temp=`pwd`
    cd $i
    R CMD BATCH --args /home/cheng/huty/softwares/script_circos/conf.merge.cazy_func.R
    # bin.cazy.conf
    cd $temp
done
# 合并

for i in Bin_all/Bin_circos/bin.*; do
    awk '{if($4 == "GH") print $0}' $i/bin.cazy.conf > $i/bin.cazy.GH
    awk '{if($4 == "GT") print $0}' $i/bin.cazy.conf > $i/bin.cazy.GT
    awk '{if($4 == "CE") print $0}' $i/bin.cazy.conf > $i/bin.cazy.CE
    awk '{if($4 == "PL") print $0}' $i/bin.cazy.conf > $i/bin.cazy.PL
    awk '{if($4 == "CBM") print $0}' $i/bin.cazy.conf > $i/bin.cazy.CBM
    awk '{if($4 == "AA") print $0}' $i/bin.cazy.conf > $i/bin.cazy.AA
done
echo "circos: GC GH CE AA PL CBM done"

echo -e "\033[32m3. GC：\033[0m"
# 3. GC
for i in Bin_all/Bin_pick/bin.*.fa; do
    base=${i##*/}
    bold=${base%.fa}
    perl /home/cheng/huty/softwares/script_circos/contig_gc.pl $i > Bin_all/Bin_circos/$bold/bin.gc  # 计算gc含量
    sed -i '1 icontig\tgc\tcompare' Bin_all/Bin_circos/$bold/bin.gc  # 加header
    
    temp=`pwd`
    cd Bin_all/Bin_circos/$bold/
    R CMD BATCH --args /home/cheng/huty/softwares/script_circos/conf.merge.gc_contig.R
    # 加start end
    sed '1d' bin.gc.conf | awk '{if($5 == "small") print $0}' | awk '{print $1,$2,$3,$4}' > bin.gc.small
    # 选择small
    sed '1d' bin.gc.conf | awk '{if($5 == "big") print $0}' | awk '{print $1,$2,$3,$4}' > bin.gc.big
    # 选择big
    cd $temp
    
done
echo "circos: GC done"

echo -e "\033[32m4. GC skew：\033[0m"
# 4. GC skew
for i in Bin_all/Bin_pick/bin.*.fa; do
    base=${i##*/}
    bold=${base%.fa}
    perl /home/cheng/huty/softwares/script_circos/contig_gc_skew.pl $i > Bin_all/Bin_circos/$bold/bin.gc.skew  # 计算gc.skew
    
    temp=`pwd`
    cd Bin_all/Bin_circos/$bold/
    sed -i '1 icontig\tvalue' bin.gc.skew  # 加header
    R CMD BATCH --args /home/cheng/huty/softwares/script_circos/conf.merge.gc.skew_contig.R
    # 加start end
    sed '1d' bin.gc.skew.conf | awk '{if($4 > 0) print $0}' > bin.gc.skew.positive
    # 选择>0
    sed '1d' bin.gc.skew.conf | awk '{if($4 < 0) print $0}' > bin.gc.skew.negative
    # 选择<0
    cd $temp

done
echo "circos: GC skew done"

echo -e "\033[32m画circos：\033[0m"
# 画circos
for i in Bin_all/Bin_circos/bin.*; do
    base=${i##*/}
    temp=`pwd`
    cd $i
    /home/cheng/huty/softwares/circos-0.69-9/bin/circos -conf /home/cheng/huty/softwares/script_circos/circos.conf --outputdir ../Figure -outputfile $base
    echo -e "\033[32m $i circos plot Done... \033[0m"
    cd $temp
done
echo "circos: plot circos done"

echo -e "\033[32m加图例：\033[0m"
# 加图例
# 做图例: R CMD BATCH --args circos.legend_luochaobing.R
temp=`pwd`
cd Bin_all/Bin_circos/Figure/
cp /home/cheng/huty/softwares/script_circos/legend.png ./
R CMD BATCH --args /home/cheng/huty/softwares/script_circos/circos.result.R
cd $temp
echo "circos: add legend done"
