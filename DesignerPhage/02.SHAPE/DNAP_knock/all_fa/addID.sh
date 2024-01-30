mkdir /hwfssz5/ST_INFECTION/Phage/P17Z10200N0246_Phage_XMF/chengli/shape/0813id_fa
for i in $(ls /hwfssz5/ST_INFECTION/Phage/P17Z10200N0246_Phage_XMF/chengli/shape/fa/*.fa)
	do
		file=`echo $(basename $i .fa)`
		echo "perl /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/cl/SHAPE/network/script/rename.pl $i $file /hwfssz5/ST_INFECTION/Phage/P17Z10200N0246_Phage_XMF/chengli/shape/id_fa/${file}_id.fa" >>addid.0813work.sh
		done
