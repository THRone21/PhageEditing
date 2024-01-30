for i in $(ls ./*fa)
	do
		file=`echo $(basename $i .fa)`
		echo "perl ID.pl internalID.txt $file" >>runID.sh
		done
