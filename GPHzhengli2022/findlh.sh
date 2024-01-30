#!/usr/bin/bash
for locs in $(cat /ldfssz1/ST_INFECTION/P17Z10200N0246_Phage_XMF/USER/thr/GPHzhengli2022/phagelocate.txt)
do
        echo -e "$locs\n$(ls -l -h $locs)" >>./lh_loc.txt
done
