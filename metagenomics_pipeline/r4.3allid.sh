cat ./*outfmt6.txt |awk '{print $1}' |sort -n |uniq >allid.txt
