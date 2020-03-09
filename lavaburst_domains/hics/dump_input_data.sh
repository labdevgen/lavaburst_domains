hic="http://genedev.bionet.nsc.ru/site/hic_out/Anopheles/hic/AcolNg/AcolNg_V3.hic";
juicer="$HOME/.local/bin/juicer_tools.1.9.8_jcuda.0.8.jar";
res="25000";
maxdist="1000000"
out=$(echo $hic | rev | cut -d "/" -f1 | rev)".$res.oe"; \
#rm $out;
#for chr in {"2L","2R","3L","3R","X"}; do
#    echo $chr;
#    java -jar $juicer dump oe KR $hic $chr $chr BP $res | \
#    awk -v chr="$chr" 'NR>1 {OFS="\t"; print chr,$0}' >> $out;
#done;
awk -v maxdist="$maxdist" '($3-$2)<=maxdist' $out > $out.${maxdist}.MB
