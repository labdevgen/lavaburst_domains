#hic="https://genedev.bionet.nsc.ru/site/hic_out/Anopheles/hic/AcolNg/AcolNg_V3.hic";
for hic in {"https://genedev.bionet.nsc.ru/site/hic_out/Anopheles/hic/AcolNg/AcolNg_V3.hic",\
"https://genedev.bionet.nsc.ru/site/hic_out/Anopheles/hic/AalbS2/AalbS2_V3.hic",\
"https://genedev.bionet.nsc.ru/site/hic_out/Anopheles/hic/AatrE3/AatrE3_V3.hic",\
"https://genedev.bionet.nsc.ru/site/hic_out/Anopheles/hic/AmerR4/AmerR4_V3.hic",\
"https://genedev.bionet.nsc.ru/site/hic_out/Anopheles/hic/AsteI2/AsteI2_V3.hic"};
do
	echo $hic;
	#juicer="$HOME/.local/bin/juicer_tools.1.9.8_jcuda.0.8.jar";
	juicer="./juicer_tools_1.11.09_jcuda.0.8.jar";
	res="25000";
	maxdist="1000000"
	out=$(echo $hic | rev | cut -d "/" -f1 | rev)".$res.oe"; \
	rm $out;
	for chr in {"2L","2R","3L","3R","X"}; do
		echo "--"$chr;
		java -jar $juicer dump oe KR $hic $chr $chr BP $res | \
		awk -v chr="$chr" 'NR>1 {OFS="\t"; print chr,$0}' >> $out;
	done;
	awk -v maxdist="$maxdist" '($3-$2)<=maxdist' $out > $out.${maxdist}.MB
done;