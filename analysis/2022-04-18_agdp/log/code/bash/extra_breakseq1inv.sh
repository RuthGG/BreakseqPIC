# Inside the individual-by indiviual breakseq results
inv=""

for ind in $(ls | grep NA); do 
	grep $inv $ind/*uni.sam | awk -v OFS="\t" -v ind=$ind '{print ind, "uni.sam", $0}' >> ${inv}.txt
	grep $inv $ind/*xun.sam | awk -v OFS="\t" -v ind=$ind '{print ind, "xun.sam", $0}' >> ${inv}.txt
	grep $inv $ind/*ini.sam | awk -v OFS="\t" -v ind=$ind '{print ind, "ini.sam", $0}' >> ${inv}.txt
	grep $inv $ind/*fil.sam | awk -v OFS="\t" -v ind=$ind '{print ind, "fil.sam", $0}' >> ${inv}.txt


done


for ind in $(ls | grep NA); do 
	cat $ind/*ref.sam | awk -v OFS="\t" -v ind=$ind '{print ind, "ref.sam", $0}' >> ${inv}_refhits.txt
done
