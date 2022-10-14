# To generate pathIndex
# from 20210325_breakseqPIC

sed "/#.*$/d"  data/raw/1KGP_30x_data/1000G_2504_high_coverage.sequence.index | awk -v OFS="\t" '{print $23, $24}' > data/use/1KGP_30x_data/indlist_all.txt
sed "/#.*$/d"  data/raw/1KGP_30x_data/1000G_698_related_high_coverage.sequence.index | awk -v OFS="\t" '{print $28, $29}' >> data/use/1KGP_30x_data/indlist_all.txt


sed "/#.*$/d"  data/raw/1KGP_30x_data/1000G_2504_high_coverage.sequence.index | awk -v OFS="\t" '{print $23, $1, $1}' > data/use/1KGP_30x_data/pathIndex.txt
sed "/#.*$/d"  data/raw/1KGP_30x_data/1000G_698_related_high_coverage.sequence.index | awk -v OFS="\t" '{print $28, $1, $1}' >> data/use/1KGP_30x_data/pathIndex.txt

