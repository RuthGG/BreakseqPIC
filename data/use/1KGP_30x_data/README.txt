# To generate pathIndex
# from 20210325_breakseqPIC

sed "/#.*$/d"  data/raw/1KGP_30x_data/1000G_2504_high_coverage.sequence.index | awk -v OFS="\t" '{print $23, $24}' > data/use/1KGP_30x_data/indlist_all.txt
sed "/#.*$/d"  data/raw/1KGP_30x_data/1000G_698_related_high_coverage.sequence.index | awk -v OFS="\t" '{print $28, $29}' >> data/use/1KGP_30x_data/indlist_all.txt


sed "/#.*$/d"  data/raw/1KGP_30x_data/1000G_2504_high_coverage.sequence.index | awk -v OFS="\t" '{print $23, $1, $1}' > data/use/1KGP_30x_data/pathIndex.txt
sed "/#.*$/d"  data/raw/1KGP_30x_data/1000G_698_related_high_coverage.sequence.index | awk -v OFS="\t" '{print $28, $1, $1}' >> data/use/1KGP_30x_data/pathIndex.txt

# To generate halfkgp without the static dataset list
$ wc -l indlist_all.txt 
3202 indlist_all.txt

$ wc -l static_dataset_list.txt 
414 static_dataset_list.txt

$ cut -f1  static_dataset_list.txt | grep -v -f - indlist_all.txt > indlist_nostatic.txt
$ wc -l indlist_nostatic.txt 
2788 indlist_nostatic.txt

$ HOWMANY=$(($(cat indlist_nostatic.txt | wc -l)/2))
$ echo $HOWMANY
1394

$ head -n$HOWMANY indlist_nostatic.txt > indlist_nostatic_half1.txt 
$ grep -v -f indlist_nostatic_half1.txt indlist_nostatic.txt > indlist_nostatic_half2.txt 

$ wc -l indlist_nostatic_half*
 1394 indlist_nostatic_haf1.txt
 1394 indlist_nostatic_haf2.txt
 2788 total
