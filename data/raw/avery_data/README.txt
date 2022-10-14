# original file was accession_list
sed "s/[abcd]//g" accession_list.txt | cut -f1   |cut -f2 -d"." | sort | uniq > indiv_list.txt

# to identify SRA for each indiv
for i in `cat indiv_list.txt` ; do 
	grep "${i}[a-z]*\s" accession_list.txt | cut -f2 | sed "s/\n/\|/g"
done


