# To generate pathIndex
# from 20210325_breakseqPIC
IDS=$(cut -f1 data/raw/1KGP_data/integrated_call_samples_v3.20130502.ALL.panel.txt | tail -n+2 )

for ID in $IDS; do

	MAIN_FILE=$(grep $ID data/raw/1KGP_data/20130502.phase3.low_coverage.alignment.index | grep '\.mapped' | cut -f1 | sed -e 's/data.*alignment\///g')
	OTHER_FILE=$(grep $ID data/raw/1KGP_data/20130502.phase3.low_coverage.alignment.index | grep '\.unmapped' | cut -f1 | sed -e 's/data.*alignment\///g')

	PREPATH="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/$ID/alignment/"

	echo -e "${ID}\t${PREPATH}${MAIN_FILE}\t${PREPATH}${OTHER_FILE}" >> data/use/1KGP_data/pathIndex.txt
 
done
