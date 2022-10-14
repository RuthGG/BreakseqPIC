#################
## FIRST ROUND ##
#################

## Pathways
data_path="../datos_librerias/"
ld_path="aux"
sv_path="aux/SV_1KGP/"


## Inversions
inversions=$( cat ../inversions_completeList.txt | sed -e 's/\n/\t/g')

# cd $ld_path/
# ld_path="./"
rm $ld_path/allinversions.vcf
for inv in $inversions; do

	if [ -d "$inv" ]; then rm -Rf $inv; fi

	echo '############'
	echo start $inv
	echo '############'

	# echo "## 1KGP samples"
	# samples_1kgp=$(cut -f1 ../../1KGP_data/integrated_call_samples_v3.20130502.ALL.panel.txt | grep -v 'sample')


	echo "## From GTypes to VCF"
	if [ $inv == "HsInv0052_region" ]; then
	chr=$(grep 'HsInv0052' $data_path/bplib.coords | cut -f2 | uniq)
	pos_left=$(grep 'HsInv0052' $data_path/bplib.coords | cut -f3,4 | sed -e 's/\t/\n/g' | sort | head -n+1)
	pos_right=$(grep 'HsInv0052' $data_path/bplib.coords | cut -f3,4 | sed -e 's/\t/\n/g' | sort | tail -n-1)
	start=$(($pos_left-500000))
	end=$(($pos_right+500000))
	elif [ $inv == "HsInv0409" ]; then
	chr=$(grep 'HsInv409' $data_path/bplib.coords | cut -f2 | uniq)
	pos_left=$(grep 'HsInv409' $data_path/bplib.coords | cut -f3,4 | sed -e 's/\t/\n/g' | sort | head -n+1)
	pos_right=$(grep 'HsInv409' $data_path/bplib.coords | cut -f3,4 | sed -e 's/\t/\n/g' | sort | tail -n-1)
	start=$(($pos_left-500000))
	end=$(($pos_right+500000))
	else
	chr=$(grep $inv $data_path/bplib.coords | cut -f2 | uniq)
	pos_left=$(grep $inv $data_path/bplib.coords | cut -f3,4 | sed -e 's/\t/\n/g' | sort | head -n+1)
	pos_right=$(grep $inv $data_path/bplib.coords | cut -f3,4 | sed -e 's/\t/\n/g' | sort | tail -n-1)
	start=$(($pos_left-500000))
	end=$(($pos_right+500000))
	fi



	bcftools view -Ov -r $chr:${start}-${end}  $sv_path/ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.vcf.gz > $ld_path/input_inversion.vcf
	grep 'INV' $ld_path/input_inversion.vcf | awk '$5=="<INV>"' | awk -v inv=$inv '$3=inv' >> $ld_path/allinversions.vcf
	grep 'INV' $ld_path/input_inversion.vcf | awk '$5=="<INV>"' | awk -v inv=$inv '$3=inv' > $ld_path/inversion.vcf


	rm $ld_path/input_inversion.vcf


	nrow=$(cat $ld_path/inversion.vcf | wc -l)

	if [ $nrow -gt 0 ]; then

		echo "## Inversion VCF into 1KGP VCF"
		if [ $chr == "X" ]; then
		vcf="../1KGP_data/vcf/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz"
		else
		vcf="../1KGP_data/vcf/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
		fi

		bcftools view -Ov -r $chr:${start}-${end} ${vcf}  > $ld_path/haplotypes.vcf

		for ids in $(seq 1 $nrow) ; do
			
			sed "${ids}q;d" $ld_path/inversion.vcf > $ld_path/lastline.txt
			cat $ld_path/haplotypes.vcf $ld_path/lastline.txt > $ld_path/haplotypes_unsorted.txt
			bcftools sort $ld_path/haplotypes_unsorted.txt > $ld_path/haplotypes_sorted.vcf
		    
			rm $ld_path/haplotypes.vcf $ld_path/haplotypes_unsorted.txt

			# echo "## Convert to plink format"
			ulimit -n 3000
			vcftools --vcf $ld_path/haplotypes.vcf --plink --chr $chr --out $ld_path/genotypes_plink_format


			# echo "## LD calculation (diferente en X??)"
			# folder_plink="/home/jon/soft/plink-1.07-x86_64"
			plink --file $ld_path/genotypes_plink_format --r2 --ld-snp $inv --ld-window-kb 1000000 --ld-window 999999 --ld-window-r2 0 --noweb

			# echo "# Save"
			mkdir $inv
			cat $ld_path/plink.ld >> $ld_path/$inv/plink.ld
			rm $ld_path/genotypes* $ld_path/haplotypes* $ld_path/plink.log $ld_path/plink.nosex


		done

	else
		echo "No 1KGP variants here"
	fi

done



