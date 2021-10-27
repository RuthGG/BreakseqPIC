
#####################
## STEP 6 - PART 2 ##
#####################
####################################################################################################################
## EXTRACT TAG SNPS FROM 1KGP ######################################################################################
####################################################################################################################

# PLEASE RUN FROM MAINDIR

## Pathways
# analysis_name="Prueba1_BS"
SCRIPTPATH=$(pwd)
data_path="${SCRIPTPATH}/data/raw/datos_librerias"
disc_path="${SCRIPTPATH}/data/use/1KGPtagSNPs"
sv_path="${SCRIPTPATH}/data/raw/SV_1KGP"

mkdir -p $disc_path
cd $disc_path

kgp_path="$disc_path"


## Inversions
# inversions=$(echo 'HsInv0003,HsInv0004,HsInv0006,HsInv0015,HsInv0016,HsInv0041,HsInv0052,HsInv0058,HsInv0059,HsInv0060,HsInv0063,HsInv0068,HsInv0081,HsInv0092,HsInv0095,HsInv0097,HsInv0098,HsInv0105,HsInv0156,HsInv0163,HsInv0164,HsInv0174,HsInv0186,HsInv0201,HsInv0260,HsInv0284,HsInv0379,HsInv0409,HsInv0960,HsInv0965,HsInv0991,HsInv1066,HsInv1075,HsInv1141,HsInv1153,HsInv1216,HsInv1269,HsInv1306,HsInv1314,HsInv1329,HsInv1376,HsInv1408,HsInv1471,HsInv1569,HsInv1614,HsInv1637,HsInv1681,HsInv1700,HsInv1751,HsInv1790' | sed -e 's/,/\t/g')
inversions=$( cat ${SCRIPTPATH}/data/raw/inversions_completeList.txt | sed -e 's/\n/\t/g')

for inv in $inversions
do

	cd $kgp_path/

	## 1KGP samples
	# samples_1kgp=$(cut -f1 ${SCRIPTPATH}/data/raw/1KGP_data/integrated_call_samples_v3.20130502.ALL.panel.txt | grep -v 'sample')

	## SMALL PANEL SAMPLSE
	samples_1kgp=$(cut -f1 ${SCRIPTPATH}/data/use/static_dataset_list.txt | grep -v 'sample')

	## From GTypes to VCF
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

	bcftools view -Ov -r $chr:${start}-${end} -s $(echo $samples_1kgp | tr " " ",") $sv_path/ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.vcf.gz  >$kgp_path/input_inversion.vcf
	grep 'INV' $kgp_path/input_inversion.vcf | awk '$5=="<INV>"' | awk -v OFS="\t" '$8="."' > $kgp_path/input_inversion_aux.vcf
	rm $kgp_path/input_inversion.vcf
	mv $kgp_path/input_inversion_aux.vcf $kgp_path/input_inversion.vcf


	inv_esp=$(cut -f3 $kgp_path/input_inversion.vcf)


	#### CHECK!
	echo '############'
	echo $inv
	echo $inv_esp
	echo '############'



	if [ "$inv_esp" != "" ]
	then

		## Inversion VCF into 1KGP VCF
		if [ $chr == "X" ]; then
			vcf="${SCRIPTPATH}/data/raw/1KGP_data/vcf/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz"
		else
			vcf="${SCRIPTPATH}/data/raw/1KGP_data/vcf/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
		fi

		bcftools view -Ov -r $chr:${start}-${end} -s $(echo $samples_1kgp | tr " " ",") ${vcf}  >$kgp_path/haplotypes.vcf


		cd $kgp_path
		cat haplotypes.vcf input_inversion.vcf  > haplotypes_with_HsInv_unsorted.vcf
		bcftools sort haplotypes_with_HsInv_unsorted.vcf > haplotypes_with_HsInv.vcf
		rm haplotypes.vcf input_inversion.vcf haplotypes_with_HsInv_unsorted.vcf
		mv haplotypes_with_HsInv.vcf haplotypes.vcf


		## Convert to plink format
		ulimit -n 3000
      	vcftools --vcf haplotypes.vcf --plink --chr $chr --out genotypes_plink_format --temp ${SCRIPTPATH}/tmp/

		number_invs=$(echo $inv_esp | tr ' ' '\n'| wc -l)

		if [ "$number_invs" == "1" ]
		then

			## LD calculation (diferente en X??)
			# folder_plink="/home/jon/soft/plink-1.07-x86_64"
			plink --file genotypes_plink_format --r2 --ld-snp $inv_esp --ld-window-kb 1000000 --ld-window 999999 --ld-window-r2 0 --noweb

			# Save
			mkdir $inv
			mv plink.ld $inv/plink.ld
			rm genotypes* haplotypes* plink.log plink.nosex
		else

			for inv_uniq in $(echo $inv_esp | tr ' ' '\n')
			do
				echo $inv_uniq

				## LD calculation (diferente en X??)
				# folder_plink="/home/jon/soft/plink-1.07-x86_64"
				plink --file genotypes_plink_format --r2 --ld-snp $inv_uniq --ld-window-kb 1000000 --ld-window 999999 --ld-window-r2 0 --noweb

				mkdir $inv
				mv plink.ld $inv/plink_$inv_uniq.ld

			done

			# Save
			rm genotypes* haplotypes* plink.log plink.nosex


		fi

	fi

done



# ###### HsInv0052 region of deletion ######
# inv=$(echo 'HsInv0052_region' | sed -e 's/,/\t/g')
# cd $kgp_path/


# ## 1KGP samples
# samples_1kgp=$(cut -f1 ${SCRIPTPATH}/data/raw/1KGP_data/integrated_call_samples_v3.20130502.ALL.panel.txt | grep -v 'sample')


# ## From GTypes to VCF
# if [ $inv == "HsInv0052_region" ]; then
# 	chr=$(grep 'HsInv0052' $data_path/bplib.coords | cut -f2 | uniq)
# 	pos_left=$(grep 'HsInv0052' $data_path/bplib.coords | cut -f3,4 | sed -e 's/\t/\n/g' | sort | head -n+1)
# 	pos_right=$(grep 'HsInv0052' $data_path/bplib.coords | cut -f3,4 | sed -e 's/\t/\n/g' | sort | tail -n-1)
# 	start=$(($pos_left-500000))
# 	end=$(($pos_right+500000))
# elif [ $inv == "HsInv0409" ]; then
# 	chr=$(grep 'HsInv409' $data_path/bplib.coords | cut -f2 | uniq)
# 	pos_left=$(grep 'HsInv409' $data_path/bplib.coords | cut -f3,4 | sed -e 's/\t/\n/g' | sort | head -n+1)
# 	pos_right=$(grep 'HsInv409' $data_path/bplib.coords | cut -f3,4 | sed -e 's/\t/\n/g' | sort | tail -n-1)
# 	start=$(($pos_left-500000))
# 	end=$(($pos_right+500000))
# else
# 	chr=$(grep $inv $data_path/bplib.coords | cut -f2 | uniq)
# 	pos_left=$(grep $inv $data_path/bplib.coords | cut -f3,4 | sed -e 's/\t/\n/g' | sort | head -n+1)
# 	pos_right=$(grep $inv $data_path/bplib.coords | cut -f3,4 | sed -e 's/\t/\n/g' | sort | tail -n-1)
# 	start=$(($pos_left-500000))
# 	end=$(($pos_right+500000))
# fi

# bcftools view -Ov -r $chr:${start}-${end} -s $(echo $samples_1kgp | tr " " ",") $sv_path/ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.vcf.gz  >$kgp_path/input_inversion.vcf
# grep 'SVTYPE=DEL' $kgp_path/input_inversion.vcf  >$kgp_path/input_inversion_aux.vcf
# rm $kgp_path/input_inversion.vcf
# mv $kgp_path/input_inversion_aux.vcf $kgp_path/input_inversion.vcf



# ## Inversion VCF into 1KGP VCF
# if [ $chr == "X" ]; then
# 	vcf="${SCRIPTPATH}/data/raw/1KGP_data/vcf/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz"
# else
# 	vcf="${SCRIPTPATH}/data/raw/1KGP_data/vcf/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
			
# fi

# bcftools view -Ov -r $chr:${start}-${end} -s $(echo $samples_1kgp | tr " " ",") ${vcf}  >$kgp_path/haplotypes.vcf

# inv_esp=$(cut -f3 $kgp_path/input_inversion.vcf)


# #### CHECK!
# echo '############'
# echo $inv
# echo $inv_esp
# echo '############'

# if [ "$inv_esp" != "" ]
# then

# 	cd $kgp_path
# 	cat haplotypes.vcf input_inversion.vcf  > haplotypes_with_HsInv_unsorted.vcf
# 	vcf-sort haplotypes_with_HsInv_unsorted.vcf > haplotypes_with_HsInv.vcf
# 	rm haplotypes.vcf input_inversion.vcf haplotypes_with_HsInv_unsorted.vcf
# 	mv haplotypes_with_HsInv.vcf haplotypes.vcf


# 	## Convert to plink format
# 	ulimit -n 3000
# 	vcftools --vcf haplotypes.vcf --plink --chr $chr --out genotypes_plink_format

# 	number_invs=$(echo $inv_esp | tr ' ' '\n'| wc -l)

# 	if [ "$number_invs" == "1" ]
# 	then

# 		## LD calculation (diferente en X??)
# 		# folder_plink="/home/jon/soft/plink-1.07-x86_64"
# 		plink --file genotypes_plink_format --r2 --ld-snp $inv_esp --ld-window-kb 1000000 --ld-window 999999 --ld-window-r2 0 --noweb

# 		# Save
# 		mkdir $inv
# 		mv plink.ld $inv/plink.ld
# 		rm genotypes* haplotypes* plink.log plink.nosex
# 	else

# 		for inv_uniq in $(echo $inv_esp | tr ' ' '\n')
# 		do
# 			echo $inv_uniq

# 			## LD calculation (diferente en X??)
# 			# folder_plink="/home/jon/soft/plink-1.07-x86_64"
# 			plink --file genotypes_plink_format --r2 --ld-snp $inv_uniq --ld-window-kb 1000000 --ld-window 999999 --ld-window-r2 0 --noweb

# 			mkdir $inv
# 			mv plink.ld $inv/plink_$inv_uniq.ld

# 		done

# 		# Save
# 		rm genotypes* haplotypes* plink.log plink.nosex


# 	fi



# fi




