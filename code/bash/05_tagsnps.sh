#!/bin/bash
# Ruth Gómez Graciani
# 02 06 2021

###############################################################################
# Description:                                                                
# Run breakseq analysis        
# This script assumes stable location of 1KGP data                    
###############################################################################

# GET VARIABLES 
# =========================================================================== #

TMPDIR=$1     
OUTDIR=$2  
REGIONS_FILE=$3   
FILTER=$4
GENOTYPES_FILE=$5
data_path=$6
SCRIPTPATH=$7

#  INITIALIZE PROCESS
# =========================================================================== #

# Make output dirs
mkdir -p $OUTDIR $TMPDIR

#  Initialize output
echo -e "INV\tTYPE\tN_SAMPLES_F${FILTER}\tmax_SNP_${FILTER}\tLD_max_SNP${FILTER}" > ${OUTDIR}/tagSNPs_max.txt
echo -e "INV\tFILTER\tCHR\tPOS\tID\tLD" > ${OUTDIR}/tagSNPs_gt08.txt


# Take regions
REGIONS=$(cut -f1 $REGIONS_FILE)


for inv in $REGIONS;   do

  echo "tagSNPs for $inv"

  #  CHECK POLYMORPHISM
  # =========================================================================== #

  # Copy inversion file of gtypes
  mkdir -p ${TMPDIR}/$inv
  grep $inv $GENOTYPES_FILE | awk '$4<=2 && $4>0{print}'  >  ${TMPDIR}/${inv}/GTypes_FinalDataSet.txt

  # Is it polymorphic?
  gtype_counter=$(cut -f3 ${TMPDIR}/${inv}/GTypes_FinalDataSet.txt | sort | uniq | wc -l)
    
  if [ $gtype_counter -eq 1 ]; then
    # Monomorphic!
    echo -e "${inv}\tMonomorphic\tNA\tNA\tNA" >> ${OUTDIR}/tagSNPs_max.txt

  else
    # Polymorphic!
  
    #  PROCESS POLYMORPHIC REGION - Include analysis regions into 1KGP VCF
    # =========================================================================== #

    # 1KGP samples - this is based on a static structure of the files
    samples_1kgp=$(cut -f1 data/raw/1KGP_data/integrated_call_samples_v3.20130502.ALL.panel.txt | grep -v 'sample')

    ##### OldCode  #####
    # if [ $inv == "HsInv0052_region" ]; then
    #   chr=$(grep 'HsInv0052' $data_path/bplib.coords | cut -f2 | uniq)
    #   pos_left=$(grep 'HsInv0052' $data_path/bplib.coords | cut -f3,4 | sed -e 's/\t/\n/g' | sort | head -n+1)
    #   pos_right=$(grep 'HsInv0052' $data_path/bplib.coords | cut -f3,4 | sed -e 's/\t/\n/g' | sort | tail -n-1)
    #   start=$(($pos_left-1500000))
    #   end=$(($pos_right+1500000))
    # elif [ $inv == "HsInv0409" ]; then
    #   chr=$(grep 'HsInv409' $data_path/bplib.coords | cut -f2 | uniq)
    #   pos_left=$(grep 'HsInv409' $data_path/bplib.coords | cut -f3,4 | sed -e 's/\t/\n/g' | sort | head -n+1)
    #   pos_right=$(grep 'HsInv409' $data_path/bplib.coords | cut -f3,4 | sed -e 's/\t/\n/g' | sort | tail -n-1)
    #   start=$(($pos_left-1500000))
    #   end=$(($pos_right+1500000))
    # else

    # Region coordinates
    CHR_REGION=$(grep "$inv" $data_path/bplib.coords | cut -f2 | uniq)
    LEFT_REGION=$(grep  "$inv" $data_path/bplib.coords | cut -f3,4 | sed -e 's/\t/\n/g' | sort | head -n+1)
    RIGHT_REGION=$(grep  "$inv" $data_path/bplib.coords | cut -f3,4 | sed -e 's/\t/\n/g' | sort | tail -n-1)

    # Analysis coordinates
    START_REGION=$(($LEFT_REGION-1500000))
    END_REGION=$(($RIGHT_REGION+1500000))

    # fi
    #### End of  OldCode #####

    # Set 1KGP VCF name; it is different depending on chromosome!
    if [ $CHR_REGION == "X" ]; then
       vcf="data/raw/1KGP_data/vcf/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz"
    else
       vcf="data/raw/1KGP_data/vcf/ALL.chr${CHR_REGION}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    fi
    
    # Select analysis coordinates from 1KGP VCF
    bcftools view -Ov -r $CHR_REGION:${START_REGION}-${END_REGION} --min-ac 2:minor -M2 -m2 ${vcf} > ${TMPDIR}/${inv}/haplotypes_unsorted.vcf
    
    # Process predicted genotypes to match 1KGP VCF format
    genotypes_std_inv=$(for sample in $samples_1kgp; do if grep -q $sample ${TMPDIR}/${inv}/GTypes_FinalDataSet.txt; then grep $sample ${TMPDIR}/${inv}/GTypes_FinalDataSet.txt | cut -f3; else echo "NA"; fi; done)
    # genotypes_GT_GL=$(echo $genotypes_std_inv | sed 's/STD/Std/g' | sed 's/REF/Std/g' | sed 's/INS/Inv/g' | sed 's/DEL/Inv/g' | sed 's/INV/Inv/g' | sed 's/HET/Het/g' | sed 's/Std/0\/0/g' | sed 's/Inv/1\/1/g' | sed 's/Het/0\/1/g' | sed 's/ND/.\/./g' | sed 's/NA/.\/./g' )C
    genotypes_GT_GL=$(echo $genotypes_std_inv | sed 's/STD\|REF/Std/g' | sed 's/INS\|INV\|ALT\|ANC\|DEL/Inv/g' | sed 's/Std/0/g' | sed 's/Inv/1/g' |  sed 's/NA/.\/./g')

    if [ $CHR_REGION == "X" ]; then
      genotypes_GT_GL_final=$(paste <(echo $genotypes_GT_GL | tr ' ' '\n') <(for sample in $samples_1kgp; do grep $sample data/raw/1KGP_data/integrated_call_samples_v3.20130502.ALL.panel.txt | cut -f4; done) | sed 's,0/0\tmale,0\tmale,g' | sed -e 's,1/1\tmale,1\tmale,g' | sed -e 's,./.\tmale,.\tmale,g' | cut -f1)
      genotypes_GT_GL=$genotypes_GT_GL_final
    fi

    # Append haplotype to 1KGP VCF file
      echo -e "$CHR_REGION\t${LEFT_REGION}\t${inv}\tA\tT\t.\tPASS\t.\tGT\t$(echo ${genotypes_GT_GL} | tr ' ' '\t')" >> ${TMPDIR}/${inv}/haplotypes_unsorted.vcf
      echo -e "$CHR_REGION\t${RIGHT_REGION}\t${inv}\tA\tT\t.\tPASS\t.\tGT\t$(echo ${genotypes_GT_GL} | tr ' ' '\t')" >> ${TMPDIR}/${inv}/haplotypes_unsorted.vcf


    # Sort combined file
    bcftools sort ${TMPDIR}/${inv}/haplotypes_unsorted.vcf > ${TMPDIR}/${inv}/haplotypes.vcf
    
    # Remove temporary files
    rm ${TMPDIR}/${inv}/haplotypes_unsorted.vcf  


    #  PROCESS POLYMORPHIC REGION - Filter samples
    # =========================================================================== #

    samples_filter=$(cat ${TMPDIR}/${inv}/GTypes_FinalDataSet.txt | awk -v filter=$FILTER '$5>filter' | cut -f2)
      
    # Check genotypes!! maybe now that is filtered it is monomorphic
    samples_filter_geno=$( cat ${TMPDIR}/${inv}/GTypes_FinalDataSet.txt | awk -v filter=$FILTER '$5>filter' | cut -f3 | sort | uniq | wc -l)
    
    if [ $samples_filter_geno -eq 1 ]; then
        samples=$( echo $samples_filter | wc -w)
        echo -e "${inv}\tPolymorphic\t${samples}\tND\tND" >> ${OUTDIR}/tagSNPs_max.txt
    else
        bcftools view -Oz -s $(echo $samples_filter | tr " " ",") --force-samples ${TMPDIR}/${inv}/haplotypes.vcf  > ${TMPDIR}/${inv}/haplotypes_filtered.vcf.gz

        #  PROCESS POLYMORPHIC REGION - LD calculation
        # =========================================================================== #

        # Convert to plink format
        ulimit -n 3000
        vcftools --gzvcf ${TMPDIR}/${inv}/haplotypes_filtered.vcf.gz --plink --chr $CHR_REGION --out ${TMPDIR}/${inv}/genotypes_plink_format_filtered

        # Run plink
        cd ${TMPDIR}/${inv}/
        plink --file genotypes_plink_format_filtered --r2 --ld-snp $inv --ld-window-kb 1000000 --ld-window 999999 --ld-window-r2 0 --noweb
        cd ${SCRIPTPATH}

    #  PROCESS POLYMORPHIC REGION - Make summary files
    # =========================================================================== #

    # A list with max LD SNP for each region, excluding itself  
    snpname=$(awk -v inv=$inv '$6!=inv' ${TMPDIR}/${inv}/plink.ld |  grep -v 'esv' | sed -e 's/^ *//g' | sed -e 's/ *$//g' | sed 's/ \+/\t/g'  | sort -g  -k7 | grep -v 'R2' | tail -n-1 | cut -f6 )
    snpld=$(awk -v inv=$inv '$6!=inv' ${TMPDIR}/${inv}/plink.ld |  grep -v 'esv' | sed -e 's/^ *//g' | sed -e 's/ *$//g' | sed 's/ \+/\t/g'  | sort -g  -k7 | grep -v 'R2' | tail -n-1 | cut -f7 )
    samples=$( echo $samples_filter | wc -w)
    echo -e "${inv}\tPolymorphic\t${samples}\t${snpname}\t${snpld}" >> ${OUTDIR}/tagSNPs_max.txt
    
    # A list with all SNPs having LD>0.8   
    awk -v inv=$inv '$6!=inv' ${TMPDIR}/${inv}/plink.ld | awk '$7>0.8' | grep -v 'esv' | sed -e 's/^ *//g' | sed -e 's/ *$//g' | sed 's/ \+/\t/g'  | sort -g  -k7 | grep -v 'R2' | awk -v OFS="\t"  -v filter=$FILTER '{print $3, filter, $6, $7}' >>  ${OUTDIR}/tagSNPs_gt08.txt

    fi
      
  fi 

    #  CLEANUP
    # =========================================================================== #

    rm -r ${TMPDIR}/$inv/

done