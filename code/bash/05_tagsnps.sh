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
REGIONS_FILE=$2  
FILTER=$3 
GENOTYPES_FILE=$4
REF_INFO=$5
LIBRARY_INFO=$6 
SCRIPTPATH=$7


#  INITIALIZE PROCESS
# =========================================================================== #

# Make output dirs
mkdir -p $TMPDIR

# Take regions
REGIONS=$REGIONS_FILE


for inv in $REGIONS;   do

  echo "tagSNPs for $inv"


  #  CHECK POLYMORPHISM
  # =========================================================================== #

  # Take those individuals with 1 or 2 reported haplotypes, and p.error <= filter
  mkdir -p ${TMPDIR}/$inv
  grep $inv $GENOTYPES_FILE | awk -v perrfil="${FILTER}" '$7 <= perrfil && $5<=2 && $5>0{print}' >  ${TMPDIR}/${inv}/GTypes_FinalDataSet.txt

  # Is it polymorphic? (count different genotypes, it can be INV, STD, INV/STD or smth like that)
  gtype_counter=$(cut -f4 ${TMPDIR}/${inv}/GTypes_FinalDataSet.txt | sort | uniq | wc -l)
    
  if [ $gtype_counter -lt 2 ]; then
    echo "Monomorphic!"
    #I want to know if I have samples at all
    sampleCount=$(wc -l <  ${TMPDIR}/${inv}/GTypes_FinalDataSet.txt  )
    echo -e "${inv}\tMonomorphic\t${sampleCount}\tNA\tNA\tNA\tNA" >  ${TMPDIR}/${inv}/tagSNPs_max.txt

  else
    echo "Polymorphic!"
  
    #  PROCESS POLYMORPHIC REGION - Include analysis regions into 1KGP VCF
    # =========================================================================== #
    echo "Taking inversion coordinates..."
    # Take inversion coordinates
      # Region coordinates
      CHR_REGION=$(grep "$inv" $LIBRARY_INFO/bplib.coords | cut -f2 | uniq)
      CHR_NUMBER=$( sed 's/chr//g' <<<  $CHR_REGION )
      LEFT_REGION=$(grep  "$inv" $LIBRARY_INFO/bplib.coords | cut -f3,4 | sed -e 's/\t/\n/g' | sort | head -n+1)
      RIGHT_REGION=$(grep  "$inv" $LIBRARY_INFO/bplib.coords | cut -f3,4 | sed -e 's/\t/\n/g' | sort | tail -n-1)

      # Analysis coordinates
      START_REGION=$(($LEFT_REGION-1500000))
      END_REGION=$(($RIGHT_REGION+1500000))

    echo "The region is $CHR_REGION:${START_REGION}-${END_REGION}"
    # Identify vcf file
      VCF_FILE=$(ls $REF_INFO | egrep "chr${CHR_NUMBER}[^0-9].*gz$")
      echo "Selecting region from VCF file in $REF_INFO/$VCF_FILE"

    # Extract samples
      samples_1kgp=$(bcftools query -l $REF_INFO/$VCF_FILE)

    # Select analysis coordinates from 1KGP VCF
      bcftools view -Ov -r $CHR_REGION:${START_REGION}-${END_REGION} --min-ac 2:minor -M2 -m2 $REF_INFO/$VCF_FILE > ${TMPDIR}/${inv}/haplotypes_unsorted.vcf
    
    # Process predicted genotypes to match 1KGP VCF format
    echo "Appending genotypes..."
    genotypes_std_inv=$(for sample in $samples_1kgp; do
      if grep -q $sample ${TMPDIR}/${inv}/GTypes_FinalDataSet.txt ; then 
        genotype=$(grep $sample ${TMPDIR}/${inv}/GTypes_FinalDataSet.txt | cut -f4)
        if [ $(echo ${#genotype}) -eq 3 ] ;  then
          echo "${genotype}/${genotype}"
        else 
          echo ${genotype}
        fi
      else
        echo "./."
      fi
    done
    )

    genotypes_GT_GL=$(echo $genotypes_std_inv | sed 's/STD\|REF/Std/g' | sed 's/INS\|INV\|ALT\|ANC\|DEL/Inv/g' | sed 's/Std/0/g' | sed 's/Inv/1/g')

    # THIS DOES NOT TAKE INTO ACCOUNT Y CHROMOSOME!
    # if [ $CHR_REGION == "X" ]; then
    #   genotypes_GT_GL_final=$(paste <(echo $genotypes_GT_GL | tr ' ' '\n') <(for sample in $samples_1kgp; do grep $sample data/raw/1KGP_data/integrated_call_samples_v3.20130502.ALL.panel.txt | cut -f4; done) | sed 's,0/0\tmale,0\tmale,g' | sed -e 's,1/1\tmale,1\tmale,g' | sed -e 's,./.\tmale,.\tmale,g' | cut -f1)
    #   genotypes_GT_GL=$genotypes_GT_GL_final
    # fi

    # Append haplotype to 1KGP VCF file
      echo -e "$CHR_REGION\t${LEFT_REGION}\t${inv}\tA\tT\t.\tPASS\t.\tGT\t$(echo ${genotypes_GT_GL} | tr ' ' '\t')" >> ${TMPDIR}/${inv}/haplotypes_unsorted.vcf
      echo -e "$CHR_REGION\t${RIGHT_REGION}\t${inv}\tA\tT\t.\tPASS\t.\tGT\t$(echo ${genotypes_GT_GL} | tr ' ' '\t')" >> ${TMPDIR}/${inv}/haplotypes_unsorted.vcf


    # Sort combined file
    echo "Sorting vcf..."
    bcftools sort ${TMPDIR}/${inv}/haplotypes_unsorted.vcf > ${TMPDIR}/${inv}/haplotypes.vcf
    
    # Remove temporary files
    rm ${TMPDIR}/${inv}/haplotypes_unsorted.vcf  

    #  PROCESS POLYMORPHIC REGION - LD calculation
    # =========================================================================== #
    # echo "Making pink format..."
    # Convert to plink format
    # ulimit -n 3000
    # vcftools --vcf ${TMPDIR}/${inv}/haplotypes.vcf --plink --chr $CHR_REGION --out ${TMPDIR}/${inv}/genotypes_plink_format_filtered
    
    echo "LD calculation..." 
    # Run plink
    cd ${TMPDIR}/${inv}/
    plink --vcf haplotypes.vcf --r2 --ld-snp $inv --ld-window-kb 1000000 --ld-window 999999 --ld-window-r2 0 --noweb
    cd ${SCRIPTPATH}

    #  PROCESS POLYMORPHIC REGION - Make summary files
    # =========================================================================== #
    echo "Summarizing..."
    # A list with max LD SNP for each region, excluding itself  
    # echo -e "INV\tTYPE\tN_SAMPLES_F${FILTER}\tchr_max_SNP_${FILTER}\tpos_max_SNP_${FILTER}\tname_max_SNP_${FILTER}\tLD_max_SNP${FILTER}" 
    snpinfo=$(awk -v inv=$inv '$6!=inv' ${TMPDIR}/${inv}/plink.ld | sort -g -k7 | tail -n-1 | awk -v OFS="\t"  '{print $4, $5, $6, $7}')
    samples=$(cat ${TMPDIR}/${inv}/GTypes_FinalDataSet.txt | wc -l )
    echo -e "${inv}\tPolymorphic\t${samples}\t${snpinfo}" >  ${TMPDIR}/${inv}/tagSNPs_max.txt
    
    # A list with all SNPs having LD>0.8       
    # echo -e "INV\tFILTER_genotypes\tFILTER_LD\tCHR\tPOS\tID\tLD"  
    awk -v OFS="\t" -v filterA=$FILTER -v filterB=0.8 -v inv=$inv  '$6!=inv && $7>=filterB {print $3, filterA, filterB, $4, $5, $6, $7}' ${TMPDIR}/${inv}/plink.ld  | sort -g -k7 | grep -v "R2" > ${TMPDIR}/${inv}/tagSNPs_gt08.txt

          
  fi 

    #  CLEANUP
    # =========================================================================== #

    rm -r ${TMPDIR}/$inv/*.vcf ${TMPDIR}/$inv/plink* ${TMPDIR}/$inv/GTypes_FinalDataSet* 

done