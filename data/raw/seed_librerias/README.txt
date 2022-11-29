CARPETA DATOS_LIBRERIAS

Generar librerías:
1) Extraemos las sondas de las páginas bplibInvs y bplibOthers del archivo excel SimpleBPLibrary_LergaJaso2021_300bp, y las ponemos todas en un archivo fasta bplib.fa
2) Extraemos las coordenadas de la página bplib_hg19coords y las ponemos en un archivo bplib.coords
Remove excluded probes marked with  sed "/#.*$/d" or similar
+ copiar ref.header-template -- cambiando el sn:1, etc por SN:chr1
+ hacer un specs.sh

OPCIONAL
7) Versión a otros assemblies
La manera más fácil para tan pocas anotaciones es https://genome.ucsc.edu/cgi-bin/hgLiftOver
awk '{print "chr"$2, $3, $4, $1}' bplib.coords > liftOver_input.bed
Pasar por el liftover, guardar como liftOver_output.bed
Re-transformar a formato orgiginal - SIN CAMBIAR EL CHR para hg38
Cambio el nombre de bplib.coords original a bplib_v37.coords
awk -v OFS="\t" '{print $4, $1, $2, $3, "300"}' liftOver_output.bed > bplib.coords


# --- automatized ------# 
3) Generar header para las sondas:

###

cd /home/jon/Desktop/BreakSeq_FunctionalPaper/Data/datos_librerias
# OJO CAMBIAR PARA LN150!!!!
grep '>' bplib.fa | sed -e 's/>//g' | sed -e 's/^/@SQ\tSN:/g' | sed -e 's/$/\tLN:300/g' > bplib.header.template.txt

echo '@HD VN:1.0 SO:unsorted' | tr ' ' '\t' > tmp
echo '@RG ID:XXX SM:XXX' | tr ' ' '\t' >> tmp
echo '@PG ID:breakseq' | tr ' ' '\t' >> tmp
cat bplib.header.template.txt >> tmp
rm bplib.header.template.txt
mv tmp bplib.header.template.txt

###
4) Generar lista de todas las inversiones:

###

cat  bplib.coords | cut -f1 | sed 's/^\w\w\w//'|sed 's/BP.*$//' | sort | uniq > inversions_completeList.txt

###

5) Copiar otros datos:
- ref.header.template
Header para el genoma hg19 (reads de 1KGP).
- [bplib.fa.fai
Fasta index. Se genera al correr el BreakSeq.]

6) Cambiar los indices a las necesidades actuales. Para ello:

####
singularity shell --bind /data/bioinfo/common/bowtie2_index:/nfs/pic.es/user/r/rgomez/20210325_breakseq/data/use/bowtie_index:rw /data/bioinfo/software/bgd-pic_breakseq.sif 

cd 20210325_breakseq/data/use/bowtie_index
cp ../../raw/datos_librerias/bplib.fa ./

bowtie2-build bplib.fa bplib
###


#### FOR CURRENT CODE: 
Step 1
Step 2
Step 5
Run 00Library commands from runCommands

