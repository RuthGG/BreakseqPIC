#!/bin/bash
#
# This will generate the following files, unless they already exist:
#
#  - $DATADIR/bplip_noA.fa (NOT overwritten)
#  - $DATADIR/bplip_noA.fa.fai (NOT overwritten)
#  - ref.header.template (NOT overwritten)
#  - $i/ref.header (NOT overwritten)
#  - $i/header (NOT overwritten)
#  - $i/uni.sorted.bam (NOT overwritten)
#  - $i/uni.sorted.bam.bai (NOT overwritten)
#  - $i/xun.sorted.bam (NOT Overwritten)
#  - $i/xun.sorted.bam.bai (NOT overwritten)
#
# The following files will be created or updated (overwritten):
#  - $i/ref.sorted.sam
#  - $i/uni.sorted.bam
#  - $i/uni.sorted.bam.bai
#  - $i/xun.sorted.bam
#  - $i/xun.sorted.bam.bai
#
# I think this is to generate bam files and filter out from
# xun the reads that map to non-bp sites of the reference.
# 

if [ $# -lt 1 ]; then
	echo "Usage: makebam.sh <dir> <ref.fasta>"
	echo "     <dir>       Current folder, e.g. '2013-01-18', where individual's directories must be."
	echo "     <ref.fasta> Path to the fasta file of the reference genome [.../data/2012-09-07/human_g1k_v37.fasta]."
	exit
fi

if [ $# -eq 1 ]; then
	REF_FASTA=/home/shareddata/Bioinformatics/BPSeq/data/2012-09-07/human_g1k_v37.fasta
else
	REF_FASTA=$2
fi

DATADIR=/gpfs/depoofs/jlerga/Imputation/Breakpoints/datos

# MIN is the minimum length of reads to be considered. It may be important to
# calculate allele frequencies later.
MIN=149

if [ ! -e $DATADIR/bplib_noA.fa ]; then
        date
        echo "Creating bplib_noA.fa."
        echo
        gawk '(/^>/){gsub(/:A$/,""); print}(/^[^>]/){print}' $DATADIR/bplib.fa > $DATADIR/bplib_noA.fa
fi

if [ ! -e $DATADIR/bplib_noA.fa.fai ]; then
        date
        echo "Indexing bplib_noA.fa."
        echo
        samtools faidx $DATADIR/bplib_noA.fa
fi

#if [ ! -e ref.header.template ]; then
#        date
#        echo "Creating the header template for alignments to the reference."
#        echo
#        echo -e "@HD\tVN:1.0\tSO:unsorted" > ref.header.template
#        echo -e "@RG\tID:XXX\tSM:XXX" >> ref.header.template
#        echo -e "@PG\tID:breakseq" >> ref.header.template
#        /gpfs/depoofs/jlerga/soft/ActivePerl/perl/bin/perl -we 'use Bio::SeqIO; $ref = Bio::SeqIO->new(-file => $ARGV[0], -format => "fasta");
#                while ($chr = $ref->next_seq) {
#                        print "\@SQ\tSN:", $chr->id, "\tLN:", $chr->length, "\n";
#                }' $REF_FASTA >> ref.header.template
#fi

for j in `cat /gpfs/depoofs/jlerga/Imputation/Breakpoints/test/ind_test.txt`; do

	i=$(echo $j | cut -c5-1000)


        if [ -d $i ] && [ ! -e $i/ref.header ]; then
                date
                echo "Creating header for $i."
                echo
                gawk -v ID=$i '{gsub(/XXX/, ID); print}' ref.header.template > $i/ref.header
        fi

        if [ -d $i ] && ls -1 $i | grep -q -P "ref\.sam$"; then
                #if [ ! -e $i/ref.sorted.sam ]; then
                        date
                        echo "Collecting and sorting alignments to the reference from individual $i."
                        echo
                        gawk '(FILENAME ~ /ref.header$/){print}((FILENAME ~ /sam$/) && (/^[^@]/)){print}' $i/ref.header $i/*.ref.sam > $i/ref.sam
                        samtools view -Sb $i/ref.sam > $i/ref.bam
                        samtools sort -n $i/ref.bam -o $i/ref.sorted.bam
                        samtools view $i/ref.sorted.bam > $i/ref.sorted.sam
                        rm $i/ref.sam $i/ref.sorted.bam $i/ref.bam
                #fi
        else
                echo "There are no hits to the human reference genome for $i."
        fi

        #
        #        date
        #        echo "Creating more headers for $i."
        #        echo
        #        echo -e "@HD\tVN:1.0\tSO:unsorted" > $i/header
        #        echo -e "@RG\tID:$i\tSM:$i" >> $i/header
        #        echo -e "@PG\tID:breakseq" >> $i/header
        #        perl -we 'use Bio::SeqIO; $bplib = Bio::SeqIO->new(-file => $ARGV[0], -format => "fasta");
        #                while ($bp = $bplib->next_seq) {
        #                        print "\@SQ\tSN:", $bp->id, "\tLN:", $bp->length, "\n";
        #                }' $DATADIR/bplib_noA.fa >> $i/header
	#
        #fi

	if [ -d $i ] && [ ! -e $i/header ]; then

                date
                echo "Creating more headers for $i."
                echo
                gawk -v ID=$i '{gsub(/XXX/, ID); print}' header > $i/header

	fi


        # It seems reasonable to eliminate *ini.sam files, but I'll wait.
        #find $i -name '*.ini.sam' -exec rm '{}' \;

        if [ -d $i ] && ls -1 $i | grep -q -P "uni\.sam$"; then
                #if [ ! -e $i/uni.sorted.bam ]; then
                        date
                        echo "Collecting and sorting unique alignments to bplib from individual $i."
                        echo
                        gawk -v RG=$i -v OFS="\t" -v MIN=$MIN '(FILENAME ~ /header$/){print}((FILENAME ~ /sam$/) && (/^[^@]/) && (length($10) >= MIN)){
                                gsub(/:A$/,"",$3)
                                print $0 "\tRG:Z:" RG}' $i/header $i/*.uni.sam > $i/uni.sam
                        samtools view -Sb $i/uni.sam > $i/uni.bam
                        samtools sort $i/uni.bam -o $i/uni.sorted.bam
                        rm $i/uni.sam $i/uni.bam
                        samtools index $i/uni.sorted.bam
                #fi
        else
                echo "There are no unique hits to breakpoints for $i."
        fi

        if [ -d $i ] && ls -1 $i | grep -q -P "xun.sam$"; then
		# Reads in xun.sam files map on breakpoints and also on the reference genome. For them
		# to be excluded as evidence of an alternative allele, they don't need to map with high
		# accuracy, nor uniquely. That's why lastly we are not using quality filters in breakseq.
		# However, for them to count as positive evidence of the reference allele, they need to
		# map to the coordinates of the reference breakpoint (not just anywhere in the genome)
		# and with some minimum quality.


                date
                echo "Collecting, filtering, and sorting non-unique alignments to bplib from individual $i."
                echo
                gawk -v RG=$i -v OFS="\t" -v MIN=$MIN '(FILENAME ~ /header$/){print}((FILENAME ~ /sam$/) && (/^[^@]/) && (length($10) >= MIN)){
                           gsub(/:A$/, "", $3)
                           print $0 "\tRG:Z:" RG}' $i/header $i/*.xun.sam > $i/xun.sam
		# Below, I am not filtering reads that happen in the ref bam file more than once, because
		# I want to respect the paired-ends, and because ambiguously mapped reads may have already
		# been filtered by the required minimum mapping quality of 15.
                perl -we 'use strict; use warnings;
                                open(OUT, ">$ARGV[0]/xun.filtered.sam") || die "I cannot open xun.filtered.sam.\n";
                                open(BP, $ARGV[1]) || die "I cannot open bplib.coords.\n";
                                my (%BPchr, %BPstart, %BPend);
				my $MinQual = 15;
                                while (<BP>) {
                                        chomp;
                                        my @line = split /\t/, $_;
                                        $BPchr{$line[0]} = $line[1];
                                        $BPstart{$line[0]} = $line[2];
                                        $BPend{$line[0]} = $line[3];
                                }
                                close BP;
                                open(REF, "$ARGV[0]/ref.sorted.sam") || die "I cannot open ref.sorted.sam.\n";
                                #my (%Freq, %Chr, %Pos, %Exclude);
				my (%Chr, %Pos, %Exclude);
                                my $num = 0;
                                while (<REF>) {
                                        $num++;
                                        my @line = split /\t/, $_;
                                        #$Freq{$line[0]}++;
                                        #if ( $Freq{$line[0]} > 1 ) { $Exclude{$line[0]} = 1 }
                                        my $good = 0;
					if ($line[4] >= $MinQual) {
	                                        for my $k (keys %BPchr) {
	                                                if (($line[3] >= $BPstart{$k} + 10) && ($line[3] <= $BPend{$k}) && ($line[2] eq $BPchr{$k})) {
	                                                        $good = 1;
	                                                }
	                                        }
					}
                                        if ($good == 0) { $Exclude{$line[0]} = 1 }
                                }
                                my $excluded = scalar keys %Exclude;
                                print "On individual $ARGV[0], out of $num reads, $excluded have been excluded.\n";
                                close REF;
                                open(XUN, "$ARGV[0]/xun.sam") || die "I cannot open xun.sam.\n";
                                while (<XUN>) {
                                        if (/^@/) {print OUT $_}
                                        else {
                                                my @line = split /\t/, $_;
                                                unless ((exists $Exclude{$line[0]}) || ($line[2] =~ /^INV/) || !($line[2] =~ /REF/)) {print OUT $_}
                                        }
                                }
                                close XUN;' $i $DATADIR/bplib.coords

                if grep -q -v -P "^@" $i/xun.filtered.sam; then
                           samtools view -Sb $i/xun.filtered.sam > $i/xun.filtered.bam
                           samtools sort $i/xun.filtered.bam -o $i/xun.sorted.bam
                           rm $i/xun.filtered.bam
                           samtools index $i/xun.sorted.bam
                else
                           echo "Indeed, it seems that all reads were excluded..."
                fi
		rm $i/xun.sam $i/xun.filtered.sam

                
        else
                echo "There are no non-unique hits to the breakpoints for $i."
        fi

done

#ME=`whoami`
#find . -user $ME ! -perm -020 -exec chmod g+w '{}' \;

