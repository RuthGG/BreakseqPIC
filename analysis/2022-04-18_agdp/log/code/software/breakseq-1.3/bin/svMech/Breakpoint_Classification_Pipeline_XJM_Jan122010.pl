#!/usr/bin/perl -w
use strict;

use FindBin '$Bin';

#ARGV[0]: deletionInfile; ARGV[1]:InsertionInfile; ARGV[2]:InversionInfile; 
#ARGV[3]: VNTRcutoff(0-1; ARGV[4]: FlankingNAHR(100-1000,even number); ARGV[5]: NAHRhomologyLength(20-400); ARGV[6]: NAHRpct(0-100); ARGV[7]: RTwindow(10-1000,even number)
#ARGV[8]: RTgap:(0-1000)
#ARGV[9] = $label
#ARGV[10]: output directory

#default parameters used in Breakseq
#del ins inv 0.5 200 50 85 200 150

unless ($ARGV[10]) {
print <<EOF;
>This program runs the SV Breakpoint classification pipeline.
>Input files should be three files:  
[params: del_file ins_file inv_file VNTRcutoff(0-1) FlankingNAHR(100-1000,even number) NAHRhomologyLength(20-400) NAHRpct(0-100) RTwindow(10-400,even number) RTgap:(0-1000) label out_dir]

>File 1: "Text_Files/CoordinatesDeletions.txt" with first two lines as fomatted as following:
>File 1: Line 1:"Chr No.	Start	End	ID"
>File 1: Line 2:"1	61886002	61892108	LIB000000001"
>...

>File 2: "Text_Files/CoordinatesInsertions.txt" with first two lines as fomatted as following:
>File 2: Line 1:"Chr No.	Start	End	InsertedSequence	ID"
>File 2: Line 2:"1       1574409 1574409 gctccaccttccgggttcacgccattctcctgcctcagcctcccaagtagctggga	LIB000000002"
>...

>File 3: "Text_Files/CoordinatesInversions.txt" with first two lines as fomatted as following:
>File 3: Line 1:"Chr No.	Start	End	ID"
>File 3: Line 2:"2	138723407	138725276	LIB000000003"
>...

>>> label will be used to name datasets...


EOF
exit;
}
my $Deletions;
my $Insertions;
my $Inversions;
my $CUTOFF = 20000000;  #20MB, events larger than $CUTOFF are not run repeatmasker on, if not NAHR or unsure 3 or 4, going to NHR
my $SIZE = 50; #SV <50bp will not be classified

my $VNTRcutoff = $ARGV[3];
my $FlankNAHR = $ARGV[4];
my $NAHRhomolen = $ARGV[5];
my $NAHRpct = $ARGV[6];
my $RTwin = $ARGV[7];
my $RTgap = $ARGV[8];
my $Window = $RTwin;
if($RTwin>$FlankNAHR){
	$Window = $NAHRhomolen;
}

my $STATUS;

$Deletions = open CoordinatesDeletions, "<$ARGV[0]" or die "cannot open $ARGV[0]\n";
$Insertions = open CoordinatesInsertions, "<$ARGV[1]" or die "cannot open $ARGV[1]\n";
$Inversions = open CoordinatesInversions, "<$ARGV[2]" or die "cannot open $ARGV[2]\n";

close CoordinatesDeletions;
close CoordinatesInsertions;
close CoordinatesInversions;

my $label = $ARGV[9];

$STATUS = (open CONFIG, "<$Bin/Config.txt");
if(!defined $STATUS){
	open CONFIG, "<$Bin/../Config.txt"
}
my $line;
my @info;
my $POS = 2;
my $out_dir = $ARGV[10];
my $genome_seq_dir;
my $RM_dir;
my $BLAST_dir;
my $RepeatMasker_dir;
while($line = <CONFIG>){
	chomp($line);
        if($line =~ /Genome_Seq_Directory/){
                @info = split/=|"/,$line;
                $genome_seq_dir = $info[$POS];
        }
        if($line =~ /RepeatMasked_Genome_Directory/){
                @info = split/=|"/,$line;
                $RM_dir = $info[$POS];
        }
        if($line =~ /BLAST_Directory/){
                @info = split/=|"/,$line;
                $BLAST_dir = $info[$POS];
        }
        if($line =~ /RepeatMasker_Program_Directory/){
                @info = split/=|"/,$line;
                $RepeatMasker_dir = $info[$POS];
        }
}
close CONFIG;

#create output directories

system "mkdir -p $out_dir";
system "mkdir -p $out_dir/RepeatMasker_Output";
system "mkdir -p $out_dir/RepeatMasker_Output/Deletions";
system "mkdir -p $out_dir/RepeatMasker_Output/Insertions";
system "mkdir -p $out_dir/RepeatMasker_Output/Inversions";
system "mkdir -p $out_dir/Output";
system "mkdir -p $out_dir/Output/Deletions";
system "mkdir -p $out_dir/Output/Insertions";
system "mkdir -p $out_dir/Output/Inversions";
system "mkdir -p $out_dir/Temp";

if($Deletions == 1){
	print "perl $Bin/Deletion_Repeatmasker_With_Software-3.pl $ARGV[0] $label $CUTOFF $out_dir $genome_seq_dir $RM_dir $SIZE ...\n";
	system "perl $Bin/Deletion_Repeatmasker_With_Software-3.pl $ARGV[0] $label $CUTOFF $out_dir $genome_seq_dir $RM_dir $SIZE";
	print "perl $Bin/SatelliteDeletion-2.pl $out_dir/Output/Deletions/$label.InputVNTRDeletion.txt $label $out_dir ...\n";
	system "perl $Bin/SatelliteDeletion-2.pl $out_dir/Output/Deletions/$label.InputVNTRDeletion.txt $label $out_dir";
}

if($Insertions == 1){
	print "perl $Bin/Insertion_Repeatmasker_With_Software.pl $ARGV[1] $label $CUTOFF $out_dir $genome_seq_dir $RepeatMasker_dir $SIZE ...\n";
	system "perl $Bin/Insertion_Repeatmasker_With_Software.pl $ARGV[1] $label $CUTOFF $out_dir $genome_seq_dir $RepeatMasker_dir $SIZE";
	print "perl $Bin/SatelliteInsertion-2.pl $out_dir/Output/Insertions/$label.InputVNTRInsertion.txt $label $out_dir ...\n";
	system "perl $Bin/SatelliteInsertion-2.pl $out_dir/Output/Insertions/$label.InputVNTRInsertion.txt $label $out_dir";
}

print "perl $Bin/combiningSatellite-2.pl $VNTRcutoff $label $out_dir ...\n";
system "perl $Bin/combiningSatellite-2.pl $VNTRcutoff $label $out_dir";

if($Deletions == 1){
 	print "perl $Bin/Flanking_Sequence_200_Deletions-3.pl $out_dir/$label.OtherVNTRDeletions-50.txt $FlankNAHR $label $out_dir ...\n";
 	system "perl $Bin/Flanking_Sequence_200_Deletions-3.pl $out_dir/$label.OtherVNTRDeletions-50.txt $FlankNAHR $label $out_dir";
	print "perl $Bin/bl2seqNew-3_Deletions.pl $out_dir/$label.FlankingSequenceDeletions.txt $label $FlankNAHR $NAHRhomolen $NAHRpct $out_dir $BLAST_dir ...\n";
	system "perl $Bin/bl2seqNew-3_Deletions.pl $out_dir/$label.FlankingSequenceDeletions.txt $label $FlankNAHR $NAHRhomolen $NAHRpct $out_dir $BLAST_dir";
	print "perl $Bin/Deletion_Repeatmasker_Step2.pl $out_dir/$label.Other_Non_NAHR_Deletions.txt $label $RTwin $RTgap $CUTOFF $out_dir ...\n";
	system "perl $Bin/Deletion_Repeatmasker_Step2.pl $out_dir/$label.Other_Non_NAHR_Deletions.txt $label $RTwin $RTgap $CUTOFF $out_dir";
	print "perl $Bin/MultiRTfragmentDeletion.pl $label $out_dir ...\n";
	system "perl $Bin/MultiRTfragmentDeletion.pl $label $out_dir";
	print "perl $Bin/PolyA_MicroHomologyFinder_NHEJ_Deletions.pl $out_dir/Output/Deletions/$label.NHEJ_Deletions.txt $label $Window $out_dir ...\n";
	system "perl $Bin/PolyA_MicroHomologyFinder_NHEJ_Deletions.pl $out_dir/Output/Deletions/$label.NHEJ_Deletions.txt $label $Window $out_dir";
	print "perl $Bin/PolyA_MicroHomologyFinderFragmentMultiRT_Deletion.pl $out_dir/Output/Deletions/$label.FragmentMultiple_Retrotransposon_Deletions.txt $label $Window $out_dir ...\n";
	system "perl $Bin/PolyA_MicroHomologyFinderFragmentMultiRT_Deletion.pl $out_dir/Output/Deletions/$label.FragmentMultiple_Retrotransposon_Deletions.txt $label $Window $out_dir";
	print "perl $Bin/PolyA_MicroHomologyFinder_SingleFragRT_Deletion.pl $out_dir/Output/Deletions/$label.Single_Retrotransposon_Deletions.txt $label $out_dir/Output/Deletions/$label.WholeMultiple_Retrotransposon_Deletions_Format3.txt $Window $out_dir ...\n";
	system "perl $Bin/PolyA_MicroHomologyFinder_SingleFragRT_Deletion.pl $out_dir/Output/Deletions/$label.Single_Retrotransposon_Deletions.txt $label $out_dir/Output/Deletions/$label.WholeMultiple_Retrotransposon_Deletions_Format3.txt $Window $out_dir";
}else{
	print "No Deletions Processed!\n";
}

if($Insertions == 1){
	print "perl $Bin/Flanking_Sequence_200_Insertions-3.pl $out_dir/$label.OtherVNTRInsertions-50.txt $FlankNAHR $label $out_dir ...\n";
	system "perl $Bin/Flanking_Sequence_200_Insertions-3.pl $out_dir/$label.OtherVNTRInsertions-50.txt $FlankNAHR $label $out_dir";
	print "perl $Bin/bl2seqNew-3_Insertions.pl $out_dir/$label.FlankingSequenceInsertions.txt $label $FlankNAHR $NAHRhomolen $NAHRpct $out_dir $BLAST_dir ...\n";
   system "perl $Bin/bl2seqNew-3_Insertions.pl $out_dir/$label.FlankingSequenceInsertions.txt $label $FlankNAHR $NAHRhomolen $NAHRpct $out_dir $BLAST_dir";
	print "perl $Bin/Insertion_Repeatmasker_Step2.pl $out_dir/$label.Other_Non_NAHR_Insertions.txt $label $RTwin $RTgap $CUTOFF $out_dir ...\n";
	system "perl $Bin/Insertion_Repeatmasker_Step2.pl $out_dir/$label.Other_Non_NAHR_Insertions.txt $label $RTwin $RTgap $CUTOFF $out_dir";
	print "perl $Bin/MultiRTfragmentInsertion.pl $label $out_dir ...\n";
	system "perl $Bin/MultiRTfragmentInsertion.pl $label $out_dir";
	print "perl $Bin/PolyA_MicroHomologyFinder_NHEJ_Insertions.pl $out_dir/Output/Insertions/$label.NHEJ_Insertions.txt $label $Window $out_dir ...\n";
	system "perl $Bin/PolyA_MicroHomologyFinder_NHEJ_Insertions.pl $out_dir/Output/Insertions/$label.NHEJ_Insertions.txt $label $Window $out_dir";
	print "perl $Bin/PolyA_MicroHomologyFinderFragmentMultiRT_Insertion.pl $out_dir/Output/Insertions/$label.FragmentMultiple_Retrotransposon_Insertions.txt $label $Window $out_dir ...\n";
	system "perl $Bin/PolyA_MicroHomologyFinderFragmentMultiRT_Insertion.pl $out_dir/Output/Insertions/$label.FragmentMultiple_Retrotransposon_Insertions.txt $label $Window $out_dir";
	print "perl $Bin/PolyA_MicroHomologyFinder_SingleFragRT_Insertion.pl $out_dir/Output/Insertions/$label.Single_Retrotransposon_Insertions.txt $label $out_dir/Output/Insertions/$label.WholeMultiple_Retrotransposon_Insertions_Format3.txt $Window $out_dir ...\n";
	system "perl $Bin/PolyA_MicroHomologyFinder_SingleFragRT_Insertion.pl $out_dir/Output/Insertions/$label.Single_Retrotransposon_Insertions.txt $label $out_dir/Output/Insertions/$label.WholeMultiple_Retrotransposon_Insertions_Format3.txt $Window $out_dir";
}else{
	print "No Insertions Processed!\n";
}

if($Inversions == 1){
	print "INVERSIONS:\nperl $Bin/Flanking_Sequence_200_Inversions-3.pl $ARGV[2] $label $FlankNAHR $out_dir $genome_seq_dir $RM_dir $SIZE ...\n";
	system "perl $Bin/Flanking_Sequence_200_Inversions-3.pl $ARGV[2] $label $FlankNAHR $out_dir $genome_seq_dir $RM_dir $SIZE";
	print "perl $Bin/bl2seqNew-3_Inversions.pl $out_dir/$label.FlankingSequenceInversions.txt $label $FlankNAHR $NAHRhomolen $NAHRpct $out_dir $BLAST_dir $genome_seq_dir...\n";
	system "perl $Bin/bl2seqNew-3_Inversions.pl $out_dir/$label.FlankingSequenceInversions.txt $label $FlankNAHR $NAHRhomolen $NAHRpct $out_dir $BLAST_dir";
	print "cp $out_dir/$label.Other_Non_NAHR_Inversions.txt $out_dir/Output/Inversions/$label.NHEJ_Inversions.txt ...\n";
	system "cp $out_dir/$label.Other_Non_NAHR_Inversions.txt $out_dir/Output/Inversions/$label.NHEJ_Inversions.txt";
	print "perl $Bin/MicroHomologyFinder_NHEJ_Inversion.pl $out_dir/Output/Inversions/$label.NHEJ_Inversions.txt $label $Window $out_dir ...\n";
	system "perl $Bin/MicroHomologyFinder_NHEJ_Inversion.pl $out_dir/Output/Inversions/$label.NHEJ_Inversions.txt $label $Window $out_dir";

}else{
	print "No Inversions Processed!\n";
}

print "perl $Bin/All_Output-2.pl $label $out_dir\n";
system "perl $Bin/All_Output-2.pl $label $out_dir";

system "perl $Bin/StandardizationMicrohomology.BreakSeq.pl $label $out_dir $genome_seq_dir";

my $text = "";

open OUT, ">$out_dir/CommandLine.$label.txt";
print OUT "delIn\tinsIn\tinvIn\tVNTRcutoff\tFlankingNAHR\tNAHRhomo\tNAHRpct\tRTwin\tRTgap\tlabel\n";
print OUT "$ARGV[0]";
foreach(1..9){
	print OUT "\t$ARGV[$_]";
}
print OUT "\n";
foreach(3..8){
	$text = $text."_".$ARGV[$_];
}
close OUT;
$text = $out_dir."/".$text;
#system "cp -r $out_dir/ $text";
#ARGV[0]: deletionInfile; ARGV[1]:InsertionInfile; ARGV[2]:InversionInfile; 
#ARGV[3]: VNTRcutoff(0-1; ARGV[4]: FlankingNAHR(100-1000,even number); ARGV[5]: NAHRhomologyLength(20-400); ARGV[6]: NAHRpct(0-100); ARGV[7]: RTwindow(10-1000,even number)
#ARGV[8]: RTgap:(0-1000)
#ARGV[9] = $label

print "DONE with \"Breakpoint_Classification_Pipeline_XJM_Jan122010.pl @ARGV\" ...\n";
print "Bye!\n";

