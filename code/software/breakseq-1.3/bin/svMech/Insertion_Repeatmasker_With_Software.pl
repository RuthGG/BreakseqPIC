#!/usr/bin/perl -w
use strict;

#This program runs repeatmask on human DNA:sequence deleted in the novel genome assembly in reference to build 36 as well as its flanking 500bp regions outside each breakpoint.
#This program is modified on Jan 11, 2010 to make it more efficient by reading the repeatmask output for Insertions from a pre-existing file downloaded from UCSC
unless ($ARGV[6]) {
	print STDERR "Exit! Insertion_Repeatmasker_With_Software.pl\n";
	exit;
}

my @fields;
my @chromosome;
my $flank;
my $n;
my $allline;
my $moreline;
my $num;
my $temp;
my $CUTOFF = $ARGV[2];  #20MB
my $CUTOFF2 = $ARGV[6];
my $FLANK_LEN = 500;

my @ref_in; #chrX: 23; chrY: 24;
my @ref_count; #chrX: 23; chrY: 24;
my $ref_line;
my @ref_info;
my $ref_chr= 0;
my $ref_start = 1;
my $ref_end = 2;
my $ref_ins = 3;
my $ref_ID = 4;
my $index;
my $count;
my $ref_rm;
my $ref_rm_count;
my $header_rm;
my $temp_line;
my @repeatmasker;
my $out_count;
my $i;
my $flank_start;
my $flank_end;
my $flank_length;
my $j;
my $out_start;
my $out_end;
my $text;
my $index_str;
my $out_dir = $ARGV[3]; #output directory
my $genome_seq_dir = $ARGV[4]; #directory for ncbi genome assembly sequences by chromosome
my $RepeatMasker_dir = $ARGV[5]; #directory for repeatmasker program installed

open COORDINATES, "<$ARGV[0]" or die;
open EDGE, ">$out_dir/Output/Insertions/$ARGV[1].UNSURE_2_500bp_Proximity_To_Either_Telemere_Insertions.txt" or die;
open ELIMINATED, ">$out_dir/Output/Insertions/$ARGV[1].UNSURE_1_Smaller_Than_Cutoff_Insertions.txt" or die;
open LARGE, ">$out_dir/Output/Insertions/$ARGV[1].LargeInsertion.txt" or die;
open OUT, ">$out_dir/Output/Insertions/$ARGV[1].InputVNTRInsertion.txt" or die;

print LARGE "chr\tstart\tend\tInsertedSeq\tID\n";
print ELIMINATED "chr\tstart\tend\tInsertedSeq\tID\n";
print OUT "chr\tstart\tend\tInsertedSeq\tID\n";
print EDGE "chr\tstart\tend\tInsertedSeq\tID\n";

foreach(1..24){
	$ref_count[$_] = 0;
}

<COORDINATES>;
while($ref_line = <COORDINATES>){
	chomp($ref_line);
	@ref_info = split/\s+/,$ref_line;
	foreach(1..22,"X","Y"){
		if($ref_info[$ref_chr] eq $_){
			if($_ eq "X"){
				$index = 23;
			}elsif($_ eq "Y"){
				$index = 24;
			}else{
				$index = $_;
			}
			$ref_count[$index]++;			
			$ref_in[$index]->[$ref_count[$index]] = [@ref_info];
		}
	}
}
close COORDINATES;

foreach $index (1..24){
	if($ref_count[$index] ==0){
		next;
	}
	$index_str = $index;
	if($index==23){
		$index_str = "X";
	}elsif($index==24){
		$index_str = "Y";
	}
	open GENOME, "<$genome_seq_dir/chr$index_str.fa";
	<GENOME>;
	$allline = "";
	while($moreline = <GENOME>){
	chomp($moreline);
		$allline = $allline.$moreline;
	}
	close GENOME;
	 
	foreach $count (1..$ref_count[$index]){
		$flank_start = ($ref_in[$index]->[$count][$ref_start]-$FLANK_LEN);
		$flank_end = ($ref_in[$index]->[$count][$ref_end]+$FLANK_LEN);
		$flank_length = length($ref_in[$index]->[$count][$ref_ins]);
		$text = "chr".$index_str."_".$ref_in[$index]->[$count][$ref_start]."-".$ref_in[$index]->[$count][$ref_end];		
		if($flank_length<$CUTOFF2){
			print ELIMINATED "$index_str\t$ref_in[$index]->[$count][$ref_start]\t$ref_in[$index]->[$count][$ref_end]\t$ref_in[$index]->[$count][$ref_ins]\t$ref_in[$index]->[$count][$ref_ID]\n";
			next;
		}	
		if(($flank_start<0) || ($flank_end>length($allline))){
			print EDGE "$index_str\t$ref_in[$index]->[$count][$ref_start]\t$ref_in[$index]->[$count][$ref_end]\t$ref_in[$index]->[$count][$ref_ins]\t$ref_in[$index]->[$count][$ref_ID]\n";
			next;	  
		}
		$flank = substr($allline,$flank_start,$FLANK_LEN).$ref_in[$index]->[$count][$ref_ins].substr($allline,$ref_in[$index]->[$count][$ref_end],$FLANK_LEN);
		$text = "chr".$index_str."_".$ref_in[$index]->[$count][$ref_start]."-".$ref_in[$index]->[$count][$ref_end].".".$ref_in[$index]->[$count][$ref_ID];
		open FASTA, ">$out_dir/RepeatMasker_Output/Insertions/Q$text.txt";
		print FASTA ">$text\t$index_str\t$ref_in[$index]->[$count][$ref_start]\t$ref_in[$index]->[$count][$ref_end]\t$ref_in[$index]->[$count][$ref_ID]\n$flank\n";
		close FASTA;
		
		if($flank_length>$CUTOFF){
			print LARGE "$index_str\t$ref_in[$index]->[$count][$ref_start]\t$ref_in[$index]->[$count][$ref_end]\t$ref_in[$index]->[$count][$ref_ins]\t$ref_in[$index]->[$count][$ref_ID]\n";
			next;
		}
		
		print OUT "$index_str\t$ref_in[$index]->[$count][$ref_start]\t$ref_in[$index]->[$count][$ref_end]\t$ref_in[$index]->[$count][$ref_ins]\t$ref_in[$index]->[$count][$ref_ID]\n";
		system "perl $RepeatMasker_dir/RepeatMasker $out_dir/RepeatMasker_Output/Insertions/Q$text.txt -s -species Human";
	}
}
close OUT;


















