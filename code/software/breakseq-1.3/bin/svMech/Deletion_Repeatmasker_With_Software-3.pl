#!/usr/bin/perl -w
use strict;

#This program runs repeatmask on human DNA:sequence deleted in the novel genome assembly in reference to build 36 as well as its flanking 500bp regions outside each breakpoint.
#This program is modified on Jan 11, 2010 to make it more efficient by reading the repeatmask output for deletions from a pre-existing file downloaded from UCSC
unless ($ARGV[6]) {
	print STDERR "Exit! Deletion_Repeatmasker_With_Software-3.pl\n";
	exit;
}

my @fields;
my @chromosome;
my $flank;
my $n;
my $line =50;
my $allline;
my $moreline;
my $numline;
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
my $ref_ID = 3;
my $index;
my $count;
my $ref_rm;
my $ref_rm_count;
my $header_rm;
my $temp_line;
my @repeatmasker;
my $out_count;
my $i;
my $RM_Start = 6;
my $RM_End = 7;
my $RM_Left = 8;
my $flank_start;
my $flank_end;
my $flank_length;
my $outStr;
my $j;
my $out_start;
my $out_end;
my $text;
my $index_str;
my $out_dir = $ARGV[3]; #output directory
my $genome_seq_dir = $ARGV[4]; #directory for ncbi genome assembly sequences by chromosome
my $RM_dir = $ARGV[5]; #directory for repeatmasked genomes by chromosome downloaded from UCSC, repeatmasker version 3.2.7

open COORDINATES, "<$ARGV[0]" or die;
open EDGE, ">$out_dir/Output/Deletions/$ARGV[1].UNSURE_2_500bp_Proximity_To_Either_Telemere_Deletions.txt" or die;
open ELIMINATED, ">$out_dir/Output/Deletions/$ARGV[1].UNSURE_1_Smaller_Than_Cutoff_Deletions.txt" or die;
open LARGE, ">$out_dir/Output/Deletions/$ARGV[1].LargeDeletion.txt" or die;
open OUT, ">$out_dir/Output/Deletions/$ARGV[1].InputVNTRDeletion.txt" or die;

print LARGE "chr\tstart\tend\tID\n";
print ELIMINATED "chr\tstart\tend\tID\n";
print OUT "chr\tstart\tend\tID\n";
print EDGE "chr\tstart\tend\tID\n";

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
	
	open RepeatMasker, "<$RM_dir/chr$index_str.RM327.fa.out";
	$header_rm = "";
	foreach(1..3){
		$temp_line = <RepeatMasker>;
		$header_rm = $header_rm.$temp_line;
	}
	$ref_rm_count = 0;
	while($temp_line = <RepeatMasker>){
		chomp($temp_line);
		@repeatmasker = split/\s+/, $temp_line;
	 	if(!($repeatmasker[0] eq "")){ 
			unshift(@repeatmasker, "");
		} 
		$ref_rm_count++;
		$ref_rm->[$ref_rm_count] = [@repeatmasker];
	}
	close RepeatMasker;
	
	foreach $count (1..$ref_count[$index]){
		if(($ref_in[$index]->[$count][$ref_end]-$ref_in[$index]->[$count][$ref_start])<$CUTOFF2){
			print ELIMINATED "$index_str\t$ref_in[$index]->[$count][$ref_start]\t$ref_in[$index]->[$count][$ref_end]\t$ref_in[$index]->[$count][$ref_ID]\n";
			next;
		}

		$flank_start = ($ref_in[$index]->[$count][$ref_start]-$FLANK_LEN);
		$flank_end = ($ref_in[$index]->[$count][$ref_end]+$FLANK_LEN);
		$flank_length = ($flank_end-$flank_start);
		
		if(($flank_start<0) || ($flank_end>length($allline))){
			print EDGE "$index_str\t$ref_in[$index]->[$count][$ref_start]\t$ref_in[$index]->[$count][$ref_end]\t$ref_in[$index]->[$count][$ref_ID]\n";
			next;	  
		}
		$flank = substr($allline,$flank_start,$flank_length);
		$text = "chr".$index_str."_".$ref_in[$index]->[$count][$ref_start]."-".$ref_in[$index]->[$count][$ref_end].".".$ref_in[$index]->[$count][$ref_ID];
		open FASTA, ">$out_dir/RepeatMasker_Output/Deletions/Q$text.txt";
		print FASTA ">$text\t$index_str\t$ref_in[$index]->[$count][$ref_start]\t$ref_in[$index]->[$count][$ref_end]\t$ref_in[$index]->[$count][$ref_ID]\n$flank\n";
		close FASTA;
		
		if(($ref_in[$index]->[$count][$ref_end]-$ref_in[$index]->[$count][$ref_start])>$CUTOFF){
			print LARGE "$index_str\t$ref_in[$index]->[$count][$ref_start]\t$ref_in[$index]->[$count][$ref_end]\t$ref_in[$index]->[$count][$ref_ID]\n";
			next;
		}
		
		print OUT "$index_str\t$ref_in[$index]->[$count][$ref_start]\t$ref_in[$index]->[$count][$ref_end]\t$ref_in[$index]->[$count][$ref_ID]\n";
		
		open OUT_RM, ">$out_dir/RepeatMasker_Output/Deletions/Q$text.txt.out" or die;
		$out_count = 0;
		# To extract information from the UCSC RM file and format it into the .out file:
		# For the start and end positions (col 6 and 7), chunk the elements at the flanking sequence start and end
		# For the start, end and left positions (last few cols) in the element, ignore it, since the only time using these information are determining whether an element is full length for STEI, and the window is 100bp. Might need to take care of this if the window changes
				
		$outStr = "";
		foreach $i (1..$ref_rm_count){
			if($ref_rm->[$i][$RM_End]<=$flank_start){
				next;
			}elsif($ref_rm->[$i][$RM_Start]>=$flank_end){
				last;
			}else{
				$out_count++;
				$out_start = ($ref_rm->[$i][$RM_Start]-$flank_start);
				$out_end = ($ref_rm->[$i][$RM_End]-$flank_start);
				if($ref_rm->[$i][$RM_Start]<=$flank_start){
					$out_start = "1";
				}
				if($ref_rm->[$i][$RM_End]>$flank_end){
					$out_end = $flank_length;
				}								
				foreach $j (1..($RM_Start-1)){
					$outStr = $outStr.$ref_rm->[$i][$j]."\t";
				}
				$outStr = $outStr.$out_start."\t".$out_end."\t(".($flank_length-$out_end).")\t";  #start, end and left
				foreach $j (($RM_Left+1)..(scalar(@repeatmasker)-2)){
					$outStr = $outStr.$ref_rm->[$i][$j]."\t";
				}
				$outStr = $outStr.$ref_rm->[$i][scalar(@repeatmasker)-1]."\n";
			}	  		  
		}
		if($out_count == 0){
			print OUT_RM "There were no repetitive sequences detected in ./Seqs/Q$text.txt";
		}else{
			print OUT_RM "$header_rm.$outStr"; 	  
		}
		close OUT_RM;	
	}
}
close OUT;


















