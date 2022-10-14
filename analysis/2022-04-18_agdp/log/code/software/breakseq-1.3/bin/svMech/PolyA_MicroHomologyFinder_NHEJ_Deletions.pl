#!/usr/bin/perl -w
use strict;
use POSIX qw(ceil floor);
#This program uses flanking sequences of NHEJs from the human Build 36 Genome and finds the poly A tail.
unless ($ARGV[3]) {  
	print STDERR "Exit!\n";
	exit;
}

my @fields;
my @leftseq;
my @rightseq;
my $mhNum;
my $mhSeq;
my $n;
my $i;
my $num2 = $ARGV[2];
my $size = 7;
my $mark1;
my $mark2;
my $window;
my $TSD = 7;
my $markTsd;
my $pA;
my $pA_text;
my $out_dir = $ARGV[3];
my $Chr = 0;
my $Start = 1;
my $End = 2;
my $Left_F = 3;
my $Right_F = 4;
my $Name = 5;
my $ID = 6;
my $SV_Size;

open FLANKINGSEQ, "<$ARGV[0]";
open MICROHOMOLOGY, ">$out_dir/Output/Deletions/$ARGV[1].MicroHomology_PolyA_NHEJ_Deletions.txt";
open UNSURE7, ">$out_dir/Output/Deletions/$ARGV[1].UNSURE_7_Potential_Processed_Pseudogene_Deletions.txt";

print MICROHOMOLOGY "ChrNo.\tStart\tEnd\tLeftFlankingSequence\tRightFlankingSequence\tName\tnmhNo.\tnmhSeq\tpolyAEixt\tTSD\tpolyA\tID\n";
print UNSURE7 "ChrNo.\tStart\tEnd\tLeftFlankingSequence\tRightFlankingSequence\tName\tnmhNo.\tnmhSeq\tpolyAEixt\tTSD\tpolyA\tID\n";

<FLANKINGSEQ>;
while(<FLANKINGSEQ>){
	chomp;
	@fields = ();
	@fields = split/\s+/, $_;
	$SV_Size = ($fields[$End]-$fields[$Start]);
	if($SV_Size>=length($fields[$Left_F])){
		@leftseq = split//,(substr($fields[$Left_F],(length($fields[$Left_F])/2-$num2/2),$num2));
		@rightseq = split//,(substr($fields[$Right_F],(length($fields[$Right_F])/2-$num2/2),$num2));
	}else{
		@leftseq = split//,(substr($fields[$Left_F],(length($fields[$Left_F])-floor($SV_Size/2)-$num2/2),$num2));
		@rightseq = split//,(substr($fields[$Right_F],(floor($SV_Size/2)-$num2/2),$num2));
	}
	$mhNum = 0;
	$mhSeq = "";
	$n = $num2/2-1;
	while($n>=0){
		if(lc($leftseq[$n]) eq lc($rightseq[$n])){
			$mhSeq = $leftseq[$n].$mhSeq;
			$mhNum++;
			$n--;
		}else{
			last;
		}
	}
	$n = $num2/2;
	while($n<=($num2-1)){
		if(lc($leftseq[$n]) eq lc($rightseq[$n])){
			$mhSeq = $mhSeq.$leftseq[$n];
			$mhNum++;
			$n++;
		}else{
			last;
		}
	}
	$mark1 = 0;
	$n = 0;
	while($n <= $num2-$size){
		$i = 0;
		$window = "";
		while($i<=$size-1){
			$window = $window.lc($leftseq[$n+$i]);
			$i++;
		}
		if($window eq "ttttttt"){
			$mark1 = 1;
		}
		$n++;
	}
	$mark2 = 0;
	$n = 0;
	while($n <= $num2-$size){
		$i = 0;
		$window = "";
		while($i<=$size-1){
			$window = $window.lc($rightseq[$n+$i]);
			$i++;
		}
		if($window eq "aaaaaaa"){
			$mark2 = 1;
		}
		$n++;
	}
	
	if($mhNum==0){
		$mhSeq = "N.A";
	}else{
		$mhSeq = uc($mhSeq);
	} 
	$markTsd = 0;
	if($mhNum >= $TSD){
		$markTsd = 1;
	}
		
	$pA=0;
	$pA_text = "N.A";
	if(($mark1 ==0)&&($mark2 ==1)){
		$pA_text = "+";
		$pA=1;
	}
	if(($mark2 ==0)&&($mark1 ==1)){
		$pA_text = "-";
		$pA=1;
	}
	if(($mark1 ==1)&&($mark2 ==1)){
		$pA_text = "+.-";
		$pA=1;
	}	
	
	if(($markTsd ==1)&&($pA==1)){
		print UNSURE7 "$fields[$Chr]\t$fields[$Start]\t$fields[$End]\t$fields[$Left_F]\t$fields[$Right_F]\t$fields[$Name]\t$mhNum\t$mhSeq\t$pA\t$markTsd\t$pA_text\t$fields[$ID]\n";
	}else{
		print MICROHOMOLOGY "$fields[$Chr]\t$fields[$Start]\t$fields[$End]\t$fields[$Left_F]\t$fields[$Right_F]\t$fields[$Name]\t$mhNum\t$mhSeq\t$pA\t$markTsd\t$pA_text\t$fields[$ID]\n";
	}
}
 

close MICROHOMOLOGY;
close FLANKINGSEQ;
close UNSURE7;


