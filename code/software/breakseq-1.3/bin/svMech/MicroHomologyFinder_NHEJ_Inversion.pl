#!/usr/bin/perl -w
use strict;
use POSIX qw(ceil floor);
#This program uses +100~-100bp flanking sequences of NHEJs from the human Build 36 Genome and finds the poly A tail.
#Given the file that holds the "ChrNo.  Start  End LeftFlankingSequence	RightFlankingSequence" coordinates and ignores the first line of title. 
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
my $out_dir = $ARGV[3];
my $Chr = 0;
my $Start = 1;
my $End = 2;
my $Left_F = 3;
my $Right_F = 4;
my $ID = 5;
my $SV_Size;

open FLANKINGSEQ, "<$ARGV[0]";
open MICROHOMOLOGY, ">$out_dir/Output/Inversions/$ARGV[1].MicroHomology_NHEJ_Inversions.txt";
print MICROHOMOLOGY "ChrNo.\tStart\tEnd\tLeftFlankingSequence\tRightFlankingSequence\tnmhNo.\tnmhSeq\tID\n";

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
	@rightseq = reverse(@rightseq);
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
	if($mhNum == 0){
		$mhSeq = "N.A";
	}else{
		$mhSeq = uc($mhSeq);
	} 
	print MICROHOMOLOGY "$fields[$Chr]\t$fields[$Start]\t$fields[$End]\t$fields[$Left_F]\t$fields[$Right_F]\t$mhNum\t$mhSeq\t$fields[$ID]\n";
}

close MICROHOMOLOGY;
close FLANKINGSEQ;


