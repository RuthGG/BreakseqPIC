#!/usr/bin/perl -w
use strict;
use List::Util qw[min max];
#This program is step 2 of program"Deletion_Repeatmasker.pl"
unless ($ARGV[5]) {  
	print STDERR "Exit!\n";
	exit;
}

my @fields;
my @repeatmasker;
my $position;
my $leftposition;
my $rightposition;
my $reference;
my $left;
my $right;
my $L;
my $R;
my $numIn = $ARGV[2]/2;
my $numOut = $ARGV[2]/2;
my $num=$ARGV[3];  #gap
my $i;
my $j;
my $mark;
my $mark2;
my $flank=500;
my $getnumline;
my @temp;
my $num_Insertion;
my $num_NHEJ=0;
my $num_NOT_SURE = 0;
my $num_SV_REPEAT = 0;
my $num_SV = 0;
my $num_Total;
my $templine;
my $CUTOFF = $ARGV[4];
my $first; 
my $last;
my $out_dir = $ARGV[5];

open SV, "<$ARGV[0]";
open SV_REPEAT,">$out_dir/Output/Deletions/$ARGV[1].Single_Retrotransposon_Deletions.txt";
open NHEJ, ">$out_dir/Output/Deletions/$ARGV[1].NHEJ_Deletions.txt";
open NOT_SURE_REPEATMASK, ">$out_dir/Output/Deletions/$ARGV[1].Multiple_Retrotransposon_Deletions.txt";
open NOT_SURE_REPEATMASK2, ">$out_dir/Output/Deletions/$ARGV[1].Multiple_Retrotransposon_Deletions_Format2.txt";

print NHEJ "ChrNo.\tStart\tEnd\tLeftFlankingSequence\tRightFlankingSequence\tName\tID\n";
print SV_REPEAT "ChrNo.\tStart\tEnd\tLeftFlankingSequence\tRightFlankingSequence\tName\tID\t\tRepeatMaskerElement\n";
print NOT_SURE_REPEATMASK "ChrNo.\tStart\tEnd\tLeftFlankingSequence\tRightFlankingSequence\tName\tID\n";
print NOT_SURE_REPEATMASK2 "ChrNo.\tStart\tEnd\tLeftFlankingSequence\tRightFlankingSequence\tName\tID\n";

<SV>;
while(<SV>){
	$num_SV++;
	chomp;
	@fields = ();
	@fields = split/\s+/, $_;
 	if(($fields[2]-$fields[1])>$CUTOFF){
		print NHEJ "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t$fields[4]\tNA_Size_Too_Large\t$fields[5]\n";
		next;  
	}
	(open REPEATMASKER,"<$out_dir/RepeatMasker_Output/Deletions/Qchr$fields[0]_$fields[1]-$fields[2].$fields[5].txt.out")or die "Cannot open repeatmasker Deletions/Qchr$fields[0]_$fields[1]-$fields[2].$fields[5].txt.out!\n";
	$templine= <REPEATMASKER>;
	@temp = ();
	@temp = split/\s+/, $templine;
 	if(($temp[0] eq "There")&&($temp[1] eq "were")&&($temp[2] eq "no")){
		$num_NHEJ++;
		print NHEJ "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t$fields[4]\tQchr$fields[0]_$fields[1]-$fields[2].$fields[5].txt.out\t$fields[5]\n";
		close REPEATMASKER;
		next;
	}
	<REPEATMASKER>;
	<REPEATMASKER>;
	$position = 0;
	$L = 0;
	$R = 0;
	while(<REPEATMASKER>){
		$mark = 0;
		chomp;
		@repeatmasker = ();
		@repeatmasker = split/\s+/, $_;
		if(!($repeatmasker[0] eq "")){ 
			unshift(@repeatmasker, "");
		} 
		if((($repeatmasker[6]>($flank-$numOut)) && ($repeatmasker[7]<($fields[2]-$fields[1]+$flank+$numOut)))&&(!(($repeatmasker[7]<=$flank)||($repeatmasker[6]>($fields[2]-$fields[1]+$flank))))){ 
			$reference->[$position] = [@repeatmasker];
			$position++;
			$left = 0;
			$right = 0;
			if(($repeatmasker[6]>($flank-$numOut))&&($repeatmasker[6]<=($flank+$numIn))){
				$left =1;
				$L = 1;
			}
			if(($repeatmasker[7]>($fields[2]-$fields[1]+$flank-$numIn))&&($repeatmasker[7]<=($fields[2]-$fields[1]+$flank+$numOut))){
				$right =1;
				$R = 1;
			}
			if (($left ==1)&&($right ==1)){
				$num_SV_REPEAT++;
				print SV_REPEAT "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t$fields[4]\tQchr$fields[0]_$fields[1]-$fields[2].$fields[5].txt.out\t$fields[5]\t\t";
				print SV_REPEAT "@repeatmasker\n";
				$mark = 1;
				last;
			}
		}
	}
	close REPEATMASKER;
 		
	if($mark == 0){
		if (($L == 0)||($R == 0)){
			$num_NHEJ++;
			print NHEJ "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t$fields[4]\tQchr$fields[0]_$fields[1]-$fields[2].$fields[5].txt.out\t$fields[5]\n";
			$mark =1;
		}
		if($mark == 0){
			$position--;	
			$first = 0;
			foreach(0..$position){
				if($reference->[$_][6]>=($flank+max($numIn,$num))){
					$first = ($_-1);
					last;
				}
			}
			$last = $position;
			for($i=$position;$i>=0;$i--){
				if($reference->[$i][7]<=($fields[2]-$fields[1]+$flank-max($numIn,$num))){
					$last = ($i+1);
					last;
				}		
			}
			$i = $first;
			while($i<=$last){
				$j = $first;
				$L = 0;
				$R = 0;
				while($j<=$last){
					if($i != $j){
	#					print "i = $i;\tj = $j;\tstart = $reference->[$i][6];\tend = $reference->[$i][7];\tjstart = $reference->[$j][6];\tjend = $reference->[$j][7]\n";
	#					print "left = ",($reference->[$i][6]-$num),";\tright = $reference->[$j][7];\n";
	#					if(!((($reference->[$i][6]-$num)>$reference->[$j][7])||(($reference->[$i][6]+$num)<$reference->[$j][6]))){
						if(!((($reference->[$i][6]-$num)>$reference->[$j][7])||(($reference->[$i][6])<$reference->[$j][6]))){
							$L = 1;
						}
	#					if(!((($reference->[$i][7]+$num)<$reference->[$j][6])||(($reference->[$i][7]-$num)>$reference->[$j][7]))){
						if(!((($reference->[$i][7]+$num)<$reference->[$j][6])||(($reference->[$i][7])>$reference->[$j][7]))){
							$R = 1;
						}
					}
					$j++;
				}
	#			print "L = $L;\tR = $R;\n";
				if(($reference->[$i][6]<($flank+$numIn))||($reference->[$i][6]<($flank+$num))){
					$L = 1;
				}
				if(($reference->[$i][7]>($fields[2]-$fields[1]+$flank-$numIn))||($reference->[$i][7]>($fields[2]-$fields[1]+$flank-$num))){
					$R = 1;
				}
				if (($L == 0)||($R == 0)){
					$num_NHEJ++;
					print NHEJ "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t$fields[4]\tQchr$fields[0]_$fields[1]-$fields[2].$fields[5].txt.out\t$fields[5]\n";
					$mark = 1;
					last;
				}
				$i++;
			}
	 
	
			if($mark ==0){
				$num_NOT_SURE++;
				print NOT_SURE_REPEATMASK "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t$fields[4]\tQchr$fields[0]_$fields[1]-$fields[2].$fields[5].txt.out\n";
				print NOT_SURE_REPEATMASK2 "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t$fields[4]\tQchr$fields[0]_$fields[1]-$fields[2].$fields[5].txt.out\t$fields[5]\t\t\t\t";
				$i = $first;
				while($i<=$last){
					$j = 0;
					while($j <=15){
						print NOT_SURE_REPEATMASK2 "$reference->[$i][$j]\t";
						print NOT_SURE_REPEATMASK "$reference->[$i][$j]\t";
					$j++;
					}
					print NOT_SURE_REPEATMASK2 "\t\t\t";
					print NOT_SURE_REPEATMASK "\n";
					$i++;
				}
				print NOT_SURE_REPEATMASK2 "\n";
			
			}
		}
	}
}

close SV;
close SV_REPEAT;
close NHEJ;
close NOT_SURE_REPEATMASK2;
close NOT_SURE_REPEATMASK;

