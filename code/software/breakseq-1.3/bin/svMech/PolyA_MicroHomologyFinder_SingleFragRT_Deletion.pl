#!/usr/bin/perl -w
use strict;
use POSIX qw(ceil floor);
#This program uses flanking sequences of NHEJs from the human Build 36 Genome and finds the poly A tail.
unless ($ARGV[4]) {  
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
my $num2 = $ARGV[3];
my $size = 7;
my $mark1;
my $mark2;
my $window;
my $position = 15;
my $position2 = 17;
my $position3 = 16;  #Tname
my $element;
my $markTsd;
my $TSD =4; 
my $ther_full = -1;
my $full;
my $Tstart;
my $Tend= 19;
my $Tleft;
my $coverage;
my $final;
my $OTHER_CUT = (280+490)*0.2;
my $CUT = 0.2;
my $LINE_CUT = 2000*$CUT;   #2000 approximately upper limite of 5'UTR length
my $Alu_CUT=280*$CUT;
my $temp_CUT;
my $Other_5pr;
my $Other_3pr;
my $pA;
my $tsd_log;
my $line;
my $Tstrand = 15;
my $Tname = 14;
my $Tclass = 13;
my $out_dir = $ARGV[4];
my $Chr = 0;
my $Start = 1;
my $End = 2;
my $Left_F = 3;
my $Right_F = 4;
my $Name = 5;
my $ID = 6;
my $SV_Size;

open FLANKINGSEQ, "<$ARGV[0]"; #Modified06072009.Single_Retrotransposon_Deletions.txt
open MICROHOMOLOGY, ">$out_dir/Output/Deletions/$ARGV[1].MicroHomologySingleFragRT_Deletions.txt";
open NOTSURE, ">$out_dir/Output/Deletions/$ARGV[1].UNSURE_5_Single_Retrotransposon_Deletions.txt";

print MICROHOMOLOGY "ChrNo.\tStart\tEnd\tLeftFlankingSequence\tRightFlankingSequence\tName\tID\tcoverage\tfull\tratio\ttherFull\tratio2\tfinal\tClass\tNames\tnmhNo.\tnmhSeq\tpolyAonCorrStrand\tStrandElement\tTSD\n";
print NOTSURE "ChrNo.\tStart\tEnd\tLeftFlankingSequence\tRightFlankingSequence\tName\tID\tcoverage\tfull\tratio\ttherFull\tratio2\tfinal\tClass\tNames\tnmhNo.\tnmhSeq\tpolyAonCorrStrand\tStrandElement\tTSD\n";
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
		$markTsd = 0;
		if($mhNum >= $TSD){
			$markTsd = 1;
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
	
	$element = "";
	$ther_full = -1;
	$coverage = 0;
	$full = 0;
	$final = "";
	$Other_5pr=0;
	$Other_3pr=0;
	$temp_CUT = 0;
	
	if($fields[$position] eq "+"){
		$Tstart = 18;
		$Tleft = 20;
	}else{
		$Tstart = 20;
		$Tleft = 18;
	}
	$fields[$Tleft] = (split/\(|\)/,$fields[$Tleft])[1];
	$full = ($fields[$Tend]+$fields[$Tleft]);
	$coverage = ($fields[$Tend]-$fields[$Tstart]+1);
	
	if($fields[$position2] =~ /LINE\/L1/){
		$element = "L1";
		$ther_full = 6000;
		$temp_CUT = $LINE_CUT;
	}
	if($fields[$position3] =~ /SVA/){
		$element = "SVA";
		$ther_full=3000;
		$temp_CUT = $OTHER_CUT;
	}
	if($fields[$position2] =~ /SINE\/Alu/){
		$element = "Alu";
		$ther_full = 280;
		$temp_CUT=$Alu_CUT;
	}
	
	if($fields[$Tstart]<$temp_CUT){
		$Other_5pr=1;
	}
	if($fields[$Tend]>($full-$temp_CUT)){
		$Other_3pr = 1;
	}
	
	if(($Other_5pr==1)&($Other_3pr==1)){
		$final = "Full";
	}elsif(($Other_5pr==0)&($Other_3pr==1)){
		$final = "5prTruncated";
	}elsif(($Other_5pr==1)&($Other_3pr==0)){
		$final = "3prTruncated";
	}elsif(($Other_5pr==0)&($Other_3pr==0)){
		$final = "bothTruncated";				
	}
	
	$pA = 0;

	if(($mark1 ==0)&&($mark2 ==1)){
		if(($fields[$position] eq "+")){
			$pA = 1;
		}
		else{
			$pA = 0;
		}
	}
	if(($mark2 ==0)&&($mark1 ==1)){
		if(($fields[$position] eq "C")){
			$pA = 1;
		}
		else{
			$pA = 0;
		}
	}
	if(($mark1 ==0)&&($mark2 ==0)){
		$pA = 0;
	}
	if(($mark1 ==1)&&($mark2 ==1)){
			$pA = 1;
	}

	if($mhNum==0){
		$mhSeq = "N.A";
	}else{
		$mhSeq = uc($mhSeq);
	}

	if(!($element eq "")){
		if(($markTsd == 1)&($pA ==1)){
			print MICROHOMOLOGY "$fields[$Chr]\t$fields[$Start]\t$fields[$End]\t$fields[$Left_F]\t$fields[$Right_F]\t$fields[$Name]\t$fields[$ID]\t$coverage\t$full\t",($coverage/$full),"\t$ther_full\t",($coverage/$ther_full),"\t$final\t$element\t$fields[$position3]\t$mhNum\t$mhSeq\t$pA\t$fields[$position]\t$markTsd\n";
		}else{
			print NOTSURE "$fields[$Chr]\t$fields[$Start]\t$fields[$End]\t$fields[$Left_F]\t$fields[$Right_F]\t$fields[$Name]\t$fields[$ID]\t$coverage\t$full\t",($coverage/$full),"\t$ther_full\t",($coverage/$ther_full),"\t$final\t$element\t$fields[$position3]\t$mhNum\t$mhSeq\t$pA\t$fields[$position]\t$markTsd\n";
		}  
	}else{
		$element = $fields[$position2];
		$final = "N.A";
		if($markTsd == 1){
			if(($element =~ /LINE/)|($element =~ /SINE/)){
				if($pA==1){
					print MICROHOMOLOGY "$fields[$Chr]\t$fields[$Start]\t$fields[$End]\t$fields[$Left_F]\t$fields[$Right_F]\t$fields[$Name]\t$fields[$ID]\t$coverage\t$full\t",($coverage/$full),"\t$ther_full\t",($coverage/$ther_full),"\t$final\t$element\t$fields[$position3]\t$mhNum\t$mhSeq\t$pA\t$fields[$position]\t$markTsd\n"; 				
				}else{
					print NOTSURE "$fields[$Chr]\t$fields[$Start]\t$fields[$End]\t$fields[$Left_F]\t$fields[$Right_F]\t$fields[$Name]\t$fields[$ID]\t$coverage\t$full\t",($coverage/$full),"\t$ther_full\t",($coverage/$ther_full),"\t$final\t$element\t$fields[$position3]\t$mhNum\t$mhSeq\t$pA\t$fields[$position]\t$markTsd\n";   	
					}
			}else{
				print MICROHOMOLOGY "$fields[$Chr]\t$fields[$Start]\t$fields[$End]\t$fields[$Left_F]\t$fields[$Right_F]\t$fields[$Name]\t$fields[$ID]\t$coverage\t$full\t",($coverage/$full),"\t$ther_full\t",($coverage/$ther_full),"\t$final\t$element\t$fields[$position3]\t$mhNum\t$mhSeq\t$pA\t$fields[$position]\t$markTsd\n";
			}
		}else{
			print NOTSURE "$fields[$Chr]\t$fields[$Start]\t$fields[$End]\t$fields[$Left_F]\t$fields[$Right_F]\t$fields[$Name]\t$fields[$ID]\t$coverage\t$full\t",($coverage/$full),"\t$ther_full\t",($coverage/$ther_full),"\t$final\t$element\t$fields[$position3]\t$mhNum\t$mhSeq\t$pA\t$fields[$position]\t$markTsd\n";   	
		}
	}
}
close FLANKINGSEQ;

open FLANKINGSEQ, "<$ARGV[2]"; #Modified06072009.WholeMultiple_Retrotransposon_Deletions_Format3.txt

<FLANKINGSEQ>;
while($line =<FLANKINGSEQ>){
	chomp($line);
	@fields = ();
	@fields = split/\s+/, $line;
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
		$markTsd = 0;
		if($mhNum >= $TSD){
			$markTsd = 1;
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

	$pA = 0;
	$tsd_log = -1;

	 if(($mark1 ==0)&&($mark2 ==1)){
		if(($fields[$Tstrand] eq "+")){
			$pA = 1;
		}
		else{
			$pA = 0;
		}
	}
	if(($mark2 ==0)&&($mark1 ==1)){
		if(($fields[$Tstrand] eq "C")){
			$pA = 1;
		}
		else{
			$pA = 0;
		}
	}
	if(($mark1 ==0)&&($mark2 ==0)){
		$pA = 0;
	}
	if(($mark1 ==1)&&($mark2 ==1)){
		$pA = 1;
	}

	if($mhNum==0){
		$mhSeq = "N.A";
	}else{
		$mhSeq = uc($mhSeq);
	}
	
	if((($fields[$Tclass] =~ /L1/)|($fields[$Tclass] =~ /Alu/)|($fields[$Tclass] =~ /SVA/))){
		if(($markTsd == 1)&($pA ==1)){
			foreach(0..$Tname){
				print MICROHOMOLOGY "$fields[$_]\t";
			}
			print MICROHOMOLOGY "$mhNum\t$mhSeq\t$pA\t$fields[$Tstrand]\t$markTsd\n";
		}else{
			foreach(0..$Tname){
				print NOTSURE "$fields[$_]\t";
			}
			print NOTSURE "$mhNum\t$mhSeq\t$pA\t$fields[$Tstrand]\t$markTsd\n";
		}  
	}else{
		if($markTsd == 1){
			if(($fields[$Tclass] =~ /LINE/)|($fields[$Tclass] =~ /SINE/)){
				if($pA==1){
					foreach(0..$Tname){
						print MICROHOMOLOGY "$fields[$_]\t";
					}
					print MICROHOMOLOGY "$mhNum\t$mhSeq\t$pA\t$fields[$Tstrand]\t$markTsd\n";
				}else{
					foreach(0..$Tname){
						print NOTSURE "$fields[$_]\t";
					}
					print NOTSURE "$mhNum\t$mhSeq\t$pA\t$fields[$Tstrand]\t$markTsd\n";
				}
			}else{
				foreach(0..$Tname){
					print MICROHOMOLOGY "$fields[$_]\t";		
				}
				print MICROHOMOLOGY "$mhNum\t$mhSeq\t$pA\t$fields[$Tstrand]\t$markTsd\n";
			}
		}else{
			foreach(0..$Tname){
				print NOTSURE "$fields[$_]\t";
			}
			print NOTSURE "$mhNum\t$mhSeq\t$pA\t$fields[$Tstrand]\t$markTsd\n";
		}
	}
}

close NOTSURE;
close MICROHOMOLOGY;
close FLANKINGSEQ;


