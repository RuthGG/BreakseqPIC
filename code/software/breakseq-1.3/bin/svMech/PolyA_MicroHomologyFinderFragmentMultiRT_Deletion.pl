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
my $mark3;
my $window;
my @elements;
my @temp;
my $position1 = 10;
my $position2 = 9;
my $position3 = 8;
my $myline;
my @elements1;
my $TSD =4; 
my $markTsd;
my $pA;
my $element;
my $element_name;
my $element_strand;
my $pA_text;
my $subOut;
my $classOut;
my @subClassOut;
my @class;
my $first;
my $repM1;
my $repM2;
my $out_dir = $ARGV[3];
my $Chr = 0;
my $Start = 1;
my $End = 2;
my $Left_F = 3;
my $Right_F = 4;
my $Name = 5;
my $ID = 6;
my $SV_Size;

open FLANKINGSEQ, "<$ARGV[0]"; #Modified06072009.FragmentMultiple_Retrotransposon_Deletions.txt
open MICROHOMOLOGY, ">$out_dir/Output/Deletions/$ARGV[1].MicroHomologyFragMultiple_Retrotransposon_Deletions.txt";
open MICROHOMOLOGY2, ">$out_dir/Output/Deletions/$ARGV[1].MicroHomologyFragMultiple_Retrotransposon_Deletions_Format2.txt";
open NOTSURE, ">$out_dir/Output/Deletions/$ARGV[1].UNSURE_6_Multiple_Retrotransposon_Deletions.txt";
open NOTSURE2, ">$out_dir/Output/Deletions/$ARGV[1].UNSURE_6_Multiple_Retrotransposon_Deletions_Format2.txt";

print MICROHOMOLOGY "ChrNo.\tStart\tEnd\tLeftFlankingSequence\tRightFlankingSequence\tName\tID\tClass\tNames\tnmhNo.\tnmhSeq\tpolyAonCorrStrand\tStrandElement\tTSD\tClassLabel\tRTnames\tstrandPolyAtail\n";
print NOTSURE "ChrNo.\tStart\tEnd\tLeftFlankingSequence\tRightFlankingSequence\tName\tID\tClass\tNames\tnmhNo.\tnmhSeq\tpolyAonCorrStrand\tStrandElement\tTSD\tClassLabel\tRTnames\tstrandPolyAtail\n";
print MICROHOMOLOGY2 "ChrNo.\tStart\tEnd\tLeftFlankingSequence\tRightFlankingSequence\tName\tID\tClass\tNames\tnmhNo.\tnmhSeq\tpolyAonCorrStrand\tStrandElement\tTSD\tClassLabel\tRTnames\tstrandPolyAtail\n";
print NOTSURE2 "ChrNo.\tStart\tEnd\tLeftFlankingSequence\tRightFlankingSequence\tName\tID\tClass\tNames\tnmhNo.\tnmhSeq\tpolyAonCorrStrand\tStrandElement\tTSD\tClassLabel\tRTnames\tstrandPolyAtail\n";
<FLANKINGSEQ>;
while($myline = <FLANKINGSEQ>){
	chomp($myline);
	@elements = ();
	@elements = split/\t\t\t/, $myline;
	@fields = ();
	@fields = split/\s+/, $elements[0];
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
	@elements1 = ();
	@elements1 = @elements;
	$element ="N.A";
	$element_name = "N.A";
	$element_strand = "N.A";
	$pA_text="N.A";
	$subOut="";
	$classOut="";
	@subClassOut=();
	@class=();
	
	if($mhNum==0){
		$mhSeq = "N.A";
	}else{
		$mhSeq = uc($mhSeq);
	}  
	shift(@elements1);
	pop(@elements1);
	if(($mark1 ==0)&&($mark2 ==1)){
		@elements1 = reverse (@elements1);
		foreach(@elements1){
			@temp = ();
			@temp = split/\s+/, $_;
			shift(@temp) if($temp[0]eq"");
		 	if(!(($temp[$position1] =~ /Low_complexity/)||($temp[$position2] =~ /Satellite/)||($temp[$position1] =~ /Simple_repeat/))){
				last;
			 }
		}
		$pA_text = "+";
		if($temp[$position1] =~ /LINE\/L1/){
			$element = "L1";
		}elsif($temp[$position2] =~ /SVA/){
			$element = "SVA";
		}elsif($temp[$position1] =~ /SINE\/Alu/){
			$element = "Alu";
		}else{
			$element = $temp[$position1];
		}
		$element_name = $temp[$position2];
		$element_strand = $temp[$position3];
		if($temp[$position3] eq "+"){
			$pA = 1;
		}else{
			$pA = 0;
		}	 
	}

	if(($mark2 ==0)&&($mark1 ==1)){
		foreach(@elements1){
			@temp = ();
			@temp = split/\s+/, $_;
				shift(@temp) if($temp[0]eq"");	  
		 	if(!(($temp[$position1] =~ /Low_complexity/)||($temp[$position2] =~ /Satellite/)||($temp[$position1] =~ /Simple_repeat/))){
				last;
			 }
		}
		if($temp[$position1] =~ /LINE\/L1/){
			$element = "L1";
		}elsif($temp[$position2] =~ /SVA/){
			$element = "SVA";
		}elsif($temp[$position1] =~ /SINE\/Alu/){
			$element = "Alu";
		}else{
			$element = $temp[$position1];
		}
		$element_name = $temp[$position2];
		$element_strand = $temp[$position3];
		if($temp[$position3] eq "C"){
			$pA = 1;
		}else{
			$pA = 0;
		}
		$pA_text = "-";
	}
	
	if(($mark1 ==0)&&($mark2 ==0)){
		$mark3 = 0;
		$pA_text = "N.A";
		foreach(@elements1){
			@temp = ();
			@temp = split/\s+/, $_;
			shift(@temp) if($temp[0]eq"");
		 	if(!(($temp[$position1] =~ /Low_complexity/)||($temp[$position2] =~ /Satellite/)||($temp[$position1] =~ /Simple_repeat/))){
				last;
			}
		}
		if((($temp[$position2] =~ /ERV/)||($temp[$position1] =~ /ERV/))&&($temp[$position3] eq "C")){
			$mark3 = 1;
			$element = $temp[$position1];
			$element_name = $temp[$position2];
			$element_strand = $temp[$position3];	
		}
		@elements1 = reverse (@elements1);   
	 	foreach(@elements1){
			@temp = ();
			@temp = split/\s+/, $_;
			shift(@temp) if($temp[0]eq"");
		 	if(!(($temp[$position1] =~ /Low_complexity/)||($temp[$position2] =~ /Satellite/)||($temp[$position1] =~ /Simple_repeat/))){
				last;
			 }
		}
		if((($temp[$position2] =~ /ERV/)||($temp[$position1] =~ /ERV/))&&($temp[$position3] eq "+")){
			$mark3 = 1;
			$element = $temp[$position1];
			$element_name = $temp[$position2];
			$element_strand = $temp[$position3];	
		}
		if($mark3 ==0){
			if($temp[$position1] =~ /LINE\/L1/){
				$element = "L1";
			}elsif($temp[$position2] =~ /SVA/){
				$element = "SVA";
			}elsif($temp[$position1] =~ /SINE\/Alu/){
				$element = "Alu";
			}else{
				$element = $temp[$position1];
			}
			$element_name = $temp[$position2];
			$element_strand = $temp[$position3];
		}
		$pA = 0;
	}

	if(($mark1 ==1)&&($mark2 ==1)){
		$pA_text = "+.-";
		foreach(@elements1){
			@temp = ();
			@temp = split/\s+/, $_;
			shift(@temp) if($temp[0]eq"");	  
			if(!(($temp[$position1] =~ /Low_complexity/)||($temp[$position2] =~ /Satellite/)||($temp[$position1] =~ /Simple_repeat/))){
				last;
			}
		}
		if($temp[$position1] =~ /LINE\/L1/){
			$element = "L1";
		}elsif($temp[$position2] =~ /SVA/){
			$element = "SVA";
		}elsif($temp[$position1] =~ /SINE\/Alu/){
			$element = "Alu";
		}else{
			$element = $temp[$position1];
		}
		$element_name = $temp[$position2];
		$element_strand = $temp[$position3];
		$pA = 1;
	}
	
	$subOut="";
	$classOut="";
	@subClassOut=();
	@class=();

	$first = 0;
	foreach(@elements1){
		@temp = ();
		@temp = split/\s+/, $_;
		shift(@temp) if($temp[0]eq"");
		if(!(($temp[$position1] =~ /Low_complexity/)||($temp[$position2] =~ /Satellite/)||($temp[$position1] =~ /Simple_repeat/))){
			if($first==0){
				$subOut = $temp[$position2];
				$classOut = $temp[$position1];
				push(@subClassOut,$temp[$position2]);
				push(@class,$temp[$position1]);
				$first=1;
			}else{
				$repM1=0;
				foreach(0..(scalar(@subClassOut)-1)){
					if($subClassOut[$_] eq $temp[$position2]){
						$repM1=1;
					}
				}
				if($repM1==0){
					$subOut = $subOut.".".$temp[$position2];
				}
				$repM2=0;       			
				foreach(0..(scalar(@class)-1)){
					if($class[$_] eq $temp[$position1]){
						$repM2=1;
					}
				}
				if($repM2==0){
					$classOut = $classOut.".".$temp[$position1];
				}
			}
		}
	}
	
	if(($element eq "L1")|($element eq "SVA")|($element eq "Alu")|($element =~ /LINE/)|($element =~ /SINE/)){
		if($pA==1){
			print MICROHOMOLOGY "$elements[0]\t$element\t$element_name\t$mhNum\t$mhSeq\t$pA\t$element_strand\t$markTsd\t$classOut\t$subOut\t$pA_text\t";
			$n = 1;
			while($n < (scalar(@elements)-1)){
				print MICROHOMOLOGY "\t\t\t$elements[$n]";
				$n++;
			}
			print MICROHOMOLOGY "\n";
			print MICROHOMOLOGY2 "$elements[0]\t$element\t$element_name\t$mhNum\t$mhSeq\t$pA\t$element_strand\t$markTsd\t$classOut\t$subOut\t$pA_text\n";
			$n = 1;
			while($n < (scalar(@elements)-1)){
				print MICROHOMOLOGY2 "$elements[$n]\n";
				$n++;
			} 
		}else{
			print NOTSURE "$elements[0]\t$element\t$element_name\t$mhNum\t$mhSeq\t$pA\t$element_strand\t$markTsd\t$classOut\t$subOut\t$pA_text\t";
			$n = 1;
			while($n < (scalar(@elements)-1)){
				print NOTSURE "\t\t\t$elements[$n]";
				$n++;
			}
			print NOTSURE "\n";
			print NOTSURE2 "$elements[0]\t$element\t$element_name\t$mhNum\t$mhSeq\t$pA\t$element_strand\t$markTsd\t$classOut\t$subOut\t$pA_text\n";
			$n = 1;
			while($n < (scalar(@elements)-1)){
				print NOTSURE2 "$elements[$n]\n";
				$n++;
			} 
		}
	}else{
		print MICROHOMOLOGY "$elements[0]\t$element\t$element_name\t$mhNum\t$mhSeq\t$pA\t$element_strand\t$markTsd\t$classOut\t$subOut\t$pA_text\t";
		$n = 1;
		while($n < (scalar(@elements)-1)){
			print MICROHOMOLOGY "\t\t\t$elements[$n]";
			$n++;
		}
		print MICROHOMOLOGY "\n";
		print MICROHOMOLOGY2 "$elements[0]\t$element\t$element_name\t$mhNum\t$mhSeq\t$pA\t$element_strand\t$markTsd\t$classOut\t$subOut\t$pA_text\n";
		$n = 1;
		while($n < (scalar(@elements)-1)){
			print MICROHOMOLOGY2 "$elements[$n]\n";
			$n++;
		} 	
	}
}
 
close NOTSURE;
close NOTSURE;
close FLANKINGSEQ;
close NOTSURE2;
close NOTSURE2;

