#!/usr/bin/perl -w
use strict;

#This program deals with VNTRs including satellite and simple_repeats and low-complexity

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

my $VNTR_len;
my $SV_Size;
my $VNTRname;
my $VNTRclass;
my $VNTRID;
my $ratio;
my $Tname = 9;
my $Tclass = 10;
my $Tstart;
my $Tend;
my $strand = 8;
my $Ql = 7;
my $Qend = 6;
my $Qstart = 5;
my @Qleft;
my $ID = 14;
my $prevQstart;
my $prevQend;
my $count = 0;
my $out_dir = $ARGV[2];

open SV, "<$ARGV[0]" or die;  #original SV file
open SV_REPEAT,">$out_dir/Output/Deletions/$ARGV[1].VNTRDeletions-2.txt" or die;  #label
open Other, ">$out_dir/Output/Deletions/$ARGV[1].NonVNTR_Deletions-2.txt" or die;

print SV_REPEAT "chr\tStart\tend\tVNTRname\tVNTRclass\tVNTRid\tVNTRlen\tSV_Size\tRatio\tID\n";
print Other "chr\tStart\tend\tID\n";
<SV>;
while(<SV>){
	$num_SV++;
	chomp;
	$count++;
	@fields = ();
	@fields = split/\s+/, $_;
	(open REPEATMASKER,"<$out_dir/RepeatMasker_Output/Deletions/Qchr$fields[0]_$fields[1]-$fields[2].$fields[3].txt.out")or die "Cannot open repeatmasker Deletions/Qchr$fields[0]_$fields[1]-$fields[2].txt.out!\n";
	$templine= <REPEATMASKER>;
	@temp = ();
	@temp = split/\s+/, $templine;
 	if(($temp[0] eq "There")&&($temp[1] eq "were")&&($temp[2] eq "no")){
		print Other "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\n";
		close REPEATMASKER;
 # system "rm $out_dir/RepeatMasker_Output/Deletions/Qchr$fields[0]_$fields[1]-$fields[2].txt*";
		next;
	}
	<REPEATMASKER>;
	<REPEATMASKER>;
	$position = 0;
	$VNTR_len = 0;
	$mark = 0;
	$mark2 = 0;
	$VNTRname = "";
	$VNTRclass = "";
	$VNTRID = "";
	$SV_Size = 0;
	$prevQstart = 0;
	$prevQend = 0;
	while(<REPEATMASKER>){		
		chomp;
		@repeatmasker = ();
		@repeatmasker = split/\s+/, $_;
		if(($repeatmasker[0] eq "")){ 
			shift(@repeatmasker);
		} 		
		if($repeatmasker[$strand] eq "+"){
			$Tstart = 11;
			$Tend = 12;
		}else{
			$Tstart = 13;
			$Tend = 12;
		}
		if($mark2==0){
			@Qleft = split/\(|\)/,$repeatmasker[$Ql];
			$SV_Size = ($Qleft[1]+$repeatmasker[$Qend]-2*$flank);
			$mark2=1;
		}  
		$L = 0;
		$R = 0;
		if(($_ =~ /Satellite/)|($_ =~ /Simple_repeat/)|($_ =~ /Low_complexity/)){
			next if(($repeatmasker[$Qend]<=$flank)|($repeatmasker[$Qstart]>=($SV_Size+$flank)));
			next if($prevQend>=$repeatmasker[$Qend]);
			if($prevQend>=$repeatmasker[$Qstart]){		#situation when repeatmasked elements overlap
				$repeatmasker[$Qstart] = ($prevQend+1);
			}
			$mark = 1;
			if($repeatmasker[$Qstart]<$flank){
				$L = ($flank-$repeatmasker[$Qstart]+1);
			}
			if($repeatmasker[$Qend]>($SV_Size+$flank)){
				$R = ($repeatmasker[$Qend]-($SV_Size+$flank));
			}
			$VNTR_len = ($VNTR_len+$repeatmasker[$Qend]-$repeatmasker[$Qstart]+1-$L-$R);
			$VNTRname =$VNTRname.".".$repeatmasker[$Tname];
			$VNTRclass = $VNTRclass.".".$repeatmasker[$Tclass];
			$VNTRID = $VNTRID.".".$repeatmasker[$ID];
			$prevQstart= $repeatmasker[$Qstart];
			$prevQend = $repeatmasker[$Qend];
		 }
	}
	close REPEATMASKER;
	if($mark==0){
		print Other "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\n";
	
	}else{
		print SV_REPEAT "$fields[0]\t$fields[1]\t$fields[2]\t$VNTRname\t$VNTRclass\t$VNTRID\t$VNTR_len\t$SV_Size\t",($VNTR_len/$SV_Size),"\t$fields[3]\n";
	}
}

close SV;
close SV_REPEAT;
close Other;
