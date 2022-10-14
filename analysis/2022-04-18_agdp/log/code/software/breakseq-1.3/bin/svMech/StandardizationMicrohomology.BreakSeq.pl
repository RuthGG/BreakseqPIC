#!/usr/bin/perl -w
use strict;

#This input file needs to be 0 based coordinates
#This program standardized the SVs to the left most position and reports the length of microhomology
#chrX and chY, "X". "Y" needs to be capitalized


unless ($ARGV[2]) {  
	print STDERR "Exit!\n";
	exit;
}

my $label = $ARGV[0];
my $out_dir = $ARGV[1];
my $genome_seq_dir = $ARGV[2];

my @ref_in;
my @ref_count;
my $ref_line;
my @ref_info;
my $index;
my $ref_header;

my $ref_chr = 0;
my $ref_start = 1;
my $ref_end = 2;
my $ref_SV_type = 3;
my $ref_inserted_seq = 23;
my $ref_mh_len = 13;
my $ref_mh_seq = 14;
my $ref_mech = 4;

my $CHR;

my $out_start;
my $out_end;
my $inserted_seq;
my $new_inserted_seq;
my $i;
my $j;
my @info;

my $mh_len;
my $mh_seq;
my $UPPER_Limit = 1000;
my $UPPER;

my $LF; #left flanking
my $RF; #right flanking
my $mark_l;
my $mark_r;
my $seq_l;
my $seq_r;
my $BREAK_end_pre = 24;
my $ADD = 5;

my $allline;
my $moreline;

my $MH_LEN = 30;
my $num_del = 0;
my $num_ins = 0;
my $num_inv = 0;
my $line;
my $P1 = 2;

system "mv $out_dir/Output/$label.All_Output.txt $out_dir/Output/$label.All_Output.v2.txt";
open MASTER, "<$out_dir/Output/$label.All_Output.v2.txt" or die;
foreach(1..24){
	$ref_count[$_] = 0;
}
$ref_header = <MASTER>;
chomp($ref_header);
while($ref_line = <MASTER>){
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
close MASTER;

open OUT, ">$out_dir/Output/$label.All_Output.txt" or die;
print OUT "$ref_header\n";

foreach $CHR (1..22,"X","Y"){
	if($CHR eq "X"){
		$index = 23;
	}elsif($CHR eq "Y"){
		$index = 24;
	}else{
		$index = $CHR;
	}	
	if($ref_count[$index] ==0){
		next;
	}
  	open GENOME, "<$genome_seq_dir/chr$CHR.fa";
  	<GENOME>;
 	 $allline = "";
 	 while($moreline = <GENOME>){
 		chomp($moreline);
  		$allline = $allline.$moreline;
 	 }
  	close GENOME;
  	
  	foreach $i (1..$ref_count[$index]){
		if($ref_in[$index]->[$i][$ref_SV_type] eq "Deletion"){
			if(($ref_in[$index]->[$i][$ref_end]-$ref_in[$index]->[$i][$ref_start])>=($UPPER_Limit/2)){
				$UPPER = $UPPER_Limit;
			}else{
				$UPPER = 2*($ref_in[$index]->[$i][$ref_end]-$ref_in[$index]->[$i][$ref_start]);
			}
 			$LF = uc(substr($allline,($ref_in[$index]->[$i][$ref_start]-$UPPER/2),$UPPER));
 			$RF = uc(substr($allline,($ref_in[$index]->[$i][$ref_end]-$UPPER/2),$UPPER));
 		 }elsif($ref_in[$index]->[$i][$ref_SV_type] eq "Insertion"){
   			$inserted_seq = $ref_in[$index]->[$i][$ref_inserted_seq];
   			if(length($inserted_seq)>=($UPPER_Limit/2)){
   				$UPPER = $UPPER_Limit;
   			}else{
   				$UPPER = 2*(length($inserted_seq));
   			}
 			$LF = uc(substr($allline,($ref_in[$index]->[$i][$ref_start]-$UPPER/2),$UPPER/2).substr($inserted_seq,0,$UPPER/2));
 			$RF = uc(substr($inserted_seq,(length($inserted_seq)-$UPPER/2),$UPPER/2).substr($allline,($ref_in[$index]->[$i][$ref_end]),$UPPER/2));
   		}elsif($ref_in[$index]->[$i][$ref_SV_type] eq "Inversion"){
			$mh_len = $ref_in[$index]->[$i][$ref_mh_len];
			$mh_seq = uc($ref_in[$index]->[$i][$ref_mh_seq]);
  		}else{
   			print "ERROR!\n$ref_in[$index]->[$i][$ref_start]\t$ref_in[$index]->[$i][$ref_end]";
   			next;
   		}
   		
		if(!($ref_in[$index]->[$i][$ref_SV_type] eq "Inversion")){
			$mark_l = $UPPER/2;
			$seq_l = substr($LF,($UPPER/2-$mark_l),$mark_l);
			for($j = 1; $j <=$UPPER/2; $j++){
				if(!(substr($LF,($UPPER/2-$j),$j) eq (substr($RF,($UPPER/2-$j),$j)))){
					$mark_l = ($j-1); #length of the microhomology to the left of the breakpoints. 
					$seq_l = substr($LF,($UPPER/2-$mark_l),$mark_l);
					last;
				}
			}
			$mark_r = $UPPER/2;
			$seq_r = substr($LF,$UPPER/2,$mark_r);
			for($j = 1; $j <=$UPPER/2; $j++){
				if(!(substr($LF,$UPPER/2,$j) eq (substr($RF,$UPPER/2,$j)))){
					$mark_r = ($j-1);        #length of the microhomology to the right of the breakpoints. 
					$seq_r = substr($LF,$UPPER/2,$mark_r);
					last;
				}
			}
			$mh_len = ($mark_l+$mark_r);
			if($mh_len == 0){
				$mh_seq = "N.A";
			}else{
				$mh_seq = $seq_l.$seq_r;
			}
		}
		if(($ref_in[$index]->[$i][$ref_mech] eq "NHR") & ($mh_len>=$MH_LEN)){
				if($ref_in[$index]->[$i][$ref_SV_type] eq "Deletion"){
					$num_del++;
				}elsif($ref_in[$index]->[$i][$ref_SV_type] eq "Insertion"){
					$num_ins++;
				}elsif($ref_in[$index]->[$i][$ref_SV_type] eq "Inversion"){
					$num_inv++;
				}
				$ref_in[$index]->[$i][$ref_mech] = "NAHR";
				foreach(0..($ref_mh_len-1)){
					print OUT "$ref_in[$index]->[$i][$_]\t";
				}
				print OUT "$mh_len\t",uc($mh_seq);
				foreach(($ref_mh_seq+1)..($BREAK_end_pre-1)){
					print OUT "\t$ref_in[$index]->[$i][$_]";
				}
				print OUT "\n";  			
		}else{
			foreach(0..($ref_mh_len-1)){
				print OUT "$ref_in[$index]->[$i][$_]\t";
			}
			print OUT "$mh_len\t",uc($mh_seq);
			foreach(($ref_mh_seq+1)..($BREAK_end_pre-1)){
				print OUT "\t$ref_in[$index]->[$i][$_]";
			}
			print OUT "\n";  			
		}
  	}
}
close OUT;

system "mv $out_dir/Output/$label.Statistics.txt $out_dir/Output/$label.Statistics.v2.txt";
open S, "<$out_dir/Output/$label.Statistics.v2.txt" or die;
open OUT, ">$out_dir/Output/$label.Statistics.txt" or die;
$line = <S>;
print OUT "$line";
while($line = <S>){
	chomp($line);
	@info = split/\s+/,$line;
	foreach $i (0..($P1-1)){
		print OUT "$info[$i]\t";
	}
	if($info[0] eq "Deletions"){
		print OUT "",($info[$P1]-$num_del),"\t",($info[$P1+1]+$num_del);
	}elsif($info[0] eq "Insertions"){
		print OUT "",($info[$P1]-$num_ins),"\t",($info[$P1+1]+$num_ins);
	}elsif($info[0] eq "Inversions"){
		print OUT "",($info[$P1]-$num_inv),"\t",($info[$P1+1]+$num_inv);
	}
	foreach $i (($P1+2)..(scalar(@info)-1)){
		print OUT "\t$info[$i]";
	}
	print OUT "\n";
}
close OUT;
close S;


system "rm $out_dir/Output/$label.Statistics.v2.txt";
system "rm $out_dir/Output/$label.All_Output.v2.txt";











