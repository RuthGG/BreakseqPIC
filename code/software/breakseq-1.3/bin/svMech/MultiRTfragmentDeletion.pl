#!/usr/bin/perl -w
use strict;

#This program tries to put together single transposable element insertion fragmented by repeatmasker from the multi transposable element insertion pool
# SVAs are fragmented into three pieces including (GGGAGA)n simple repeats or (TCTCCC)n
my $line;
my @info;
my $reference;
my $ref_count = 0;
#my $temp = 0;
my @fields;
my $Tname = 9;
my $Tclass = 10;
my $Tstart;
my $Tend = 12;
my $Tleft;
my $strand = 8;
my $Qleft = 7;
my $Qend = 6;
my $Qstart = 5;
my $mark = 0;
my $mark2 =0;
my $mark3 = 0;
my $mark4 = 0;
my $flag = 0;
my $element = "";
my $full;
my $coverage = 0;
my $cov_start= 0;
my $cov_end = 0;
my $flank = 500;
my $prev_Qend;
my $ther_full;
my $final;
my $element_strand;
my $Other_count;
my $Other_count_start;
my $Other_count_end;
my $otherM;
my $otherM2;
my $OTHER_CUT = (280+490)*0.2;
my $CUT = 0.2;
my $LINE_CUT = 2000*$CUT;   #2000 approximately upper limite of 5'UTR length
my $LINE_OVER = 0.5;  #largest over lap of LINE1 fragments
my $overlap;
my $distance;
my $AluFull = 280;
my $Other_5pr;
my $Other_3pr;
my $Other_between;
my $i;
my $j; 
my $classOut;
my @subClassOut = ();
my $subOut="";
my $mark_rep = 0;
my $out_dir = $ARGV[1];
	
open MULTI, "<$out_dir/Output/Deletions/$ARGV[0].Multiple_Retrotransposon_Deletions_Format2.txt";
open OUT_WH, ">$out_dir/Output/Deletions/$ARGV[0].WholeMultiple_Retrotransposon_Deletions.txt"; #frangments can be put together
open OUT_WH2, ">$out_dir/Output/Deletions/$ARGV[0].WholeMultiple_Retrotransposon_Deletions_Format2.txt";
open OUT_WH3, ">$out_dir/Output/Deletions/$ARGV[0].WholeMultiple_Retrotransposon_Deletions_Format3.txt";
open OUT_FR, ">$out_dir/Output/Deletions/$ARGV[0].FragmentMultiple_Retrotransposon_Deletions.txt"; #frangments cannot be put together
open OUT_FR2, ">$out_dir/Output/Deletions/$ARGV[0].FragmentMultiple_Retrotransposon_Deletions_Format2.txt";
open OUT_FR3, ">$out_dir/Output/Deletions/$ARGV[0].FragmentMultiple_Retrotransposon_Deletions_Format3.txt";
print OUT_WH "chr\tStart\tEnd\tleftF\trightF\tQname\tID\tcoverage\tfull\tratio\ttherFull\tratio2\tfinal\tClass\tNames\tstrandElement\n";
print OUT_WH2 "chr\tStart\tEnd\tleftF\trightF\tQname\tID\tcoverage\tfull\tratio\ttherFull\tratio2\tfinal\tClass\tNames\tstrandElement\n";
print OUT_WH3 "chr\tStart\tEnd\tleftF\trightF\tQname\tID\tcoverage\tfull\tratio\ttherFull\tratio2\tfinal\tClass\tNames\tstrandElement\n";
print OUT_FR "chr\tStart\tEnd\tleftF\trightF\tQname\tID\n";
print OUT_FR2 "chr\tStart\tEnd\tleftF\trightF\tQname\tID\n";
print OUT_FR3 "chr\tStart\tEnd\tleftF\trightF\tQname\tID\n";
<MULTI>;
while($line = <MULTI>){
	chomp($line);
	@info = ();
	@info = split/\t\t\t/,$line;
	$ref_count = (scalar(@info)-2);
	$mark = 0;
	$mark2 = 0; #0 means only one class of RT; 1 means more than one class
	$mark3 = 0;
	$mark4 = 0; #0 means that strands are same;
	$coverage = 0;
	$cov_start = 0;
	$cov_end = 0;
	$full = 0;
	$Other_count_start = -1;
	$Other_count_end = -1;
	$ther_full = -1;
	$element = "";
	$element_strand= "";
	$prev_Qend = -1;
	$classOut = "";
	@subClassOut = ();
	foreach(0..($ref_count-1)){
		@fields = (split/\s+/,$info[$_+1]);
		if($fields[0]eq""){
			shift(@fields);
		}
		if($fields[$strand] eq "+"){
			$Tstart = 11;
			$Tleft = 13;
		}else{
			$Tstart = 13;
			$Tleft = 11;
		}
		$fields[$Tleft] = (split/\(|\)/,$fields[$Tleft])[1];
		$fields[$Qleft] = (split/\(|\)/,$fields[$Qleft])[1];
		next if(($fields[$Qend]<=$flank)|(($fields[$Qend]+$fields[$Qleft]-$fields[$Qstart])<=$flank));
		$reference->[$_]=[@fields];
		if($mark == 0){
			next if(($fields[$Tclass] =~/Low_complexity/)|($fields[$Tclass] =~/Satellite/)|($fields[$Tclass] =~/Simple_repeat/));
			$mark = 1;
			$element = $fields[$Tclass];
			$full = ($fields[$Tend]+$fields[$Tleft]);
			$cov_start = $fields[$Qstart];
			$element_strand = $fields[$strand];
			$Other_count_start = $_;
			$classOut = $fields[$Tclass];
			push(@subClassOut,$fields[$Tname]);
			$subOut = $fields[$Tname];
		}else{
			if($mark2==0){
				if(!($fields[$Tclass] eq $element)){
					$mark2 = 1;
					if(!(($fields[$Tclass] =~/Low_complexity/)|($fields[$Tclass] =~/Satellite/)|($fields[$Tclass] =~/Simple_repeat/))){
						$mark3 = 1;
						$cov_end = $prev_Qend;
						$Other_count_end = ($_-1);
					}else{
						if($element =~/Other/){						
							$mark2 = 0;
						}else{
							$cov_end = $prev_Qend;   #Tend of the previous element
						}					
					}
				}else{
					$mark4 = 1 if(!($fields[$strand] eq $element_strand));
				}
			}else{
				$mark3 = 1 if(!(($fields[$Tclass] =~/Low_complexity/)|($fields[$Tclass] =~/Satellite/)|($fields[$Tclass] =~/Simple_repeat/)));
			}
		}
		$prev_Qend = $fields[$Qend];
	}	
	if($classOut eq "Other"){
		$classOut = "SVA";
	}elsif($classOut=~/LINE\/L1/){
		$classOut = "L1";
	}
	if(($mark==1)&($mark3==0)&(!($element eq ""))&($mark4==0)){
		if($cov_end<=0){
			if($fields[$Qleft]<$flank){
				$cov_end = ($fields[$Qend]+$fields[$Qleft]-$flank);
			}else{
				$cov_end = $fields[$Qend];
			}
		}
		if($cov_start<=$flank){
			$cov_start = ($flank+1);
		}	
		$coverage = ($cov_end-$cov_start+1);
		if($element =~ /L1/){		
			$ther_full = 6000;
		}elsif($element =~ /Other/){
			$ther_full = 3000;
		}elsif($element =~ /ERVK/){
			$ther_full = 10000;
		}elsif($element =~ /Alu/){
			$ther_full = 280;
		}
		$final= "N.A";		 	
		if($Other_count_end==-1){
			for($i = ($ref_count-1);$i>=0;$i--){
				if(!(($reference->[$i][$Tclass] =~/Low_complexity/)|($reference->[$i][$Tclass] =~/Satellite/)|($reference->[$i][$Tclass] =~/Simple_repeat/))){
					$Other_count_end = $i;
					last;
				}
			}
		}		
		foreach($Other_count_start..$Other_count_end){
			$i = $_;
			$mark_rep = 0;
			next if((($reference->[$i][$Tclass] =~/Low_complexity/)|($reference->[$i][$Tclass] =~/Satellite/)|($reference->[$i][$Tclass] =~/Simple_repeat/)));
			foreach(0..(scalar(@subClassOut)-1)){
				$j = $_;
				if($subClassOut[$j] eq $reference->[$i][$Tname]){
					$mark_rep =1;
				}
			}
			if($mark_rep==0){
				push(@subClassOut,$reference->[$i][$Tname]);
				$subOut = $subOut.".".$reference->[$i][$Tname];
			}
		}		
		
		$flag = 1;
		
		if($element eq "Other"){
			$Other_5pr = 0;
			$Other_3pr = 0;
			if($element_strand eq "+"){
				$Tstart = 11;
				$Tleft = 13;
				foreach($Other_count_start..$Other_count_end){
					if($_==$Other_count_start){
						if($reference->[$_][$Tstart]<$OTHER_CUT){
							$Other_5pr = 1;
						}
					}
					if($_==$Other_count_end){
						if($reference->[$_][$Tend]>($full-$OTHER_CUT)){
							$Other_3pr = 1;
						}
					}
				} 
			}else{
				$Tstart = 13;
				$Tleft = 11;
				for($i = $Other_count_end;$i>=$Other_count_start;$i--){
					if($i==$Other_count_end){
						if($reference->[$i][$Tstart]<$OTHER_CUT){
							$Other_5pr = 1;
						}
					}
					if($i==$Other_count_start){
						if($reference->[$i][$Tend]>($full-$OTHER_CUT)){
							$Other_3pr = 1;
						}
					}				
				}
			}			
			if($flag==0){
				print OUT_FR "$line\n";
				print OUT_FR3 "$info[0]\n";
				foreach(0..(scalar(@info)-2)){
					print OUT_FR2 "$info[$_]\n";
				}
				next;
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
		}   #end of other/SVA
		
		if($element =~ /LINE/){
			 $Other_5pr = 0;
			 $Other_3pr = 0;
			 $Other_between = 1;
			 $overlap = 0;
			 $distance = 0;
			if($element_strand eq "+"){
				$Tstart = 11;
				$Tleft = 13;
				foreach($Other_count_start..$Other_count_end){
					if($_==$Other_count_start){
						if($reference->[$_][$Tstart]<$LINE_CUT){
							$Other_5pr = 1;
						}
					}
					if($_==$Other_count_end){
						if($reference->[$_][$Tend]>($full-$LINE_CUT)){
							$Other_3pr = 1;
						}
					}
					if(!($_==$Other_count_start)){
						if(($reference->[$_][$Tstart]<$reference->[$_-1][$Tstart])|($reference->[$_][$Tend]<$reference->[$_-1][$Tend])){
							$flag = 0;
						}else{
							if($reference->[$_][$Tstart]<=$reference->[$_-1][$Tend]){
								$overlap=($reference->[$_-1][$Tend]-$reference->[$_][$Tstart]+1);
								if($reference->[$_][$Qstart]<=$reference->[$_-1][$Qend]){
									$overlap = ($overlap-($reference->[$_-1][$Qend]-$reference->[$_][$Qstart]+1));  #For cases where repeatmasker output has "*"
								}
								if(($overlap>$LINE_OVER*($reference->[$_-1][$Tend]-$reference->[$_-1][$Tstart]+1))|($overlap>$LINE_OVER*($reference->[$_][$Tend]-$reference->[$_][$Tstart]+1))){
									$flag = 0;
								}
							}else{
								$distance = ($reference->[$_][$Tstart]-$reference->[$_-1][$Tend]);
								if($distance >$LINE_CUT){
									$flag= 0;
								}
							}
						}
					}
				}
			}else{
				$Tstart = 13;
				$Tleft = 11;
				for($i = $Other_count_end;$i>=$Other_count_start;$i--){
					if($i==$Other_count_end){
						if($reference->[$i][$Tstart]<$LINE_CUT){
							$Other_5pr = 1;
						}
					}
					if($i==$Other_count_start){
						if($reference->[$i][$Tend]>($full-$LINE_CUT)){
							$Other_3pr = 1;
						}
					}
					if(!($i==$Other_count_end)){
						if(($reference->[$i][$Tstart]<$reference->[$i+1][$Tstart])|($reference->[$i][$Tend]<$reference->[$i+1][$Tend])){
							$flag = 0;
						}else{
							if($reference->[$i][$Tstart]<=$reference->[$i+1][$Tend]){
								$overlap=($reference->[$i+1][$Tend]-$reference->[$i][$Tstart]+1);
								if($reference->[$i][$Qstart]<=$reference->[$i+1][$Qend]){
									$overlap = ($overlap-($reference->[$i+1][$Qend]-$reference->[$i][$Qstart]+1));  #For cases where repeatmasker output has "*"
								}
								if(($overlap>$LINE_OVER*($reference->[$i+1][$Tend]-$reference->[$i+1][$Tstart]+1))|($overlap>$LINE_OVER*($reference->[$i][$Tend]-$reference->[$i][$Tstart]+1))){
									$flag = 0;
								}
							}else{
								$distance = ($reference->[$i][$Tstart]-$reference->[$i+1][$Tend]);
								if($distance >$LINE_CUT){
									$flag= 0;
								}
							}
						}
					}
				}
			}
			if($flag==0){
				print OUT_FR "$line\n";
				print OUT_FR3 "$info[0]\n";
				foreach(0..(scalar(@info)-2)){
					print OUT_FR2 "$info[$_]\n";
				}
				next;
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
		}	#End of LINE
			
		print OUT_WH "$info[0]\t$coverage\t$full\t",($coverage/$full),"\t$ther_full\t",($coverage/$ther_full),"\t$final\t$classOut\t$subOut\t$element_strand";		
		print OUT_WH2 "$info[0]\t$coverage\t$full\t",($coverage/$full),"\t$ther_full\t",($coverage/$ther_full),"\t$final\t$classOut\t$subOut\t$element_strand\n";
		print OUT_WH3 "$info[0]\t$coverage\t$full\t",($coverage/$full),"\t$ther_full\t",($coverage/$ther_full),"\t$final\t$classOut\t$subOut\t$element_strand\n";
		foreach(1..(scalar(@info)-2)){
			print OUT_WH2 "$info[$_]\n";
			print OUT_WH "\t\t\t$info[$_]";
		}
		print OUT_WH "\n";
	}else{    #This modifications deal with SVA+Alu (On the opposite strand) +...
		$flag = 0;
		if(($element eq "Other")|($element=~/Alu/)){
			$otherM = 0;
			$otherM2 = 0;
			$Other_count = 0;
			$Other_count_start = -1;
			$Other_count_end = -1;
			$Other_between = 1;
			$Other_5pr = 0;
			$Other_3pr = 0;
			$final = "N.A";
			foreach(0..($ref_count-1)){
				$Other_count++;
				$i = $_;
				if(!(($reference->[$i][$Tclass] =~/Low_complexity/)|($reference->[$i][$Tclass] =~/Satellite/)|($reference->[$i][$Tclass] =~/Simple_repeat/))){
					if($otherM==0){
						$Other_count_start=$i;
						$otherM2 = 1;
					}
					$otherM = 1;
				}else{
					if($otherM2==1){
						$Other_count_end = ($i-1);
						last;
					}
				}
			}
			if($Other_count==$ref_count){
				if($Other_count_end==-1){
					$Other_count_end = ($ref_count-1);
				}
			}
			if($reference->[$Other_count_start][$Tclass] =~/Alu/){
				$flag = 1;
				$Tstart = 11;
				$Tleft = 13;
				$full = ($reference->[$Other_count_start+1][$Tend]+$reference->[$Other_count_start+1][$Tleft]);
				$classOut = "SVA";
				$element_strand = $reference->[$Other_count_start+1][$strand];
				@subClassOut = ();
				push(@subClassOut,$reference->[$Other_count_start+1][$Tname]);
				$subOut = $reference->[$Other_count_start+1][$Tname];
				foreach(($Other_count_start+1)..($Other_count_end)){
					$i = $_;
					if(!($reference->[$i][$Tclass] eq "Other")){
						$flag = 0;
						last;
					}
					if(($reference->[$i][$strand] eq $reference->[$Other_count_start][$strand])){
						$flag = 0;
						last;
					}
					if($i==$Other_count_end){
						if($reference->[$i][$Tend]>($full-$OTHER_CUT)){
							$Other_3pr = 1;
							$final = "Full";
						}else{
							$final = "3prTruncated";
						}
					}
					$mark_rep = 0;
					foreach(0..(scalar(@subClassOut)-1)){
						next if((($reference->[$i][$Tclass] =~/Low_complexity/)|($reference->[$i][$Tclass] =~/Satellite/)|($reference->[$i][$Tclass] =~/Simple_repeat/)));
						$j = $_;
						if($subClassOut[$j] eq $reference->[$i][$Tname]){
							$mark_rep =1;
						}
					}
					if($mark_rep==0){
						push(@subClassOut,$reference->[$i][$Tname]);
						$subOut = $subOut.".".$reference->[$i][$Tname];
					}
				}
				$flag = 0 if($Other_count_end==$Other_count_start);   #Modified June 19 2009
			}elsif($reference->[$Other_count_end][$Tclass] =~/Alu/){
				$flag = 1;
				$Tstart = 13;
				foreach(($Other_count_start)..($Other_count_end-1)){
					$i = $_;
					if(!($reference->[$i][$Tclass] eq "Other")){
						$flag = 0;
						last;
					}
					if(($reference->[$i][$strand] eq $reference->[$Other_count_end][$strand])){
						$flag = 0;
						last;
					}
					if($i==$Other_count_start){
						if($reference->[$i][$Tend]>($full-$OTHER_CUT)){
							$Other_3pr = 1;
							$final = "Full";
						}else{
							$final = "3prTruncated";
						}
					}
					$mark_rep = 0;
					foreach(0..(scalar(@subClassOut)-1)){
						$j = $_;
						next if((($reference->[$i][$Tclass] =~/Low_complexity/)|($reference->[$i][$Tclass] =~/Satellite/)|($reference->[$i][$Tclass] =~/Simple_repeat/)));
						if($subClassOut[$j] eq $reference->[$i][$Tname]){
							$mark_rep =1;
						}
					}
					if($mark_rep==0){
						push(@subClassOut,$reference->[$i][$Tname]);
						$subOut = $subOut.".".$reference->[$i][$Tname];
					}
				}		
				$flag = 0 if($Other_count_end==$Other_count_start); #Modified June 19 2009
			}
		}
		if($flag==0){
			print OUT_FR "$line\n";
			print OUT_FR3 "$info[0]\n";
			foreach(0..(scalar(@info)-2)){
				print OUT_FR2 "$info[$_]\n";
			}
		}else{
			$coverage = ($reference->[$Other_count_end][$Qend]-$reference->[$Other_count_start][$Qstart]+1);
			$full = -1;
			$ther_full = 3000;
			print OUT_WH "$info[0]\t$coverage\t$full\t",($coverage/$full),"\t$ther_full\t",($coverage/$ther_full),"\t$final\t$classOut\t$subOut\t$element_strand";		
			print OUT_WH2 "$info[0]\t$coverage\t$full\t",($coverage/$full),"\t$ther_full\t",($coverage/$ther_full),"\t$final\t$classOut\t$subOut\t$element_strand\n";
			print OUT_WH3 "$info[0]\t$coverage\t$full\t",($coverage/$full),"\t$ther_full\t",($coverage/$ther_full),"\t$final\t$classOut\t$subOut\t$element_strand\n";
			foreach(1..(scalar(@info)-2)){
				print OUT_WH2 "$info[$_]\n";
				print OUT_WH "\t\t\t$info[$_]";
			}
			print OUT_WH "\n";						
		}
	}
}
close MULTI;






