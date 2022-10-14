#!/usr/bin/perl -w
use strict;
use POSIX qw(ceil floor);
#This program uses the -100~+100 bp sequences of the left & right breakpoint of the Inversions in the Novel Genome
#This program excutes "blseq" and prints the calls with >50bp sequence, >85%identity
#Check whether homologous regions spans the breakpoints
#Check whether the homologous region doesn't have more than $OFFSET of offset at the two breakpoints
unless ($ARGV[6]) {  
	print STDERR "Exit!\n";
	exit;
}

my @fields;
my @bl2seq;
my $identity;
my $Strand;
my $Query;
my $Sbjct;
my @strand;
my @query;
my @sbjct;
my @homo;
my $homology;
my @percent;
my $percentage;
my $n;
my @elements;
my $line;
my $mark;
my $line1;
my @offsets;
my $mark2;
my @numbers;
my $i;
my @all_line_query;
my @all_line_sbjct;
my $j;
my $k;
my @temp;
my $query_lines;
my $sbjct_lines;
my $length_f = $ARGV[2]/2;
my $start_Query;
my $start_Sbjct;
my $end_Sbjct;
my $end_Query;
my $position_query;
my $position_sbjct;
my $mark3;
my $offset;
my $minoffset;
my $OFFSET = 20;
my $minoffset2;
my $lambda;
my $HOMO = $ARGV[3];
my $PCT = $ARGV[4];
my $out_dir = $ARGV[5];
my $BLAST_dir = $ARGV[6];
my $SV_Size; #Feb 2010
my $Chr = 0;
my $Start = 1;
my $End = 2;
my $Left_F = 3;
my $Right_F = 4;
my $ID = 5;

open FLANKINGSEQ, "<$ARGV[0]";
open NAHR, ">$out_dir/Output/Inversions/$ARGV[1].NAHR_Inversions.txt";
open OTHER, ">$out_dir/$ARGV[1].Other_Non_NAHR_Inversions.txt";
open NOTSURE,">$out_dir/Output/Inversions/$ARGV[1].UNSURE_3_OffsetOrNOtSpanBP_In_NAHR_Inversions.txt";
open NOTSURE2,">$out_dir/Output/Inversions/$ARGV[1].UNSURE_3_OffsetOrNOtSpanBP_In_NAHR_Inversions_bl2seqDetail.txt";
open NOTSURE3,">$out_dir/Output/Inversions/$ARGV[1].UNSURE_4_Potential_NAHR_Inversions.txt";
open NOTSURE4,">$out_dir/Output/Inversions/$ARGV[1].UNSURE_4_Potential_NAHR_Inversions_bl2seqDetail.txt";

print NAHR "ChrNo.\tStart\tEnd\tLeftFlankingSequence\tRightFlankingSequence\tID\tIdentity\tHomology\tPercentage\n";
print OTHER "ChrNo.\tStart\tEnd\tLeftFlankingSequence\tRightFlankingSequence\tID\n";
print NOTSURE "ChrNo.\tStart\tEnd\tLeftFlankingSequence\tRightFlankingSequence\tID\n";
print NOTSURE2 "ChrNo.\tStart\tEnd\tLeftFlankingSequence\tRightFlankingSequence\tID\n";
print NOTSURE3 "ChrNo.\tStart\tEnd\tLeftFlankingSequence\tRightFlankingSequence\tID\n";
print NOTSURE4 "ChrNo.\tStart\tEnd\tLeftFlankingSequence\tRightFlankingSequence\tID\n";

<FLANKINGSEQ>;
while($line1=<FLANKINGSEQ>){
	chomp($line1);
	@fields = ();
	@fields = split/\s+/, $line1;
	$SV_Size = ($fields[$End]-$fields[$Start]); #Feb 2010
	open INPUT1, ">$out_dir/Temp/Input1.$fields[$ID].$ARGV[1].txt";
	open INPUT2, ">$out_dir/Temp/Input2.$fields[$ID].$ARGV[1].txt";
	print INPUT1 "$fields[$Left_F]\n";
	print INPUT2 "$fields[$Right_F]\n";
	close INPUT1;
	close INPUT2;
	system "$BLAST_dir/bl2seq -i $out_dir/Temp/Input1.$fields[$ID].$ARGV[1].txt -j $out_dir/Temp/Input2.$fields[$ID].$ARGV[1].txt -p blastn -o $out_dir/Temp/bl2seq.$fields[$ID].$ARGV[1].txt -F F";
	system "rm -rf $out_dir/Temp/Input1.$fields[$ID].$ARGV[1].txt $out_dir/Temp/Input2.$fields[$ID].$ARGV[1].txt";
	open BL2SEQ, "<$out_dir/Temp/bl2seq.$fields[$ID].$ARGV[1].txt";

	@numbers = ();
	$numbers[0] = "";

	$n = 1;
	$lambda = 0;
	while(<BL2SEQ>){
		@elements = ();
		@elements = split/\s+/, $_;
		foreach (@elements){
			if($_ eq "Score"){
				if($numbers[0] eq ""){
					$numbers[0] = $n;
				}else{
					push(@numbers, $n);
				}    		  	
			}
			if($_ eq "Lambda"){
				if(!($numbers[0] eq "")){
					push(@numbers, ($n+1));
				}
				$lambda = 1;
			}
		}
		if($lambda == 1){
			last;
		}
		$n++;
	}
	close BL2SEQ;
	
	if($numbers[0] eq ""){
		 print OTHER "$fields[$Chr]\t$fields[$Start]\t$fields[$End]\t$fields[$Left_F]\t$fields[$Right_F]\t$fields[$ID]\tNo_hits\n";
	}else{								
		open BL2SEQ, "<$out_dir/Temp/bl2seq.$fields[$ID].$ARGV[1].txt";
		$i = 0;
		$n = 1;
		$mark = 0;
		$mark2 = 0;
		$mark3 = 0;
		$minoffset = ($length_f*2+1);
		$minoffset2 = ($length_f*2+1);
		while($i < @numbers-1){
			$query_lines = "";
			$sbjct_lines = "";
		
			while($n<$numbers[$i+1]){
				$line =<BL2SEQ>;
				chomp($line);
				if($n == ($numbers[$i]+1)){
					$identity = $line;
					@bl2seq = ();
					@bl2seq = split/\s+/, $identity;
					@homo = ();
					@homo = split//, $bl2seq[3];
					$homology = "";
					foreach(@homo){
						if($_ eq "/"){
							last;
						}
						$homology= $homology.$_;
					}
					@percent = ();
					@percent = split//, $bl2seq[4];
					$percentage = "";
					foreach(@percent){
						if(!(($_ eq "%")||($_ eq "(")||($_ eq ")")||($_ eq ","))){
							$percentage= $percentage.$_;
						}
					}
				}
				if($n == ($numbers[$i]+2)){
					$Strand = $line;
					@strand = ();
					@strand = split/\s+/, $Strand;
						}
				if($n == ($numbers[$i]+5)){   		
					@temp = ();
					@temp = split/\s+/, $line;
					$start_Query = $temp[1];
				}   
				if($n == ($numbers[$i]+7)){   		
					@temp = ();
					@temp = split/\s+/, $line;
					$start_Sbjct = $temp[1];
				}   
				if($n == ($numbers[$i+1]-6)){   		
					@temp = ();
					@temp = split/\s+/, $line;
					$end_Query = $temp[3];
				}   
				if($n == ($numbers[$i+1]-4)){   		
					@temp = ();
					@temp = split/\s+/, $line;
					$end_Sbjct = $temp[3];
				}
				if($n>$numbers[$i]){
					if(($n- $numbers[$i])%5 ==0){
						@temp = ();
						@temp = split/\s+/, $line;
						if(!($line eq "")){
							$query_lines = $query_lines.$temp[2];
						}
					}		
					if(($n- ($numbers[$i]+7))%5 ==0){
						@temp = ();
						@temp = split/\s+/, $line;
						$sbjct_lines = $sbjct_lines.$temp[2];
					}
				}
			
				if($n == ($numbers[$i+1]-1)){

					if (($homology >=$HOMO) && ($percentage >= $PCT)&&(!($strand[3] eq $strand[5]))){
						$mark = 1;
						if($SV_Size >= 2*$length_f){   #Feb 2010		 
							if(($start_Query <= ($length_f+1))&&($start_Sbjct >= ($length_f))&&($end_Query >= ($length_f))&&($end_Sbjct <= ($length_f+1))){   #check whether spanning breakpoint
								$k = 0;
								@temp = ();
								@temp = split//, $query_lines;
								$j = 0;
								while($k <= ($length_f-$start_Query+1)){
									if($start_Query == ($length_f+1)){    #Modified June 2009
										$position_query = 0;
										last;
									}
									if($k == ($length_f-$start_Query+1)){
										$position_query = $j;
										last;
									}                               
									if(!($temp[$j] eq "-")){
										$k++;
									}
									$j++;
								}
								$k = 0;
								@temp = ();
								@temp = split//, $sbjct_lines;
								$j = 0;
								while($k <= ($start_Sbjct-$length_f)){
									if($start_Sbjct == ($length_f)){      #Modified June 2009
										$position_sbjct = 0;
										last;
									}
									if($k == ($start_Sbjct-$length_f)){
										$position_sbjct = $j;
										last;
									}                             
									if(!($temp[$j] eq "-")){
										$k++;
									}
									$j++;
								}
								if($position_sbjct>$position_query){
									$offset = ($position_sbjct - $position_query);
								}else{
									$offset = ($position_query-$position_sbjct);
								}
								if($offset < $minoffset){
									$minoffset = $offset;
								}
					
								if($minoffset<=$OFFSET){
									if($mark2 == 0){
										print NAHR "$fields[$Chr]\t$fields[$Start]\t$fields[$End]\t$fields[$Left_F]\t$fields[$Right_F]\t$fields[$ID]\t$identity\t$homology\t$percentage\n";	
										$mark2 = 1;
									}
								}
							}
						}else{
							if(($start_Query <= ($length_f+1))&&($start_Sbjct >= floor($SV_Size/2))&&($end_Query >= ($length_f))&&($end_Sbjct <= (floor($SV_Size/2)+1))){   #check whether spanning breakpoint
								$k = 0;
								@temp = ();
								@temp = split//, $query_lines;
								$j = 0;
								while($k <= ($length_f-$start_Query+1)){
									if($start_Query == ($length_f+1)){    #Modified June 2009
										$position_query = 0;
										last;
									}
									if($k == ($length_f-$start_Query+1)){
										$position_query = $j;
										last;
									}                               
									if(!($temp[$j] eq "-")){
										$k++;
									}
									$j++;
								}
								$k = 0;
								@temp = ();
								@temp = split//, $sbjct_lines;
								$j = 0;
								while($k <= ($start_Sbjct-floor($SV_Size/2))){
									if($start_Sbjct == floor($SV_Size/2)){      #Modified June 2009
										$position_sbjct = 0;
										last;
									}
									if($k == ($start_Sbjct-floor($SV_Size/2))){
										$position_sbjct = $j;
										last;
									}                             
									if(!($temp[$j] eq "-")){
										$k++;
									}
									$j++;
								}
								if($position_sbjct>$position_query){
									$offset = ($position_sbjct - $position_query);
								}else{
									$offset = ($position_query-$position_sbjct);
								}
								if($offset < $minoffset){
									$minoffset = $offset;
								}
					
								if($minoffset<=$OFFSET){
									if($mark2 == 0){
										print NAHR "$fields[$Chr]\t$fields[$Start]\t$fields[$End]\t$fields[$Left_F]\t$fields[$Right_F]\t$fields[$ID]\t$identity\t$homology\t$percentage\n";	
										$mark2 = 1;
									}
								}
							}						
						}
					}
				
					if((($homology <$HOMO) &&($homology >30)&& ($percentage >= $PCT)&&(!($strand[3] eq $strand[5])))||(($homology >=$HOMO) && ($percentage < $PCT)&&(!($strand[3] eq $strand[5])))){
						if($SV_Size >= 2*$length_f){
							if(($start_Query <= ($length_f+1))&&($start_Sbjct >= ($length_f-1))&&($end_Query >= ($length_f-1))&&($end_Sbjct <= ($length_f+1))){
								$k = 0;
								@temp = ();
								@temp = split//, $query_lines;
								$j = 0;
								while($k <= ($length_f-$start_Query+1)){
									if($start_Query == ($length_f+1)){    #Modified June 2009
										$position_query = 0;
										last;
									}
									if($k == ($length_f-$start_Query+1)){
										$position_query = $j;
										last;
									}                               
									if(!($temp[$j] eq "-")){
										$k++;
									}
									$j++;
								}
								$k = 0;
								@temp = ();
								@temp = split//, $sbjct_lines;
								$j = 0;
								while($k <= ($start_Sbjct-$length_f)){
									if($start_Sbjct == ($length_f)){      #Modified June 2009
										$position_sbjct = 0;
										last;
									}
									if($k == ($start_Sbjct-$length_f)){
										$position_sbjct = $j;
										last;
									}                             
									if(!($temp[$j] eq "-")){
										$k++;
									}
									$j++;
								}
								if($position_sbjct>$position_query){
									$offset = ($position_sbjct - $position_query);
								}else{
									$offset = ($position_query-$position_sbjct);
								}
								if($offset < $minoffset2){
									$minoffset2 = $offset;
								}
					
								if($minoffset2<=$OFFSET){
									if($mark3 == 0){
										$mark3 = 1;
									}
								}
							}
						}else{
							if(($start_Query <= ($length_f+1))&&($start_Sbjct >= floor($SV_Size/2))&&($end_Query >= ($length_f))&&($end_Sbjct <= (floor($SV_Size/2)+1))){   #check whether spanning breakpoint
								$k = 0;
								@temp = ();
								@temp = split//, $query_lines;
								$j = 0;
								while($k <= ($length_f-$start_Query+1)){
									if($start_Query == ($length_f+1)){    #Modified June 2009
										$position_query = 0;
										last;
									}
									if($k == ($length_f-$start_Query+1)){
										$position_query = $j;
										last;
									}                               
									if(!($temp[$j] eq "-")){
										$k++;
									}
									$j++;
								}
								$k = 0;
								@temp = ();
								@temp = split//, $sbjct_lines;
								$j = 0;
								while($k <= ($start_Sbjct-floor($SV_Size/2))){
									if($start_Sbjct == (floor($SV_Size/2))){      #Modified June 2009
										$position_sbjct = 0;
										last;
									}
									if($k == ($start_Sbjct-floor($SV_Size/2))){
										$position_sbjct = $j;
										last;
									}                             
									if(!($temp[$j] eq "-")){
										$k++;
									}
									$j++;
								}
								if($position_sbjct>$position_query){
									$offset = ($position_sbjct - $position_query);
								}else{
									$offset = ($position_query-$position_sbjct);
								}
								if($offset < $minoffset2){
									$minoffset2 = $offset;
								}
					
								if($minoffset2<=$OFFSET){
									if($mark3 == 0){
										$mark3 = 1;
									}
								}
							}
						
						}
					}
				}
				$n++;
			}
			$i++;
		} 
				#mark: see whether >=$HOMO bp & >=85% identity
				#mark2: See whether have large offset or do not span BP
				#mark3: See whether have (1)30bp< length <50bp, >=85%, no large offset, spaning BP;
																#or (2) length >= 50bp, <85%, no large offset, spaning BP;
		if(($mark == 0)&&($mark3 == 0)){
			print OTHER "$fields[$Chr]\t$fields[$Start]\t$fields[$End]\t$fields[$Left_F]\t$fields[$Right_F]\t$fields[$ID]\n";
		}
		#NOTSURE holds">=50bp & >=85% identity", but have large offset or do not span BP: might be mistake in the novel genome assembly
		#NOTSURE2 holds all NOTSURE Cases & their detail bl2sseq results
		if(($mark ==1)&&($mark2 == 0)){
			print NOTSURE "$fields[$Chr]\t$fields[$Start]\t$fields[$End]\t$fields[$Left_F]\t$fields[$Right_F]\t$fields[$ID]\noffset=$minoffset\n";
			print NOTSURE2 "$fields[$Chr]\t$fields[$Start]\t$fields[$End]\t$fields[$Left_F]\t$fields[$Right_F]\t$fields[$ID]\noffset=$minoffset\n";
			open BL2SEQ, "<$out_dir/Temp/bl2seq.$fields[$ID].$ARGV[1].txt";
			while(<BL2SEQ>){
				print NOTSURE2 "$_";
			}
			close BL2SEQ;
		}
		#NOTSURE3 holds calls with either of the two mark3 conditions: might be potential NAHRs
		#NOTSURE4 holds all NOTSURE3 Cases & their detail bl2sseq results
		if(($mark ==0)&&($mark3 == 1)){
			print NOTSURE3 "$fields[$Chr]\t$fields[$Start]\t$fields[$End]\t$fields[$Left_F]\t$fields[$Right_F]\t$fields[$ID]\n";
			print NOTSURE4 "$fields[$Chr]\t$fields[$Start]\t$fields[$End]\t$fields[$Left_F]\t$fields[$Right_F]\t$fields[$ID]\n";
			open BL2SEQ, "<$out_dir/Temp/bl2seq.$fields[$ID].$ARGV[1].txt";
			while(<BL2SEQ>){
				print NOTSURE4 "$_";
			}
			close BL2SEQ;
		}       	
	}
	system "rm -rf $out_dir/Temp/bl2seq.$fields[$ID].$ARGV[1].txt"; #clean
}
						
						


close FLANKINGSEQ;
close NAHR;
close OTHER;
close NOTSURE;
close NOTSURE2;
close NOTSURE3;
close NOTSURE4;
