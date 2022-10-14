#!/usr/bin/perl -w
use strict;

my $N1=0;
my $N2=0;
my $N3=0;
my $N4=0;
my $N5=0;
my $N6=0;
my $N7=0;
my $N8=0;
my $N9=0;
my $N10=0;
my $N11=0;
my $N12=0;
my $N13=0;
my $N14=0;
my $N15=0;
my $N16=0;
my $N17=0;
my $N18=0;
my $N19=0;
my $N20=0;
my $N21=0;
my $N22=0;
my $N23=0;
my $N24=0;
my $N25=0;
my $N26=0;
my $N27=0;
my $N28=0;
my $N29=0;
my $N30=0;
my $m1=0;
my $m2=0;
my $m3=0;
my $m4=0;
my $m5=0;
my $m6=0;
my $m7=0;
my $m8=0;
my $m9=0;
my $m10=0;
my $m11=0;
my $m12=0;
my $m13=0;
my $m14=0;
my $m15=0;
my $m16=0;
my $m17=0;
my $m18=0;
my $m19=0;
my $m20=0;
my $m21=0;
my $m22=0;
my $m23=0;
my $m24=0;
my $m25=0;
my $m26=0;
my $m27=0;
my $m28=0;
my $m29=0;
my $m30=0;
my @fields;

my $label = $ARGV[0];
my $out_dir = $ARGV[1];

$N1 = (open NAHR_Deletions, "<$out_dir/Output/Deletions/$label.NAHR_Deletions.txt");
$N2 = (open SINGLE_RT_Deletions,"<$out_dir/Output/Deletions/$label.MicroHomologySingleFragRT_Deletions.txt");
$N3 = (open NHEJ_Deletions, "<$out_dir/Output/Deletions/$label.MicroHomology_PolyA_NHEJ_Deletions.txt");
$N4 = (open MULTI_RT_Deletions, "<$out_dir/Output/Deletions/$label.MicroHomologyFragMultiple_Retrotransposon_Deletions.txt");
$N5 = (open NAHR_Insertions, "<$out_dir/Output/Insertions/$label.NAHR_Insertions.txt");
$N6 = (open SINGLE_RT_Insertions,"<$out_dir/Output/Insertions/$label.MicroHomologySingleFragRT_Insertions.txt");
$N7 = (open NHEJ_Insertions, "<$out_dir/Output/Insertions/$label.MicroHomology_PolyA_NHEJ_Insertions.txt");
$N8 = (open MULTI_RT_Insertions, "<$out_dir/Output/Insertions/$label.MicroHomologyFragMultiple_Retrotransposon_Insertions.txt");
$N9 = (open NAHR_Inversions, "<$out_dir/Output/Inversions/$label.NAHR_Inversions.txt");
$N10 = (open NHEJ_Inversions, "<$out_dir/Output/Inversions/$label.MicroHomology_NHEJ_Inversions.txt");
$N11 = (open UNSURE1_Deletions, "<$out_dir/Output/Deletions/$label.UNSURE_1_Smaller_Than_Cutoff_Deletions.txt");
$N12 = (open UNSURE2_Deletions, "<$out_dir/Output/Deletions/$label.UNSURE_2_500bp_Proximity_To_Either_Telemere_Deletions.txt");
$N13 = (open UNSURE3_Deletions, "<$out_dir/Output/Deletions/$label.UNSURE_3_OffsetOrNOtSpanBP_In_NAHR_Deletions.txt");
$N14 = (open UNSURE4_Deletions, "<$out_dir/Output/Deletions/$label.UNSURE_4_Potential_NAHR_Deletions.txt");
$N15 = (open UNSURE1_Insertions, "<$out_dir/Output/Insertions/$label.UNSURE_1_Smaller_Than_Cutoff_Insertions.txt");
$N16 = (open UNSURE2_Insertions, "<$out_dir/Output/Insertions/$label.UNSURE_2_500bp_Proximity_To_Either_Telemere_Insertions.txt");
$N17 = (open UNSURE3_Insertions, "<$out_dir/Output/Insertions/$label.UNSURE_3_OffsetOrNOtSpanBP_In_NAHR_Insertions.txt");
$N18 = (open UNSURE4_Insertions, "<$out_dir/Output/Insertions/$label.UNSURE_4_Potential_NAHR_Insertions.txt");
$N19 = (open UNSURE1_Inversions, "<$out_dir/Output/Inversions/$label.UNSURE_1_Smaller_Than_Cutoff_Inversions.txt");
$N20 = (open UNSURE2_Inversions, "<$out_dir/Output/Inversions/$label.UNSURE_2_500bp_Proximity_To_Either_Telemere_Inversions.txt");
$N21 = (open UNSURE3_Inversions, "<$out_dir/Output/Inversions/$label.UNSURE_3_OffsetOrNOtSpanBP_In_NAHR_Inversions.txt");
$N22 = (open UNSURE4_Inversions, "<$out_dir/Output/Inversions/$label.UNSURE_4_Potential_NAHR_Inversions.txt");
$N23 = (open UNSURE5_Deletions, "<$out_dir/Output/Deletions/$label.UNSURE_5_Single_Retrotransposon_Deletions.txt");
$N24 = (open UNSURE6_Deletions, "<$out_dir/Output/Deletions/$label.UNSURE_6_Multiple_Retrotransposon_Deletions.txt");
$N25 = (open UNSURE7_Deletions, "<$out_dir/Output/Deletions/$label.UNSURE_7_Potential_Processed_Pseudogene_Deletions.txt");
$N26 = (open UNSURE5_Insertions, "<$out_dir/Output/Insertions/$label.UNSURE_5_Single_Retrotransposon_Insertions.txt");
$N27 = (open UNSURE6_Insertions, "<$out_dir/Output/Insertions/$label.UNSURE_6_Multiple_Retrotransposon_Insertions.txt");
$N28 = (open UNSURE7_Insertions, "<$out_dir/Output/Insertions/$label.UNSURE_7_Potential_Processed_Pseudogene_Insertions.txt");
$N29 = (open VNTR_Deletions, "$out_dir/Output/Deletions/$ARGV[0].VNTRDeletions-50.txt");
$N30 = (open VNTR_Insertions, "<$out_dir/Output/Insertions/$ARGV[0].VNTRInsertions-50.txt");

my $TOTAL = 21;
my $coverage = 5;
my $class = 11;
my $nmhNo = 13;
my $polyAonCorrStrand = 15;
my $StrandElement= 16;
my $TSD = 17;
my $ClassLabel = 18;
my $strandPolyAtail = 20;
my $UNSURE = 21;

my $NUM_NAHR = $coverage;
my $NUM_singleRT1 = ($coverage+1);
my $NUM_singleRT2 = ($ClassLabel+1);
my $NUM_NHR = $coverage;
my $NUM_NHR2 = $nmhNo;  
my $NUM_NHR3 = $ClassLabel;
my $NUM_NHR4 = $polyAonCorrStrand;
my $NUM_MultiRT0 = $coverage;
my $NUM_MultiRT1 = $class;
my $NUM_MultiRT2 = $UNSURE;
my $NUM_UNSURE1 = $coverage;
my $NUM_UNSURE2 = $UNSURE;
my $count;

open STATISTIC, ">$out_dir/Output/$label.Statistics.txt";

open COMBINE, ">$out_dir/Output/$label.All_Output.txt";
print COMBINE "chrNo.\tStart\tEnd\tSV_Class\tMechanism\tcoverage\tfull\tratio\ttherFull\tratio2\tfinal\tClass\tNames\tnmhNo.\tnmhSeq\tpolyAonCorrStrand\tStrandElement\tTSD\tClassLabel\tRTnames\tstrandPolyAtail\tUNSUREclass\tID\tInsertedSeq\n";
if($N1 == 1){
	<NAHR_Deletions>;
	while(<NAHR_Deletions>){
		$m1++;
		chomp;
		@fields = ();
		@fields = split/\s+/, $_;
		print COMBINE "$fields[0]\t$fields[1]\t$fields[2]\tDeletion\tNAHR";
		foreach($NUM_NAHR..$TOTAL){
			print COMBINE "\tN.A";
		}
		print COMBINE "\t$fields[5]\tN.A\n";
	}
}


if($N2 == 1){
	<SINGLE_RT_Deletions>;
	while(<SINGLE_RT_Deletions>){
		$m2++;
		chomp;
		@fields = ();
		@fields = split/\s+/, $_;
		print COMBINE "$fields[0]\t$fields[1]\t$fields[2]\tDeletion\tSTEI";
		foreach($NUM_singleRT1..$TOTAL){
			if($_<$NUM_singleRT2){
				print COMBINE "\t$fields[$_+1]";
			}else{
				print COMBINE "\tN.A";
			}
		}
		print COMBINE "\tN.A\t$fields[6]\tN.A\n";
	}
}

if($N3 == 1){
	<NHEJ_Deletions>;
	while(<NHEJ_Deletions>){
		$m3++;
		chomp;
		@fields = ();
		@fields = split/\s+/, $_;
		print COMBINE "$fields[0]\t$fields[1]\t$fields[2]\tDeletion\tNHR";	
		foreach($NUM_NHR..$TOTAL){
			if($_<$NUM_NHR2){
				print COMBINE "\tN.A";
			}elsif($_<=$NUM_NHR4){
				print COMBINE "\t$fields[$_-7]";
			}
		}
		print COMBINE "\tN.A\t$fields[9]\tN.A\tN.A\t$fields[10]\tN.A\t$fields[11]\tN.A\n";			
	}
}

if($N4 == 1){
	<MULTI_RT_Deletions>;
	while(<MULTI_RT_Deletions>){
		$m4++;
		chomp;
		@fields = ();
		@fields = split/\s+/, $_;
		print COMBINE "$fields[0]\t$fields[1]\t$fields[2]\tDeletion\tMTEI";
		foreach($NUM_MultiRT0..$TOTAL){
			if($_<$NUM_MultiRT1){
				print COMBINE "\tN.A";
			}elsif($_<$NUM_MultiRT2){
				print COMBINE "\t$fields[$_-4]";
			}
		}
		print COMBINE "\tN.A\t$fields[6]\tN.A\n";			
	}
}

if($N5 == 1){
	<NAHR_Insertions>;
	while(<NAHR_Insertions>){
		$m5++;
		chomp;
		@fields = ();
		@fields = split/\s+/, $_;
		print COMBINE "$fields[0]\t$fields[1]\t$fields[2]\tInsertion\tNAHR";
		foreach($NUM_NAHR..$TOTAL){
			print COMBINE "\tN.A";
		}
		print COMBINE "\t$fields[6]\t$fields[5]\n";
	}
}


if($N6 == 1){
	<SINGLE_RT_Insertions>;
	while(<SINGLE_RT_Insertions>){
		$m6++;
		chomp;
		@fields = ();
		@fields = split/\s+/, $_;
		print COMBINE "$fields[0]\t$fields[1]\t$fields[2]\tInsertion\tSTEI";
		foreach($NUM_singleRT1..$TOTAL){
			if($_<$NUM_singleRT2){
				print COMBINE "\t$fields[$_+1]";
			}else{
				print COMBINE "\tN.A";
			}
		}
		print COMBINE "\tN.A\t$fields[7]\t$fields[6]\n";
	}
}

if($N7 == 1){
	<NHEJ_Insertions>;
	while(<NHEJ_Insertions>){
		$m7++;
		chomp;
		@fields = ();
		@fields = split/\s+/, $_;
		print COMBINE "$fields[0]\t$fields[1]\t$fields[2]\tInsertion\tNHR";	
		foreach($NUM_NHR..$TOTAL){
			if($_<$NUM_NHR2){
				print COMBINE "\tN.A";
			}elsif($_<=$NUM_NHR4){
				print COMBINE "\t$fields[$_-7]";
			}
		}
		print COMBINE "\tN.A\t$fields[9]\tN.A\tN.A\t$fields[10]\tN.A\t$fields[12]\t$fields[11]\n";			
	}
}

if($N8 == 1){
	<MULTI_RT_Insertions>;
	while(<MULTI_RT_Insertions>){
		$m8++;
		chomp;
		@fields = ();
		@fields = split/\s+/, $_;
		print COMBINE "$fields[0]\t$fields[1]\t$fields[2]\tInsertion\tMTEI";
			foreach($NUM_MultiRT0..$TOTAL){
			if($_<$NUM_MultiRT1){
				print COMBINE "\tN.A";
			}elsif($_<$NUM_MultiRT2){
				print COMBINE "\t$fields[$_-3]";
			}
		}
		print COMBINE "\tN.A\t$fields[7]\t$fields[6]\n";			
	}
}

if($N9 == 1){
	<NAHR_Inversions>;
	while(<NAHR_Inversions>){
		$m9++;
		chomp;
		@fields = ();
		@fields = split/\s+/, $_;
		print COMBINE "$fields[0]\t$fields[1]\t$fields[2]\tInversion\tNAHR";
		foreach($NUM_NAHR..$TOTAL){
			print COMBINE "\tN.A";
		}
		print COMBINE "\t$fields[5]\tN.A\n";
	}
}


if($N10 == 1){
	<NHEJ_Inversions>;
	while(<NHEJ_Inversions>){
		$m10++;
		chomp;
		@fields = ();
		@fields = split/\s+/, $_;
		print COMBINE "$fields[0]\t$fields[1]\t$fields[2]\tInversion\tNHR";	
		foreach($NUM_NHR..$TOTAL){
			if($_<$NUM_NHR2){
				print COMBINE "\tN.A";
			}elsif($_<$NUM_NHR4){
				print COMBINE "\t$fields[$_-8]";
			}else{
				print COMBINE "\tN.A";
			}
		}
		print COMBINE "\t$fields[7]\tN.A\n";
	}
}

if($N11 == 1){
	<UNSURE1_Deletions>;
	while(<UNSURE1_Deletions>){
		$m11++;
		chomp;
		@fields = ();
		@fields = split/\s+/, $_;
		print COMBINE "$fields[0]\t$fields[1]\t$fields[2]\tDeletion\tUNSURE";
		foreach($NUM_UNSURE1..($NUM_UNSURE2-1)){
			print COMBINE "\tN.A";
		}
		print COMBINE "\t1\t$fields[3]\tN.A\n";
	}
}

if($N12 == 1){
	<UNSURE2_Deletions>;
	while(<UNSURE2_Deletions>){
		$m12++;
		chomp;
		@fields = ();
		@fields = split/\s+/, $_;
		print COMBINE "$fields[0]\t$fields[1]\t$fields[2]\tDeletion\tUNSURE";
		foreach($NUM_UNSURE1..($NUM_UNSURE2-1)){
			print COMBINE "\tN.A";
		}
		print COMBINE "\t2\t$fields[3]\tN.A\n";
	}
}

if($N13 == 1){
	<UNSURE3_Deletions>;
	while(<UNSURE3_Deletions>){
		$m13++;
		chomp;
		@fields = ();
		@fields = split/\s+/, $_;
		if(scalar(@fields)>2){
		print COMBINE "$fields[0]\t$fields[1]\t$fields[2]\tDeletion\tUNSURE";
		foreach($NUM_UNSURE1..($NUM_UNSURE2-1)){
			print COMBINE "\tN.A";
		}
		print COMBINE "\t3\t$fields[5]\tN.A\n";
		}
	}
}

if($N14 == 1){
	<UNSURE4_Deletions>;
	while(<UNSURE4_Deletions>){
		$m14++;
		chomp;
		@fields = ();
		@fields = split/\s+/, $_;
		print COMBINE "$fields[0]\t$fields[1]\t$fields[2]\tDeletion\tUNSURE";
		foreach($NUM_UNSURE1..($NUM_UNSURE2-1)){
			print COMBINE "\tN.A";
		}
		print COMBINE "\t4\t$fields[5]\tN.A\n";
	}
}

if($N15 == 1){
	<UNSURE1_Insertions>;
	while(<UNSURE1_Insertions>){
		$m15++;
		chomp;
		@fields = ();
		@fields = split/\s+/, $_;
		print COMBINE "$fields[0]\t$fields[1]\t$fields[2]\tInsertion\tUNSURE";
		foreach($NUM_UNSURE1..($NUM_UNSURE2-1)){
			print COMBINE "\tN.A";
		}
		print COMBINE "\t1\t$fields[4]\t$fields[3]\n";
	}
}

if($N16 == 1){
	<UNSURE2_Insertions>;
	while(<UNSURE2_Insertions>){
		$m16++;
		chomp;
		@fields = ();
		@fields = split/\s+/, $_;
		print COMBINE "$fields[0]\t$fields[1]\t$fields[2]\tInsertion\tUNSURE";
		foreach($NUM_UNSURE1..($NUM_UNSURE2-1)){
			print COMBINE "\tN.A";
		}
		print COMBINE "\t2\t$fields[4]\t$fields[3]\n";
	}
}

if($N17 == 1){
	<UNSURE3_Insertions>;
	while(<UNSURE3_Insertions>){
		$m17++;
		chomp;
		@fields = ();
		@fields = split/\s+/, $_;
		if(scalar(@fields)>2){
			print COMBINE "$fields[0]\t$fields[1]\t$fields[2]\tInsertion\tUNSURE";
			foreach($NUM_UNSURE1..($NUM_UNSURE2-1)){
				print COMBINE "\tN.A";
			}
		print COMBINE "\t3\t$fields[6]\t$fields[5]\n";
		}
	}
}

if($N18 == 1){
	<UNSURE4_Insertions>;
	while(<UNSURE4_Insertions>){
		$m18++;
		chomp;
		@fields = ();
		@fields = split/\s+/, $_;
		print COMBINE "$fields[0]\t$fields[1]\t$fields[2]\tInsertion\tUNSURE";
		foreach($NUM_UNSURE1..($NUM_UNSURE2-1)){
			print COMBINE "\tN.A";
		}
		print COMBINE "\t4\t$fields[6]\t$fields[5]\n";
	}
}


if($N19 == 1){
	<UNSURE1_Inversions>;
	while(<UNSURE1_Inversions>){
		$m19++;
		chomp;
		@fields = ();
			@fields = split/\s+/, $_;
		print COMBINE "$fields[0]\t$fields[1]\t$fields[2]\tInversion\tUNSURE";
		foreach($NUM_UNSURE1..($NUM_UNSURE2-1)){
			print COMBINE "\tN.A";
		}
		print COMBINE "\t1\t$fields[3]\tN.A\n";
	}
}

if($N20 == 1){
	<UNSURE2_Inversions>;
	while(<UNSURE2_Inversions>){
		$m20++;
		chomp;
		@fields = ();
		@fields = split/\s+/, $_;
		print COMBINE "$fields[0]\t$fields[1]\t$fields[2]\tInversion\tUNSURE";
		foreach($NUM_UNSURE1..($NUM_UNSURE2-1)){
			print COMBINE "\tN.A";
		}
		print COMBINE "\t2\t$fields[3]\tN.A\n";
	}
}

if($N21 == 1){
	<UNSURE3_Inversions>;
	while(<UNSURE3_Inversions>){
		$m21++;
		chomp;
		@fields = ();
		@fields = split/\s+/, $_;
		if(scalar(@fields)>2){
			print COMBINE "$fields[0]\t$fields[1]\t$fields[2]\tInversion\tUNSURE";
			foreach($NUM_UNSURE1..($NUM_UNSURE2-1)){
				print COMBINE "\tN.A";
			}
			print COMBINE "\t3\t$fields[5]\tN.A\n";
		}
	}
}

if($N22 == 1){
	<UNSURE4_Inversions>;
	while(<UNSURE4_Inversions>){
		$m22++;
		chomp;
		@fields = ();
		@fields = split/\s+/, $_;
		print COMBINE "$fields[0]\t$fields[1]\t$fields[2]\tInversion\tUNSURE";
		foreach($NUM_UNSURE1..($NUM_UNSURE2-1)){
			print COMBINE "\tN.A";
		}
		print COMBINE "\t4\t$fields[5]\tN.A\n";
	}
}

if($N23 == 1){
	<UNSURE5_Deletions>;
	while(<UNSURE5_Deletions>){
		$m23++;
		chomp;
		@fields = ();
		@fields = split/\s+/, $_;
		print COMBINE "$fields[0]\t$fields[1]\t$fields[2]\tDeletion\tUNSURE";
		foreach($NUM_singleRT1..($UNSURE-1)){
			if($_<$NUM_singleRT2){
				print COMBINE "\t$fields[$_+1]";
			}else{
				print COMBINE "\tN.A";
			}
		}
		print COMBINE "\tN.A\t5\t$fields[6]\tN.A\n";
	}
}

if($N24 == 1){
	<UNSURE6_Deletions>;
	while(<UNSURE6_Deletions>){
		$m24++;
		chomp;
		@fields = ();
		@fields = split/\s+/, $_;
		print COMBINE "$fields[0]\t$fields[1]\t$fields[2]\tDeletion\tUNSURE";
			foreach($NUM_MultiRT0..$TOTAL){
			if($_<$NUM_MultiRT1){
				print COMBINE "\tN.A";
			}elsif($_<$NUM_MultiRT2){
				print COMBINE "\t$fields[$_-5]";
			}
		}
		print COMBINE "\t6\t$fields[6]\tN.A\n";			
	}
}

if($N25 == 1){
	<UNSURE7_Deletions>;
	while(<UNSURE7_Deletions>){
		$m25++;
		chomp;
		@fields = ();
		@fields = split/\s+/, $_;
		print COMBINE "$fields[0]\t$fields[1]\t$fields[2]\tDeletion\tUNSURE";
		foreach($coverage..$TOTAL){
			if($_<$nmhNo){
				print COMBINE "\tN.A";
			}elsif($_<$StrandElement){
				print COMBINE "\t$fields[$_-7]";
			}
		}
		print COMBINE "\tN.A\t$fields[9]\tN.A\tN.A\t$fields[10]\t7\t$fields[11]\tN.A\n";	
	}
}

if($N26 == 1){
	<UNSURE5_Insertions>;
	while(<UNSURE5_Insertions>){
		$m26++;
		chomp;
		@fields = ();
		@fields = split/\s+/, $_;
		print COMBINE "$fields[0]\t$fields[1]\t$fields[2]\tInsertion\tUNSURE";
		foreach($NUM_singleRT1..($UNSURE-1)){
			if($_<$NUM_singleRT2){
				print COMBINE "\t$fields[$_+1]";
			}else{
				print COMBINE "\tN.A";
			}
		}
		print COMBINE "\tN.A\t5\t$fields[7]\t$fields[6]\n";
	}
}

if($N27 == 1){
	<UNSURE6_Insertions>;
	while(<UNSURE6_Insertions>){
		$m27++;
		chomp;
		@fields = ();
		@fields = split/\s+/, $_;
		print COMBINE "$fields[0]\t$fields[1]\t$fields[2]\tInsertion\tUNSURE";
			foreach($NUM_MultiRT0..$TOTAL){
			if($_<$NUM_MultiRT1){
				print COMBINE "\tN.A";
			}elsif($_<$NUM_MultiRT2){
				print COMBINE "\t$fields[$_-5]";
			}
		}
		print COMBINE "\t6\t$fields[7]\t$fields[6]\n";			
	}
}

if($N28 == 1){
	<UNSURE7_Insertions>;
	while(<UNSURE7_Insertions>){
		$m28++;
		chomp;
		@fields = ();
		@fields = split/\s+/, $_;
		print COMBINE "$fields[0]\t$fields[1]\t$fields[2]\tInsertion\tUNSURE";
		foreach($coverage..$TOTAL){
			if($_<$nmhNo){
				print COMBINE "\tN.A";
			}elsif($_<$StrandElement){
				print COMBINE "\t$fields[$_-7]";
			}
		}
		print COMBINE "\tN.A\t$fields[9]\tN.A\tN.A\t$fields[10]\t7\t$fields[12]\t$fields[11]\n";	
	}
}

if($N29 == 1){
	<VNTR_Deletions>;
	while(<VNTR_Deletions>){
		$m29++;
		chomp;
		@fields = ();
		@fields = split/\s+/, $_;
		print COMBINE "$fields[0]\t$fields[1]\t$fields[2]\tDeletion\tVNTR";
		foreach($coverage..$TOTAL){
			print COMBINE "\tN.A";
		}
		print COMBINE "\t$fields[9]\tN.A\n";	
	}
}

if($N30 == 1){
	<VNTR_Insertions>;
	while(<VNTR_Insertions>){
		$m30++;
		chomp;
		@fields = ();
		@fields = split/\s+/, $_;
		print COMBINE "$fields[0]\t$fields[1]\t$fields[2]\tInsertion\tVNTR";
		foreach($coverage..$TOTAL){
			print COMBINE "\tN.A";
		}
		print COMBINE "\t$fields[10]\t$fields[3]\n";	
	}
}

#print STATISTIC "Notes On Unsure Cases:\n";
#print STATISTIC "Unsure1:Smaller_Than_200bp;\n";
#print STATISTIC "Unsure2:500bp_Proximity_To_Either_Telemere, Flanking sequence not available;\n";
#print STATISTIC "Unsure3:Large_Offset_Or_DoNOt_Span_BP when blasting two flanking sequences;\n";
#print STATISTIC "Unsure4:Potential NAHR, LongHomology_Or_MediumLengthWithHighPercentage, EITHER(1)30bp< homology length <50bp, >=85%, no large offset, spaning BP; OR (2) homology length >= 50bp, <85%, no large offset, spaning BP;\n";
#print STATISTIC "Unsure5:Single_Transposon, No TSD, OR No PolyA tail on the correct strand;\n";
#print STATISTIC "Unsure6:Multiple_Transposon, No PolyA tail on the correct strand;\n";
#print STATISTIC "Unsure7:Potential_Processed_Pseudogene, Both TSD AND PolyA tail;\n\n\n";

print STATISTIC "SV\tVNTR\tNHR\tNAHR\tSTEI\tMTEI\tUnsure1\tUnsure2\tUnsure3\tUnsure4\tUnsure5\tUnsure6\tUnsure7\tUnsureTotal\tTotal\n";

print STATISTIC "Deletions\t$m29\t$m3\t$m1\t$m2\t$m4\t$m11\t$m12\t",($m13/2),"\t$m14\t$m23\t$m24\t$m25\t",($m11+$m12+$m13/2+$m14+$m23+$m24+$m25),"\t",($m3+$m1+$m2+$m4+$m11+$m12+$m13/2+$m14+$m23+$m24+$m25+$m29),"\n";

print STATISTIC "Insertions\t$m30\t$m7\t$m5\t$m6\t$m8\t$m15\t$m16\t",($m17/2),"\t$m18\t$m26\t$m27\t$m28\t",($m15+$m16+$m17/2+$m18+$m26+$m27+$m28),"\t",($m7+$m5+$m6+$m8+$m15+$m16+$m17/2+$m18+$m26+$m27+$m28+$m30),"\n";

print STATISTIC "Inversions\t0\t$m10\t$m9\t0\t0\t$m19\t$m20\t",($m21/2),"\t$m22\t0\t0\t0\t",($m19+$m20+$m21/2+$m22),"\t",($m10+$m9+$m19+$m20+$m21/2+$m22),"\n";



close NAHR_Deletions;
close SINGLE_RT_Deletions;
close NHEJ_Deletions;
close MULTI_RT_Deletions;
close NAHR_Insertions;
close SINGLE_RT_Insertions;
close NHEJ_Insertions;
close MULTI_RT_Insertions;
close NAHR_Inversions;
close NHEJ_Inversions;
close COMBINE;
close STATISTIC;
close UNSURE1_Deletions;
close UNSURE2_Deletions;
close UNSURE3_Deletions;
close UNSURE4_Deletions;
close UNSURE1_Insertions;
close UNSURE2_Insertions;
close UNSURE3_Insertions;
close UNSURE4_Insertions;
close UNSURE1_Inversions;
close UNSURE2_Inversions;
close UNSURE3_Inversions;
close UNSURE4_Inversions;
close UNSURE5_Deletions;
close UNSURE6_Deletions;
close UNSURE7_Deletions;
close UNSURE5_Insertions;
close UNSURE6_Insertions;
close UNSURE7_Insertions;
close VNTR_Deletions;
close VNTR_Insertions;
