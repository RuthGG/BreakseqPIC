#!/usr/bin/perl -w
use strict;
use Bio::Tools::GFF;
# Modified on Feb 19 to take gff as input and output gff to STDOUT
# Only look at features with gene_type eq "protein_coding" or IG_gene or "polymorphic_pseudogene"
# and transcrip_type eq "protein_coding"or IG_gene or "polymorphic_pseudogene"
# Genecode v3b features:
# CDS
# UTR
# exon
# gene
# start_codon
# stop_codon
# transcript

unless ($ARGV[2]) {  
	print STDERR "Args: InFile OutFile GeneAnnotationGFF \nExit!\n";
	exit;
}


my @all;
my @tag = ("gene_id","transcript_id","gene_type","transcript_type");
my @type = ("protein_coding","IG_C_gene","IG_D_gene","IG_J_gene","IG_V_gene", "polymorphic_pseudogene");
my $temp;
my $feature;
my $count = 1;
my $gffstr;
my $index;
my @feats = ("gene","exon","transcript");

my $SV_Class = 2;
my $Chr = 0;
my $Start = 3;
my $End = 4;
my $header_line;
my $line;
my @info;
my $i;
my $j;
my $key_gene;
my $value_gene;
my $key;
my $value;
my $ref;
my $ref_trx;
my $ref_exon;
my $Num_Gene_Affected;
my $Dist_Closest_Gene;
my $Single_Gene_Class;
my $Gene_ID;
my $dist1;
my $dist2;
my $CUT = 200000000;
my %trx;
my $Num;
my $all_keys;

my $Num_No_Gene=0;
my $Num_Multi_Gene=0;
my $Num_Sing_Gene_Removal=0;
my $Num_Sing_Gene_Intron=0;
my $Num_Sing_Gene_AllTrx=0;
my $Num_Sing_Gene_ParTrx=0;
my $Num_Sing_Gene_Other=0;
my $Num_All_SV=0;
my $m1;
my $m2;

#my $gffio = Bio::Tools::GFF->new(-file => "/Users/Jasmine/Human_Genome_Data/Gencode_v3b/gencode.v3b.annotation.NCBI36.gtf" , -gff_version => 2.5);
my $gffio = Bio::Tools::GFF->new(-file => "$ARGV[2]" , -gff_version => 2.5);
while($feature = $gffio->next_feature()) {
	last if($count>10000);
	next if(!($feature->has_tag($tag[0])));
	next if(($feature->primary_tag eq "CDS") || ($feature->primary_tag eq "UTR")|| ($feature->primary_tag eq "start_codon")|| ($feature->primary_tag eq "stop_codon"));    #Requires too much memory if included
	next if((!($feature->has_tag($tag[2])))||(!($feature->has_tag($tag[3]))));
	$m1 = 0;
	$m2 = 0;
	foreach $i (@type){
		$m1 = 1 if(sprintf("%s",$feature->get_tag_values($tag[2])) eq $i);
		$m2 = 1 if(sprintf("%s",$feature->get_tag_values($tag[3])) eq $i);
	}
	next if($m1*$m2==0);
	$temp = (split/chr/,($feature->seq_id))[1];
	$index = "";	
	foreach(1..22,"X","Y"){
		if($temp eq $_){
			if($temp eq "X"){
				$index = 23;
			}elsif($temp eq "Y"){
				$index = 24;
			}else{
				$index = $temp;
			}
		}
	}
	next if($index eq "");
	if(not exists $all[$index]->{sprintf("%s",$feature->get_tag_values($tag[0]))}->{$feature->primary_tag}){
		$all[$index]->{sprintf("%s",$feature->get_tag_values($tag[0]))}->{$feature->primary_tag}->[0] = $feature;
	}else{
		$temp = $all[$index]->{sprintf("%s",$feature->get_tag_values($tag[0]))}->{$feature->primary_tag};
		push(@$temp,$feature);
	}
	$count++;
}
$gffio->close();

open IN, "<$ARGV[0]" or die;
open OUT, ">$ARGV[1]" or die;
while($line = <IN>){
	$Num_Gene_Affected = 0;
	$Dist_Closest_Gene = $CUT;
	$Single_Gene_Class = "N.A";
	$Gene_ID = "N.A";

	chomp($line);
	@info = split/\s+/,$line;
	$index = (split/chr/,($info[$Chr]))[1];
	if(uc($index) eq "X"){
		$index = 23;
	}elsif(uc($index) eq "Y"){
		$index = 24;
	}
	$ref = $all[$index];
	if($info[$SV_Class] eq "Deletion"){
		while(($key_gene,$value_gene) = each %$ref){
			$dist1 = ($value_gene->{$feats[0]}->[0]->start-$info[$End]);
			$dist2 = ($info[$Start]-$value_gene->{$feats[0]}->[0]->end);
			if(($dist1>0)||($dist2>0)){
				if($Dist_Closest_Gene!=0){
					if($dist1>0){
						if((abs($Dist_Closest_Gene) >$dist1)){
							$Dist_Closest_Gene=$dist1;
							$Gene_ID = $key_gene;
						}
					}elsif($dist2>0){
						if((abs($Dist_Closest_Gene) >$dist2)){
							$Dist_Closest_Gene=((-1)*$dist2);
							$Gene_ID = $key_gene;
						}
					}
				}
			}else{
				$Num_Gene_Affected++;
				$Dist_Closest_Gene = 0;
				if($Num_Gene_Affected==1){
					$Gene_ID = $key_gene;
					$Single_Gene_Class = $key_gene;
				}else{
					$Gene_ID = $Gene_ID.",".$key_gene;
					$Single_Gene_Class = "N.A";
				}
			}
		}
		if(!($Single_Gene_Class eq "N.A")){
			$temp = $ref->{$Single_Gene_Class}->{$feats[0]}->[0];
			if(($temp->start>=$info[$Start]) &&($temp->end<=$info[$End])){
				$Single_Gene_Class = "Gene_Removal";
			}else{
				if(not exists $ref->{$Single_Gene_Class}->{$feats[2]}){
					$Single_Gene_Class = "Other";
				}else{
					$ref_trx = $ref->{$Single_Gene_Class}->{$feats[2]};
					%trx = ();
					foreach $i (0..(scalar(@$ref_trx)-1)){
						if($ref_trx->[$i]->has_tag($tag[1])){
							$trx{sprintf("%s",$ref_trx->[$i]->get_tag_values($tag[1]))} = 0;					
						}
					}
					$ref_exon = $ref->{$Single_Gene_Class}->{$feats[1]};
					foreach $i (0..(scalar(@$ref_exon)-1)){
						if(!(($ref_exon->[$i]->start>$info[$End])||($ref_exon->[$i]->end<$info[$Start]))){
							if($ref_exon->[$i]->has_tag($tag[1])){
								$trx{sprintf("%s",$ref_exon->[$i]->get_tag_values($tag[1]))}++;					
							}
						}
					}
					$Num = 0;
					while(($key,$value) = each %trx){
						if($value==0){
							$Num++;
						}
					}
					$all_keys = keys %trx;
					if($Num == $all_keys){
						$Single_Gene_Class = "Intron_Only";
					}elsif($Num == 0){
						$Single_Gene_Class = "All_Trx";
					}else{
						$Single_Gene_Class = "Partial_Trx";
					}
				}
			}
		}
	}elsif($info[$SV_Class] eq "Insertion"){
			while(($key_gene,$value_gene) = each %$ref){
			$dist1 = ($value_gene->{$feats[0]}->[0]->start-$info[$End]+1);
			$dist2 = ($info[$Start]+1-$value_gene->{$feats[0]}->[0]->end);
			if(($dist1>0)||($dist2>0)){
				if($Dist_Closest_Gene!=0){
					if($dist1>0){
						if((abs($Dist_Closest_Gene) >$dist1)){
							$Dist_Closest_Gene=$dist1;
							$Gene_ID = $key_gene;
						}
					}elsif($dist2>0){
						if((abs($Dist_Closest_Gene) >$dist2)){
							$Dist_Closest_Gene=((-1)*$dist2);
							$Gene_ID = $key_gene;
						}
					}
				}
			}else{
				$Num_Gene_Affected++;
				$Dist_Closest_Gene = 0;
				if($Num_Gene_Affected==1){
					$Gene_ID = $key_gene;
					$Single_Gene_Class = $key_gene;
				}else{
					$Gene_ID = $Gene_ID.",".$key_gene;
					$Single_Gene_Class = "N.A";
				}
			}
		}
		if(!($Single_Gene_Class eq "N.A")){
			$temp = $ref->{$Single_Gene_Class}->{$feats[0]}->[0];
			if(($temp->start>$info[$Start]) &&($temp->end<$info[$End])){
				$Single_Gene_Class = "Gene_Removal";
			}else{
				if(not exists $ref->{$Single_Gene_Class}->{$feats[2]}){
					$Single_Gene_Class = "Other";
				}else{
					$ref_trx = $ref->{$Single_Gene_Class}->{$feats[2]};
					%trx = ();
					foreach $i (0..(scalar(@$ref_trx)-1)){
						if($ref_trx->[$i]->has_tag($tag[1])){
							$trx{sprintf("%s",$ref_trx->[$i]->get_tag_values($tag[1]))} = 0;					
						}
					}
					$ref_exon = $ref->{$Single_Gene_Class}->{$feats[1]};
					foreach $i (0..(scalar(@$ref_exon)-1)){
						if(!(($ref_exon->[$i]->start>=$info[$End])||($ref_exon->[$i]->end<=$info[$Start]))){
							if($ref_exon->[$i]->has_tag($tag[1])){
								$trx{sprintf("%s",$ref_exon->[$i]->get_tag_values($tag[1]))}++;					
							}
						}
					}
					$Num = 0;
					while(($key,$value) = each %trx){
						if($value==0){
							$Num++;
						}
					}
					$all_keys = keys %trx;
					if($Num == $all_keys){
						$Single_Gene_Class = "Intron_Only";
					}elsif($Num == 0){
						$Single_Gene_Class = "All_Trx";
					}else{
						$Single_Gene_Class = "Partial_Trx";
					}
				}
			}
		}
	}
	if($Dist_Closest_Gene==$CUT){
		$Dist_Closest_Gene = "N.A";
	}
	print OUT "$line; Num_Gene \"$Num_Gene_Affected\"; Dist_Closet_Gene \"$Dist_Closest_Gene\"; ";
	if(!($Single_Gene_Class eq "N.A")){
		print OUT "Single_Gene_Type \"$Single_Gene_Class\"; ";
	}
	if($Dist_Closest_Gene eq "0"){
		print OUT "Gene \"$Gene_ID\"";
	}
	print OUT "\n";
	$Num_All_SV++;
	if($Num_Gene_Affected==0){
		$Num_No_Gene++;
	}elsif($Num_Gene_Affected>1){
		$Num_Multi_Gene++;
	}elsif($Num_Gene_Affected==1){
		if($Single_Gene_Class eq "Gene_Removal"){
			$Num_Sing_Gene_Removal++;
		}elsif($Single_Gene_Class eq "Intron_Only"){
			$Num_Sing_Gene_Intron++;
		}elsif($Single_Gene_Class eq "All_Trx"){
			$Num_Sing_Gene_AllTrx++;
		}elsif($Single_Gene_Class eq "Partial_Trx"){
			$Num_Sing_Gene_ParTrx++;
		}elsif($Single_Gene_Class eq "Other"){
			$Num_Sing_Gene_Other++;
		}
	}
}
close IN;
close OUT;


















