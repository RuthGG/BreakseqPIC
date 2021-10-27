#!/usr/bin/perl -w
use strict;
#this program combines VNTR with coverage >=50% ($CUT). 
unless ($ARGV[2]) {  
	print STDERR "Exit!\n";
	exit;
}

my $line;
my @fields;
my $CUT = $ARGV[0];
my $out_dir = $ARGV[2];
my $ID =9;
my $PCT = 8;

open Del50, ">$out_dir/Output/Deletions/$ARGV[1].VNTRDeletions-50.txt" or die;
open DelOther, ">$out_dir/$ARGV[1].OtherVNTRDeletions-50.txt" or die;
open Del, "<$out_dir/Output/Deletions/$ARGV[1].VNTRDeletions-2.txt" or die;  #$label
$line = <Del>;
print Del50 "$line";
print DelOther "chr\tstart\tend\tID\n";
while($line = <Del>){
	chomp($line);
	@fields = split/\s+/,$line;
	if($fields[$PCT]>=$CUT){
		print Del50 "$line\n";
	}else{
		print DelOther "$fields[0]\t$fields[1]\t$fields[2]\t$fields[$ID]\n";
	}
}
open Del2, "<$out_dir/Output/Deletions/$ARGV[1].NonVNTR_Deletions-2.txt" or die;  #$label
<Del2>;
while($line = <Del2>){
	chomp($line);
	@fields = split/\s+/,$line;
	print DelOther "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\n";
}
close Del2;
close Del;
close DelOther;
close Del50;


open Ins50, ">$out_dir/Output/Insertions/$ARGV[1].VNTRInsertions-50.txt" or die;
open InsOther, ">$out_dir/$ARGV[1].OtherVNTRInsertions-50.txt" or die;
open Ins, "<$out_dir/Output/Insertions/$ARGV[1].VNTRInsertions-2.txt" or die;  #$label
$line = <Ins>;
print Ins50 "$line";
print InsOther "chr\tstart\tend\tInserted_sequence\tID\n";
while($line = <Ins>){
	chomp($line);
	@fields = split/\s+/,$line;
	if($fields[$PCT+1]>=$CUT){
		print Ins50 "$line\n";
	}else{
		print InsOther "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t$fields[$ID+1]\n";
	}
}
open Ins2, "<$out_dir/Output/Insertions/$ARGV[1].NonVNTR_Insertions-2.txt" or die;  #$label
<Ins2>;
while($line = <Ins2>){
	chomp($line);
	@fields = split/\s+/,$line;
	print InsOther "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t$fields[4]\n";
}
close Ins2;
close Ins;
close InsOther;
close Ins50;


