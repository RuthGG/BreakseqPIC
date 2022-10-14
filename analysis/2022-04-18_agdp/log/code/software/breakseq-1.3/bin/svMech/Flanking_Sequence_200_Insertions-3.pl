#!/usr/bin/perl -w
use strict;
use POSIX qw(ceil floor);
#This program grabs -$NAHRwin/2~+$NAHRwin/2bp flanking sequences from left and right breakpoints
unless ($ARGV[3]) {  
	print STDERR "Exit!\n";
	exit;
}


my $temp;
my @fields;
my $RPflank = 500;
my $NAHRwin = $ARGV[1];  #1000
my $leftflank;
my $rightflank;
my $seq;
my $out_dir = $ARGV[3];
my $SV_Size; #Feb 2010

open COORDINATES, "<$ARGV[0]"; #output from non-VNTR-50
open FLANKINGSEQ, ">$out_dir/$ARGV[2].FlankingSequenceInsertions.txt";

print FLANKINGSEQ "Chr\tStart\tEnd\tLeftFlankingSequence\tRightFlankingSequence\tInserted_seq\tID\n";

<COORDINATES>;
while($temp = <COORDINATES>){
	chomp($temp);
	@fields = ();
	@fields = split/\s+/, $temp;
	open IN, "<$out_dir/RepeatMasker_Output/Insertions/Qchr$fields[0]_$fields[1]-$fields[2].$fields[4].txt" or die;
	<IN>;
	$seq = <IN>;
	chomp($seq);
	close IN;
	$SV_Size = (length($seq)-2*$RPflank);  #Feb 2010
 		$leftflank = substr($seq,($RPflank-$NAHRwin/2),$NAHRwin);
 		$rightflank = substr($seq,(length($fields[3])+$RPflank-$NAHRwin/2),$NAHRwin);
	print FLANKINGSEQ "$fields[0]\t$fields[1]\t$fields[2]\t$leftflank\t$rightflank\t$fields[3]\t$fields[4]\n";
}
close COORDINATES;

open  COORDINATES, "<$out_dir/Output/Insertions/$ARGV[2].LargeInsertion.txt";  #Large file
<COORDINATES>;
while($temp = <COORDINATES>){
	chomp($temp);
	@fields = ();
	@fields = split/\s+/, $temp;
	open IN, "<$out_dir/RepeatMasker_Output/Insertions/Qchr$fields[0]_$fields[1]-$fields[2].$fields[4].txt" or die;
	<IN>;
	$seq = <IN>;
	chomp($seq);
	close IN;
	$SV_Size = (length($seq)-2*$RPflank);  #Feb 2010
 		$leftflank = substr($seq,($RPflank-$NAHRwin/2),$NAHRwin);
 		$rightflank = substr($seq,(length($fields[3])+$RPflank-$NAHRwin/2),$NAHRwin);
	print FLANKINGSEQ "$fields[0]\t$fields[1]\t$fields[2]\t$leftflank\t$rightflank\t$fields[3]\t$fields[4]\n";
}
close COORDINATES;
close FLANKINGSEQ;
