#!/usr/bin/perl


use warnings;
use Getopt::Std;

my %opts = ('o' => '-'); # direct to standard output by default.
getopts('s:i:o:h', \%opts);

sub usage {
	print "result_table.2.pl -i <input_prefix> -o <output> -s <sample_id>\n";
	print "The input files are assumed to be <input_prefix>.xun and <input_prefix>.uni.\n";
	print "This version will only analyse one individual at a time.\n";
	exit;
}

if (exists($opts{'h'}) || !exists($opts{'i'}) || !exists($opts{'s'})) {
	usage();
}


#############################
#### Open output file
#############################

open(FILE2, ">$opts{'o'}");

#############################
#### Print header
#############################

#print FILE2 "#Inversion\tAllele\tBreakpoint\tIndividual\tReads\n";
# I do not print header, to be able to sort the result later.

#################################################
#### Open files .uni and .xun for each individual
#################################################

if (-e $opts{'i'} . ".uni") {
	open(UNI, ($opts{'i'} . ".uni")) or die "I cannot open $opts{'i'}.uni\n";
	while (<UNI>) {
	        my @line = split /\t/, $_;
	        #if (($line[0] !~ /^ST.+REF$/) && ($line[0] =~ /^(ST|INV)(\d+)BP(1|2)(.+)$/)) {
	        if ($line[0] =~ /^([A-Z]{3})(\D*0*\d+)BP(0|1|2)(.+)$/) {
	                #print FILE2 $2, "\t", $1, "_", $4, "\tBP", $3, "\t", $opts{'s'}, "\t", $line[1], "\n";
	                print FILE2 $2, "\t", $1, "\tBP", $3, "\t", $opts{'s'}, "\t", $line[1], "\n";
	        }
	}
	close(UNI);
}

if (-e $opts{'i'} . ".xun") {
	open(XUN, ($opts{'i'} . ".xun")) or die "I cannot open $opts{'i'}.xun.\n";
	while (<XUN>) {
	        my @line = split /\t/, $_;
	        #if ($line[0] =~	/^ST(\d+)BP(1|2)REF$/) {
	        if ($line[0] =~ /^([A-Z]{3})(\D*0*\d+)BP(0|1|2)(.+)$/) {
	                #print FILE2 $1,	"\tST_REF\tBP", $2, "\t", $opts{'s'}, "\t", $line[1], "\n";
	                print FILE2 $2, "\t", $1, "\tBP", $3, "\t", $opts{'s'}, "\t", $line[1], "\n";
	        }
	        #if ($line[0] =~ /^([A-Z]{3})(\w+\d+)BP(1|2)(.+)$/) {
	        #    print FILE2 $2, "\t", $1, "\t", $4, "\tBP", $3, "\t", $opts{'s'}, "\t", $line[1], "\n";
	        #}
	}
	close(XUN);
}
close(FILE2);
