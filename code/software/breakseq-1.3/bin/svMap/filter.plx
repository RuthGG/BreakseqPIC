#!/usr/bin/perl
#
# This is a substitution of the script filter.py, to take advantage
# of bowtie2's ability of mapping longer reads. Because the queries
# downloaded include reads of many different lengths, it is no longer
# acceptable to set their length a priori, but each read needs to be
# checked for its coverage of the breakpoint individually.
#
# I respect the original usage scheme:
#
# perl $BIN/filter.plx $LIB_LEN $SEED_LEN $MIN_COVER $OUTPUT.ini.sam > $OUTPUT.fil.sam

use strict;
use warnings;

my $lib_len = $ARGV[0];
my $seed_len = $ARGV[1];
my $min_cover = $ARGV[2];

open(SAM, $ARGV[3]) || die "Failed to open $ARGV[3].\n";
while (<SAM>) {
	my @sam = split /\t/, $_;
	$seed_len = length($sam[9]);
	if (($sam[3] >= $lib_len / 2 + $min_cover - $seed_len + 1) && ($sam[3] <= $lib_len / 2 - $min_cover + 1)) {
		print $_;
	}
}

close SAM;
