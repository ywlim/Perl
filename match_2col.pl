#!/usr/bin/perl
# match_2col.pl by Yoong Wearn Lim
# compare 2 files
# based on two identifiers in file 1 and file 2
# print file 1 and file 2
# user should use awk to select for desired columns in final output file

use strict; use warnings;

die "usage: match_2col.pl <file 1> <file 1 identifier col 1> <file 1 identifier col 2> <file 2> <file 2 identifier col 1> <file 1 identifier col 2>\n" unless @ARGV == 6;

open (IN1, $ARGV[0]) or die "error opening $ARGV[0]\n";
open (IN2, $ARGV[3]) or die "error opening $ARGV[3]\n";

my @col; my @col2;
my %stuff;

# minus one to identifiers because arrays are 0-based
my ($id1) = $ARGV[1] - 1;
my ($id2) = $ARGV[4] - 1;
my ($id1b) = $ARGV[2] - 1;
my ($id2b) = $ARGV[5] - 1;

# file 2
while (my $line2 = <IN2>)
{
	chomp $line2;
	@col2 = split("\t", $line2);
	$stuff{$col2[$id2]}{$col2[$id2b]} = $line2;
#	print "$id2	$col2[$id2]	$stuff{$col2[$id2]}\n";
}

while (my $line = <IN1>)
{
	chomp $line;
	@col = split("\t", $line);
	if (exists $stuff{$col[$id1]}{$col[$id1b]})
	{
		print "$line	$stuff{$col[$id1]}{$col[$id1b]}\n";
	}
}

close IN1;
close IN2;
