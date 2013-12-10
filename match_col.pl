#!/usr/bin/perl
# match_col.pl by Yoong Wearn Lim
# compare 2 files
# based on an identifier in file 1 and file 2
# print file 1 and file 2
# user should use awk to select for desired columns in final output file

use strict; use warnings;

die "usage: match_col.pl <file 1> <file 1 identifier col> <file 2> <file 2 identifier col>\n" unless @ARGV == 4;

open (IN1, $ARGV[0]) or die "error opening $ARGV[0]\n";
open (IN2, $ARGV[2]) or die "error opening $ARGV[2]\n";

my @col; my @col2;
my %stuff;

# minus one to identifiers because arrays are 0-based
my ($id1) = $ARGV[1] - 1;
my ($id2) = $ARGV[3] - 1;

# file 2
while (my $line2 = <IN2>)
{
	chomp $line2;
	@col2 = split("\t", $line2);
	$stuff{$col2[$id2]} = $line2;
#	print "$id2	$col2[$id2]	$stuff{$col2[$id2]}\n";
}

while (my $line = <IN1>)
{
	chomp $line;
	@col = split("\t", $line);
	if (exists $stuff{$col[$id1]})
	{
		print "$line	$stuff{$col[$id1]}\n";
	}
}

close IN1;
close IN2;
