#!/usr/bin/perl
# bed_to_wig_var.pl by Yoong Wearn Lim
# convert bed format to wig variable step format
# my attempt to improve the weird BedToWig.sh that always assumes span=12 or something

use strict; use warnings;

die "usage: bed_to_wig_var.pl <bed file> \n" unless @ARGV == 1;
open (IN, $ARGV[0]) or die "error opening $ARGV[0]\n";

my $oldspan;
my $oldchr = "";
while (my $line = <IN>)
{
	chomp $line;
	my ($chr, $start, $end, $value) = split("\t", $line);
	my $span = $end - $start;
	print "variableStep chrom=$chr span=$span\n" unless (($chr eq $oldchr) && ($span == $oldspan));
	print "$start	$value\n";
	$oldspan = $span;
	$oldchr = $chr;
}

close IN;
