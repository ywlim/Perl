#!/usr/bin/perl
# unique_bg.pl by Yoong Wearn Lim
# to clean up entries that are duplicates in a bedgraph file

use strict; use warnings;

my %seen;
while (my $line = <>)
{
	chomp $line;
	my ($chr, $start, $end) = $line =~ m/(\w+)\t(\d+)\t(\d+)/;
	if (!exists $seen{$chr}{$start})
	{
		print "$line\n";
	}
	$seen{$chr}{$start} = 1;
}
