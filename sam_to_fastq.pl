#!/usr/bin/perl
# sam_to_fastq.pl by Yoong Wearn Lim
# convert sam file to fastq file 

use strict; use warnings;

die "usage: sam_to_fastq.pl <sam file>\n" unless @ARGV == 1;
open (IN, $ARGV[0]) or die "error opening $ARGV[0]\n";

while (my $line = <IN>)
{
	chomp $line;
	next if ($line =~ m/^@/);	# skip header
	my @stuff = split("\t", $line);
	my ($id) = $stuff[0];
	my ($read) = $stuff[9];
	my ($quality) = $stuff[10];
	print "\@$id\n$read\n+\n$quality\n";
}

close IN;

