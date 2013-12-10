#!/usr/bin/perl
# extract_wig_chr.pl by Yoong Wearn Lim
# extract certain chr from wig file

use strict; use warnings;

die "usage: extract_wig_chr.pl <wig file> <chr number>\n" unless @ARGV == 2;
open (IN, $ARGV[0]) or die "error opening $ARGV[0]\n";

my $chr_wanted = $ARGV[1];
my $trigger = 0;

while (my $line = <IN>)
{
	chomp $line;
	if ($line =~ m/^track/)	{print "$line\n"}	# print header
	if ($line =~ m/variableStep\schrom=chr\w+\sspan/)
	{
		my ($chr) = $line =~ m/variableStep\schrom=chr(\w+)\sspan/;
		if ($chr eq $chr_wanted)
		{
			print "$line\n";
			$trigger = 1;
		}
		else
		{
			$trigger = 0;
		}
	}

	elsif ($line =~ m/^\d+/)
	{
		print "$line\n" if ($trigger == 1);
	}
}

close IN;