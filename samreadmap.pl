#!/usr/bin/perl
# samreadmap.pl by Yoong Wearn Lim
# a script to find the number of mapped reads in a sam file

use strict; use warnings;

die "usage: samreadmap.pl <sam file>\n" unless @ARGV == 1;
my ($filename) = $ARGV[0];

open(IN, $filename) or die "error reading $filename";

my $nomap = 0;
my $map = 0;
my $unknown = 0;
my $count = 0;
my $chr;
my $start;
my %exist;
my %exist2;
my $unique = 0;
my $notunique = 0;
my $read;
my $uniqueread = 0;
my $dupread = 0;

while (my $line = <IN>)
{	
	chomp $line;
	if ($line =~ m/^@/)
	{
	}
	
	elsif ($line =~ m/\S+\t\d+\tchr(\w*)/)
	{
		$map++;
		$count++;
	
		($chr, $start) = $line =~ m/\S+\t\d+\tchr(\w*)\t(\d+)/;
		
		if (!defined $exist{$chr}{$start})
		{
			$unique++;
		}
		
		else 
		{
			$notunique++;
		}
		
		$exist{$chr}{$start} = $start;
	}
	
	elsif ($line =~ m/\S+\t\d+\t\*/)
	{
		$nomap++;
		$count++;
		
		($read) = $line =~ m/\S+\t\d+\t\*\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t(\w+)\t/;
		
		if (!defined $exist2{$read})
		{
			$uniqueread++;
		}
		
		else 
		{
			$dupread++;
		}
		
		$exist2{$read} = $read;
		
	}
}

close IN;

my $percentmap = $map / $count * 100;
my $percentnomap = $nomap / $count * 100;
my $percentunique = $unique / $count *100;
my $percentdup = $notunique / $count *100;
my $percentunique2 = $uniqueread / $count *100;
my $percentdup2 = $dupread / $count *100;

print "File: $filename\n";
print "\ntotal reads:	$count\n\n";
print "total mapped reads:	$map	$percentmap%\n";
print "	unique mapped reads:	$unique	$percentunique%\n";
print "	duplicate mapped reads:	$notunique	$percentdup%\n\n";
print "total unmapped reads:	$nomap	$percentnomap%\n";
print "	unique unmapped reads:	$uniqueread	$percentunique2%\n";
print "	duplicate unmapped reads:	$dupread	$percentdup2%\n\n";
