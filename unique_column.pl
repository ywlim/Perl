#!/usr/bin/perl
# unique_column.pl by Yoong Wearn Lim
# to clean up entries that are duplicates in a specific column in a tab deliminated bed file

use strict; use warnings;

die "usage: unique_column.pl <file> <column>\n" unless @ARGV == 2;
my ($file) = $ARGV[0];
my $k = $ARGV[1] - 1;	# minus 1 because array is zero based
open(IN, $file) or die "error opening $file\n";
my %seen;
while (my $line = <IN>)
{
        chomp $line;
        my (@column) = split("\t", $line);
	my $value = $column[$k];
        if (!exists $seen{$value})
        {
                print "$line\n";
        }
        $seen{$value} = 1;
}
close IN;

