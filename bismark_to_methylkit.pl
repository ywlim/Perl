#!/usr/bin/perl
# bismark_to_methylkit.pl by YWL
# bismark file format has to be like this:
# chr1    10496   10496   98.5981308411215        422     6
# and the input file names has to be something like AGS2_CpG_OT_meth.txt and AGS2_CpG_OB_meth.txt
# methylkit format is:
# $chr\.$start\t$chr\t$start\t$strand_meth\t$covs\t$freqC\t$freqT

use strict; use warnings;

my ($input) = $ARGV[0];
die "usage: bismark_to_methylkit.pl <bismark file prefix>\n" unless @ARGV == 1;
my $top = "$ARGV[0]" . "_CpG_OT_meth.txt";
my $bottom = "$ARGV[0]" . "_CpG_OB_meth.txt";

open (my $in, "<", $top) or die "Cannot read from $top\n";

while (my $line = <$in>) 
{
	chomp($line);
	my ($chr, $start, $end, $junk, $C, $T) = split("\t", $line);
        my $freqC = sprintf '%.2f', $C / ($C + $T);      # Not Converted = Methylated
        my $freqT = sprintf '%.2f', $T / ($C + $T); # Converted = Not Methylated
	my $covs = $T + $C;
	print "$chr\.$start\t$chr\t$start\tF\t$covs\t$freqC\t$freqT\n";
}
close $in;

# lazy to use loop, just copy and paste code for bottom strand
open (my $in, "<", $bottom) or die "Cannot read from $bottom\n";

while (my $line = <$in>) 
{
        chomp($line);
        my ($chr, $start, $end, $junk, $C, $T) = split("\t", $line);
        my $freqC = sprintf '%.2f', $C / ($C + $T);      # Not Converted = Methylated
        my $freqT = sprintf '%.2f', $T / ($C + $T); # Converted = Not Methylated
        my $covs = $T + $C;
        print "$chr\.$start\t$chr\t$start\tR\t$covs\t$freqC\t$freqT\n";
}
close $in;


