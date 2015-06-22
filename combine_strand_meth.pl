#!/usr/bin/perl
# combine_strand_meth.pl by Yoong Wearn Lim
# this script combines methylation info on both plus and minus strands
# by processing output files of Bismark's genome_methylation_bismark2bedGraph_v4.pl

# Example of how to go from a fastq file to methylation info:

# bismark /data/Bismark_indexes/hg19 AGS1.fastq
# methylation_extractor --single-end AGS1.fastq_bismark.sam
# genome_methylation_bismark2bedGraph_v4.pl --counts --split_by_chromosome CpG_OT_AGS1.fastq_bismark.sam.txt > AGS1_CpG_OT_meth.txt
# genome_methylation_bismark2bedGraph_v4.pl --counts --split_by_chromosome CpG_OB_AGS1.fastq_bismark.sam.txt > AGS1_CpG_OB_meth.txt
# awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5+$6 "\t+"}' AGS1_CpG_OT_meth.txt > temp1
# awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5+$6 "\t-"}' AGS1_CpG_OB_meth.txt > temp2
# cat temp1 temp2 | sort -k1,1 -k2,2n > temp3
# combine_strand_meth_v2.pl temp3 > AGS1_meth.bed

use strict; use warnings;

die "usage: combind_strand_meth.pl <bismark_OT_OB_file>\n" unless @ARGV == 1;
open (IN, $ARGV[0]) or die "error opening $ARGV[0]\n";

my ($oldchr, $oldstart, $oldend, $oldmeth, $oldcount, $oldstrand);

# initialize previous strand to be -
# so that if the file starts with -
# it won't compare with previous strand
my $previous = "-";

while (my $line = <IN>)
{
    chomp $line;
    my ($chr, $start, $end, $meth, $count, $strand) = split("\t", $line);

    # if see plus, don't print yet
    # store it in memory
    if ($strand eq "+")
    {
        if ($previous eq "+")
        {
            # print previously saved plus strand
            print "$oldchr\t$oldstart\t$oldend\t$oldmeth\t$oldcount\n";
        }

        ($oldchr, $oldstart, $oldend, $oldmeth, $oldcount, $oldstrand) = ($chr, $start, $end, $meth, $count, $strand);
        $previous = "+";
    }

    elsif ($strand eq "-")
    {
        if ($previous eq "+")
        {
            # minus strand is complimentary of previous plus strand
            # calculate new meth
            # print combined info (use plus strand coordinates)
            if ($start == ($oldstart + 1))
            {
                my $newcount = $count + $oldcount;
                my $newmeth = (($meth * $count) + ($oldmeth * $oldcount)) / $newcount;
                print "$oldchr\t$oldstart\t$oldend\t$newmeth\t$newcount\n";
            }

            # minus strand is not complimentary of previous plus strand
            else
            {
                # print previously saved plus strand
                print "$oldchr\t$oldstart\t$oldend\t$oldmeth\t$oldcount\n";
                # print current minus strand
                print "$chr\t$start\t$end\t$meth\t$count\n";
            }
        }
        elsif ($previous eq "-")
        {
            print "$chr\t$start\t$end\t$meth\t$count\n";
        }
        $previous = "-";
    }
}

close IN;
