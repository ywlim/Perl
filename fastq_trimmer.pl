#!/usr/bin/perl
# fastq_trimmer.pl by Yoong Wearn Lim
# to detect and trim of illumina adapter sequence from a fastq file
# illumina adapters:
# >Solexa_forward_contam
# AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAA
# >Solexa_reverse_contam
# AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTGAAAAA
# usually for single end reads we should only get forward contamination
# at the 3' end
# but strangely we have some reverse contamination
# that starts at 5' end
# that most trimmer program won't trim
# so my script checks for that too
# only keep reads >= 15 bp after trimming

use strict; use warnings;

die "usage: fastq_trimmer <fastq>\n" unless @ARGV == 1;

my ($stat) = "trim_stats.txt";

if ($ARGV[0] =~ /\.gz$/)
{
	open(IN, "gunzip -c $ARGV[0] |") or die "can't open pipe to $ARGV[0]";
}
else
{
	open (IN, $ARGV[0]) or die "error opening $ARGV[0]\n";
}
open (STAT, ">>$stat") or die "error opening $stat\n";

my $count1; my $count2; my $count3; my $count4; my $count5; my $count6;
my $total;
my $remain1; my $remain2; my $remain3; my $remain4; my $remain5;
my $totremain1; my $totremain2; my $totremain3; my $totremain4; my $totremain5;
my $kept; my $keptlength;

while (my $text = <IN>)
{
	my $ID = $text;
	chomp $ID;

	my $line = <IN>;	# the read
	chomp $line;

	<IN>;	# third line, don't care
	my $quality = <IN>;	# forth line, quality
	chomp $quality;

	if ($line =~ m/AGATCGGAAGAGCGGT.*/)
	{
		$count1++;
		my ($match) = $line =~ m/(AGATCGGAAGAGCGGT.*)/;
		$remain1 = length($line) - length($match);
		$totremain1 += $remain1;

		if ($remain1 >= 15)	# only keep the read if it's >= 15 bp long
		{
			$line =~ s/AGATCGGAAGAGCGGT.*//;
			$quality = substr($quality, 0, $remain1);

			print "$ID\n";
			print "$line\n";
			print "+\n";
			print "$quality\n";
			$kept++;
			$keptlength += length($line);
		}
	}

	elsif ($line =~ m/AGATCGGAAGAGCGG$/)
	{
		$count2++;
		my ($match) = $line =~ m/(AGATCGGAAGAGCGG$)/;
		$remain2 = length($line) - length($match);
		$totremain2 += $remain2;

		if ($remain2 >= 15)	# only keep the read if it's >= 15 bp long
		{
			$line =~ s/AGATCGGAAGAGCGG$//;
			$quality = substr($quality, 0, $remain2);

			print "$ID\n";
			print "$line\n";
			print "+\n";
			print "$quality\n";
			$kept++;
			$keptlength += length($line);
		}
	}

	elsif ($line =~ m/AGATCGGAAGAGCGTC.*/)
	{
		$count3++;
		my ($match) = $line =~ m/(AGATCGGAAGAGCGTC.*)/;
		$remain3 = length($line) - length($match);
		$totremain3 += $remain3;

		if ($remain3 >= 15)	# only keep the read if it's >= 15 bp long
		{
			$line =~ s/AGATCGGAAGAGCGTC.*//;
			$quality = substr($quality, 0, $remain3);

			print "$ID\n";
			print "$line\n";
			print "+\n";
			print "$quality\n";
			$kept++;
			$keptlength += length($line);
		}
	}

	elsif ($line =~ m/AGATCGGAAGAGCGT$/)
	{
		$count4++;
		my ($match) = $line =~ m/(AGATCGGAAGAGCGT$)/;
		$remain4 = length($line) - length($match);
		$totremain4 += $remain4;

		if ($remain4 >= 15)	# only keep the read if it's >= 15 bp long
		{
			$line =~ s/AGATCGGAAGAGCGT$//;
			$quality = substr($quality, 0, $remain4);

			print "$ID\n";
			print "$line\n";
			print "+\n";
			print "$quality\n";
			$kept++;
			$keptlength += length($line);
		}
	}


	# bizzarre case where the we have reverse contamination
	# although it's single end read
	# after checking for quality
	# i realize that the many reads start with this
	# and the entire read is junk
	# with bad quality

	elsif ($line =~ m/GATCGGAAGAGCGGT.*/)
	{
		$count5++;
		my ($match) = $line =~ m/(GATCGGAAGAGCGGT.*)/;
		$remain5 = length($line) - length($match);
		$totremain5 += $remain5;

		if ($remain5 >= 15)	# only keep the read if it's >= 15 bp long
		{
			$line =~ s/GATCGGAAGAGCGGT.*//;
			$quality = substr($quality, 0, $remain5);
			print "$ID\n";
			print "$line\n";
			print "+\n";
			print "$quality\n";
			$kept++;
			$keptlength += length($line);
		}

	}

	# no adapter contamination
	else
	{
		$count6++;

		print "$ID\n";
		print "$line\n";
		print "+\n";
		print "$quality\n";
		$kept++;
		$keptlength += length($line);
	}

	$total++;
}

close IN;

print STAT "File: $ARGV[0]\n";
print STAT "Total reads: $total\n";
print STAT "Reads kept: $kept\n";
print STAT "Average kept read length:	", $keptlength / $kept, "\n\n";
print STAT "Read adapter contamination: \n\n";
print STAT "No contamination:	", $count6 / $total * 100, " %\n\n";
print STAT "AGATCGGAAGAGCGGT.*	", $count1 / $total * 100, " %\n";
print STAT "Average remaining read length:	", $totremain1 / $count1, "\n\n";
print STAT "AGATCGGAAGAGCGG\$	", $count2 / $total * 100, " %\n";
print STAT "Average remaining read length:	", $totremain2 / $count2, "\n\n";
print STAT "AGATCGGAAGAGCGTC.*	", $count3 / $total * 100, " %\n";
print STAT "Average remaining read length:	", $totremain3 / $count3, "\n\n";
print STAT "AGATCGGAAGAGCGT\$	", $count4 / $total * 100, " %\n";
print STAT "Average remaining read length:	", $totremain4 / $count4, "\n\n";
print STAT "GATCGGAAGAGCGGT.*	", $count5 / $total * 100, " %\n";
print STAT "Average remaining read length:	", $totremain5 / $count5, "\n\n";
__END__

### NGB00362.1:1-61 Illumina Paired End PCR Primer 2.0
perl -p -i -e 's/AGATCGGAAGAGCGGT.*//' z-trim-test.trim
perl -p -i -e 's/AGATCGGAAGAGCGG$//' z-trim-test.trim

### NGB00361.1:1-92 Illumina PCR Primer
perl -p -i -e 's/AGATCGGAAGAGCGTC.*//' z-trim-test.trim
perl -p -i -e 's/AGATCGGAAGAGCGT$//' z-trim-test.trim

### Common region for Illumina PCR Primer and Paired End PCR Primer 2.0
perl -p -i -e 's/AGATCGGAAGAGCG$//' z-trim-test.trim
perl -p -i -e 's/AGATCGGAAGAGC$//' z-trim-test.trim
perl -p -i -e 's/AGATCGGAAGAG$//' z-trim-test.trim
perl -p -i -e 's/AGATCGGAAGA$//' z-trim-test.trim
perl -p -i -e 's/AGATCGGAAG$//' z-trim-test.trim
perl -p -i -e 's/AGATCGGAA$//' z-trim-test.trim
perl -p -i -e 's/AGATCGGA$//' z-trim-test.trim

Removing of homopolymer tails:

perl -p -i -e 's/^A{8,}//' z-trim-test.trim
perl -p -i -e 's/^T{8,}//' z-trim-test.trim
perl -p -i -e 's/^C{8,}//' z-trim-test.trim
perl -p -i -e 's/^G{8,}//' z-trim-test.trim

perl -p -i -e 's/A{8,}$//' z-trim-test.trim
perl -p -i -e 's/T{8,}$//' z-trim-test.trim
perl -p -i -e 's/G{8,}$//' z-trim-test.trim
perl -p -i -e 's/C{8,}$//' z-trim-test.trim
