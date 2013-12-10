#!/usr/bin/perl
# wig_to_metaplot.pl by Yoong Wearn Lim on 2/3/13
# input wig file, and some bed files, and generate metaplot for the bed region
# my attempt to rewrite Paul's avg_read_v4.pl


use strict; use warnings;
use Getopt::Std;	# first time using get opt! excited!
our ($opt_h, $opt_w, $opt_b, $opt_n);
getopts('hw:b:n:');	# w and b take arguments, thus the : after them

my $usage = "usage: wig_to_metaplot.pl -w wigfile.wig -b bedfile1_plus.bed,bedfile1_minus.bed,bedfile2_plus.bed,bedfile2_minus.bed... -n name1,name2\n";

if ($opt_h)
{
	die $usage;
}

if (!$opt_b)	{die "-b not given\n"}
if (!$opt_w)	{die "-w not given\n"}
if (!$opt_n)	{die "-n not given\n"}

my @bed = split(",", $opt_b);
my @name = split(",", $opt_n);
my $num_bed = @bed;
my $num_name = @name;

#print "number of bed entries: ", $num_bed, "\n";
#print "number of name entries: ", $num_name, "\n";

die "Numbers of entry of bed and name don't match! Check that number of bed files = number of names * 2\n" unless ($num_name == $num_bed / 2);

#print "	bed is @bed\n
#		name is @name\n
#		wig is $opt_w\n";

# check bed files before spending time processing wig file

for (my $f = 0; $f < @bed; $f++)	# loop each bed file
{
	open (BED, $bed[$f]) or die "can't open $bed[$f] bed file, make sure to use full path\n";
	close BED;
}

#### WIG FILE ############

open (WIG, $opt_w) or die "can't open $opt_w wig file\n";

my $chr; my $span;
my %val;

print "Processing wig file\n";
while (my $line = <WIG>)
{
	chomp $line;
	next if (($line !~ m/^variableStep/) and ($line !~ /^\d+/));
	die "sorry, fixedStep wig file not supported\n" if ($line =~ m/^fixedStep/);

	if ($line =~ m/^variableStep/)
	{
		($chr, $span) = $line =~ m/^variableStep\schrom=chr(\w+)\sspan=(\d+)/;
		# print "chr is $chr and span is $span\n";
		print "Now processing chr$chr...\n";
	}

	elsif ($line =~ m/\d+\t\d+/)
	{
		my ($position, $value) = $line =~ m/(\d+)\t(\S+)/;
		# print "$position	$value\n";
		for (my $i = $position; $i <= ($position + $span); $i++)
		{
			if (!defined $val{$chr}{$i})
			{
				$val{$chr}{$i} = 0;
			}
			$val{$chr}{$i} += $value;
			# print "$chr	$i	$val{$chr}{$i}\n";
		}
	}
}

print "Done processing wig file\n";
close WIG;

#### BED FILES #############

my @depth; my @count;
my $k;
for (my $f = 0; $f < @bed; $f++)	# loop each bed file
{
	open (BED, $bed[$f]) or die "can't open $bed[$f] bed file\n";

	print "Processing bed file: $bed[$f]\n";
	while (my $line = <BED>)
	{
		chomp $line;
		my ($chro, $start, $end) = $line =~ m/^chr(\w+)\t(\d+)\t(\d+)/;
		# print "$chro	$start	$end\n";

		$k = 1;	# position within the bed window
		for (my $j = $start; $j <= $end; $j++)
		{
			if (exists $val{$chro}{$j})
			{
				$depth[$f][$k] += $val{$chro}{$j};
				$count[$f][$k]++;
			#	print "chr$chro	$j	$val{$chro}{$j}	$k	$depth[$k]	$count[$k]\n";
			}
			# debug when no wig value at the bed position
			#else
			#{
			#	print "chr$chro	$j	novalue	$k	$depth[$k]	$count[$k]\n";
			#}

			$k++;
		}
	}
}

print "Done processing bed files\n";
close BED;

#### RESULTS ############
# print result (average depth at each position h)
# assume that the bed coordinates are centered (eg. +- 500 TSS)

open (OUT, ">metaplot.txt") or die "can't write to metaplot.txt\n";
print "Now printing result\n";
my $header = join("\t", @name);
print OUT "bp\t$header\n";

my $h_adjusted = 0 - ($k / 2);
for (my $h = 1; $h < $k; $h++)
{
	print OUT $h_adjusted + $h, "\t";
	for (my $f = 0; $f < @bed; $f+=2)	# loop each bed set (plus and minus)
	{
		#print "$depth[$f][$h]	$depth[$f+1][-$h]\t";

		# in case no value was found at that position
		if ((!defined $count[$f][$h]) and (!defined $count[$f+1][-$h]))
		{
			print OUT "NA\t";
		}
		# only plus strand data available
		elsif (!defined $count[$f+1][-$h])
		{
			my $avg_depth = $depth[$f][$h] / $count[$f][$h];
			print OUT "$avg_depth\t";
		}
		# only minus strand data available
		elsif (!defined $count[$f][$h])
		{
			my $avg_depth = $depth[$f+1][-$h] / $count[$f+1][-$h];
			print OUT "$avg_depth\t";
		}
		else
		{
			my $avg_depth = ($depth[$f][$h] + $depth[$f+1][-$h]) / ($count[$f][$h] + $count[$f+1][-$h]);	# get average for plus and minus strand, flipping minus strand coordinate with minus array
			print OUT "$avg_depth\t";
		}
	}
	print OUT "\n";
}

close OUT;

############ make R script for graphing result ##############

open (R, ">metaplot.R") or die "can't open metaplot.R\n";
print R "library(ggplot2)\n";
print R "library(reshape)\n";
print R "pdf(file=\"metaplot.pdf\", family=\"Helvetica\", width=12, height=8)\n";
print R "plot<-read.table(\"metaplot.txt\", header=T)\n";
print R "plot.melt <- melt(plot[,c('bp', ";

for (my $w = 0; $w < @name; $w++)
{
	print R "'$name[$w]'";
	print R ", " unless ($w == $num_name - 1);
}

print R ")], id.vars=1)\n";
print R "ggplot(plot.melt, aes(x=bp, y=value, colour=variable, group=variable)) + geom_smooth() + theme_bw() + opts(title=\"$opt_w\", panel.grid.minor=theme_blank()) + scale_colour_brewer(palette=\"Set1\", name=\"Bed\") + ylim(0,100)\n";

close R;

################ run that R script! ##############

`R --vanilla < metaplot.R`
