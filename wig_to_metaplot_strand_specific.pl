#!/usr/bin/perl
# wig_to_metaplot_strand_specific.pl by Yoong Wearn Lim on 11/19/13
# modified from wig_to_metaplot.pl written on 2/3/13
# input wig file, and some bed files, and generate metaplot for the bed region


use strict; use warnings;
use Getopt::Std;	# first time using get opt! excited!
our ($opt_h, $opt_w, $opt_b, $opt_n);
getopts('hw:b:n:');	# w and b take arguments, thus the : after them

my $usage = "usage: wig_to_metaplot.pl -w plus.wig,minus.wig -b bed1_plus,bed1_minus,bed2_plus,bed2_minus... -n name1, name2...\n";

if ($opt_h)
{
	die $usage;
}

if (!$opt_b)	{die "-b not given\n"}
if (!$opt_w)	{die "-w not given\n"}
if (!$opt_n)	{die "-n not given\n"}

my @wig = split(",", $opt_w);
my @bed = split(",", $opt_b);
my @name = split(",", $opt_n);
#my $num_bed = @bed;
#my $num_name = @name;

#print "number of bed entries: ", $num_bed, "\n";
#print "number of name entries: ", $num_name, "\n";

#die "Numbers of entry of bed and name don't match! Check that number of bed files = number of names * 2\n" unless ($num_name == $num_bed / 2);

#print "	bed is @bed\n
#		name is @name\n
#		wig is $opt_w\n";

# check bed files before spending time processing wig file

#for (my $f = 0; $f < @bed; $f++)	# loop each bed file
#{
#	open (BED, $bed[$f]) or die "can't open $bed[$f] bed file, make sure to use full path\n";
#	close BED;
#}
my %val;
#### WIG FILE ##########
# q = 0 is plus wig file
# q = 1 is minus wig file
print "Processing wig files...\n";
for (my $q = 0; $q < 2; $q++)
{
	open (WIG, $wig[$q]) or die "can't open $wig[$q] wig file\n";

	my $chr; my $span;



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
				if ($q == 0)
				{				
					if (!defined $val{plus}{$chr}{$i})
					{
						$val{plus}{$chr}{$i} = 0;
					}
					$val{plus}{$chr}{$i} += $value;
				}
				elsif ($q == 1)
				{				
					if (!defined $val{minus}{$chr}{$i})
					{
						$val{minus}{$chr}{$i} = 0;
					}
					$val{minus}{$chr}{$i} += $value;
				}
			}
		}
	}

	print "Done processing $wig[$q]\n";
	close WIG;
}

#### BED FILES #############

my %depth; my %count;
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
			# get value from plus wig file			
			if (exists $val{plus}{$chro}{$j})
			{
				$depth{plus}[$f][$k] += $val{plus}{$chro}{$j};
				$count{plus}[$f][$k]++;
			}
			
			# get value from minus wig file
			if (exists $val{minus}{$chro}{$j})
			{
				$depth{minus}[$f][$k] += $val{minus}{$chro}{$j};
				$count{minus}[$f][$k]++;
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

# output file header
my @longname;
for (my $w = 0; $w < @name; $w++)
{ 
	push (@longname, ($name[$w] . "_sense"));	
	push (@longname, ($name[$w] . "_antisense"));	
}

my $num_longname = @longname;

open (OUT, ">metaplot.txt") or die "can't write to metaplot.txt\n";
print "Now printing result\n";
my $header = join("\t", @longname);
print OUT "bp\t$header\n";

my %avg_depth;

my $h_adjusted = 0 - ($k / 2);
for (my $h = 1; $h < $k; $h++)
{
	print OUT $h_adjusted + $h, "\t";
	for (my $f = 0; $f < @bed; $f+=2)	# loop each bed set (plus and minus)
	{
		# get average for plus and minus strand, flipping minus strand coordinate with minus array [-$h]
		# sense: plus ($f) gene with plus wig signal; minus ($f+1) gene with minus wig signal

		# sense
		if ((!defined $count{plus}[$f][$h]) and (!defined $count{minus}[$f+1][-$h]))	
		{
			$avg_depth{sense} = "NA";
		}
		elsif (!defined $count{plus}[$f][$h])
		{
			$avg_depth{sense} = $depth{minus}[$f+1][-$h] / $count{minus}[$f+1][-$h];
		}
		elsif (!defined $count{minus}[$f+1][-$h])
		{
			$avg_depth{sense} = $depth{plus}[$f][$h] / $count{plus}[$f][$h];
		}
		else 
		{
			$avg_depth{sense} = ($depth{plus}[$f][$h] + $depth{minus}[$f+1][-$h]) / ($count{plus}[$f][$h] + $count{minus}[$f+1][-$h]);	
		}

		# antisense (* -1 to get negative values)
		if ((!defined $count{minus}[$f][$h]) and (!defined $count{plus}[$f+1][-$h]))	
		{
			$avg_depth{antisense} = "NA";
		}
		elsif (!defined $count{plus}[$f+1][-$h])
		{
			$avg_depth{antisense} = ($depth{minus}[$f][$h] / $count{minus}[$f][$h]) * -1;
		}
		elsif (!defined $count{minus}[$f][$h])
		{
			$avg_depth{antisense} = ($depth{plus}[$f+1][-$h] / $count{plus}[$f+1][-$h]) * -1;
		}
		else
		{
			$avg_depth{antisense} = (($depth{minus}[$f][$h] + $depth{plus}[$f+1][-$h]) / ($count{minus}[$f][$h] + $count{plus}[$f+1][-$h])) * -1;	
		}

		print OUT "$avg_depth{sense}\t$avg_depth{antisense}\t";
		
	}
	print OUT "\n";
}

close OUT;

############ make R script for graphing result ##############
### R script needs fixing: need to combine bedname with sense or antisense



open (R, ">metaplot.R") or die "can't open metaplot.R\n";
print R "library(ggplot2)\n";
print R "library(reshape)\n";
print R "pdf(file=\"metaplot.pdf\", family=\"Helvetica\", width=12, height=8)\n";
print R "plot<-read.table(\"metaplot.txt\", header=T)\n";
print R "plot.melt <- melt(plot[,c('bp', ";

for (my $w = 0; $w < @longname; $w++)
{
	print R "'$longname[$w]'";
	print R ", " unless ($w == $num_longname - 1);
}

print R ")], id.vars=1)\n";
print R "ggplot(plot.melt, aes(x=bp, y=value, colour=variable, group=variable)) + geom_smooth() + theme_bw() + opts(panel.grid.minor=theme_blank()) + scale_colour_brewer(palette=\"Set1\", name=\"Bed\")\n";

close R;

################ run that R script! ##############

`R --vanilla < metaplot.R`
