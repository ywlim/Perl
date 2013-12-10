#!/usr/bin/perl
# wig_to_metaplot_strand_specific.pl by Yoong Wearn Lim on 11/19/13
# modified from wig_to_metaplot.pl written on 2/3/13
# input wig file, and some bed files, and generate metaplot for the bed region


use strict; use warnings;
use Getopt::Std;	# first time using get opt! excited!
use FileHandle;

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
my $span;
#### WIG FILE ##########
# q = 0 is plus wig file
# q = 1 is minus wig file
# breaking wig file up
my $trigger = 0;
print "Pre-processing wig files...\n";
for (my $q = 0; $q < 2; $q++)
{
	open (WIG, $wig[$q]) or die "can't open $wig[$q] wig file\n";

	my $chr;

	while (my $line = <WIG>)
	{
		chomp $line;
		next if (($line !~ m/^variableStep/) and ($line !~ /^\d+/));
		die "sorry, fixedStep wig file not supported\n" if ($line =~ m/^fixedStep/);

		if ($line =~ m/^variableStep/)
		{
			close TEMP if ($trigger == 1);
			($chr, $span) = $line =~ m/^variableStep\schrom=chr(\w+)\sspan=(\d+)/;
			# print "chr is $chr and span is $span\n";
			print "Breaking chr$chr...\n";
			my $wig_out = $chr . "_" . $q . "_wig.temp";
			open (TEMP, ">$wig_out") or die "error writing to $wig_out\n";
		}
		elsif ($line =~ m/\d+\t\d+/)
		{
			print TEMP "$line\n";
			$trigger = 1;
		}
	}
}

# breaking bed files up
my %FH;
my @chromosome;
my %seen;
for (my $f = 0; $f < @bed; $f++)	# loop each bed file
{
	open (BED, $bed[$f]) or die "can't open $bed[$f] bed file\n";

	print "Pre-processing bed file: $bed[$f]\n";
	while (my $line = <BED>)
	{
		chomp $line;
		my ($chrom) = $line =~ m/^chr(\w+)/;
        if (!defined $seen{$chrom})
        {
            push (@chromosome, $chrom); # these are the chromosomes with bed values
            $seen{$chrom} = $chrom;
        }

		if (!exists $FH{$chrom}{$f})	# see chrom for the first time, open a new file
		{

			$FH{$chrom}{$f} = new FileHandle;
			$FH{$chrom}{$f}->open(">$chrom\_$f\_bed.temp");
			$FH{$chrom}{$f}->print("$line\n");
		}

		else
		{
			$FH{$chrom}{$f}->print("$line\n");
		}
	}
	close BED;
}


my %depth; my %count;
my $k;
for (my $c = 0; $c < @chromosome; $c++)
{
	print "Loading chr$chromosome[$c] miniwig into memory...\n";
	# load miniwig into memory
	for (my $q = 0; $q < 2; $q++)
	{
		open (MINIWIG, "$chromosome[$c]\_$q\_wig.temp") or last;    # last because that wig for this chro does't exist;

		while (my $line = <MINIWIG>)
		{
			chomp $line;
			my ($position, $value) = $line =~ m/(\d+)\t(\S+)/;

			for (my $i = $position; $i <= ($position + $span); $i++)
			{
				if ($q == 0)
				{
					if (!defined $val{plus}{$chromosome[$c]}{$i})
					{
						$val{plus}{$chromosome[$c]}{$i} = 0;
					}
					$val{plus}{$chromosome[$c]}{$i} += $value;
				}
				elsif ($q == 1)
				{
					if (!defined $val{minus}{$chromosome[$c]}{$i})
					{
						$val{minus}{$chromosome[$c]}{$i} = 0;
					}
					$val{minus}{$chromosome[$c]}{$i} += $value;
				}
			}
		}
		close MINIWIG;
	}

	# process the same chrom bed file


	for (my $f = 0; $f < @bed; $f++)	# loop each bed file
	{
		# need to close all bed filehandles
		if (exists $FH{$chromosome[$c]}{$f})
		{
    	    $FH{$chromosome[$c]}{$f}->close;
    	}

		print "Working on chr$chromosome[$c] of $bed[$f]...\n";
		open (MINIBED, "$chromosome[$c]\_$f\_bed.temp") or last;

		while (my $line2 = <MINIBED>)
		{
			chomp $line2;
			my ($chro, $start, $end) = $line2 =~ m/^chr(\w+)\t(\d+)\t(\d+)/;
			#print "$chro	$start	$end\n";

			$k = 1;	# position within the bed window
			for (my $j = $start; $j <= $end; $j++)
			{
				# get value from plus wig file
				if (exists $val{plus}{$chromosome[$c]}{$j})
				{
					$depth{plus}[$f][$k] += $val{plus}{$chromosome[$c]}{$j};
					$count{plus}[$f][$k]++;
				}

				# get value from minus wig file
				if (exists $val{minus}{$chromosome[$c]}{$j})
				{
					$depth{minus}[$f][$k] += $val{minus}{$chromosome[$c]}{$j};
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
		close MINIBED;
	}

	# remove miniwig from memory
	undef %val;
}



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

`R --vanilla < metaplot.R`;
`rm *.temp`;
