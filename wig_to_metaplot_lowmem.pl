#!/usr/bin/perl
# wig_to_metaplot_lowmem.pl by Yoong Wearn Lim
# input wig file, and some bed files, and generate metaplot for the bed region

# for instance, if the wig file contains methylation info
# while bed files contains TSS coordinates (one file for each strand)
# the script will generate metaplot of methylation around TSS

# this script takes multiple sets of bed files
# so we can use, for example, bed coordinates for TSS and TTS in the same run

# note that bed coordinates have to be centered (eg. +- 500 TSS)

use strict; use warnings;
use Getopt::Std;	
use FileHandle;

our ($opt_h, $opt_w, $opt_b, $opt_n);
getopts('hw:b:n:');	# w and b take arguments, thus the : after them

my $usage = "usage: wig_to_metaplot.pl -w wigfile.wig -b bedfile1_plus.bed,bedfile1_minus.bed,bedfile2_plus.bed,bedfile2_minus.bed.. -n name1,name2..

For example:

wig_to_metaplot.pl -w AGS1_meth.wig -b tss_plus.bed,tss_minus.bed,tts_plus.bed,tts_minus.bed -n tss,tts

example of tss_plus.bed (+-2kb around tss for genes on the plus strand):

chr1	9873	13873
chr1	67090	71090
chr1	321891	325891
chr1	365658	369658

example of tss_minus.bed (+-2kb around tss for genes on the minus strand):

chr1	27370	31370
chr1	34081	38081
chr1	138566	142566
chr1	620034	624034\n";

if (($opt_h) or (!$opt_b) or (!$opt_w) or (!$opt_n))
{
	die $usage;
}

my $wig = $opt_w;
my @bed = split(",", $opt_b);
my @name = split(",", $opt_n);
my $num_bed = @bed;
my $num_name = @name;

die "Numbers of entry of bed and name don't match! Check that number of bed files = number of names * 2\n" unless ($num_name == $num_bed / 2);

# check bed files before spending time processing wig file

for (my $f = 0; $f < @bed; $f++)	# loop each bed file
{
	open (BED, $bed[$f]) or die "can't open $bed[$f] bed file\n";
	close BED;
}

#### WIG FILE ############

# breaking wig file up by chromosome
# as a way to save memory

my $trigger = 0;
print "Pre-processing wig files...\n";

if ($wig =~ /\.gz$/)
{
    open (WIG, "gunzip -c $wig |") || die "can't open pipe to $wig\n";
}
else
{
    open (WIG, $wig) || die "can't open $wig\n";
}

my $chr; my $span; my %val;

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
        my $wig_out = $chr . "_wig.temp";
        open (TEMP, ">$wig_out") or die "error writing to $wig_out\n";
    }
    elsif ($line =~ m/\d+\t\d+/)
    {
        print TEMP "$line\n";
        $trigger = 1;
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

my @depth; my @count;
my $k;
for (my $c = 0; $c < @chromosome; $c++)
{
	print "Loading chr$chromosome[$c] miniwig into memory...\n";
	# load miniwig into memory

    open (MINIWIG, "$chromosome[$c]\_wig.temp") or last;    # last because that wig for this chro does't exist;

    while (my $line = <MINIWIG>)
    {
        chomp $line;
        my ($position, $value) = $line =~ m/(\d+)\t(\S+)/;

        for (my $i = $position; $i <= ($position + $span); $i++)
        {
            if (!defined $val{$chromosome[$c]}{$i})
            {
                $val{$chromosome[$c]}{$i} = 0;
            }
            $val{$chromosome[$c]}{$i} += $value;
        }
    }
    close MINIWIG;

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
				# get value from wig file
				if (exists $val{$chromosome[$c]}{$j})
				{
					$depth[$f][$k] += $val{$chromosome[$c]}{$j};
					$count[$f][$k]++;
				}

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

`rm *.temp`;

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
