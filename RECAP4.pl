#!/usr/bin/perl
# ==================================================================
# RECAP Re-Mix
# RECAP is a wrapper algorithm that resamples ChIP-seq and control
# data to estimate and control for biases built into peak calling
# algorithms.
# The purpose of this script is to recalibrate p-values based on
# the original ChIP-seq/control and re-mixed ChIP-seq/control
# data. A detailed explanation of the algorithm can be found
# under PUBLICATIONS.
#
# HISTORY:
#   02/01/2018 - v1.0.0 - First Creation
#   02/08/2018 - v1.0.1 - Outputs header and new column for RECAP p-values
#   02/12/2018 - v1.0.2 - Binary search instead of first index search
#   03/11/2018 - v1.0.3 - Filters off downregulated diffReps p-values
#   03/14/2018 - v1.0.4 - Calculates BH values for RECAP r-values
#   03/22/2018 - v1.0.5 - Calculates the LFDR binned by order of magitude
#   05/12/2018 - v1.0.6 - LFDR binned by half decade
#
# CREDITS:
# RECAP was developed by Justin G. Chitpin and Theodore J. Perkins.
# Development of RECAP was carried out at the Ottawa Hospital
# Research Institute in the Perkins Lab.
#
# PUBLICATIONS:
# If you use RECAP, please cite the following paper:
# <INSERT PUBLICATION HERE>
#
# QUESTIONS:
# Please contact tperkins@ohri.ca
# ==================================================================


# ==================================================================
# Version
our $version = "1.0.6";

# Modules
use warnings;
use strict;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use List::BinarySearch qw(binsearch_pos);
use List::Util qw(min max);
use Math::Utils qw(ceil floor);
use Statistics::KernelEstimation;
use autodie;
$| = 1;
# ==================================================================

######################## Begin Script Here #########################
####################################################################

# Start job time
BEGIN { our $start_run = time(); }

# Variables
my @fileOriginal; # Holds contents of original file
my @fileRemix;    # Holds contents of re-mixed file
my @columnsOriginal; # Holds contents of original file split by delimiter
my @columnsRemix;    # Holds contents of original file split by delimiter
my @pvalOrig;    # Original p-values
my @pvalRemix;   # Re-mixed p-values
my $matching;    # Number of re-mixed p-values equal or less than original p
my @RECAP;       # Re-calibrated r-values
my @fileHeader;  # File header for output file
my $delimOutput; # Delimiter type passed
my @indexFDR;    # Index of FDR values corresponding to original file
my @joinFDR;     # FDR values indexed by $indexFDR
my @RECAPhist;   # Re-calibrated r-values with 1.0 converted to 0.99
                 # Done to correctly implement histograms
my @binsLog;     # -log histogram bin edges. ex. 0:1 -> 1:0.1
my @bins;        # raw histogram bin edges. ex. 0:1 -> 1:0.1
my @histCount;   # Histogram count of the LFDR values by @bins
my @theoreticalBins; # Theoretical histogram count of LFDR values by @bins
my @LFDR;        # LFDR values calculated by histogram based estimation
my $KDE;         # Calculator object
my $LFDRKDE;     # Smoothed LFDR values by KDE method
my $wBandwidth;  # Bandwidth of KDE
my $minKDE;      # Minimum LFDR value for KDE
my $maxKDE;      # Maximum LFDR value for KDE


# Arguments
GetOptions(
  'dirOrig=s'    => \my $dirOriginal,
  'nameOrig=s'   => \my $nameOriginal,
  'dirRemix=s'   => \my $dirRemix,
  'nameRemix=s'  => \my $nameRemix,
  'dirOutput=s'  => \my $dirOutput,
  'nameOutput=s' => \my $nameOutput,
  'header=i'     => \my $header,
  'pvalCol=i'    => \my $pvalCol,
  'delim=s'      => \my $delim,
  'software=s'   => \my $software,
  'help'         => \my $help,
) or die "ERROR: Invalid options passed to $0\n";

# List parameters and check if all arguments filled
if ( $help || ! defined $dirOriginal
           || ! defined $nameOriginal
           || ! defined $dirRemix
           || ! defined $dirOutput
           || ! defined $nameOutput
           || ! defined $nameRemix
           || ! defined $header
           || ! defined $pvalCol
           || ! defined $delim
           || ! defined $software) {
  print BOLD, "  USAGE:\n", RESET;
  print "  --dirOrig\tOriginal file directory\n";
  print "  --nameOrig\tOriginal file name\n";
  print "  --dirRemix\tRe-mixed file directory\n";
  print "  --nameRemix\tRe-mixed file name\n";
  print "  --dirOutput\tOutput file directory\n";
  print "  --nameOutput\tOutput file name\n";
  print "  --header\tHeader lines in files\n";
  print "  --pvalCol\tColumn containing p-values\n";
  print "  --delim\tDelimiter type of file (t)ab or (c)omma\n";
  print "  --software\tPeak caller (M)ACS2 or (D)iffReps or (O)ther\n";
  print "  --help\tDisplay this help and exit\n";
  exit;
}

# Standard output
print "##################################################\n";
print "############       RECAP v$version       ############\n";
print "##################################################\n";
print "\n";
print "Input directory (Original): $dirOriginal\n";
print "Input file name (Original): $nameOriginal\n";
print "Input directory (Re-Mixed): $dirRemix\n";
print "Input file name (Re-Mixed): $nameRemix\n";
print "Output directory: $dirOutput\n";
print "Output file name: $nameOutput\n";
print "Header number   : $header\n";
print "p-value column  : $pvalCol\n";
print "Delimiter type  : $delim\n";
print "Peak caller type: $software\n";

# Input validation
if ( ! -d $dirOriginal ) {
  die "ERROR: Directory of original file does not exist\n";
}
chdir $dirOriginal;
if ( ! -e $nameOriginal ) {
  die "ERROR: Original file does not exist\n";
}
if ( ! -d $dirRemix ) {
  die "ERROR: Directory of re-mixed file does not exist\n";
}
chdir $dirRemix;
if ( ! -e $nameRemix ) {
  die "ERROR: Re-mixed file does not exist\n";
}
if ( ! -d $dirOutput ) {
  die "ERROR: Directory for output does not exist\n";
}

if ( ! $header == [0-9]*([0-9]) ) {
  die "ERROR: Specify a whole number for the header\n";
}
if ( ! $pvalCol == [0-9]*([0-9]) ) {
  die "ERROR: Specify a whole number for the p-value column\n";
} else {
  $pvalCol--; # p-value column must be n-1 for split() to work
}

$delim = lc $delim;
if ( $delim =~ /[^ct]/ ) {
  die "ERROR: Specify 'c' for comma or 't' for tab\n";
} else {
  if ( $delim eq "c" ) {
    $delim = qr/,/;
    $delimOutput = ",";
  } elsif ( $delim eq "t" ) {
    $delim = qr/\t/;
    $delimOutput = "\t";
  }
}

$software = uc $software;
if ( $software =~ /[^MDO]/ ) {
  die "ERROR: Specify 'M' for MACS, 'D' for diffReps, 'O' for other\n";
}

# Reading original summary file
print "\n";
print "Reading original summary file\n";
chdir($dirOriginal);
open(my $importOrig, '<', $nameOriginal);
while(<$importOrig>) {
  next unless $. > $header;
  push @fileOriginal, $_;
}
close $importOrig;
print "Check!\n";

# Reading re-mixed summary file
print "Reading re-mixed summary file\n";
chdir($dirRemix);
open(my $importRemix, $nameRemix);
while(<$importRemix>) {
  next unless $. > $header;
  push @fileRemix, $_;
}
close $importRemix;
print "Check!\n";

# Removing downregulated p-values if diffReps summary file
if ( $software eq "D" ) {
  print "Filtering off downregulated diffReps p-values from original/re-mixed files\n";
  @fileOriginal = grep( /Up/, @fileOriginal );
  @fileRemix = grep( /Up/, @fileRemix );
}

# Extracting p-value columns from summary files
# Exponentiate p-values if MACS summary files
print "Scanning p-value column of original file\n";
foreach my $line (@fileOriginal) {
  @columnsOriginal = split(/$delim/, "$line");
  if ( $software eq "M" ) {
    push @pvalOrig, 10**-($columnsOriginal[$pvalCol]);
  } else {
    push @pvalOrig, $columnsOriginal[$pvalCol];
  }
}
print "Check!\n";
print "Scanning p-value column of re-mixed file\n";
foreach my $line (@fileRemix) {
  @columnsRemix = split(/$delim/, "$line");
  if ( $software eq "M" ) {
    push @pvalRemix, 10**-($columnsRemix[$pvalCol]);
  } else {
    push @pvalRemix, $columnsRemix[$pvalCol];
  }
}
print "Check!\n";

# Sort re-mixed p-values
print "Sorting re-mixed p-values\n";
my @sorted_pvalRemix = sort { $a <=> $b } @pvalRemix;
print "Check!\n";

## RECAP algorithm
print "Recalibrating via RECAP procedure\n";

# Implementing RECAP algorithm with binary search
for my $idx ( 0 .. $#pvalOrig ) {
  $matching = binsearch_pos { $a <=> $b } $pvalOrig[$idx], @sorted_pvalRemix;
  $RECAP[$idx] = $matching / scalar @sorted_pvalRemix;
}

## Benjamini-Hochberg adjustment of RECAP r-values
print "FDR-adjusting RECAP r-values\n";

# Join RECAP values with its original index
@indexFDR = map{ $_ } 1 .. scalar @RECAP;
@joinFDR = map[ $indexFDR[$_], $RECAP[$_] ], 0 .. $#indexFDR;

# Sort 2D array by RECAP r-value
@joinFDR = sort{ $a->[1] <=> $b->[1] } @joinFDR;

# Implementing BH, allowing r-values of the same rank
{
  my $rank=0;
  for ( my $i=0; $i<scalar @joinFDR; $i++ ) {
    unless ( $i == scalar @joinFDR && $joinFDR[$i] == $joinFDR[$i+1] ) {
      $rank++;
    }
    $joinFDR[$i]->[1] = $joinFDR[$rank-1][1] / ($rank / scalar @joinFDR);
  }
}

# Adjusting FDR for the next higher raw r-value if it is smaller than r-value times m/i
for ( my $i=$#joinFDR; $i>0; $i-- ) {
  if ( $joinFDR[$i][1] < $joinFDR[$i-1][1] ) {
    $joinFDR[$i-1]->[1] = $joinFDR[$i][1];
  }
}

# Re-sort 2D array of FDR r-values by original index
@joinFDR = sort { $a->[0] <=> $b->[0] } @joinFDR;
# Only retain column of FDR r-values
@joinFDR = map{ $$_[1] } @joinFDR;

## Calculating the LFDR using a simple histogram-based estimation approach
print "Calculating the LFDR\n";

# Converting all 1.0 r-values to 0.99. This is a workaround to make
# the histogram count calculations accurate during the LFDR calculation.
# The actual r-values of 1.0 should be printed in the output column
foreach ( @RECAP ) {
  if ( $_ == 1 ) {
    $_ = 0.9999999999;
  }
  # Removing all r-values equal to 0. This is because the LFDR cannot
  # be calculated from r-values equal to 0.
  # The actual r-values of 0 should be printed in the output column
  push @RECAPhist, $_ if ($_) != 0;
}

# Calculate optimal histogram bin ranges
{
  my $max = max(@RECAPhist);
  my $min = min(@RECAPhist);
  $max = floor(-(log($max) / log(10)));
  $min = ceil(-(log($min) / log(10)));

  for ( my $i=$max; $i<=$min; $i=$i+0.5 ) { # CHANGED to half decade
    push @binsLog, $i;
  }
  # Raw half decade bin edges (ex. 1.0, 0.316, 0.1, 0.0316 ...)
  @bins = map { 10**(-$_) } @binsLog;
}

# Histogram count by binary searching index
{
  # Sort r-values for histogram counting
  my @RECAPsort = sort { $a <=> $b } @RECAPhist;
  my $index;
  #for my $i ( 0 .. ($#bins-1) ) {
  for ( my $i=0; $i<(scalar @bins - 1); $i++ ) {
    $index = (binsearch_pos {$a <=> $b} $bins[$i], @RECAPsort) - (binsearch_pos {$a <=> $b} $bins[$i+1], @RECAPsort);
    $histCount[$i] = $index / (scalar @RECAPsort);
  }
}

# Calculating theoretical bins under uniform distribution
{
  for ( my $i=0; $i<(scalar @bins - 1); $i++ ) {
    $theoreticalBins[$i] = $bins[$i] - $bins[$i+1];
  }
}

# Generate theoretical histogram under uniform distribution
# Calculate the local false discovery rate. Uniform(p) / pdf(RECAP)
# Histogram-based estimation
{
  #for $i ( 0 .. ($#bins-1) ) {
  for ( my $i=0; $i<(scalar @theoreticalBins); $i++ ) {
    unless ( $histCount[$i] == 0 ) {
      $LFDR[$i] = $theoreticalBins[$i] / $histCount[$i];
    } else {
      $LFDR[$i] = 0;
    }
    if ( $LFDR[$i] > 1 ) {
      $LFDR[$i] = 1;
    }
  }
}

# Assigning LFDR value to each row of RECAP r-values
my @LFDRorder;
for my $element (0 .. $#RECAP) {
  my $j=0;
  while ( $RECAP[$element] != 0 && $bins[$j] > $RECAP[$element] && $RECAP[$element] <= $bins[$j+1] ) {
    $j++;
  }
  unless ( $RECAP[$element] == 0 ) {
    push ( @LFDRorder, $LFDR[$j] );
  } else {
    push ( @LFDRorder, 1 );
  }
}

=pod
# Smoothing empirical LFDR by kernel density estimation
print "KDE smoothing of empirical LFDR histogram\n";

# Instantiate Gaussian kernel
$KDE = Statistics::KernelEstimation->new();

#my @LFDRorder2= grep { $_ != 0 } @LFDRorder;

# Add data
for $LFDRKDE ( @LFDRorder ) {
  $KDE->add_data( $LFDRKDE );
}

# Select method for optimizing bandwidth
$wBandwidth = $KDE->default_bandwidth();

# Calculating KDE for histogram bin estimated LFDR values
{ 
  #( $minKDE, $maxKDE ) = $KDE->extended_range();
  ( $minKDE, $maxKDE ) = $KDE->range();
  for( $LFDRKDE=$minKDE; $LFDRKDE<=$maxKDE; $LFDRKDE+=($maxKDE-$minKDE)/20 ) {
    #print $LFDRKDE, "\t", $KDE->pdf( $LFDRKDE, $wBandwidth ), "\t", $KDE->cdf( $LFDRKDE, $wBandwidth ), "\n";
    print $LFDRKDE, "\t", $KDE->pdf( $LFDRKDE ), "\t", $KDE->cdf( $LFDRKDE ), "\n";
  }
}
=cut
=pod
# Saving empirical LFDR distribution to text file
print "Saving histogram-based estimation of the LFDR to output\n";
chdir($dirOutput);
open ( my $exportLFDR, ">", "TESThbeLFDR.txt" );
# Printing header
print $exportLFDR "-log10(Bin)\tbin(Edge)\tbin(Width)\thist(RECAP)\thist(LFDR)\n";
# Printing contents
{
  for my $idx ( 0 .. $#histCount ) {
    #print "$bins[$idx]\n";
    push my @contents, $binsLog[$idx], $bins[$idx], $theoreticalBins[$idx], $histCount[$idx], $LFDR[$idx], "\n";
    print $exportLFDR join "\t", @contents; 
  }
}
close $exportLFDR;
print "Check!\n";
=cut

# Save recalibrated p-values to file
print "Saving RECAP summary file to output\n";
chdir($dirOutput);
open ( my $exportRECAP, ">", $nameOutput );

# Saving header from original file and adding RECAP BH(RECAP) column headers
if ( $header != 0 ) {
  chdir($dirOriginal);
  open($importOrig, '<', $nameOriginal);
  while(<$importOrig>) {
    next unless $. <= $header;
    push @fileHeader, $_;
  }
chomp($fileHeader[$header-1]);
push @fileHeader, "${delimOutput}RECAP${delimOutput}BH(RECAP)${delimOutput}LFDR\n";
print $exportRECAP join "", @fileHeader;
}

# Add column of r-values to input, original summary file
# Add column of Benjamini-Hochberg corrected r-value
# Add column of LFDR values
for my $idx ( 0 .. $#fileOriginal ) {
  @columnsOriginal = split(/$delim/, "$fileOriginal[$idx]");
  chomp(@columnsOriginal);
  push @columnsOriginal, $RECAP[$idx], $joinFDR[$idx], $LFDRorder[$idx], "\n";
  print $exportRECAP join $delimOutput, @columnsOriginal;
}
close $exportRECAP;
print "Check!\n";

# Job time
my $end_run = time();
my $run_time = $end_run - our $start_run;
print "Job took $run_time seconds\n";
