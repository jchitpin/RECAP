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
our $version = "1.0.1";

# Modules
use warnings;
use strict;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use List::BinarySearch qw(binsearch_pos);
use autodie;
$| = 1;
# ==================================================================

######################## Begin Script Here #########################
####################################################################

# Start job time
BEGIN { our $start_run = time(); }

# Variables
my @fileOriginal;
my @fileRemix;
my @columnsOriginal;
my @columnsRemix;
my @pvalOrig;
my @pvalRemix;
my $matching;
my @RECAP;
my @fileHeader;

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
  'MACS=s'       => \my $MACS,
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
           || ! defined $delim) {
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
  print "  --MACS\tMACS summary file? (y)es or (n)o\n";
  print "  --help\tDisplay this help and exit\n";
  exit;
}

print "##################################################\n";
print "######       RECAP RECALIBRATE v$version       ######\n";
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
print "MACS file type? : $MACS\n";

# Input validation
if ( ! -d $dirOriginal ) {
  die "ERROR: Directory of original file does not exist\n";
} elsif ( ! -e $nameOriginal ) {
  die "ERROR: Original file does not exist\n";
}
if ( ! -d $dirRemix ) {
  die "ERROR: Directory of re-mixed file does not exist\n";
} elsif ( ! -e $nameRemix ) {
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
  } elsif ( $delim eq "t" ) {
    $delim = qr/\t/;
  }
}

$MACS = lc $MACS;
if ( $MACS =~ /[^yn]/ ) {
  die "ERROR: Specify 'y' for MACS or 'n' otherwise\n";
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


# Extracting p-value columns from summary files
# Exponentiate p-values if MACS summary files
print "Scanning p-value column of original file\n";
foreach my $line (@fileOriginal) {
  @columnsOriginal = split(/$delim/, "$line");
  if ( $MACS eq "y" ) {
    push @pvalOrig, 10**-($columnsOriginal[$pvalCol]);
  } else {
    push @pvalOrig, $columnsOriginal[$pvalCol];
  }
}
print "Check!\n";
print "Scanning p-value column of re-mixed file\n";
foreach my $line (@fileRemix) {
  @columnsRemix = split(/$delim/, "$line");
  if ( $MACS eq "y" ) {
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

# RECAP algorithm
print "Beginning RECAP procedure\n";
print "This may take some time...\n";

# Implementing RECAP algorithm with binary search
for my $idx ( 0 .. $#pvalOrig ) {
  $matching = binsearch_pos { $a <=> $b} $pvalOrig[$idx], @sorted_pvalRemix;
  $RECAP[$idx] = $matching / scalar @sorted_pvalRemix;
}

# Save recalibrated p-values to file
print "Saving RECAP summary file to output\n";
open ( my $export, ">", $nameOutput );

# Saving header from original file
if ($header != 0) {
  chdir($dirOriginal);
  open($importOrig, '<', $nameOriginal);
  while(<$importOrig>) {
    next unless $. <= $header;
    push @fileHeader, $_;
  }
chomp(@fileHeader[$header-1]);
push @fileHeader, "\tRECAP\n";
print $export join "", @fileHeader;
}

# Add column of recalibrated p-values to input, original summary file
chdir($dirOutput);
for my $idx ( 0 .. $#fileOriginal ) {
  @columnsOriginal = split(/$delim/, "$fileOriginal[$idx]");
  chomp(@columnsOriginal);
  push @columnsOriginal, $RECAP[$idx], "\n";
  print $export join "\t", @columnsOriginal;
}
close $export;
print "Check!\n";

# Job time
my $end_run = time();
my $run_time = $end_run - our $start_run;
print "Job took $run_time seconds\n";
