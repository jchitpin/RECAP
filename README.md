# Introduction

ChIP-seq is used extensively to identify sites of transcription factor binding or regions of epigenetic modifications to the genome. While many programs have been designed to analyze ChIP-seq data, nearly all fall into the statistical trap of using the data twice--once to determine candidate enriched regions, and a second time to assess enrichment by methods of classical statistical hypothesis testing. This double use of the data has the potential to invalidate the statistical significance assigned to enriched regions, or "peaks", and as a consequence, to invalidate false discovery rate estimates. Thus, the true significance or reliability of peak calls remains unknown. We propose a new wrapper algorithm called RECAP, that uses resampling of ChIP-seq and control data to estimate and correct for biases built into peak calling algorithms. RECAP is a powerful new tool for assessing the true statistical significance of ChIP-seq peak calls.

# Installation

Download RECAP and extract to your desired directory. There should be two scripts:
1.  RECAP_Re-Mix.sh
1.  RECAP.pl

The shell and Perl script should be runnable on any Unix system with Bash and Perl installed. 

# Usage

> bash RECAP_Re-Mix.sh [-i] [-t] [-c] [-o] [-m] [-b]  
> perl RECAP.pl --dirOrig --nameOrig --dirRemix --nameRemix --dirOutput --nameOutput --header --pvalCol --delim --MACS 


## Parameters

## Example

