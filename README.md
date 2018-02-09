## Introduction

ChIP-seq is used extensively to identify sites of transcription factor binding or regions of epigenetic modifications to the genome. While many programs have been designed to analyze ChIP-seq data, nearly all fall into the statistical trap of using the data twice--once to determine candidate enriched regions, and a second time to assess enrichment by methods of classical statistical hypothesis testing. This double use of the data has the potential to invalidate the statistical significance assigned to enriched regions, or "peaks", and as a consequence, to invalidate false discovery rate estimates. Thus, the true significance or reliability of peak calls remains unknown. We propose a new wrapper algorithm called RECAP, that uses resampling of ChIP-seq and control data to estimate and correct for biases built into peak calling algorithms. RECAP is a powerful new tool for assessing the true statistical significance of ChIP-seq peak calls.

## Installation

Download RECAP and extract to your desired directory. There should be two scripts:
1.  RECAP_Re-Mix.sh
1.  RECAP.pl

Both scripts should be runnable on any system with Bash and Perl installed with no extra CPAN modules required.

## Usage

```
bash RECAP_Re-Mix.sh [-i] [-t] [-c] [-o] [-m] [-b]  
```

Argument | Description
------------ | -------------
-i, --input | Input treatment/control BED file directory
-t, --treatment | Treatment BED file
-c, --control | Control BED file
-o, --output | Output file directory
-m, --method | Method of re-mixing*
-b, --bootstrap | Number of re-mixes*

### Options(\*)

**-m, --method**  
Choose either *equal* or *unequal*. *Equal* distributes the treatment and control file reads into two equally sized re-mixed BED files. *Unequal* distributes the treatment and control file reads into two re-mixed BED files with the same read counts as the input files. *Unequal* is the recommended parameter.

**-b, --bootstrap**  
*Bootstrap* is the number of times the treatment and control files are re-mixed to generate re-mixed BED files.

```
perl RECAP.pl [--dirOrig] [--nameOrig] [--dirRemix] [--nameRemix] 
              [--dirOutput] [--nameOutput] [--header] [--pvalCol] 
              [--delim] [--MACS] 
```

Argument | Description
------------ | -------------
--dirOrig | Input original peak calling summary file directory
--nameOrig | Original peak calling summary file
--dirRemix | Input re-mixed peak calling summary file directory
--nameRemix | Re-mixed peak calling summary file
--dirOutput | Output directory
--nameOutput | Original peak calling summary file with RECAP 
--header | Number of header lines in peak calling summary file
--pvalCol | Column number containing *p*-values in summary file
--delim | Delimiter type*
--MACS | Whether peak calling summary file belongs to MACS*

### Options(\*)

**--delim**
Choose either *(c)omma* or *(t)ab* delimiters depending on the output of your peak caller. **NOTE: Version 1.0.1 cannot handle .xls files. Please convert them to .txt for RECAP.pl to work.**

**--MACS**
Choose either *(y)es* or *(n)o* if the peak caller used is MACS. Choosing *yes* antilogs the *p*-values, a necessary step during *p*-value recalibration with RECAP.

### Example Workflow

Suppose we are interested in analyzing a treatment and control file with MACS and recalibrating the resulting *p*-values.

1.  Re-mix treatment and control BED files: 

    ```bash RECAP_Re-Mix.sh -i ~/ChIP-Seq/files -t Treatment.bed -c Control.bed -o ~/ChIP-Seq/files/ -m unequal -b 1```
    
    This will create a new directory `~/ChIP-Seq/files/re-mix` with files `Treatment.bootstrap_1.bed` and `Control.bootstrap_1.bed`
    
1.  Analyze "original" files with MACS:

    ```macs2 callpeak -t Treatment.bed -c Control.bed --nomodel -p 0.1 -n Analysis --outdir ~/ChIP-Seq/analysis/```
    
    This will create several files including ```Analysis_peaks.xls```
    
1.  Analyze "re-mixed" files with MACS:

    ```macs2 callpeak -t Treatment.bootstrap_1.bed -c Control.bootstrap_1.bed --nomodel -p 0.00001 -n Analysis.bootstrap_1 -- outdir ~/ChIP-Seq/analysis/```
    
    This will create several files including ```Analysis.bootstrap_1_peaks.xls```
    
1.  Convert the MACS summary files to .txt extension. One way of doing that on the command line is:

    ```for file in *.xls; do```
    
      ```mv "$file" "`basename "$file" .xls`.txt"```
      
    ```done```

1.  Recalibrate the *p*-values:

    ```perl RECAP.pl --dirOrig ~/ChIP-Seq/analysis/ --nameOrig Analysis_peaks.txt --dirRemix ~/ChIP-Seq/analysis --nameRemix Analysis.bootstrap_1_peaks.txt --dirOutput ~/ChIP-Seq/analysis/ Analysis.RECAP.bootstrap_1.txt --header 28 --pvalCol 7 --delim t --MACS y```
    
    There are generally 28 header lines in the MACS summary file. The 7th column contains the *p*-values. The output file `Analysis.RECAP.bootstrap_1.txt` will retain the same header as the original summary file but contain a new column of recalibrated *p*-values.
