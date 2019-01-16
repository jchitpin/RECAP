Please see <https://github.com/theodorejperkins/RECAP> for the latest stable build.  
Bleeding edge builds will be posted here (none currently).

## Introduction

ChIP-seq is used extensively to identify sites of transcription factor binding or regions of epigenetic modifications to the genome. The fundamental bioinformatics problem is to take ChIP-se read data and data representing some kind of control, and determine genomic regions that are enriched in the ChIP-se versus the control, also called "peak calling." While many programs have been designed to solve this task, nearly all fall into the statistical trap of using the data twice--once to determine candidate enriched regions, and a second time to assess enrichment by methods of classical statistical hypothesis testing. This double use of the data has the potential to invalidate the statistical significance assigned to enriched regions, or "peaks", and as a consequence, to invalidate false discovery rate estimates. Thus, the true significance or reliability of peak calls remains unknown. We propose a new wrapper algorithm, RECAP, that uses resampling of ChIP-seq and control data to estimate and correct for biases built into peak calling algorithms. RECAP is a powerful new tool for assessing the true statistical significance of ChIP-seq peak calls.

## Installation

Download RECAP and extract to your desired directory. There should be five scripts:
1.  RECAP_Re-Mix.sh
1.  RECAP.pl
1. RECAP_MACS.sh
1. RECAP_SICER.sh
1. RECAP_diffReps.sh

The first two scripts should be runnable on any system with Bash and Perl installed. Several CPAN modules are required. To install them:
* cpan
* install List::BinarySearch
* install List::Util
* install List::MoreUtils
* install Math::Utils

## Usage

### Example Workflow using Wrapper Scripts

Argument | Description
------------ | -------------
-i, --input | Input treatment/control BED file directory (absolute path)
-t, --treatment | Treatment BED file
-c, --control | Control BED file
-o, --output | Output file directory (absolute path and must exist)
-b, --bootstrap | Number of re-mixes
-e, --header | Header number of peak calling output files
-h, --help | Display this help and exit


#### MACS

Suppose we are interested in analyzing a treatment and control file with MACS.

1. Open RECAP_MACS.sh and modify any peak-calling preferences in line 101 and 109 *except* the p-value threshold. Possible options are listed on the MACS Github page https://github.com/taoliu/MACS. Currently, the recommended default MACS settings are used for regular peak calling. 

2. Run the script with the arguments below. Absolute directory paths must be specified and the output directory must already exist. A bootstrap of 1 is recommended and the header should be set to 29 (28 if using MACS --nomodel parameter).    
   ```bash RECAP_MACS.sh -i ~/ -t treatment_file.bed -c control_file.bed -o ~/output_directory -b 1 -e 29 header_lines_in_output_file   ```
   
3. Check the output directory containing re-mixed bed files in ```re-mix```, the original peak calling output files in ```MACS_original```, re-mixed peak calling output files in ```MACS_re-mix```, and the final recalibrated output files in ```MACS_RECAP```.

```
bash RECAP_Re-Mix.sh [-i] [-t] [-c] [-o] [-m] [-b] [-h]  
```

Argument | Description
------------ | -------------
-i, --input | Input treatment/control BED file directory (absolute path)
-t, --treatment | Treatment BED file
-c, --control | Control BED file
-o, --output | Output file directory (absolute path and must exist)
-m, --method | Method of re-mixing*
-b, --bootstrap | Number of re-mixes*
-h, --help | Display this help and exit

### Options(\*)

**-m, --method**  
Choose either *equal* or *unequal*. *Equal* distributes the treatment and control file reads into two equally sized re-mixed BED files. *Unequal* distributes the treatment and control file reads into two re-mixed BED files with the same read counts as the input files. *Unequal* is the recommended parameter.

**-b, --bootstrap**  
*Bootstrap* is the number of times the treatment and control files are re-mixed to generate re-mixed BED files.

```
perl RECAP.pl [--dirOrig]   [--nameOrig]   [--dirRemix] [--nameRemix] 
              [--dirOutput] [--nameOutput] [bootstrap]  [--header] 
              [--pvalCol]   [--delim]      [--software] [--help] 
```

Argument | Description
------------ | -------------
--dirOrig | Input original peak calling summary file directory (absolute path)
--nameOrig | Original peak calling summary file
--dirRemix | Input re-mixed peak calling summary file directory (absolute path)
--nameRemix | Re-mixed peak calling summary file ending in '.bootstrap_#.bed
--dirOutput | Output directory (absolute path and must exist)
--nameOutput | Original peak calling summary file with RECAP 
--bootstrap | Number of re-mixing procedures*
--header | Number of header lines in peak calling summary file
--pvalCol | Column number containing *p*-values in summary file
--delim | Delimiter type*
--software | Type of peak caller used*
-- help | Display this help and exit

### Options(\*)

**--bootstrap**
Ensure that the re-mixed peak calling file ends in '.boostrap_#.bed. Replace '#' with the bootstrap number.

**--delim**
Choose either *(c)omma* or *(t)ab* delimiters depending on the output of your peak caller.

**--software**
Choose either *(M)ACS* for MACS2, or *(D)iffReps* for diffReps, or *(O)ther* for another type of peak caller. Choosing *M* negative antilogs the *p*-values, a necessary step during *p*-value recalibration with RECAP. Choosing 'D' filters off any downregulated p-values, a special feature of diffReps that must be removed during *p*-value recalibration with RECAP.

### Notes

The re-mixing process takes minutes to perform. Recalibrating the *p*-values with the Perl script should take seconds to minutes.

### Example (Manual) Workflow

Suppose we are interested in analyzing a treatment and control file with MACS and recalibrating the resulting *p*-values.

1.  Re-mix treatment and control BED files: 

    ```bash RECAP_Re-Mix.sh -i ~/ChIP-Seq/files -t Treatment.bed -c Control.bed -o ~/ChIP-Seq/files/ -m unequal -b 1```
    
    This will create a new directory `~/ChIP-Seq/files/re-mix` with files `Treatment.bootstrap_1.bed` and `Control.bootstrap_1.bed`
    
1.  Analyze "original" files with MACS:

    ```macs2 callpeak -t Treatment.bed -c Control.bed --nomodel -p 0.1 -n Analysis --outdir ~/ChIP-Seq/analysis/```
    
    This will create several files including ```Treatment_peaks.xls```. **NOTE:** Please retain only the ```_peaks.xls``` file which is to be recalibrated.
    
1.  Analyze "re-mixed" files with MACS (please use *p*=0.1 for MACS and *p*=0.99 for SICER/diffReps):

    ```macs2 callpeak -t Treatment.bootstrap_1.bed -c Control --nomodel -p 0.1 -n Treatment.bootstrap_1 --outdir ~/ChIP-Seq/analysis/```
    
    This will create several files including ```Treatment.bootstrap_1_peaks.xls```. Please retain only the ```_peaks.xls``` file which is to be recalibrated.
    
1.  Recalibrate the *p*-values:

    ```perl RECAP.pl --dirOrig ~/ChIP-Seq/analysis/ --nameOrig Treatment_peaks.xls --dirRemix ~/ChIP-Seq/analysis --nameRemix Treatment --dirOutput ~/ChIP-Seq/analysis/ --nameOutput Treatment.RECAP.bootstrap_1.txt --bootstrap 1 --header 28 --pvalCol 7 --delim t --software M```
    
    **NOTE:** There are generally 29 header lines in the MACS summary file (28 if using --nomodel). The 7th column contains the *p*-values. The output file `Analysis.RECAP.bootstrap_1.txt` will retain the same header as the original summary file but contain a new column of recalibrated *p*-values and FDR-adjusted recalibrated *p*-values.
