#!/bin/bash
# ===============================================================
# RECAP Wrapper
# RECAP is a wrapper algorithm that resamples ChIP-seq and control
# data to estimate and control for biases built into peak calling
# algorithms. 
# This wrapper script is designed to conveniently peak call 
# ChIP-seq data using diffReps and recalibrate peak p-values using 
# RECAP  in one convenient package.
# A detailed explanation of the algorithm can be found under 
# PUBLICATIONS. 
#
# HISTORY:
#   29/08/2017 - v1.0.0 - First Creation
#   16/11/2018 - v1.0.1 - Support for SLURM batch processing
#
# CREDITS:
# RECAP was developed by Justin G. Chitpin, Aseel Awdeh, 
# and Theodore J. Perkins.
# Development of RECAP was carried out at the Ottawa Hospital
# Research Institute in the Perkins Lab.
#
# PUBLICATIONS:
# If you use RECAP, please cite the following paper:
# <INSERT PUBLICATION HERE>
#
# QUESTIONS:
# Please contact tperkins@ohri.ca
# ===============================================================


## OPTIONAL SLURM DIRECTIVES FOR BATCH PROCESSING WITH A CLUSTER
## Double comment if not using batch processing
## Else, fill in the two parameters located below
# ===============================================================
#SBATCH --job-name=JobRECAP_diffReps
#SBATCH -c 1
#SBATCH --mem=5000
#SBATCH -t 16:00:00
#SBATCH --array=0-9:1
#SBATCH -o /global/home/hpc3862/Projects/LFDR_Bias_ChIP-Seq/SIMULATIONS/Simulation_01/log/JobD_RECAP_%A_%a.out
#SBATCH -e /global/home/hpc3862/Projects/LFDR_Bias_ChIP-Seq/SIMULATIONS/Simulation_01/log/JobD_RECAP_%A_%a.err
# ===============================================================
# Location and name file containing all ChIPs to be analyzed
CHIP_LOC="/global/home/hpc3862/Projects/LFDR_Bias_ChIP-Seq/SIMULATIONS/Simulation_01/files/0.MatchingPairs/ChIP_Names.1.txt"
# Location of all matching control files to be analyzed
CONT_LOC="/global/home/hpc3862/Projects/LFDR_Bias_ChIP-Seq/SIMULATIONS/Simulation_01/files/0.MatchingPairs/ChIP_Names.2.txt"
ARRAY_CHIP=($(awk -F, 'FNR > 0 {print $1}' $CHIP_LOC ))
ARRAY_CONT=($(awk -F, 'FNR > 0 {print $1}' $CONT_LOC ))
# ===============================================================


## PLEASE FILL THE FOLLOWING PARAMETERS
# ===============================================================
# ChIP/Control directory
INPUT_DIR="/global/home/hpc3862/Projects/LFDR_Bias_ChIP-Seq/SIMULATIONS/Simulation_01/files/1.Raw/"
# ChIP name
#CHIP_NAME="Simulation_01.Sim.1.1.bed"
CHIP_NAME=${ARRAY_CHIP[${SLURM_ARRAY_TASK_ID}]}
# Control name
#CONTROL_NAME="Simulation_01.Sim.1.2.bed"
CONTROL_NAME=${ARRAY_CONT[${SLURM_ARRAY_TASK_ID}]}
# Output directory for subsequent diffReps and RECAP analyses
OUTPUT_DIR="/global/home/hpc3862/Projects/LFDR_Bias_ChIP-Seq/SIMULATIONS/Simulation_01/files/2.Results/diffReps/"
# Number of remixes for RECAP recalibration (default=1)
BOOTSTRAP=10
# Number of header lines in diffReps summary file (default=33)
HEADER=33
# ===============================================================


# ===============================================================
# Provide a variable for the location of this and other scripts
SCRIPT_PATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REMIX_PATH=$(find ~/ -type f -name "RECAP_Re-Mix.sh" | head -n 1)
PERL_PATH=$(find ~/ -type f -name "RECAP.pl" | head -n 1)
# ===============================================================


## PLEASE EDIT AND ADD YOUR DESIRED diffReps PARAMETERS IN 2) AND 3)
# ===============================================================
# 0) Timer function for entire ChIP-seq analysis and RECAP process
start=`date +%s`

# 1) Re-mix ChIP and control bed files
#bash $REMIX_PATH -i $INPUT_DIR -t $CHIP_NAME -c $CONTROL_NAME -o $OUTPUT_DIR -m unequal -b $BOOTSTRAP

# 2) Call original peaks using diffReps
# Please specify your own diffReps parameters!
# NOTE: p-value threshold must be set to 1.0 for diffReps
#cd $OUTPUT_DIR
#mkdir -p diffReps_original
#diffReps.pl --treatment "$INPUT_DIR/$CHIP_NAME" --control "$INPUT_DIR/$CONTROL_NAME" -meth gt -gname hg19 --pval 1 --mode n --report "$OUTPUT_DIR/diffReps_original/${CHIP_NAME%.bed}.txt" --nohs --noanno --nsd 20

# 3) Call re-mixed peaks using diffReps specifying desired parameters
# Please specify your own diffReps parameters!
# NOTE: p-value threshold must be set to 1.0 for diffReps
#cd $OUTPUT_DIR
#mkdir -p diffReps_re-mix
#for (( i=1; i<=$BOOTSTRAP; i++ ))
#do
	#diffReps.pl --treatment "$OUTPUT_DIR/re-mix/${CHIP_NAME%.bed}.bootstrap_$i.bed" --control "$OUTPUT_DIR/re-mix/${CONTROL_NAME%.bed}.bootstrap_$i.bed" -meth gt -gname hg19 --pval 1 --mode n --report "$OUTPUT_DIR/diffReps_re-mix/${CHIP_NAME%.bed}.bootstrap_$i.txt" --nohs --noanno --nsd 20
#done

# 4) Recalibrate original peak p-values using RECAP
# NOTE: Check for correct header and p-value column if you obtain any errors here
cd $OUTPUT_DIR
mkdir -p diffReps_RECAP

if [ $BOOTSTRAP -eq 1 ]
then
	perl $PERL_PATH --dirOrig "$OUTPUT_DIR/diffReps_original" --nameOrig "${CHIP_NAME%.bed}.txt" --dirRemix "$OUTPUT_DIR/diffReps_re-mix" --nameRemix "${CHIP_NAME%.bed}.bootstrap_$i.txt" --dirOutput "$OUTPUT_DIR/diffReps_RECAP" --nameOutput "${CHIP_NAME%.bed}.RECAP.bootstrap_1.txt" --bootstrap 1 --header $HEADER --pvalCol 13 --delim t --software D
else
	for i in {1,$BOOTSTRAP}
	do
		perl $PERL_PATH --dirOrig "$OUTPUT_DIR/diffReps_original" --nameOrig "${CHIP_NAME%.bed}.txt" --dirRemix "$OUTPUT_DIR/diffReps_re-mix" --nameRemix "${CHIP_NAME%.bed}" --dirOutput "$OUTPUT_DIR/diffReps_RECAP" --nameOutput "${CHIP_NAME%.bed}.RECAP.bootstrap_$i.txt" --bootstrap $i --header $HEADER --pvalCol 13 --delim t --software D
	done
fi

# 5) Finish timer and display runtime
end=`date +%s`
runtime=$((end-start))
echo "Runtime was $runtime seconds"
# ===============================================================
