#!/bin/bash
# ===============================================================
# RECAP Wrapper
# RECAP is a wrapper algorithm that resamples ChIP-seq and control
# data to estimate and control for biases built into peak calling
# algorithms. 
# This wrapper script is designed to conveniently peak call 
# ChIP-seq data using SICER and recalibrate peak p-values using 
# RECAP in one convenient package.
# A detailed explanation of the algorithm can be found under 
# PUBLICATIONS. 
#
# HISTORY:
#   29/08/2017 - v1.0.0 - First Creation
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
#SBATCH --job-name=JobRECAP_SICER
#SBATCH -c 1
#SBATCH --mem=10000
#SBATCH -t 48:00:00
#SBATCH --array=0-1:1
#SBATCH -o /global/home/hpc3862/Projects/LFDR_Bias_ChIP-Seq/SIMULATIONS/Simulation_01/log/JobS_RECAP_%A_%a.out
#SBATCH -e /global/home/hpc3862/Projects/LFDR_Bias_ChIP-Seq/SIMULATIONS/Simulation_01/log/JobS_RECAP_%A_%a.err
# ===============================================================
# Load Python modules for SICER
module load anaconda/2.7.13
#module load anaconda/3.5.3
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
# Window size
WINDOW=100
# Gap size
GAP=200
# Output directory for subsequent SICER and RECAP analyses
OUTPUT_DIR="/global/home/hpc3862/Projects/LFDR_Bias_ChIP-Seq/SIMULATIONS/Simulation_01/files/2.Results/SICER"
# Number of remixes for RECAP recalibration (default=1)
BOOTSTRAP=10
# Number of header lines in SICER summary file (default=0)
HEADER=0
# ===============================================================


# ===============================================================
# Provide a variable for the location of this and other scripts
SCRIPT_PATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REMIX_PATH=$(find ~/ -type f -name "RECAP_Re-Mix.sh" | head -n 1)
PERL_PATH=$(find ~/ -type f -name "RECAP.pl" | head -n 1)
SICER_PATH=$(find ~/ -type f -name "SICER.sh" | head -n 1)
# ===============================================================


## PLEASE EDIT AND ADD YOUR DESIRED SICER PARAMETERS IN 2) AND 3)
# ===============================================================
# 0) Timer function for entire ChIP-seq analysis and RECAP process
start=`date +%s`

# 1) Re-mix ChIP and control bed files
#bash $REMIX_PATH -i $INPUT_DIR -t $CHIP_NAME -c $CONTROL_NAME -o $OUTPUT_DIR -m unequal -b $BOOTSTRAP

# 2) Call original peaks using SICER
# Please specify your own SICER parameters!
# NOTE: p-value threshold must be set to 1.0 for SICER

cd $OUTPUT_DIR

if [ -n "${SLURM_ARRAY_TASK_ID}" ]
then

	# Switching parent directory for SICER temporary files generated during each analysis
	mkdir -p SICER_original
	mkdir -p SICER_Temp
	cd SICER_Temp
	rm -rf ${SLURM_ARRAY_TASK_ID}
	mkdir ${SLURM_ARRAY_TASK_ID}
	cd ${SLURM_ARRAY_TASK_ID}

	bash $SICER_PATH  $INPUT_DIR $CHIP_NAME $CONTROL_NAME "$OUTPUT_DIR/SICER_original" hg38 1 $WINDOW 100 1 $GAP 1
else
	mkdir -p SICER_original
	cd "$OUTPUT_DIR/SICER_original"
	bash $SICER_PATH  $INPUT_DIR $CHIP_NAME $CONTROL_NAME "$OUTPUT_DIR/SICER_original" hg38 1 $WINDOW 100 1 $GAP 1
fi

# 3) Call re-mixed peaks using SICER specifying desired parameters
# Please specify your own SICER parameters!
# NOTE: p-value threshold must be set to 1.0 for SICER
cd $OUTPUT_DIR

mkdir -p SICER_re-mix
mkdir -p SICER_Temp_re-mix
cd SICER_Temp_re-mix
rm -rf ${SLURM_ARRAY_TASK_ID}
mkdir ${SLURM_ARRAY_TASK_ID}
cd ${SLURM_ARRAY_TASK_ID}

for (( i=1; i<=$BOOTSTRAP; i++ ))
do
	bash $SICER_PATH  "$OUTPUT_DIR/re-mix" "${CHIP_NAME%.bed}.bootstrap_$i.bed" "${CONTROL_NAME%.bed}.bootstrap_$i.bed" "$OUTPUT_DIR/SICER_re-mix" hg38 1 $WINDOW 100 1 $GAP 1
done

# SICER generates those silly temporary .chr files in the folder containing the script (I think; need to check)
# I probably need to generate a temporary file, cd into it, then run the bash script.
# Afterwards, I need to move the contents back into $OUTPUT_DIR/SICER_original

# All non-SICER summary files in SICER_re-mix must be deleted if $BOOTSTRAP > 1
cd "$OUTPUT_DIR/SICER_re-mix"
find . -type f ! -name '*-islands-summary' -delete 

# 4) Recalibrate original peak p-values using RECAP
# NOTE: Check for correct header and p-value column if you obtain any errors here
cd $OUTPUT_DIR
mkdir -p SICER_RECAP

if [ $BOOTSTRAP -eq 1 ]
then
	perl $PERL_PATH --dirOrig "$OUTPUT_DIR/SICER_original" --nameOrig "${CHIP_NAME%.bed}-W${WINDOW}-G${GAP}-islands-summary" --dirRemix "$OUTPUT_DIR/SICER_re-mix" --nameRemix "${CHIP_NAME%.bed}" --dirOutput "$OUTPUT_DIR/SICER_RECAP" --nameOutput "${CHIP_NAME%.bed}.RECAP.bootstrap_1-W${WINDOW}-G${GAP}-islands-summary" --bootstrap 1 --header $HEADER --pvalCol 6 --delim t --software O
else
	for i in {1,$BOOTSTRAP}
	do
		perl $PERL_PATH --dirOrig "$OUTPUT_DIR/SICER_original" --nameOrig "${CHIP_NAME%.bed}-W${WINDOW}-G${GAP}-islands-summary" --dirRemix "$OUTPUT_DIR/SICER_re-mix" --nameRemix "${CHIP_NAME%.bed}" --dirOutput "$OUTPUT_DIR/SICER_RECAP" --nameOutput "${CHIP_NAME%.bed}.RECAP.bootstrap_${i}-W${WINDOW}-G${GAP}-islands-summary" --bootstrap $i --header $HEADER --pvalCol 6 --delim t --software O
	done
fi
rm chr.list

# 5) Finish timer and display runtime
end=`date +%s`
runtime=$((end-start))
echo "Runtime was $runtime seconds"
# ===============================================================
