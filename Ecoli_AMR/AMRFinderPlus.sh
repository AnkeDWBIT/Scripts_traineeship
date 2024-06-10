#!/bin/bash
# Script to loop through an input directory and run AMRFinderPlus serially for all fasta files

# STEP 1 : VALIDATE INPUT PARAMETERS
################################################################################################################
function usage(){
	errorString="This AMR script requires 1 parameter:\n
		1. Full path to input directory with fasta-files.";
	echo -e ${errorString};
	exit 1;
}
if [ "$#" -ne 1 ]; then
	usage
fi

# STEP 2 : GO THROUGH THE INPUT DIRECTORY, MAKE AN OUTPUT DIRECTORY AND RUN AMRFINDERPLUS ON ALL .fasta FILES  
################################################################################################################
inputFolder=$1;
# Remove trailing slash for the input parameter if this is last char
len=${#inputFolder};
lastPos=$(expr $len - 1);
lastChar=${inputFolder:$lastPos:1};
if [[ $lastChar == '/' ]]; then
	inputFolder=${inputFolder:0:$lastPos};
fi

# Check if the directory exists and make an output directory
outputFolder="${inputFolder}/AMR_output"
if [ ! -d "$outputFolder" ]; then
    # If it doesn't exist, create it
    mkdir -p "$outputFolder"
    echo "Output folder created."
else
    echo "Output folder already exists."
fi


# Change the file extension from .FASTA to .fasta
for isolate in $(ls ${inputFolder}/*.FASTA);
    do
        mv $isolate ${isolate/.FASTA/.fasta}
    done 

# Go through the fasta-files in the input directory, run AMRFinderPlus on them and save the output in the output directory as tsv-file
for isolate in $(ls ${inputFolder}/*.fasta); 
    do
        # Remove .fasta from the basename
        base=$(basename $isolate)
        base="${base%.*}"
        # Construct the output filename
        output_file="${outputFolder}/${base}_output.tsv"
        # Run AMFinderPlus
        amrfinder -n $isolate --organism Escherichia -o $output_file
    done