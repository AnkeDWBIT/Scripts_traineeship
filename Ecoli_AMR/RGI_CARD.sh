#!/bin/bash
# Script to loop through an input directory and run RGI CARD serially for all fasta files

# STEP 1 : VALIDATE INPUT PARAMETERS
################################################################################################################
function usage(){
	errorString="This RGI CARD script requires 1 parameter:\n
		1. Full path to input directory with fasta-files.";
	echo -e ${errorString};
	exit 1;
}
if [ "$#" -ne 1 ]; then
	usage
fi

# STEP 2 : GO THROUGH THE INPUT DIRECTORY, MAKE AN OUTPUT DIRECTORY AND RUN RGI CARD ON ALL .fasta FILES  
################################################################################################################
inputFolder=$1;
# Remove trailing slash for the input parameter if this is last char
len=${#inputFolder};
lastPos=$(expr $len - 1);
lastChar=${inputFolder:$lastPos:1};
if [[ $lastChar == '/' ]]; then
	inputFolder=${inputFolder:0:$lastPos};
fi

echo inputFolder: $inputFolder

# Check if the directory exists and make an output directory
outputFolder="${inputFolder}/RGI_output"
if [ ! -d "$outputFolder" ]; then
    # If it doesn't exist, create it
    mkdir -p "$outputFolder"
    echo "Output folder $outputFolder created."
else
    # If the folder already exists, stop the script
    echo "Output folder already exists. Please remove the folder or change the output folder name to avoid overwriting."
    #exit 1
fi

# Change the file extension from .FASTA to .fasta
for isolate in $(ls ${inputFolder}/*.FASTA);
    do
        mv $isolate ${isolate/.FASTA/.fasta}
    done 

# Go through the fasta-files in the input directory, run RGI CARD on them and save the output in the output directory
for isolate in $(ls ${inputFolder}/*.fasta); 
    do
        # Remove .fasta from the basename
        base=$(basename $isolate)
        base="${base%.*}"
        # Run RGI CARD (installed as Docker)
        docker run -v "$inputFolder":/data quay.io/biocontainers/rgi:6.0.3--pyha8f3691_0 rgi main -i "/data/${base}.fasta" -o "/data/RGI_output/${base}_output" --clean
        # Print message when RGI is run for a file
        echo "RGI CARD run for $isolate"
    done
    


