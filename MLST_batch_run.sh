#!/bin/bash
# This script is used to run classical MLST for multiple genomes at once


# VALIDATE INPUT
function usage(){
	errorString="This MLST script requires 2 parameters:\n
		1. Path of the folder with input files (fasta) to run MLST on (files can't be in subdirectories).\n
        2. Specify the species (e.g. paeruginosa).";

	echo -e ${errorString};
	exit 1;
}
if [ "$#" -ne 2 ]; then
	usage
fi

# 1. INPUT FOLDER CONTAINING .fna files
inputFolder=$1;
# Remove trailing slash if this is last char
len=${#inputFolder};
lastPos=$(expr $len - 1);
lastChar=${inputFolder:$lastPos:1};
if [[ $lastChar == '/' ]]; then
	inputFolder=${inputFolder:0:$lastPos};
fi

# 2. SPECIES
species=$2

# 3. RUN MLST FOR ALL GENOMES
# Concatenate filenames
for i in $(ls ${inputFolder}/*.fna); do
	inputFile="${i}";
    # Remove .fna
    posKeep=$(expr ${#i} - 4);
	baseNameTmp=${i:0:$posKeep};
    # Remove path to get sample name for output folder
	baseName=${baseNameTmp/"$inputFolder"/""};
    # Remove leading slash, if present
    baseName="${baseName#/}"
    # Show inputfiles
    echo "### Inputfile: $inputFile (basename: $baseName) ###";
    echo "### Running MLST ###";
    # Compose command
    MLSTCommand="python mlst.py -i ${inputFile} -s ${species} -o . -x --tmp_dir tmp_MLST --database database";

    # Show command
    echo -e "$MLSTCommand";

    # Execute
   output=$(eval $MLSTCommand);
    # Show output
   echo -e "$output";
   # Rename output file to basename.json
    mv "data.json" "MLST_${baseName}.json"
done