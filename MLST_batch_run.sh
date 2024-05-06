#!/bin/bash
# This script is used to run classical MLST for multiple genomes at once


# VALIDATE INPUT
function usage(){
	errorString="This MLST script requires 2 parameters:\n
		1. Full path to folder with input files (fasta) to run MLST on (files can't be in subdirectories).\n
        2. Specify the species (e.g. paeruginosa).\n
        3. Full path to output folder (will create if it doens't exist yet).\n
        4. Create extended output (y/n)?";

	echo -e ${errorString};
	exit 1;
}
if [ "$#" -ne 4 ]; then
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

# 3. OUTPUT FOLDER
outputFolder=$3
# Check if the directory exists
if [ ! -d "$outputFolder" ]; then
    # If it doesn't exist, create it
    mkdir -p "$outputFolder"
    echo "Output folder created."
else
    echo "Output folder already exists."
fi

# 4. EXTENDED OUPUT
# Extended output creates files results.txt, results_tab.tsv, Hit_in_genome_seq.fsa, MLST_allele_seq.fsa
# If user specifies "y"", the argument -x will be added to the MLSTCommand in section 5 of this code.


# 5. RUN MLST FOR ALL GENOMES
# Concatenate filenames
for i in $(ls ${inputFolder}/*.fna); do
	inputFile="${i}";
    # Remove .fna
    posKeep=$(expr ${#i} - 4);
	baseNameTmp=${i:0:$posKeep};
    # Remove path to get genome assembly for output folder
	baseName=${baseNameTmp/"$inputFolder"/""};
    # Remove leading slash, if present
    baseName="${baseName#/}";
    # Shorten the basename (only retain the primary identifier e.g. GCF_000168335.1_ASM16833v1_genomic -> GCF_000168335.1)
    baseName=$(echo "$baseName" | cut -d '_' -f 1,2)
    
    # Show inputfiles
    echo "### Inputfile: $inputFile (basename: $baseName) ###";
    echo "### Running MLST ###";

    # Compose command
    MLSTCommand="python mlst.py -i ${inputFile} -s ${species} -o ${outputFolder} --tmp_dir tmp_MLST --database database";
    if [ "$4" = "y" ]; then
    MLSTCommand="$MLSTCommand -x"
fi

    # Show command
    echo -e "$MLSTCommand";

    # Execute
   output=$(eval $MLSTCommand);

    # Show output
   echo -e "$output";

   # Rename output files
   # Rename standard output file (data.json)
    mv "${outputFolder}/data.json" "${outputFolder}/MLST_${baseName}.json"
    # Rename extended output files
    if [ "$4" = "y" ]; then
    mv "${outputFolder}/results.txt" "${outputFolder}/results_${baseName}.txt"
    mv "${outputFolder}/results_tab.tsv" "${outputFolder}/results_tab_${baseName}.tsv"
    mv "${outputFolder}/Hit_in_genome_seq.fsa" "${outputFolder}/Hit_in_genome_seq_${baseName}.fsa"
    mv "${outputFolder}/MLST_allele_seq.fsa" "${outputFolder}/MLST_allele_seq_${baseName}.fsa"
    fi

done