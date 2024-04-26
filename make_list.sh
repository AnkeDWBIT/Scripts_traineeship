#!/bin/bash

#Script to automatically make a list of the paths to all files in a directory, 1 per line.
#script requires 2 parameters: full path to inputdirectory with all fastafiles and paht to an outputfile
#outputfile will be created. Don't use the same directory for inputfiles and outputfile.
#Usage: make_ani_list.sh inputdirectory outputfile

inputdirectory=$1;
outputfile=$2


# Check if the correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 inputdirectory outputfile"
    exit 1
fi

inputdirectory="$1"
outputfile="$2"

# Check if input directory exists and is accessible
if [ ! -d "$inputdirectory" ]; then
    echo "Error: Input directory '$inputdirectory' does not exist or is not accessible."
    exit 1
fi

# Create or overwrite the output file
if ! touch "$outputfile"; then
    echo "Error: Unable to create or overwrite the output file '$outputfile'."
    exit 1
fi

# Use find to list all files in the input directory, exclude the script itself, and write their full paths to the output file
find "$inputdirectory" -type f ! -name "$(basename "$0")" -exec readlink -f {} + > "$outputfile"

echo "File paths written to $outputfile successfully."
