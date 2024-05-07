#!/usr/bin/env python
# Script to write a short summary of BUSCO output to an text file
import sys
import os
import json
import re

# STEP 1 : checking the command line arguments
############################################################################################################
# Check if any command-line arguments have been provided
if len(sys.argv) < 3:
	print("Error: Need to provide the correct command-line arguments.")
	print("Usage: python scriptname.py [1] [2] [3]")
	print("\t[1] = Full path to directory with BUSCO output files")
	print("\t[2] = Choose an identifier to add to the file-name (e.g. Paeurginosa_1 -> busco_summary_Paeurginosa_1.txt)")
	print("\t[3] = Full path to output directory)")
	sys.exit(1)

# Store command-line argument(s) in objects & construct the excel file name
input_dir = sys.argv[1]
file_id = sys.argv[2]
output_file = "busco_summary_" + file_id + ".txt"
output_dir = sys.argv[3]

# STEP 2 : Making a new text file to store RefSeq-IDs & short summary of BUSCO output int
############################################################################################################
# Path to the text file within the output directory
file_path = os.path.join(output_dir, output_file)
# Open the text file in append mode
with open(file_path, 'a') as file:
# Write headers
    file.write("RefSeq Identifier\tBusco Summary\n")

    # STEP 3 : Looking for short_summary.json files in the BUSCO output directory and reading the data
    ############################################################################################################
    # Iterate over subdirectories
    for subdir in os.listdir(input_dir):
        subdir_path = os.path.join(input_dir, subdir)
        if os.path.isdir(subdir_path):
            # Iterate over files in each subdirectory
            for filename in os.listdir(subdir_path):
                if "short_summary" in filename and filename.endswith(".json"):
                    # Parse out the RefSeq identifier for later use
                    refseq_id = filename.split('.')[-4]
                    print(refseq_id)
                    # Contstruct paths to the JSON short summary files
                    short_summary_path = os.path.join(subdir_path, filename)
                    print(short_summary_path)
                    # Open the JSON files and read the data
                    with open(short_summary_path) as f:
                        data = json.load(f)
                        # Retrieve one_line_summary from the JSON file
                        onelinesum = data['results']['one_line_summary']
                        print(onelinesum)
                        # Write data to the text file
                        file.write(f"{refseq_id}\t{onelinesum}\n")
                        # Close the JSON file
                        f.close()

# Close the text file
file.close()

print("Output written to file:" + file_path)