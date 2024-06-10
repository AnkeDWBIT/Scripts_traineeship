#!/usr/bin/env python
# Script to convert .tsv file to .csv file
import sys, os , re

# STEP 1 : CHECK & PROCESS COMMAND LINE ARGUMENTS
###########################################################################################################
# Check if any command-line arguments have been provided
if len(sys.argv) < 2:
	print("Error: Need to provide command-line arguments.")
	print("Usage: python scriptname.py [1]")
	print("\t[1] = Full path to input directory (=folder containing only tsv-files)")
	sys.exit(1)

# Store the command-line argument(s) in an object
input_dir = sys.argv[1]
# If input directory does not end with a slash, add it
if not input_dir.endswith("/"):
    input_dir += "/"

# STEP 2 : LOOK FOR TSV FILES AND CONVERT FORMAT TO CSV     
###########################################################################################################
# Iterate over files in the input directory
for file in os.listdir(input_dir):
    # Construct the full path to the file
    file_path = os.path.join(input_dir, file)
    # Open the tsv file and convert it to csv
    with open(file_path, 'r') as tsv_file:
        # Copy the file path and add a csv identifier
        csv_file_path = file_path + "_csv"
        # Open the new csv file in write mode
        with open(csv_file_path, 'w') as csv_file:
            for line in tsv_file:
                # Replace every tab with comma
                fileContent = re.sub("\t", ",", line)
                # Write the content to the csv file
                csv_file.write(fileContent)
                        

# Print message to indicate that the script has finished
print("Successfully made csv-files.")
