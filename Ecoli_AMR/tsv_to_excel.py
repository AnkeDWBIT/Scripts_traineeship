#!/usr/bin/env python
# Script to convert .tsv file to Excel-file
import sys, os , csv
from pathlib import Path
from xlsxwriter.workbook import Workbook

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

# STEP 2 : MAKE AN OUTPUT FOLDER TO STORE ALL THE EXCEL FILES
###########################################################################################################
# If input directory does not end with a slash, add it
if not input_dir.endswith("/"):
    input_dir += "/"
    # Create the path to the output directory
    output_dir = os.path.dirname(input_dir) + "/AMR_output_excel/"
else:
    output_dir = input_dir + "AMR_output_excel/"

# Check if the output directory exists, if not create it
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# STEP 3 : LOOK FOR TSV FILES AND CONVERT FORMAT TO EXCEL    
###########################################################################################################
# Iterate over files in the input directory
for file in os.listdir(input_dir):
    if file.endswith(".tsv"):
        # Retrieve the isolate name from the file name
        isolate_name = file.split(".")[0]
        # Create an XlsxWriter workbook object and add a worksheet.
        workbook = Workbook(output_dir + isolate_name + ".xlsx")
        worksheet = workbook.add_worksheet()
        # Create path to file
        file_path = input_dir + file
        # Create a TSV file reader.
        tsv_reader = csv.reader(open(file_path, 'r'), delimiter='\t')
        # Read the row data from the TSV file and write it to the XLSX file.
        for row, data in enumerate(tsv_reader):
            worksheet.write_row(row, 0, data)
        # Close the XLSX file.
        workbook.close()                   

# Print message to indicate that the script has finished
print("Successfully made Excel-files.")
