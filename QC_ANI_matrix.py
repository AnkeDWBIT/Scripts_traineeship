#!/usr/bin/env python
# Script for quality control of ANI matrix
# Look for genomes that consistenly have a percent identity lower than 95 with all other genomes of the dataset (indicative of a diferent species)

import sys
from openpyxl import Workbook
from openpyxl import load_workbook
import itertools

# PROCESS INPUT ARGUMENTS
######################################################################################################################################
# Check if any command-line arguments have been provided
if len(sys.argv) < 2:
    print("Error: Need to provide the correct command-line arguments.")
    print("Usage: python scriptname.py [1] [2]")
    print("\t[1] = Full path to Excel file with ANI matrix and MLST results")
    print("\t[2] = Worksheet with ANI matrix in Excel file [1] e.g. Pseudomonas_aeruginosa")
    sys.exit(1)

# Store command-line argumentS
excel_file_path = sys.argv[1]
ANI_sheet_name = sys.argv[2]

# LOAD THE ANI MATRIX FROM THE EXCEL FILE  
######################################################################################################################################
# Load the workbook
#excel_file_path = "/home/guest/BIT11_Traineeship/Scripts_traineeship/FastANI_matrix_Pseudomonas_aeruginosa_4thcopy.xlsx"
wb = load_workbook(excel_file_path)
# Select the worksheet named "MLST" & "Pseudomonas_aeruginosa"
#ws_MLST=wb["MLST"]
ws_ANI=wb[ANI_sheet_name]

# LOOK FOR GENOMES WITH AN AVERAGE ANI VALUE BELOW 95
######################################################################################################################################
# Iterate over each row in the matrix except row 1 (header row) and the first values in each row (header values)
genomes_to_remove = {}
for row in range(2, ws_ANI.max_row + 1):
    # Look for rows where the average ANI value is below 0.95
    genome_name = ws_ANI.cell(row=row, column=1).value
    ani_values = [cell.value for cell in ws_ANI[row][2:]]  # Exclude the first cell
    avg_ANI = sum(ani_values) / len(ani_values)
    if avg_ANI < 95:
        # Add the genome name and average ANI value to a dictionary (key = genome name, value = average ANI value) if avg_ANI < 95
        genomes_to_remove[genome_name] = avg_ANI

# MAKE A NEW FILE WITH QC RESULTS OF THE GENOMES
######################################################################################################################################
# If the dictionary is empty, add a short statement to the file
# If the dictionary is not empty, write the genome names and average ANI values to the file
if not genomes_to_remove:
    with open("genomes_to_remove.txt", "w") as file:
         file.write("No genomes were found with an average ANI value below 95.")
else:
        with open("genomes_to_remove.txt", "w") as file:
            file.write("Genomes to remove \n(average ANI <95 is indicative of a different species than other genomes in the dataset):\n\n")
            for genome, avg_ANI in genomes_to_remove.items():
                 file.write(f"{genome}\t{avg_ANI}\n")