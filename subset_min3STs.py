#!/usr/bin/env python
# Script to make new ANI & MLST excel worksheets that only contain genomes having and ST that occurs at least 3 times in the dataset
from openpyxl import load_workbook
import sys
import itertools

# STEP 1 : HANDLE COMMAND-LINE ARGUMENTS
############################################################################################################
# Check if any command-line arguments have been provided
if len(sys.argv) < 3:
    print("Error: Need to provide the correct command-line arguments.")
    print("Usage: python scriptname.py [1] [2] [3]")
    print("\t[1] = Full path to Excel file with ANI matrix and MLST results")
    print("\t[2] = Worksheet with ANI matrix in Excel file [1] e.g. Pseudomonas_aeruginosa or ANI_matrix_reviewed")
    print("\t[3] = Worksheet with MLST results in Excel file [1] e.g. MLST or MLST_reviewed")
    sys.exit(1)

# Store command-line arguments
excel_file_path = sys.argv[1]
ANI_sheet_name = sys.argv[2]
MLST_sheet_name = sys.argv[3]

# STEP 2 : PROCESS EXCEL FILE
############################################################################################################
# Load the workbook & some of it's worksheets
#excel_file_path = "/home/guest/BIT11_Traineeship/Scripts_traineeship/FastANI_matrix_Pseudomonas_aeruginosa_results.xlsx"
wb = load_workbook(excel_file_path)
ws_MLST=wb[MLST_sheet_name]
ws_ANI=wb[ANI_sheet_name]

# Make a dictionary with the STs as key and genomes as values
ST_dict = {}
for row in range(2, ws_MLST.max_row + 1):
    ST = ws_MLST.cell(row=row, column=2).value
    genome = ws_MLST.cell(row=row, column=1).value
    if ST in ST_dict:
        ST_dict[ST].append(genome)
    else:
        ST_dict[ST] = [genome]

# Reduce the dictionary to only the genomes with STs that occur at least 3 times in the dataset
ST_dict = {key: value for key, value in ST_dict.items() if len(value) >= 3} 

# Make a list of all genomes (values) in the dictionary
genomes = []
for value in ST_dict.values():
    genomes.extend(value)

# STEP 3 : MAKE NEW WORKSHEETS IN THE EXCEL FILE WITH REDUCED DATA
############################################################################################################
# Make a new worksheet for the subset ANI matrix & MLST results
ws_ANI_subset = wb.create_sheet("Subset_ANI_matrix_3STs")
ws_MLST_subset = wb.create_sheet("Subset_MLST_3STs")
# Write column headers to the new MLST worksheet
ws_MLST_subset.cell(row=1, column=1, value="RefSeq ID")
ws_MLST_subset.cell(row=1, column=2, value="ST")

# Write RefSeq identifiers from the genomes list to the new worksheets & write the STs to the ws_MLST_subset worksheet
GCF_index = 0
index = 2
for genome in genomes:
    # Write the RefSeq ID as column and row headers to the ws_ANI_subset worksheet
    ws_ANI_subset.cell(row=index, column=1, value=genome)
    ws_ANI_subset.cell(row=1, column=index, value=genome)
    # Write the RefSeq ID as row header to the ws_MLST_subset worksheet
    ws_MLST_subset.cell(row=index, column=1, value=genome)
    # Write the STs of the genomes in the GCF list to the ws_MLST_subset worksheet
    for row in range(1, ws_MLST.max_row + 1):
        genome_MLST = ws_MLST.cell(row=row, column=1).value
        if genome_MLST == genome:
            ST = ws_MLST.cell(row=row, column=2).value
            if genome_MLST in genomes:
                ws_MLST_subset.cell(row=index, column=2, value=ST)
                GCF_index += 1
                index += 1

# Make all possible combinations of the genomes in the list
for genome_pair in itertools.product(genomes, genomes):
    # Extract the column- & row index of the genomes in the ANI matrix
    for row in range(2, ws_ANI.max_row+1):
            if (ws_ANI.cell(row, 1).value) == genome_pair[0]:
                row_index = row
    for col in range(2, ws_ANI.max_row+1):
            if (ws_ANI.cell(1, col).value) == genome_pair[1]:
                col_index = col           
    # Find the ANI-value of the genome pair in the Excel file
    ANI_value = ws_ANI.cell(row_index, col_index).value
    #print(genome_pair[0], genome_pair[1], ANI_value)
    # Extract the column - & row index of the genomes in the subset ANI matrix
    for row in range(2, ws_ANI_subset.max_row+1):
        if ws_ANI_subset.cell(row, 1).value == genome_pair[0]:
            row_index_subset = row
    for col in range(2, ws_ANI_subset.max_row+1):
        if ws_ANI_subset.cell(1, col).value == genome_pair[1]:
            col_index_subset = col
    # Write the ANI-value to the new worksheet
    ws_ANI_subset.cell(row_index_subset, col_index_subset, value=ANI_value)
    

# Save & close the workbook
wb.save(excel_file_path)
wb.close()

