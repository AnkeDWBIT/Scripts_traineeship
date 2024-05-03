#!/usr/bin/env python
# Script to write ANI & MLST results to excel
# First worksheet is a full ANI matrix, second worksheet contains the strain types assigned by MLST
import xlsxwriter
import os
import json
import re

#STEP 1 : adding the original FastANI matrix to an Excel file
############################################################################################################

# Opening & reading the FastANI matrix from the output file
#file_path = input("Location to the FastANI matrix: ")
file_path = "/home/guest/BIT11_Traineeship/03_FastANI/fastANI_subset5.matrix"

with open(file_path) as file:
    lines = file.readlines()
    
# Extract the matrix values
fastani_matrix = []
for line in lines:
    row = line.strip().split()
    row.append("100")
    fastani_matrix.append(row)

# Remove the first row (contains amount of sequences)
del(fastani_matrix[0])

# Replace path to the genomic-fasta-file to the RefSeq-identifier (e.g. GCF_000168335.1)
for list in fastani_matrix:
    list[0] = list[0].split('/')[-2]
            
# Initialize Excel Workbook, & add a worksheet
#species = input("What is the target species (_ as delimiter, for example Pseudomonas_aeruginosa): ")
species = "Pseudomonas_aeruginosa"
file_name = "FastANI_matrix_"+ species +".xlsx"
wb = xlsxwriter.Workbook(file_name)
ws = wb.add_worksheet(species)

# Add each item of the matrix in a cell of the excel worksheet
for row_num, row_data in enumerate(fastani_matrix):
    for col_num, col_data in enumerate(row_data):
        ws.write(row_num+1, col_num, col_data)



#STEP 2 : mirroring the martrix diagonally to construc a full matrix
############################################################################################################

# Only retain ANI-values in the lists (= remove RefSeq-ID and value 100)
ani_values = [list[1:-1] for list in fastani_matrix]

# Add the mirrored ANI-values to the matrix
startcol = 1
for col_num, col_data in enumerate(ani_values):
    for row_num, row_data in enumerate(col_data):
        ws.write(row_num+1, col_num + startcol, row_data)

 # Add the RefSeq-IDs as column names in the first row
RefSeq_IDs = []
for list in fastani_matrix:
    RefSeq_IDs.append(list[0])  
for col_num, refseq_id in enumerate(RefSeq_IDs, start=1):
    ws.write(0, col_num, refseq_id)
    

# STEP 3 : add strain types (ST) assigned by the MLST tool to a new worksheet in the excel file
############################################################################################################
ws2 = wb.add_worksheet("MLST")  

#mlst_results_dir = input("Location to the JSON files from MLST analysis: ")
mlst_results_dir = "/home/guest/BIT11_Traineeship/04_MLST/mlst/"

# Write column headers
ws2.write(0, 1, "ST")  # Write ST at row 0, column 1

# Track the current row for writing data
current_row = 1

for file in os.listdir(mlst_results_dir):
    if "MLST_GCF" in file:
        result_file_path = os.path.join(mlst_results_dir, file)
        # Regular expression pattern to match the RefSeq identifier
        pattern = r'GCF_\d+\.\d+'
         # Search for RefSeq identifier in the filename
        genome_refseq_match = re.search(pattern, file)
        # Extract genome_refseq if found
        if genome_refseq_match:
            genome_refseq = genome_refseq_match.group()
        # Opening JSON file
        with open(result_file_path) as f:
            data = json.load(f)
        # Retrieve the sequence type (ST)
        ST = data['mlst']['results']['sequence_type']
        # Write data to worksheet
        ws2.write(current_row, 0, genome_refseq)
        ws2.write(current_row, 1, ST)
        # Increment row counter for the next entry    
        current_row += 1
        # Close the JSON file  
        f.close()


# Auto adjust column-width
ws.autofit()

# Save the workbook
wb.close()
