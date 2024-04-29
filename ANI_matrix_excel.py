#!/usr/bin/env python
# Script to construct a full matrix from the FastANI matrix and writing it to an excel file
import xlsxwriter

#STEP 1 : adding the original FastANI matrix to an Excel file
############################################################################################################

# Opening & reading the FastANI matrix from the output file
#file_path = input("Location to the FastANI matrix: ")
file_path = "/home/guest/BIT11_Traineeship/FastANI/fastANI_output/fastANI_subset5.matrix"

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
file_name = "FastANI_matrix_2_"+ species +".xlsx"
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

# Auto adjust column-width
ws.autofit()

# Save the workbook
wb.close()
