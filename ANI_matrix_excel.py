#!/usr/bin/env python
# Script to construct a full matrix from the FastANI matrix and writing it to an excel file
import xlsxwriter

# Opening & reading the FastANI matrix from the output file
file_path = input("Location to the FastANI matrix: ")
#file_path = '/home/guest/BIT11_Traineeship/FastANI/fastANI_output/fastANI_subset5.matrix'

with open(file_path) as file:
    lines = file.readlines()
    
# Extract the matrix values
fastani_matrix = []
for line in lines:
    row = line.strip().split()
    row.append('100')
    fastani_matrix.append(row)

# Remove the first row (contains amount of sequences)
del(fastani_matrix[0])

# Replace path to the genomic-fasta-file to the RefSeq-identifier (e.g. GCF_000168335.1)
for list in fastani_matrix:
    list[0] = list[0].split('/')[-2]
    print(list[0])
            
# Initialize Excel Workbook, & add a worksheet
species = input("What is the target species (_ as delimiter, for example Pseudomonas_aeruginosa): ")
file_name = 'FastANI_matrix_'+ species +'.xlsx'
wb = xlsxwriter.Workbook(file_name)
ws = wb.add_worksheet(species)

# Add each item of the matrix in a cell of the excel worksheet
for row_num, row_data in enumerate(fastani_matrix):
    for col_num, col_data in enumerate(row_data):
        ws.write(row_num, col_num, col_data)

# Save the workbook
wb.close()