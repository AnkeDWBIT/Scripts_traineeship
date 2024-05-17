import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from openpyxl import load_workbook

# AGGLOMERATIVE CLUSTERING OF ANI VALUES
# WRITE RESULTS TO EXCEL AND COMPARE PER GENOME IF ANI CLUSTER & MLST TYPE CORRESPOND
######################################################################################################################################
# Load the percent identity matrix from an Excel file into a dataframe
file_path = "/home/guest/BIT11_Traineeship/Scripts_traineeship/FastANI_matrix_Pseudomonas_aeruginosa (another copy).xlsx"
sheet_name = "Subset_ANI_matrix"
percent_identity_df = pd.read_excel(file_path, sheet_name=sheet_name, index_col=0)

# Convert the DataFrame to a Numpy array & Convert percent identity to distance (assuming values are in percentage form)
percent_identity_matrix = percent_identity_df.to_numpy()
distance_matrix = 100 - percent_identity_matrix

# Load the workbook where to find MLST values and store ANI cluster labels
excel_file_path = "/home/guest/BIT11_Traineeship/Scripts_traineeship/FastANI_matrix_Pseudomonas_aeruginosa (another copy).xlsx"
wb = load_workbook(excel_file_path)
# Select the worksheet named with MLST and ANI results
ws_MLST=wb["Subset_MLST_results"]
ws_ANI=wb[sheet_name]

# Perform agglomerative clustering
# Find the best value for the distance threshold by testing a range of values
col = 3
"""
for value in range(5, 0, -1):
    value = value/10
"""
values_range = [0.5, 0.4, 0.3, 0.2, 0.1, 0.05]
for value in values_range:
    clustering = AgglomerativeClustering(metric='precomputed', linkage='average', distance_threshold=value, n_clusters=None)
    clustering.fit(distance_matrix)
    # Get the cluster labels
    labels = clustering.labels_
    print("ANI cluster labels (using distance_threshold = ", value, ")", labels)
    # Add the labels from the Agglomerative clustering of ANI values to the ws_MLST worksheet staring form row 2, column 3
    header = "ANI cluster (distance_threshold = " + str(value) + ")"
    ws_MLST.cell(row=1, column=col, value=header)
    index = 2
    for label in labels:
        ws_MLST.cell(row=index, column=col, value=label)
        index += 1
    col += 1
 
# Save & close the workbook
wb.save(excel_file_path)
wb.close()