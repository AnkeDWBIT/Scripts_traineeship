from shiny import App, render, ui, reactive
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from openpyxl import load_workbook
from collections import OrderedDict
from sklearn.metrics import adjusted_rand_score

# Shiny for python script that asks user for input of an excel file and select the worksheets with ANI & MLST results
# User can specify a value for distance_threshold via the slider
# Single-linkage agglomeratie clustering is performed
# Adjusted randis computed to test the congruence between the clustering and the MLST results

# DEFINE THE UI OF THE APP
################################################################################################################
app_ui = ui.page_sidebar(
    # Sidebar for user input (upload Excel file & select worksheets for input data)
    ui.sidebar(
        ui.h2("Input file"),
        ui.input_file("file", "Upload Excel File", accept=".xlsx"),
        ui.h4("Select worksheets"),
        ui.input_select("select_ANI", "ANI matrix", choices=[],),
        ui.input_select("select_MLST", "MLST (metadata)", choices=[]),
    ),
    # Main content area
    ui.card(
        ui.card(
            ui.h2("Optimizing distance_threshold in Agglomerative Clustering of ANI values"),
            ui.input_slider("n", "Distance_threshold", min=0, max=1, value=0.3),
            ui.input_action_button("go", "Perform Clustering", class_="btn-success"),
            ),
        # Output area where the user can see the results of the clustering
        ui.card(
            ui.h2("Congruence of partitioning methods (MLST & ANI clustering)"),
            ui.output_text("txt"),
            ),
        ui.card(
            ui.h2("Results table: Ordered MLST results with ANI cluster labels"),
            ui.output_table("ordered_MLST_table"),
        )))

# DEFINE THE SERVER LOGIC OF THE APP
################################################################################################################
def server(input, output, session):

    @reactive.Effect
    @reactive.event(input.file) # Event that triggers the function when the file is uploaded
    # Function that opens the Excel file and finds the worksheets in it, used for the select input
    def update_sheet_choices():
        file_info = input.file()
        if file_info:
            # Extract the actual file path
            file_data = file_info[0]['datapath']
            # Read the Excel file and get the sheet names
            xls = pd.ExcelFile(file_data)
            worksheets = xls.sheet_names
            # Update the choices of the select input for sheets
            ui.update_select("select_ANI", choices=worksheets)
            ui.update_select("select_MLST", choices=worksheets)

    # Process data from the excel file (agglomerative clustering, combine results with metadata)
    @reactive.Effect
    @reactive.event(input.go) # Event that triggers the function when the button is clicked
    # Function that performs agglomerative clustering on the ANI values
    def process_data():
        file_info = input.file()
        ANI_sheet_name = input.select_ANI()
        MLST_sheet_name = input.select_MLST()
        
        if file_info and ANI_sheet_name and MLST_sheet_name:
            # Load the percent identity matrix from an Excel file into a dataframe
            df = pd.read_excel(file_info[0]['datapath'], sheet_name=ANI_sheet_name, engine='openpyxl', header=0, index_col=0)
            # Convert the DataFrame to a Numpy array & Convert percent identity to distance (assuming values are in percentage form)
            percent_identity_matrix = df.to_numpy()
            distance_matrix = 100 - percent_identity_matrix

            # Perform agglomerative clustering
            dist_thres = input.n()
            clustering = AgglomerativeClustering(metric='precomputed', linkage='average', distance_threshold=dist_thres, n_clusters=None)
            clustering.fit(distance_matrix)
            
            # Get the cluster labels
            labels = clustering.labels_

        # Make an output table with the MLST results and ANI cluster labels, ordered as the genomes in the ANI matrix
        @output
        @render.table
        def ordered_MLST_table():
            file_info = input.file()
            ANI_sheet_name = input.select_ANI()
            MLST_sheet_name = input.select_MLST()
            
            if file_info and ANI_sheet_name and MLST_sheet_name:
                # Load the data from the Excel file
                wb = load_workbook(file_info[0]['datapath'])
                ws_MLST=wb[MLST_sheet_name]
                ws_ANI=wb[ANI_sheet_name]
                # Store the genomes of the ANI sheet in a list
                ANI_genomes = []
                for row in range(2, ws_ANI.max_row + 1):
                    ANI_genomes.append(ws_ANI.cell(row=row, column=1).value)
                # Store the genomes and STs of the MLST sheet in a dictionary
                MLST_dict = {}
                for row in range(2, ws_MLST.max_row + 1):
                    MLST_dict[ws_MLST.cell(row=row, column=1).value] = ws_MLST.cell(row=row, column=2).value
                # Create an ordered dictionary (MLST dictionary in order of ANI genomes list)
                ordered_MLST = OrderedDict((key, MLST_dict[key]) for key in ANI_genomes)
                # Convert ordered dictionary to a DataFrame for display
                df_MLST = pd.DataFrame(list(ordered_MLST.items()), columns=['Genome', 'ST'])
                # Add the cluster labels to the DataFrame
                df_MLST['Cluster'] = labels
                return df_MLST

            # Add a text output with the adjusted rand value to the output area               
            @output
            @render.text
            def txt():  
                adjusted_rand = adjusted_rand_score(df_MLST['ST'], labels) # Calculate the adjusted rand value (congruence between two typin methods)
                return f"Agglomerative clustering of ANI values with distance_threshold = {input.n()} \n Data retrieved from the Excel file on worksheets '{input.select_ANI()}' & '{input.select_MLST()} and clustering results summarized in the table \n Clustering congruence with MLST: {adjusted_rand}'"

# Create the Shiny app
app = App(app_ui, server)