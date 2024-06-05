from shiny import App, render, ui, reactive
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from openpyxl import load_workbook

# Shiny for python script that asks user for input of an excel file and select the worksheets with ANI & MLST results
# User can specify a value for distance_threshold via the slider

# DEFINE THE UI OF THE APP
################################################################################################################
app_ui = ui.page_sidebar(
    # Sidebar where used can upload an Excel file and select a worksheet to use an input data
    ui.sidebar(
        ui.h2("Input file"),
        ui.input_file("file", "Upload Excel File", accept=".xlsx"),
        ui.h4("Select worksheets"),
        ui.input_select(  "select_ANI", "ANI matrix", choices=[],),
        ui.input_select(  "select_MLST", "MLST (metadata)", choices=[])
    ),
    # Main content area where the user can adjust the distance_threshold parameter & where agglomerative clustering is performed
    ui.card(
        ui.h2("Optimizing distance_threshold in AgglomerativeClustering of ANI values"),
        ui.input_slider("n", "Distance_threshold", min=0, max=1, value=0.3),
        ui.output_text_verbatim("txt")
))

# DEFINE THE SERVER LOGIC OF THE APP
################################################################################################################
def server(input, output, session):
    @reactive.Effect
    @reactive.event(input.file)

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

    @output
    # Adds text to the page that shows the distance_threshold value and the selected worksheets
    @render.text
    def txt():
        return f"Agglomerative clustering of ANI values with distance_threshold = {input.n()} \n Data retrieved from the Excel file on worksheets '{input.select_ANI()}' & '{input.select_MLST()}'"

# Create the Shiny app
app = App(app_ui, server)