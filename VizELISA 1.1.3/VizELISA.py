"""
For consideration:
    List of LM antibodies
    https://plantcellwalls.leeds.ac.uk/wp-content/uploads/sites/103/2021/11/JPKab2021.pdf

"""





import customtkinter as ctk
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import random
from scipy import stats
from scipy.stats import f_oneway
import seaborn as sns
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import subprocess
import sys
import time
import tkinter

##########################    O U T P U T    F O L D E R    ###################

current_timestamp = datetime.now().strftime('%Y-%m-%d_%Hh%Mm%Ss')
out_folder_name = current_timestamp



########################## G L O B A L    V A R I A B L E S ###################

A1 = None       #
final_DPI = 600 
heatmap_dimensions = 17,10   #
outlier_label = "outlier" 
extraction_count = 2 
extract1, extract2, extract3 = "", "", ""  
extraction_labels = []  
percent_blank_cut_off_threshold = 20 # only values that are x% higher than blank are kept
P_01, P_05 = "**", "*"      #
graph_title, x_title, y_title = "", "", ""  #
fontsize_title, fontsize_axis = 17, 11      #
min_cutoff = 5
x_rotation, y_rotation = 0, 0
colour_scheme = "Greens"
image_save_as = "ELISA heatmap"
auto_outlier_removal = True
use_median = False
export_stats = False    #
generate_barplot = False    #
generate_linegraph = False  #
use_antibody_colours = False
bold_graph_title = False
bold_axis_titles = False
DEBUGGING = False

barplot_bar_count = 2






###################### DATA PROCESSING PREFERENCES ############################
per_extraction_normalisation = False     # Set this value to True to normalise the heatmaps for each extraction individualy
normalise = True # Impacts bar graphs. Set this value to False to plot unnormalised values on bar graphs
multiple_templates = True
discriminate_by_extraction = True
remove_extraction_label = True      # This will remove the extraction label prefix as given by user under the variable extraction_labels
lower_cut_off_threshold = 5 # Eliminates heatmap values below this value after normalisation


############################ GLOBAL VARIABLES #################################
global_heatmap = 0 # This is reserved for a fully formatted heatmap before exporting it as an image
subplots_processed = 0 # This is a place holder for a dictionary (each key is the extraction label corresponding to a heatmap matrix)
heatmap = pd.DataFrame()
heatmap_ordered = pd.DataFrame()
stdev = {} # The dictionary below will keep track of the stdev between indicated replicates
antibodies = [
    "Partially methyl esterified HG (JIM5)",
    "Methyl esterified HG (JIM7)",
    "Partially methyl esterified HG (LM18)",
    "Unesterified HG (LM19)",
    "Methyl esterified HG (LM20)",
    "Partially methyl esterified HG (LM7)",
    "Backbone of RG-I (RU2)",
    "Backbone of RG-I (RU1)",
    "Beta-1,4-galactan (LM5)",
    "Alpha-1,5-arabinan (LM6)",
    "Galactosylated RG-I (LM16)",
    "Beta-1,4-mannan (LM21)",
    "Xyloglucan (XXXG motif) (LM15)",
    "Xyloglucan (XXLG & XLLG motifs) (LM24)",
    "Xyloglucan (XXLG & XLLG motifs) (LM25)",
    "Beta-1,4-xylan (LM11)",
    "Xylosyl residues (LM23)",
    "Arabinogalactan-protein (AGP) (JIM8)",
    "Arabinogalactan-protein (AGP) (JIM13)",
    "Arabinogalactan-protein (AGP) (JIM15)",
    "Galactosyl residue of AGP (JIM16)",
    "Glucuronosyl residue of AGP (LM2)",
    "(1-3)-beta-D-glucan (BS-400-2)",
    "(1-3; 1-4)-beta-glucans (BS-400-3)",
    "Cellulose (crystalline) (CBM3a)"
    ]


polymer_classes = ["Pectin", "AGP", "Hemicellulose", "Cellulose and ß-glucan", "Cellulose and Hemicellulose"]




mapping_file_name = "_mapping_template"









class MyApp(ctk.CTk):
    def __init__(self):
        super().__init__()
        
        
        
        self.geometry("800x600")
        self.title("VizELISA 1.1.3")

        self.rowconfigure(0, weight=1)  # Top section housing Data Input, Quick Parameters and Graph Specifications
        self.rowconfigure(1, weight=1)  # Bottom section (Python's standard output will be redirected here)

        # Top half divided into 3 sections
        self.top_frame = ctk.CTkFrame(self, corner_radius=0)
        self.top_frame.pack(side="top", fill="both", expand=True)

        
        # TEXT INPUT DICTIONARY - houses the text inputs for both DATA INPUT and GRAPH SPECS
        self.text_inputs = {}



        ###############################################
        ###############################################
        ###         D A T A    I N P U T            ###
        ###############################################
        ###############################################

        # Section 1: Scrollable frame with buttons and checkboxes - DATA INPUT
        self.button_frame = ctk.CTkScrollableFrame(self.top_frame)
        self.button_frame.pack(side="left", fill="both", expand=True)
        ctk.CTkLabel(self.button_frame, text="Data Input", font=("Arial", 20)).pack(pady=10)
        # Add buttons and checkboxes here
        
        data_input_buttons = ["Open data mapping file", 
                         "Check data validity", 
                         "Antibody list", 
                         "Polymer class groupings",
                         "Exclude antibodies",
                         "Barplot sample groups",
                         "Bulk file renaming",
                         "List of valid colours",
                         "INSTRUCTIONS"]
        
        for i in range(len(data_input_buttons)):
            btn = ctk.CTkButton(self.button_frame, text=data_input_buttons[i], command=lambda i=i: self.button_click(i, data_input_buttons[i]))
            btn.pack(fill="x")
            

        # Text input checkboxes under DATA INPUT section            
        data_input_checkboxes = ["Global normalisation", 
                                 "Local normalisation", 
                                 "Remove extraction labels"]
        self.checkbox_states = [False, False, False]



        
        for i in range(len(data_input_checkboxes)):
            ctk.CTkCheckBox(self.button_frame, text=data_input_checkboxes[i], command=lambda i=i: self.checkbox_click_data_input(i, data_input_checkboxes[i])).pack(fill="x")

        # Data input text fields
        data_inputs = {"A1 position on 96-well plate ": "A1",
                       "Plate outlier label ": "outlier",
                       "Number of extractions ": "2", 
                       "Extraction label 1 ": "C", 
                       "Extraction label 2 ": "N", 
                       "Extraction label 3 ": "",
                       "% of blank cut-off ": 20,
                       "Remove values below ": 5,
                       "P<0.05 label ": "*",
                       "P<0.01 label ": "**"}
        
        for field in data_inputs:
            frame = ctk.CTkFrame(self.button_frame)
            frame.pack(fill="x")
            ctk.CTkLabel(frame, text=field).pack(side="left")
            entry = ctk.CTkEntry(frame)
            entry.insert(0, data_inputs[field]) ###################################################### default values
            entry.pack(side="right", fill="x")
            self.text_inputs[field] = entry






        ###############################################
        ###############################################
        ###        G R A P H    S P E C S           ###
        ###############################################
        ###############################################
        

        # Section 2: Scrollable frame with text input fields - GRAPH SPECIFICATIONS
        self.textinput_frame = ctk.CTkScrollableFrame(self.top_frame)
        self.textinput_frame.pack(side="right", fill="both", expand=True)
        ctk.CTkLabel(self.textinput_frame, text="Graph Specifications", font=("Arial", 20)).pack(pady=10)

        # Add text input fields here
        graph_inputs = {"Graph title ": "",
                        "Fig size (x,y) ": "(17,10)",
                        "X-axis title ": "Monoclonal antibodies", 
                        "Y-axis title ": "", 
                        "Title font size ": "17",
                        "Axis font size ": "11", 
                        "X-label rotation (°) ": 30,
                        "Y-label rotation (°) ": 0,
                        "Colour scheme ": "Greens",
                        "BARPLOT #pairs": 2,
                        "Save image as ": "ELISA heatmap", 
                        "Image resolution ": 500
                        }
        
        for field in graph_inputs:
            frame = ctk.CTkFrame(self.textinput_frame)
            frame.pack(fill="x")
            ctk.CTkLabel(frame, text=field).pack(side="left")
            entry = ctk.CTkEntry(frame)
            entry.insert(0, graph_inputs[field]) ###################################################### default values
            entry.pack(side="right", fill="x")
            self.text_inputs[field] = entry


        # generate_low_res_img button
        self.confirm_button = ctk.CTkButton(self.textinput_frame, text="Confirm Parameters", command=self.confirm_parameters)
        self.confirm_button.pack(fill="x")

        # Generate a button for bar plot only
        self.confirm_button = ctk.CTkButton(self.textinput_frame, text="Generate barplots", command=self.barplot)
        self.confirm_button.pack(fill="x")







        ###############################################
        ###############################################
        ###    Q U I C K    P A R A M E T E R S     ###
        ###############################################
        ###############################################
        
        quick_checkboxes = [
            "Automatic outlier removal",
            "Use median instead of mean",
            "Export descriptive statistics",
            "Generate bar plot",
            "Generate correlation matrix", #TODO
            "Generate line graph",
            "Colour code antibodies by class",
            "Bold graph title",
            "Bold axis titles",
            "DEBUGGING MODE"
            
            ]
        
        
        # Section 3: Scrollable frame with checkboxes - QUICK PARAMETERS
        self.checkbox_frame = ctk.CTkScrollableFrame(self.top_frame)
        self.checkbox_frame.pack(side="right", fill="both", expand=True)
        ctk.CTkLabel(self.checkbox_frame, text="Quick Parameters", font=("Arial", 20)).pack(pady=10)

        # Add checkboxes here
        self.checkbox_states_QP = [False] * 10
        for i in range(len(quick_checkboxes)):
            btn = ctk.CTkCheckBox(self.checkbox_frame, text=quick_checkboxes[i], command=lambda i=i: self.checkbox_click_QP(i, quick_checkboxes[i]))
            btn.pack(fill="x")
            
        
        
        # Confirm button
        self.confirm_button = ctk.CTkButton(self.checkbox_frame, text="Load Files", command=self.load_files)
        self.confirm_button.pack(fill="x")





        ###############################################
        ###############################################
        ###          S T D     O U T P U T          ###
        ###############################################
        ###############################################

        # Redirect stdout to the stdout frame
        self.stdout_frame = ctk.CTkScrollableFrame(self, corner_radius=0, height = 10)
        self.stdout_frame.config(bg="#333333")
        self.stdout_frame.pack(side="bottom", fill="x")
        self.stdout_text = tkinter.Text(self.stdout_frame, font=("Arial", 12), height=10, fg="white", bg="#333333")
        self.stdout_text.pack(fill="both", expand=True)
        sys.stdout = self

        # QUIT button =========================================================
        self.quit_button = ctk.CTkButton(self.stdout_frame, text="GENERATE FINAL FIGURE", command=self.generate_final_figure)
        self.quit_button.pack(fill="x", anchor="s")
        self.quit_button.focus_set()
        
        
        
        
        
        
        

        


    def button_click(self, i, button_label):
        
        if button_label == "Open data mapping file":
            print("Opening mapping file...")
            os.system("_mapping_template.xlsx")
            
        elif button_label == "Check data validity":
            self.load_files()
            
        elif button_label == "Antibody list":
            os.system("antibody_list.txt")
            
        elif button_label == "Polymer class groupings":
            print("polly classes")
            os.startfile("barplot_groups.xlsx")
            
        elif button_label == "Exclude antibodies":
            print("Opening exclusion list...")
            
        elif button_label == "Barplot sample groups":
            os.system("barplot samples.xlsx") #TODO
            
        elif button_label == "Samples":
            print(6969)
            
        elif button_label == "Bulk file renaming":
            print(420420)
            
        elif button_label == "INSTRUCTIONS":
            print("Opening instructions...")
            os.system("Instructions.txt")
            
        elif button_label == "List of valid colours":
            print("Opening list of colours...")
            os.system("colours.txt")
        

        





    def checkbox_click_data_input(self, i, label_DI):

        global global_normalisation, local_normalisation, remove_extraction_label



        self.checkbox_states[i] = not self.checkbox_states[i]

        if self.checkbox_states[i] == True:        
            print(f"{label_DI}:    enabled")
        elif self.checkbox_states[i] == False:
            print(f"{label_DI}:    disabled")

        #print(label_DI, "=", self.checkbox_states[i])
        
        
        if label_DI == "Global normalisation":
            global_normalisation = truefalse(self.checkbox_states[i], True, False)

        elif label_DI == "Local normalisation":
            local_normalisation = truefalse(self.checkbox_states[i], True, False)
        elif label_DI == "Remove extraction labels":
            remove_extraction_label = truefalse(self.checkbox_states[i], True, False)




    def checkbox_click_QP(self, i, checkbox_name_QP):
        """
        Parameters
        ----------
        i : int
            Number corresponding to position of checkbox.
        checkbox_name_QP : string
            Name of checkbox found under Quick Parameters.

        Returns
        -------
        None.

        """
        
        
        self.checkbox_states_QP[i] = not self.checkbox_states_QP[i]
        # print(f"Checkbox {i} toggled", checkbox_name, "ENABLED")
        if self.checkbox_states_QP[i]:
            print(checkbox_name_QP + ":    enabled")
        elif not self.checkbox_states_QP[i]:
            print(checkbox_name_QP + ":    disabled")
        
        global auto_outlier_removal, use_median, export_stats, generate_barplot, generate_linegraph, use_antibody_colours, bold_graph_title, bold_axis_titles, DEBUGGING
        
  
        if checkbox_name_QP == "Automatic outlier removal":
            auto_outlier_removal = truefalse(self.checkbox_states_QP[i], True, False)
        if checkbox_name_QP == "Use median instead of mean":
            use_median = truefalse(self.checkbox_states_QP[i], True, False)
        elif checkbox_name_QP == "Export descriptive statistics":
            export_stats = truefalse(self.checkbox_states_QP[i], True, False)
        elif checkbox_name_QP == "Generate bar plot":
            generate_barplot = truefalse(self.checkbox_states_QP[i], True, False)
        elif checkbox_name_QP == "Generate line graph":
            generate_linegraph = truefalse(self.checkbox_states_QP[i], True, False)
        elif checkbox_name_QP == "Colour code antibodies by class":
            use_antibody_colours = truefalse(self.checkbox_states_QP[i], True, False)
        elif checkbox_name_QP == "Bold graph title":
            bold_graph_title = truefalse(self.checkbox_states_QP[i], True, False)
        elif checkbox_name_QP == "Bold axis titles":
            bold_axis_titles = truefalse(self.checkbox_states_QP[i], True, False)
        elif checkbox_name_QP == "DEBUGGING MODE":
            DEBUGGING = truefalse(self.checkbox_states_QP[i], True, False)
        
        
        
        

    def generate_final_figure(self):
        # Top part is only for testing purposes
        global global_heatmap
        if not global_heatmap:
            print("ERROR - Please confirm parameters before final figure can be generated!")
        
        if global_heatmap:
            global_heatmap.figure.savefig(image_save_as+".png", dpi=final_DPI, bbox_inches = "tight")
            os.system("TEMP.png")
            
    def barplot(self):
        self.update_global_input_variables()
        initiate_barplots()
        
    def confirm_parameters(self):
        
        # text_inputs is a dictionary in the format {text_field: user_provided_parameter}
        text_inputs = self.get_text_inputs()
        self.update_global_input_variables()
        

        
        if DEBUGGING:
            for key in text_inputs:
                print(key, "-", text_inputs[key])
            if text_inputs["Graph title "]:
                print("Override default graph title!")
        
        if not subplots_processed:
            print(":(")
            print("Error - Please load data files")
        else:
            create_heatmap(subplots)
            print("\n\n")
            print(">>> DONE <<<\nAdjust settings under Graph Specifications and confirm parameters again\n\nor\n\nGenerate final high-resolution image with Generate Final Figure")

        
    def load_files(self):  
        
        # Here we create our output folder for this analysis
        global out_folder_name
        os.mkdir(out_folder_name)
        
        
        
        self.update_global_input_variables()
        
        template_list = []
        for file in os.listdir():
            if file[:14] == "ELISA_template":
                template_list.append(file)
                
        
        global number_of_templates
        number_of_templates = 2 #TODO: get input from user
        for sheetnumber in range(number_of_templates):
            template = "_mapping_template.xlsx"
            coords = get_coordinates(template, sheetnumber)
            extract_data(coords, sheetnumber+1)
            
            
        
        
        """
        
        The code below allows us to separate our data by extraction in addition
        to ELISA plate number. The user indicates extraction labels. For each
        extraction label, the samples are scanned for the extraction prefix and
        is counted. Each extraction label is given a count number, which will
        be used to separate rows into subplots according to extraction,
        UNDER THE ASSUMPTION THAT:
            a) The samples are in order of extractions as indicated by the user
            b) The samples from each extraction are grouped together
        
        """

        if discriminate_by_extraction:
            row_labels = heatmap_ordered.index
            global extract_count
            extract_count = {}
            global subplots
            subplots = {}
            running_row_count = 0
            for extract in extraction_labels:
                counter = 0
                extract_name_length = len(extract)
                
                for label in row_labels:
                    if label[:extract_name_length] == extract:
                        counter += 1
                extract_count[extract] = counter
                
                
                subplots[extract] = heatmap_ordered.iloc[running_row_count : running_row_count + counter]
                running_row_count += counter
            
            
            
            # Here we remove the extraction label if the user indicates so
            for plot in subplots:
                for column in subplots[plot]:
                    for cell_name, value in subplots[plot][column].items():
                        for extract in extraction_labels:
                            if cell_name[:len(extract)] == extract and remove_extraction_label:
                                subplots[plot] = subplots[plot].rename(index = {cell_name: cell_name[len(extract):]})
            
            
            
            
                if per_extraction_normalisation and normalise: # The code below will allow for the normalisation of each subplot individually if the user sets global_normalisation to False
                        for column in subplots[plot]:
                            print(subplots[plot])
                            col_max = -1
                            for row_name, cell_value in subplots[plot][column].items():
                                debugging = True
                                if debugging:
                                    print("cell value: ", cell_value)
                                    print("row name: ", row_name)
                                if cell_value > col_max:
                                    col_max = cell_value
                            if col_max > 0:
                                multiplier = 100 / col_max
                            else:
                                multiplier = 0
                            
                                
                            for row_name, cell_value in subplots[plot][column].items():
                                print(row_name)
                                subplots[plot][column][row_name] *= multiplier
            
            
            
            ##################################################################################################################################
            global subplots_processed
            subplots_processed = subplots
            #create_heatmap(subplots)
            
            
            print("\n\n\n>>> DONE <<<\nScroll up for details")
            print("Proceed to Graph Specifications and press Confirm Parameters when done")
            
            
            
            
            
            # Here we update the variables supplied by the user via text fields
            self.update_global_input_variables()
            
            
            
    def get_text_inputs(self):
        return {field: entry.get() for field, entry in self.text_inputs.items()}

    def on_closing(self):
        self.destroy()

    
    def open_data_mapping_file(self):
        os.system("mapping_file_name")

    def write(self, text):
        self.stdout_text.insert("end", text + "\n")
        self.stdout_text.see("end")

    def flush(self):
        pass
    
    def run_file(self):
        for i in range(3):
            print(random.randint(1,1000))
        subprocess.run(["python", "test_file.py"])
        
        
    def update_global_input_variables(self):
        
        
        global A1, outlier_label, extraction_count, extract1, extract2, extract3, extraction_labels, graph_title, x_title, y_title, fontsize_title, fontsize_axis, min_cutoff, x_rotation, y_rotation, colour_scheme, image_save_as, P_01, P_05, percent_blank_cut_off_threshold, final_DPI, heatmap_dimensions, barplot_bar_count
        
        
        text_inputs = self.get_text_inputs()



        A1 = text_inputs["A1 position on 96-well plate "]
        outlier_label = text_inputs["Plate outlier label "]
        heatmap_dimensions = text_inputs["Fig size (x,y) "]
        heatmap_dimensions = eval(heatmap_dimensions)
        
        
        extraction_count = int(text_inputs["Number of extractions "])
        extract1 = text_inputs["Extraction label 1 "]
        extract2 = text_inputs["Extraction label 2 "]
        extract3 = text_inputs["Extraction label 3 "]
        extraction_labels = [extract1]
        if len(extract2) > 0:  extraction_labels += [extract2]
        if len(extract3) > 0:  extraction_labels += [extract3]
        
        percent_blank_cut_off_threshold = int(text_inputs["% of blank cut-off "])
        P_05 = text_inputs["P<0.05 label "]
        P_01 = text_inputs["P<0.01 label "]
        graph_title = text_inputs["Graph title "]
        x_title = text_inputs["X-axis title "]
        y_title = text_inputs["Y-axis title "]
        fontsize_title = int(text_inputs["Title font size "])
        fontsize_axis = int(text_inputs["Axis font size "])
        min_cutoff = int(text_inputs["Remove values below "])
        x_rotation = int(text_inputs["X-label rotation (°) "])
        y_rotation = int(text_inputs["Y-label rotation (°) "])
        colour_scheme = text_inputs["Colour scheme "]
        final_DPI = int(text_inputs["Image resolution "])
        image_save_as = text_inputs["Save image as "]
        
        barplot_bar_count = int(text_inputs["BARPLOT #pairs"])


logo = """
.::         .::            .::::::::.::      .::  .:: ::        .:       
 .::       .::  .:         .::      .::      .::.::    .::     .: ::     
  .::     .::     .:::: .::.::      .::      .:: .::          .:  .::    
   .::   .::   .::     .:: .::::::  .::      .::   .::       .::   .::   
    .:: .::    .::   .::   .::      .::      .::      .::   .:::::: .::  
     .::::     .::  .::    .::      .::      .::.::    .:: .::       .:: 
      .::      .::.::::::::.::::::::.::::::::.::  .:: ::  .::         .::    version 1.1.3"""



def bar_plot_3_groups_extractions_combined(sample_class, polymer_class, dataframe_subplots, group_ranges:list):
    save_barplots = False
    
    combine_extractions_barplot = True
    
    display_full_probe_name = True
    display_additional_probe_code = True # This will display the antibody code underneath the target name on the x-axis

    if True: # Indicates which polymers to plot on bar plots based on polymer class supplied to this function
        if polymer_class.lower() == "pectin":
            polymers = [ 
                "JIM5",
                "JIM7",
                "LM18",
                "LM19",
                "LM20",
                "LM7",
                "RU2",
                "RU1",
                "LM5",
                "LM6",
                "LM16"]
            
        elif polymer_class.lower() == "pectin_for_viz_figure":
            polymers = [ 
                
                "JIM7",
                "LM18",
                "LM20",
                
                
                "RU1",
                "LM5",
                "LM6",
                ]
            
        elif polymer_class.lower() == "agp":
            polymers = [
                "JIM8",
                "JIM13",
                "JIM15",
                "JIM16",
                "LM2"]
        elif polymer_class.lower() == "hemicellulose":
            polymers = [
                "LM21",
                "LM15",
                "LM24",
                "LM25",
                "LM11",
                "LM23"]
        elif polymer_class.lower() == "cellulose and ß-glucan":
            polymers = [
                "BS-400-2",
                "BS-400-3",
                "CBM3a"]
        elif polymer_class.lower() == "cellulose and hemicellulose":
            polymers = [
                "CBM3a",
                "BS-400-2",
                "BS-400-3",           
                "LM21",
                "LM15",
                "LM24",
                "LM25",
                "LM11",
                "LM23"]
        elif polymer_class.lower() == "hg & rg-i":
            polymers = ["JIM5", "LM16"]
        elif polymer_class.lower() == "lm25+lm15":
            polymers = ["LM15", "LM25"]
        else:
            print("INVALID POLYMER CLASS!")
            print(polymer_class, "not found")
            sys.exit()

    
    ##########################################################################
    #                   D A T A     S E T T I N G S                          #
    ##########################################################################


    # This will be used to plot a varable number of probes on a single plot without needing manual adjustment
    index = np.arange(len(polymers))
    
    group1_range = group_ranges[0]
    group2_range = group_ranges[1]
    group3_range = group_ranges[2]
    
    
    group1_data = []
    group1_error = []
    group1_negative_data = []
    group1_negative_error = []
    
    group2_data = []
    group2_error = []
    group2_negative_data = []
    group2_negative_error = []
    
    
    group3_data = []
    group3_error = []
    group3_negative_data = []
    group3_negative_error = []



    # TUKEY'S HSD TESTING
    CDTA_tukey_groups = []
    NaOH_tukey_groups = []
    
    
    
    
    for df_no, dataframe in enumerate(dataframe_subplots):
        
        if display_full_probe_name:
            probe_names = []
        else:
            probe_names = polymers[:]
        
        
        
        for column in polymers:  # Use the alternative version of this loop-conditional pair if anything goes wrong with sequence and labels on bar plots
            if column in dataframe_subplots[dataframe]:
                
        # for column in dataframe: ###### REUSE THIS if anything is broken with bar plots and order
        #     if column in polymers:
                
            
                data_1 = []
                data_2 = []
                data_3 = []
                error_1 = []
                error_2 = []
                error_3 = []
                
                ### Here we get out full epitope name out of the list of antibodies defined at the start of the program
                if display_full_probe_name:
                    name_check = column + ")"
                    for name in antibodies:
                        if name_check in name:
                            trimmed_name = name[:-len(name_check)-2]
                            probe_names.append(trimmed_name)

                
                
                
                #Group1
                group1 = dataframe_subplots[dataframe][column].iloc[group1_range[0]: group1_range[1]]
                group1_list = list(group1)
                
                # Here we account for outliers that are marked as "-1"
                if -1 in group1_list:
                    #input("STOP STOP STOP STOP")
                    while -1 in group1_list:
                        group1_list.remove(-1)
                    print("group 1:", group1_list)
                    #input("STOP STOP STOP STOP")
                
                standard_error_group1 = stats.sem(group1_list)
                group1_average = np.mean(group1_list)
                data_1.append(group1_average)
                error_1.append(standard_error_group1)
                #print(column, group1_name+" average:", group1_average)
                
                
                #Group2
                group2 = dataframe_subplots[dataframe][column].iloc[group2_range[0]: group2_range[1]]
                group2_list = list(group2)
                
                # Here we account for outliers that are marked as "-1"
                if -1 in group2_list:
                    #input("STOP STOP STOP STOP")
                    while -1 in group2_list:
                        group2_list.remove(-1)
                    print("group 2:", group2_list)
                    #input("STOP STOP STOP STOP")
                    
                group2_average = np.mean(group2_list)
                standard_error_group2 = stats.sem(group2_list)
                data_2.append(group2_average)
                error_2.append(standard_error_group2)
                #print(column, group2_name+" average:", group2_average)
                
                
                
                #Group3
                group3 = dataframe_subplots[dataframe][column].iloc[group3_range[0]: group3_range[1]]
                group3_list = list(group3)
                
                # Here we account for outliers that are marked as "-1"
                if -1 in group3_list:
                    #input("STOP STOP STOP STOP")
                    while -1 in group3_list:
                        group3_list.remove(-1)
                    print("group 3:", group3_list)
                    #input("STOP STOP STOP STOP")
                
                standard_error_group3 = stats.sem(group3_list)
                group3_average = np.mean(group3_list)
                data_3.append(group3_average)
                error_3.append(standard_error_group3)
                #print(column, group1_name+" average:", group3_average)
                

                
                # This removes really small values from significance indication on bar plots
                for i in range(len(group1_list)):
                    if group1_list[i] < 5:
                        group1_list[i] = 0
                for i in range(len(group2_list)):
                    if group2_list[i] < 5:
                        group2_list[i] = 0
                for i in range(len(group3_list)):
                    if group3_list[i] < 5:
                        group3_list[i] = 0
                
                
                if df_no == 0:
                    group1_data += data_1
                    group1_error += error_1
                    group2_data += data_2
                    group2_error += error_2
                    group3_data += data_3
                    group3_error += error_3
                    
                    tukey_codes = tukey_HSD([group1_list, group2_list, group3_list])
                    CDTA_tukey_groups.append(tukey_codes)
                    
                    
                elif df_no == 1:
                    group1_negative_data += data_1
                    group1_negative_error += error_1
                    group2_negative_data += data_2
                    group2_negative_error += error_2
                    group3_negative_data += data_3
                    group3_negative_error += error_3
                    
                    tukey_codes = tukey_HSD([group1_list, group2_list, group3_list])
                    NaOH_tukey_groups.append(tukey_codes)
            else:
                print(column, " ERROR (*o*)")
           
    

    # The following will take each group (the CDTA and NaOH data for root, stem, and leaf for a SINGLE CULTIVAR and normalise it to 100)
    if normalise:
        for i in range(len(group1_data)):
            multiplier = 100 / (max([group1_data[i], group2_data[i], group3_data[i], group1_negative_data[i], group2_negative_data[i], group3_negative_data[i]]))
            group1_data[i] *= multiplier
            group2_data[i] *= multiplier
            group3_data[i] *= multiplier
            group1_negative_data[i] *= multiplier 
            group2_negative_data[i] *= multiplier
            group3_negative_data[i] *= multiplier
    
            group1_error[i] *= multiplier
            group2_error[i] *= multiplier
            group3_error[i] *= multiplier
            group1_negative_error[i] *= multiplier 
            group2_negative_error[i] *= multiplier
            group3_negative_error[i] *= multiplier
    
    
    
    # Here all NaOH values are multiplied by (-1) so that we can plot it as negative values on our bar plots
    for data_list in [group1_negative_data, group2_negative_data, group3_negative_data]:
        for i in range(len(data_list)):
            data_list[i] *= -1
    
    # Here we find the upper and lower limits of the bar plot figure by determining the max height of the data and error bars combines (allowing for +10 for significance codes to be added)
    y_upper = 120
    y_lower = -120
    for i in range(len(group1_data)):
        if group1_data[i] + group1_error[i] + 10 > y_upper:
            y_upper = group1_data[i] + group1_error[i] + 10
        if group2_data[i] + group2_error[i] + 10 > y_upper:
            y_upper = group2_data[i] + group2_error[i] + 10
        if group3_data[i] + group3_error[i] + 10 > y_upper:
            y_upper = group3_data[i] + group3_error[i] + 10

    for i in range(len(group1_data)):
        if group1_negative_data[i] - group1_negative_error[i] - 10 < y_lower:
            y_upper = group1_negative_data[i] - group1_negative_error[i] - 10
        if group2_negative_data[i] - group2_negative_error[i] - 10 < y_lower:
            y_upper = group2_negative_data[i] - group2_negative_error[i] - 10
        if group3_negative_data[i] - group3_negative_error[i] - 10 < y_lower:
            y_upper = group3_negative_data[i] - group3_negative_error[i] - 10


        
            
    ##########################################################################
    #               F I G U R E     S E T T I N G S                          #
    ##########################################################################
            
    plt.figure()  
    fig, ax = plt.subplots(figsize = (len(polymers)*1.5, 5)) # This will auto-adjust the figure to the right size. Play around with this if you want other figure dimentions

    #plt.style.use("default")
    #fig.patch.set_facecolor("white")
    
    #print("group 1 error", group1_error)
    #print("group 1 data", group1_data)


    name1 = "Leaf" 
    name2 = "Stem"
    name3 = "Root"
    
    bar_colours = ["darkseagreen",
                   "orange",
                   "cornflowerblue"]
    
    # cute pastel purple and blue
    bar_colours = ["#534B62",
                   "#A499B3",
                   "#D0BCD5"]
    
    bar_colours = ["#06A77D",
                   "#005377",
                   "#052F5F"]
    
       
    bar_colours = ["#5C80BC",
                   "#30323D",
                   "#4D5061"]
    
    # works well with #DBA159 (earth yellow - BEST SO FAR) as NaOH background
    # and UCLA blue #2274A5
    # and sepia #6D3D14
    # SECOND BEST - #8B9EB7 cadet grey
    bar_colours = ["#81C14B", 
                   "#2E933C",
                   "#297045"]


    bar_width = 0.22            
    ax.bar(index - bar_width,  group1_data,  bar_width, yerr = group1_error, color = bar_colours[0], label = name1)
    ax.bar(index,  group2_data,  bar_width, yerr = group2_error, color = bar_colours[1], label = name2)          
    ax.bar(index + bar_width,  group3_data,  bar_width, yerr = group3_error, color = bar_colours[2], label = name3)
    ax.axhline(y=0, color="black", linestyle="-", linewidth = 0.4)
    
    
    # EXPERIMENTAL - adding negative equivalents
    if combine_extractions_barplot:
        bar_width = 0.22
        ax.bar(index - bar_width,  group1_negative_data,  bar_width, yerr = group1_negative_error, color = bar_colours[0], label = name1 + " negative")
        ax.bar(index,  group2_negative_data,  bar_width, yerr = group2_negative_error, color = bar_colours[1], label = name2 + " negative")
        ax.bar(index + bar_width,  group3_negative_data,  bar_width, yerr = group3_negative_error, color = bar_colours[2], label = name3 + " negative")
    
    
    # Background colour
    y_lower, y_upper = ax.get_ylim()
    y_lower = -200
    y_upper = 200
    ax.axhspan(y_lower, 0, facecolor="#5EB1BF", alpha=0.5, zorder=0)
    ax.axhspan(0, y_upper, facecolor="white", alpha=0.5, zorder=0)

    
    title = polymer_class + " in " + sample_class + " organs"
    ax.set_title(title, fontdict={'fontsize': 16, 'fontweight': 'bold', 'fontname': 'Cambria'})
    ax.set_ylabel("Absorbance at 450 nm", fontdict={'fontsize': 12, 'fontweight': 'bold', 'fontname': 'Cambria'})
    
    
    ### Here we add the significance codes on the bar plots
    for i in range(len(polymers)):
        text_CDTA = str(CDTA_tukey_groups[i][0]) + " " + str(CDTA_tukey_groups[i][1]) + " " + str(CDTA_tukey_groups[i][2])
        text_NaOH = str(NaOH_tukey_groups[i][0]) + " " + str(NaOH_tukey_groups[i][1]) + " " + str(NaOH_tukey_groups[i][2])
        #plt.text(index[i], 20, text_CDTA, ha = "center", va = "bottom", fontdict = {'fontsize': 20, 'fontweight': 'bold', 'fontname': 'Cambria'})
        #plt.text(index[i], -20, text_NaOH, ha = "center", va = "bottom", fontdict = {'fontsize': 20, 'fontweight': 'bold', 'fontname': 'Cambria'})
        
        

        ### Here we plot the relevant letter codes on the bar plots to indicate significance
        
        #group1
        code = CDTA_tukey_groups[i][0]
        plt.text(index[i] - bar_width, group1_data[i]+group1_error[i], code, ha = "center", va = "bottom", fontdict = {'fontsize': 10, 'fontweight': 'bold', 'fontname': 'Cambria'})
        #group1_negative
        code = NaOH_tukey_groups[i][0]
        plt.text(index[i] - bar_width, group1_negative_data[i]-10 - group1_negative_error[i]*1.5, code, ha = "center", va = "bottom", fontdict = {'fontsize': 10, 'fontweight': 'bold', 'fontname': 'Cambria'})
        #group2
        code = CDTA_tukey_groups[i][1]
        plt.text(index[i], group2_data[i]+group2_error[i], code, ha = "center", va = "bottom", fontdict = {'fontsize': 10, 'fontweight': 'bold', 'fontname': 'Cambria'})
        #group2_negative
        code = NaOH_tukey_groups[i][1]
        plt.text(index[i], group2_negative_data[i]-10 - group2_negative_error[i]*1.5, code, ha = "center", va = "bottom", fontdict = {'fontsize': 10, 'fontweight': 'bold', 'fontname': 'Cambria'})
        #group3
        code = CDTA_tukey_groups[i][2]
        plt.text(index[i] + bar_width, group3_data[i]+group3_error[i], code, ha = "center", va = "bottom", fontdict = {'fontsize': 10, 'fontweight': 'bold', 'fontname': 'Cambria'})
        #group3_negative
        code = NaOH_tukey_groups[i][2]
        plt.text(index[i] + bar_width, group3_negative_data[i]-10 - group3_negative_error[i]*1.5, code, ha = "center", va = "bottom", fontdict = {'fontsize': 10, 'fontweight': 'bold', 'fontname': 'Cambria'})
        
    
    # The code below adds the antibody abbreviation alongside the target
    if display_additional_probe_code:
        for i in range(len(probe_names)):
            spaces = " " * round(len(probe_names[i])/2)
            #print(spaces, probe_names[i])
            temp_name = probe_names[i] + "\n" + spaces + "(" + polymers[i] + ")"
            probe_names[i] = temp_name
    
    if not len(probe_names) < len(index):
        plt.xticks(index, probe_names)
        if display_full_probe_name:
            ax.set_xticklabels(probe_names, rotation = 50, ha = "right")
        else:
            ax.set_xticklabels(probe_names)
    
    
    """
    ADD functionality of limiting y-axis to height of graph + error
    
    ADD Tukey significance codes here
    
    Old format:
    plt.ylim(0, max(max(group1_data), max(group2_data), max(group3_data), max(group4_data))*1.4)
    """
    
    
        
    ax.legend(loc = "upper left", bbox_to_anchor = (1,1)) #USE bbox_to_anchor = (1,1) for legend outside of plot frame; USE loc = "upper/lower left/right" to change location
    #ax.set_ylim(y_lower, y_upper)
    
    if save_barplots:
        settings = "__normalised__"

        fig_name = sample_class + " " + polymer_class + " - " + settings + ".png"    
        
        create_folder_if_not_exists(sample_class)
        folder = os.getcwd() + "\\" + sample_class + "\\"
        name = folder + fig_name

        plt.savefig(name, dpi = 300, bbox_inches = "tight")    

def bar_plot_3_groups(sample_class, polymer_class, dataframe, extract_used, group_ranges:list):
    save_barplots = False
    
    combine_extractions_barplot = True
    
    display_full_probe_name = True
    display_additional_probe_code = True # This will display the antibody code underneath the target name on the x-axis

    if True: # Indicates which polymers to plot on bar plots based on polymer class supplied to this function
        if polymer_class.lower() == "pectin":
            polymers = [
                "JIM5",
                "JIM7",
                "LM18",
                "LM19",
                "LM20",
                "LM7",
                "LM8",
                "RU2",
                "RU1",
                "LM5",
                "LM6",
                "LM13",
                "LM16"]
        elif polymer_class.lower() == "pectic polymers":
            polymers = [
                ]
        elif polymer_class.lower() == "agp":
            polymers = [
                "JIM8",
                "JIM13",
                "JIM15",
                "JIM16",
                "LM2"]
        elif polymer_class.lower() == "hemicellulose":
            polymers = [
                "LM21",
                "LM15",
                "LM24",
                "LM25",
                "LM11",
                "LM23"]
        elif polymer_class.lower() == "cellulose and ß-glucan":
            polymers = [
                "BS-400-2",
                "BS-400-3",
                "CBM3a"]
        elif polymer_class.lower() == "cellulose and hemicellulose":
            polymers = [
                "CBM3a",
                "BS-400-2",
                "BS-400-3",           
                "LM21",
                "LM15",
                "LM24",
                "LM25",
                "LM11",
                "LM23"]
        elif polymer_class.lower() == "hg & rg-i":
            polymers = ["JIM5", "LM16"]
        elif polymer_class.lower() == "lm25+lm15":
            polymers = ["LM15", "LM25"]
        else:
            print("INVALID POLYMER CLASS!")
            print(polymer_class, "not found")
            sys.exit()


    significance = [] #################### USE BRUTE FORCE CODE TO GET CODES FROM TUKEY TEST
    
    # This will be used to plot a varable number of probes on a single plot without needing manual adjustment
    index = np.arange(len(polymers))
    
    
    group1_range = group_ranges[0]
    group2_range = group_ranges[1]
    group3_range = group_ranges[2]

    
    
    
    group1_data = []
    group1_error = []
    group2_data = []
    group2_error = []
    group3_data = []
    group3_error = []

    
    if display_full_probe_name:
        probe_names = []
    else:
        probe_names = polymers[:]
        
        
    for column in polymers:  # Use the alternative version of this loop-conditional pair if anything goes wrong with sequence and labels on bar plots
        if column in dataframe:
            
    # for column in dataframe: ###### REUSE THIS if anything is broken with bar plots and order
    #     if column in polymers:
            
            ### Here we get out full epitope name out of the list of antibodies defined at the start of the program
            if display_full_probe_name:
                name_check = column + ")"
                for name in antibodies:
                    if name_check in name:
                        #print(name)
                        trimmed_name = name[:-len(name_check)-2]
                        probe_names.append(trimmed_name)
            
            
            
            #Group1
            group1 = dataframe[column].iloc[group1_range[0]: group1_range[1]]
            group1_list = list(group1)
            
            # Here we account for outliers that are marked as "-1"
            if -1 in group1_list:
                #input("STOP STOP STOP STOP")
                while -1 in group1_list:
                    group1_list.remove(-1)
                print("group 1:", group1_list)
                #input("STOP STOP STOP STOP")
            
            standard_error_group1 = stats.sem(group1_list)
            group1_average = np.mean(group1_list)
            group1_data.append(group1_average)
            group1_error.append(standard_error_group1)
            #print(column, group1_name+" average:", group1_average)
            
            
            #Group2
            group2 = dataframe[column].iloc[group2_range[0]: group2_range[1]]
            group2_list = list(group2)
            
            # Here we account for outliers that are marked as "-1"
            if -1 in group2_list:
                #input("STOP STOP STOP STOP")
                while -1 in group2_list:
                    group2_list.remove(-1)
                print("group 2:", group2_list)
                #input("STOP STOP STOP STOP")
                
            group2_average = np.mean(group2_list)
            standard_error_group2 = stats.sem(group2_list)
            group2_data.append(group2_average)
            group2_error.append(standard_error_group2)
            #print(column, group2_name+" average:", group2_average)
            
            
            
            #Group3
            group3 = dataframe[column].iloc[group3_range[0]: group3_range[1]]
            group3_list = list(group3)
            
            # Here we account for outliers that are marked as "-1"
            if -1 in group3_list:
                #input("STOP STOP STOP STOP")
                while -1 in group3_list:
                    group3_list.remove(-1)
                print("group 3:", group3_list)
                #input("STOP STOP STOP STOP")
            
            standard_error_group3 = stats.sem(group3_list)
            group3_average = np.mean(group3_list)
            group3_data.append(group3_average)
            group3_error.append(standard_error_group3)
            #print(column, group1_name+" average:", group3_average)
            
            
            
            
            # This removes really small values from significance indication on bar plots
            for i in range(len(group1_list)):
                if group1_list[i] < 5:
                    group1_list[i] = 0
            for i in range(len(group2_list)):
                if group2_list[i] < 5:
                    group2_list[i] = 0
            for i in range(len(group3_list)):
                if group3_list[i] < 5:
                    group3_list[i] = 0

    
    """
    ***************
    ***************
    Tukey's HSD test needs to happen here.
    Required:
        list of bool values indicating significant differences, e.g., [True, False, False, False, True]
    ***************
    ***************
    ***************
    """
            

            
    plt.figure()    
    fig, ax = plt.subplots(figsize = (len(polymers)*1.5, 5)) # This will auto-adjust the figure to the right size. Play around with this if you want other figure dimentions
    #plt.style.use("default")
    #fig.patch.set_facecolor("white")
    
    #print("group 1 error", group1_error)
    #print("group 1 data", group1_data)
    
    name1 = "Leaf" 
    name2 = "Stem"
    name3 = "Root"

    

    group1_negative = group1_data[:]
    group2_negative = group2_data[:]
    group3_negative = group3_data[:]
    
    for i in range(len(group1_negative)):
        group1_negative[i]*=-1
        group2_negative[i]*=-1
        group3_negative[i]*=-1
    
    bar_width = 0.2            
    ax.bar(index - bar_width,  group1_data,  bar_width, yerr = group1_error, color = "darkseagreen", label = name1)
    ax.bar(index,  group2_data,  bar_width, yerr = group2_error, color = "orange", label = name2)          
    ax.bar(index + bar_width,  group3_data,  bar_width, yerr = group3_error, color = "cornflowerblue", label = name3)
    
    # EXPERIMENTAL - adding negative equivalents
    if combine_extractions_barplot:
        ax.bar(index - bar_width,  group1_negative,  bar_width, yerr = group1_error, color = "darkslategrey", label = name1 + " negative")
        ax.bar(index,  group2_negative,  bar_width, yerr = group2_error, color = "darkgoldenrod", label = name2 + " negative")
        ax.bar(index + bar_width,  group3_negative,  bar_width, yerr = group3_error, color = "blue", label = name3 + " negative")
    
    
    title = sample_class + " " + polymer_class + " (" + extract_used + ")"
    ax.set_title(title, fontdict={'fontsize': 16, 'fontweight': 'bold', 'fontname': 'Cambria'})
    ax.set_ylabel("Absorbance at 450 nm", fontdict={'fontsize': 12, 'fontweight': 'bold', 'fontname': 'Cambria'})
    
    
    
    # The code below adds the antibody abbreviation alongside the target
    if display_additional_probe_code:
        for i in range(len(probe_names)):
            spaces = " " * round(len(probe_names[i])/2)
            #print(spaces, probe_names[i])
            temp_name = probe_names[i] + "\n" + spaces + "(" + polymers[i] + ")"
            probe_names[i] = temp_name
    
    if not len(probe_names) < len(index):
        plt.xticks(index, probe_names)
        if display_full_probe_name:
            ax.set_xticklabels(probe_names, rotation = 50, ha = "right")
        else:
            ax.set_xticklabels(probe_names)
    
    
    """
    ADD functionality of limiting y-axis to height of graph + error
    
    ADD Tukey significance codes here
    
    Old format:
    plt.ylim(0, max(max(group1_data), max(group2_data), max(group3_data), max(group4_data))*1.4)
    """
    
    
        
    ax.legend(loc = "upper left", bbox_to_anchor = (1,1)) #USE bbox_to_anchor = (1,1) for legend outside of plot frame; USE loc = "upper/lower left/right" to change location
    
    if save_barplots:
        settings = "__normalised__"

        fig_name = sample_class + " " + polymer_class + " - " + extract_used + " " + settings + ".png"    
        
        create_folder_if_not_exists(sample_class)
        folder = os.getcwd() + "\\" + sample_class + "\\"
        name = folder + fig_name

        plt.savefig(name, dpi = 300, bbox_inches = "tight")
        
def get_coordinates(template_name: str, sheet_number: int):
    
    # This line gets an absolute path to your template file
    template_path = os.path.abspath(template_name)
    
    template_df = pd.read_excel(template_path, sheet_name = sheet_number, index_col = 0, engine = "openpyxl")
    
    
    # Here we read the user's preferred order of sample names before reducing the size of the dataframe to the standard 96-well layout
    # NOTE - order for all samples is read from the FIRST TEMPLATE FILE ONLY
    if sheet_number == 0:
        global order
        order = list(template_df.iloc[9:, 13].dropna().reset_index(drop = True))
    
    
    # the following trims the template to only contain the cells of interest [row_start:row_end, col_start:col_end] in the standard 96-well format
    template_df = template_df.iloc[9:17, 0:12].reset_index(drop = True)
    # The line below resets the column indexes to start at 0
    template_df.columns = range(len(template_df.columns))
    
    
    # The dictionary below will store the coordinates of the sample replicates for each unique sample name (including blanks)
    sample_coordinates = {}
    
    
    # Here we iterate through each row and column. 
    # If the cell is not empty, the coordinates and sample name is stored
    # This accounts for some samples not having the same number of replicates 
    
    for row_index, row in template_df.iterrows():
        for column_index, cell_value in row.items():
            not_empty = pd.notnull(cell_value)
            if not_empty:
                sample_name = cell_value
                coordinates = [row_index, column_index]
                if sample_name in sample_coordinates:
                    sample_coordinates[sample_name].append(coordinates)
                else:
                    sample_coordinates[sample_name] = [coordinates]
        

    return sample_coordinates

def extract_data(sample_coordinates, template_number):
    
    # Here we get the relative file path to the data folder
    data_folder = os.path.abspath("ELISA_template_1.xlsx")[:-21] + "data"
    
    temp_df = pd.DataFrame()
    
    with open("_outliers.txt", "a") as f: #TODO eliminate
        f.write("\n\n\n")
        time = datetime.now()
        time = time.strftime("OUTPUT: %Y-%m-%d | %H:%M:%S ===================================================================================\n")
        message = "Group " + str(template_number) + " - " + time
        f.write(message)
    
    for file in os.listdir(data_folder):
        
        if file[-5:] == ".xlsx" and file[:17] != "_mapping_template" and int(file[-6]) == template_number:
            
            print("Opening:", file)
            probe_name = file
            
            # This line gets an absolute path to your data file
            probe_path = data_folder + "\\" + probe_name
            
            probe_temp_df = pd.read_excel(probe_path, sheet_name = 0, index_col = 0, engine = "openpyxl")
            
            # the following trims the template to only contain the cells of interest [row_start:row_end, col_start:col_end]
            probe_temp_df = probe_temp_df.iloc[9:17, 0:12].reset_index(drop = True)
            # The line below resets the column indexes to start at 0
            probe_temp_df.columns = range(len(probe_temp_df.columns))
            
            
            # Now we cross-reference the template coordinates with the actual data file for the current probe
            probe_values = {}
            for sample in sample_coordinates:
                probe_values[sample] = []
                for coord in sample_coordinates[sample]:
                    probe_values[sample] += [probe_temp_df.iloc[coord[0], coord[1]]]
                    
            stdev_temp = {}
            
            for sample in probe_values:
                print("#", sample)
                
                if not outlier_label in probe_values[sample]:
                    stdev_temp[sample] = np.std(probe_values[sample])
                    if use_median:
                        probe_values[sample] = sorted(list(probe_values[sample]))[len(probe_values[sample])//2]
                    else:
                        probe_values[sample] = sum (probe_values[sample])/len(probe_values[sample]) #TODO Add ability to indicate median instead of mean if MORE THAN 3 REPS
                
                elif outlier_label in probe_values[sample]:
                    with open("_outliers.txt", "a") as f:
                        f.write("\n")
                        f.write(file)
                        f.write(": ")
                        f.write(sample)
                        f.write(" - ")
                        count_outliers = probe_values[sample].count(outlier_label)
                        if count_outliers == 3:
                            f.write("Biological replicate removed...")
                        elif count_outliers < 3:
                            f.write("Number of technical replicates removed: ")
                            f.write(str(count_outliers))
                        f.write("\n")
                    print(sample, probe_values[sample])
                    
                    for i in range(probe_values[sample].count("outlier")):
                        probe_values[sample].remove("outlier")
                    
                    if probe_values[sample] == []:
                        probe_values[sample] = [-1, -1, -1]
                    print(sample, probe_values[sample])
                    #input("ENTER")
                    stdev_temp[sample] = np.std(probe_values[sample])
                    if use_median:
                        probe_values[sample] = sorted(list(probe_values[sample]))[len(probe_values[sample])//2]
                    else:
                        probe_values[sample] = sum(probe_values[sample])/len(probe_values[sample])
                    print(probe_values[sample])
            stdev[file] = stdev_temp
            # Here the blank value is subtracted from each sample average
            for sample in probe_values:
                if sample != "blank" and probe_values[sample] != -1:
                    probe_values[sample] -= probe_values["blank"]
                    
                    # Here set values smaller than an indicated % of the blank to 0
                    if probe_values[sample] <= probe_values["blank"]*(percent_blank_cut_off_threshold/100):
                        if probe_values[sample] == -1:
                            pass
                        else:
                            probe_values[sample] = 0
            # Here the blank value is removed from the dictionary
            probe_values.pop("blank")
            # Here we append the final heatmap dataframe with the current probe
            sample_names = list(probe_values.keys())
            sample_values = list(probe_values.values())
            
            if multiple_templates:
                probe_name = probe_name[:-7] 
            
            temp_df[probe_name] = pd.DataFrame(sample_values, index = sample_names, columns = [probe_name])
            
            
            """
                stdev_temp[sample] = np.std(probe_values[sample])
                probe_values[sample] = sum(probe_values[sample])/len(probe_values[sample]) #TODO median function
            stdev[file] = stdev_temp
                
            # Here the blank value is subtracted from each sample average
            for sample in probe_values:
                if sample != "blank":
                    probe_values[sample] -= probe_values["blank"]
                    
                    # Here set values smaller than an indicated % of the blank to 0
                    if probe_values[sample] <= probe_values["blank"]*(percent_blank_cut_off_threshold/100):
                        probe_values[sample] = 0
            # Here the blank value is removed from the dictionary
            probe_values.pop("blank")
                
                
                
            # Here we append the final heatmap dataframe with the current probe
            sample_names = list(probe_values.keys())
            sample_values = list(probe_values.values())
            
            if multiple_templates:
                probe_name = probe_name[:-7] 
            
            temp_df[probe_name] = pd.DataFrame(sample_values, index = sample_names, columns = [probe_name])
            """
            
            print("done \n") # This is included for troubleshooting purposes
            
    global heatmap
    if template_number == 1:
        heatmap = temp_df
    elif template_number > 1:
        heatmap = pd.concat([heatmap, temp_df], axis=0)
        
        
    # Now we use the indicated preferred order of sample names
    global heatmap_ordered, heatmap_ordered_unnormalised
    heatmap_ordered = heatmap.reindex(order)
    heatmap_ordered_unnormalised = heatmap_ordered[:]
    
    
    # Here we reorder the probes on the x-axis of the heatmap
    ##### The user either has to manually indicate the order in a python list OR BAKE INTO CODE FOR AUTOMATIC DETECTION >>>>>>>
    #custom_order = ["LM5", "LM6", "LM7", "LM16", "LM18", "LM19", "LM20", "JIM5", "JIM7", "RU_1", "RU2", "LM21", "LM11", "LM15", "LM24", "LM25", "JIM15", "JIM16", "JIM8", "JIM13", "LM2", "LM23", "BS-400-2", "BS-400-3", "CBM3a"]
    custom_order = ["JIM5", "JIM7", "LM18", "LM19", "LM20", "LM7", "RU2", "RU1", "LM5", "LM6", "LM16", "LM21", "LM15", "LM24", "LM25", "LM11", "LM23", "JIM8", "JIM13", "JIM15", "JIM16", "LM2", "BS-400-2", "BS-400-3", "CBM3a"]
    print("\n\n\n\nNumber of probes detected:", len(custom_order))
    try:
        heatmap_ordered = heatmap_ordered[custom_order]
    except KeyError:
        print("ERROR: Please check custom probe order for omissions or spelling errors")
        
    
    
    # Now the data is normalised to 100
    # The max value per column is set to 100 and the rest relative to that
    for column in heatmap_ordered:
        #print(column)
        col_max = -1
        for row_name, cell_value in heatmap_ordered[column].items():
            debugging = False
            if debugging:
                print("cell value: ", cell_value)
                print("row name: ", row_name)
            if cell_value > col_max:
                col_max = cell_value
        if not col_max == 0:
            multiplier = 100 / col_max
        
        if not normalise:
            multiplier = 1
        
        for row_name, cell_value in heatmap_ordered[column].items():
            heatmap_ordered[column][row_name] *= multiplier
        
        
        for row_name, cell_value in heatmap_ordered[column].items():
            if cell_value < min_cutoff:
                heatmap_ordered[column][row_name] *= 0
 
def create_heatmap(heatmap_dataframe):
    
    """
    This function will either take a dataframe from which to make a heatmap
    OR a dictionary containing multiple heatmaps, depending on how many
    extractions the user indicated. The subplots can either be below or beside
    each other.
    
    Make sure custom_order is set on heatmap and that no error occurs
    Otherwise the heatmap will not correspond to the heatmap labels.
    """
    
    # Change save_image to True if you want to save the plot
    save_image = True
    # Increase/decrease this value for higher/lower resolution image output
    image_dpi = 50
    
    coloured_columns = use_antibody_colours
    try_cmap(colour_scheme)
    
    # Change the values below to fine-tune the final heatmap
    decimal_places = 0
    x_axis_font_scale = fontsize_axis
    y_axis_font_scale = fontsize_axis
    title_font_scale = fontsize_title
    figure_dimensions = heatmap_dimensions
    x_label_rotation = x_rotation
    y_label_rotation = y_rotation
    numbers_visible = True # change to False for heatmap without values
    
    column_colours = ['blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'red', 'red', 'red', 'red', 'red', 'red', 'green', 'green', 'green', 'green', 'green', 'purple', 'purple', 'orange']
    column_colours = ['navy', 'navy', 'navy', 'navy', 'navy', 'navy', 'navy', 'navy', 'navy', 'navy', 'navy', 'maroon', 'maroon', 'maroon', 'maroon', 'maroon', 'maroon', 'black', 'black', 'black', 'black', 'black', 'chocolate', 'chocolate', 'purple']
    column_colours = []
    for i in range(11):
        column_colours.append("navy")
    for i in range(6):
        column_colours.append("maroon")
    for i in range(5):
        column_colours.append("black")
    for i in range(2):
        column_colours.append("chocolate")
    for i in range(1):
        column_colours.append("purple")


    if DEBUGGING:
        plt.figure()
        hmap = sns.heatmap(heatmap_dataframe["C"], cmap = "Greens")
        # ["JIM7", "LM18", "LM20", "RU1", "LM5", "LM6"]
        hmap.set_title("CDTA")
        plt.show()
        plt.figure()
        hmap = sns.heatmap(heatmap_dataframe["N"], cmap = "Greens")
        hmap.set_title("NaOH")
        plt.show()




    # Check whether a single heatmap or a dictionary with subplots was passed to create_heatmap()
    if type(heatmap_dataframe) == dict:
        fig, axes = plt.subplots(extraction_count, 1, figsize = figure_dimensions) # N rows (extraction plots); 1 column (underneath each other)
        # This configuration will place subplots underneath each other
        

        
        # heatmap_dataframe = OrderedDict([(), ()])
        
        for i, plot in enumerate(heatmap_dataframe):

            subplot = heatmap_dataframe[plot].round(decimals = decimal_places)
            sns.set(font_scale=0.8)      
            
            # Here a custom colour palette is indicated using hex codes
            # Built-in colour schemes can also be used (see bottom of script)
            custom_colours = ["#FFFFFF", "#EADAF0", "#E7C8F2", "#E49EDF", "#D77AD1", "#BF5DB8" ,"#A4459D" ,"#812A7C", "#690F63"]
            custom_palette = sns.color_palette(custom_colours)


            if bold_axis_titles:
                font_weight = "bold"
            else:
                font_weight = "light"            

            
            if extraction_count > 1:
                axx = axes[i]
            elif extraction_count == 1:
                axx = axes
            
            heatmap_plot = sns.heatmap(subplot, ax = axx, linewidths = 0, cmap = colour_scheme, annot=numbers_visible, fmt = "1g", cbar_kws = {"orientation": "vertical"}, vmin=0, vmax=100) # Use argument annot=True for numbers on heatmap
            heatmap_plot.set_yticklabels(heatmap_plot.get_yticklabels(), rotation = y_label_rotation, fontsize = y_axis_font_scale, fontdict = {"family": "serif", "size": y_axis_font_scale, "weight": font_weight}) ### Here we get the y labels from the existing dataframe and set the rotation of the labels to 0
            #heatmap_plot.set_xticklabels(heatmap_plot.get_xticklabels(), rotation = x_label_rotation, fontsize = x_axis_font_scale, fontdict = {"family": "serif", "size": x_axis_font_scale}) ## Here we do the same for x but set the labels at an angle for easier readability
            boi = heatmap_plot.get_xticklabels()
            
            
            

            
            # Here the bottom x-labels are removed
            heatmap_plot.set_xticklabels(antibodies, rotation = x_label_rotation, fontsize = x_axis_font_scale, fontdict = {"family": "serif", "size": x_axis_font_scale, "weight": font_weight}, ha = "left")
            heatmap_plot.xaxis.tick_top()
            
            y_title = plot
            if i == 0:
                y_title = "CDTA"
            elif i == 1:
                y_title = "NaOH"
            
            
            if coloured_columns:
                for j, tick in enumerate(heatmap_plot.get_xticklabels()):
                    tick.set_color(column_colours[j])
            
            
            heatmap_plot.set_ylabel(y_title, fontdict = {"family": "serif", "size": fontsize_axis, "weight": font_weight})
            
            
            #Here we add a title if the user indicates so
            if len(graph_title) > 0:
                heatmap_plot.set_title(graph_title, fontdict = {"family": "serif", "size": fontsize_title, "weight": "bold"})

        plt.tight_layout()


        heatmap_plot.figure.savefig("TEMP"+".png", dpi=image_dpi, bbox_inches = "tight")
        os.system("TEMP.png")
        
        
        global global_heatmap
        global_heatmap = heatmap_plot
            
    else:
        pass
        """
        This area is reserved for the case where no subplots will be made.
        In other words, either a single extraction is visualised, or
        multiple extractions will be visualised on the same heatmap.
        """

def create_grouped_bar_plot(probe, probe_index, to_be_removed):
    
    
    save_barplot = False
    barplot_dpi = 250
    
    
    global barplot_frame
    barplot_frame = pd.DataFrame()
    try:
        barplot_frame["AU"] = heatmap_ordered_unnormalised[probe]
    except:
        print("\nXXXXXX\n\nInvalid probe indicated for grouped bar plot!\n\nXXXXXX")
    
    barplot_frame = barplot_frame.drop(to_be_removed)
    barplot_frame["Extraction"] = ["C", "C", "C", "C", "C", "N", "N", "N", "N", "N" ]
    
    no_extraction_names = []
    for name in barplot_frame.index:
        no_extraction_names.append(name[2:])
    
    barplot_frame["Sample"] = no_extraction_names
    
    probe_CDTA = stdev[probe + "_1.xlsx"]
    probe_NaOH = stdev[probe + "_2.xlsx"]
    global probe_error_combined
    probe_error_combined = {**probe_CDTA, **probe_NaOH}
    
    # Here we remove irrelevant samples from the std dev dictionary. Always remove "blank"
    for sample in to_be_removed:
        probe_error_combined.pop(sample)
    probe_error_combined.pop("blank")
    
    barplot_frame["stdev"] = probe_error_combined.values()
    barplot_frame["check std matches w/ sample"] = probe_error_combined.keys()

    
    global leaf_order, cdta_au, cdta_error, naoh_au, naoh_error

    mid = int(len(barplot_frame)/2)
    
    leaf_order = barplot_frame["Sample"].iloc[:mid]
    leaf_order_no_BB = []
    for name in leaf_order:
        leaf_order_no_BB.append("leaf " + name[-1])
    cdta_au, cdta_error = barplot_frame["AU"].iloc[:mid], barplot_frame["stdev"].iloc[:mid]
    naoh_au, naoh_error = barplot_frame["AU"].iloc[mid:], barplot_frame["stdev"].iloc[mid:]
    
    
    
    # Here we create the physical barplot
    plt.style.use("default")
    plt.figure(figsize=(10,2))
    x = np.arange(len(leaf_order))  # the label locations
    width = 0.35  # the width of the bars
    
    fig, ax = plt.subplots(figsize = (10,5))
    fig.patch.set_facecolor("white")
    
    ax.bar(x - width/2, cdta_au, width, label='CDTA', color='darkseagreen', yerr=cdta_error, capsize=0)
    ax.bar(x + width/2, naoh_au, width, label='NaOH', color='slategray', yerr=naoh_error, capsize=0)
    
    
    # Below we customise the figure axes and titles
    ax.set_xticks(x)
    ax.set_xticklabels(leaf_order_no_BB)
    ax.legend(loc = "upper left", bbox_to_anchor = (1, 1))
    ax.set_xlabel("Leaf order", fontdict={'fontsize': 12, 'fontweight': 'bold', 'fontname': 'Cambria'})
    ax.set_ylabel("Absorbance at 450 nm", fontdict={'fontsize': 12, 'fontweight': 'bold', 'fontname': 'Cambria'})
    ax.set_title(antibodies[probe_index], fontdict={'fontsize': 16, 'fontweight': 'bold', 'fontname': 'Cambria'})
    

    if save_barplot: # I had to change the "/" in LM25's label so avoid it being seen as a directory
        plt.savefig(antibodies[probe_index] + ".png", dpi = barplot_dpi, bbox_inches = "tight")

    plt.show()

def create_correlation_table(probe_dataframe, samples_to_remove, extraction: str):
    """
    This function is used to 1) Create a correlation matrix from a dataframe,
    and 2) to create a heatmap using said matrix.
    
    This function needs to be called for each fraction of your cell wall material.
    The samples_to_remove parameter asks for a set of samples to be excluded from the matrix.
    
    e.g., for extract in ["CDTA_data", "NaOH_data"]:
                create_correlation_table(extract, ["bad_sample_1", "bad_sample_2"])
    
    
    
    
    TAKE SPECIAL NOTE:
        * If there is no variation in a numerical column, the .corr() method
        will replace the column with "NaN" (Not a Number). Indicate the name of 
        these probes by providing a list to the samples_to_remove argument.
    """
    
    full_correlation_table = False
    
    
    x_axis_font_scale = 16
    y_axis_font_scale = 16
    title_font_scale = 25
    fig_size = (25, 25)
    colour_palette = "RdYlGn"
    
    
    global corr_matrix
    
    if samples_to_remove:
        corr_dataframe = probe_dataframe.drop(samples_to_remove)
    else:
        corr_dataframe = probe_dataframe
        
    # Check whether "spearman" or "pearson" correlation is applicable
    corr_matrix = corr_dataframe.corr(method = "spearman").round(2)
    
    plt.figure(figsize = fig_size)
    corr_matrix_rounded = corr_matrix.round(decimals = 1)
    
    # Now we use the seaborn library to create a correlation table from the generated matrix
    sns.set(font_scale=1.0)
    
    
    if not full_correlation_table:
        sns.set_style("white")
        mask = np.triu(np.ones_like(corr_matrix_rounded, dtype = bool), k=1)
        corr_heatmap = sns.heatmap(corr_matrix_rounded, mask = mask, annot = True, linewidths = 0, fmt = "1g", cmap = colour_palette, cbar_kws = {"orientation": "vertical", "shrink": 0.5, "aspect": 35}, square = True)
    else:
        corr_heatmap = sns.heatmap(corr_matrix_rounded, annot = True, linewidths = 0, fmt = "1g", cmap = colour_palette)
    
    corr_heatmap.set_yticklabels(corr_heatmap.get_yticklabels(), rotation = 0, fontsize = y_axis_font_scale, fontdict = {"family": "serif", "size": y_axis_font_scale}) ### Here we get the y labels from the existing dataframe and set the rotation of the labels to 0
    corr_heatmap.set_xticklabels(corr_heatmap.get_xticklabels(), rotation = 50, fontsize = x_axis_font_scale, fontdict = {"family": "serif", "size": x_axis_font_scale}, ha = "right") ## Here we do the same for x but set the labels at an angle for easier readability
    #corr_heatmap.set_title("CORRELATION TABLE: " + extraction, fontsize = title_font_scale, fontdict = {"family": "serif", "size": title_font_scale})
    
    boi = random.randint(1,10000)
    #plt.savefig("Maked corr table - spearman" + str(boi), dpi=100, bbox_inches = "tight")

def create_folder_if_not_exists(folder_name):
    abs_path = os.path.join(os.getcwd(), folder_name)
    
    path_exists = os.path.exists(abs_path)
    
    if not path_exists:
        os.makedirs(abs_path)
        print(f"Folder '{abs_path}' created")
    else:
        print(f"Folder '{abs_path}' already exists")
    
def tukey_HSD(data_list_of_lists) -> list :
    """
    This function takes a 2-dimensional array of datapoints,
    e.g., [group1, group2, group3] where each group refers to a list of data points.
    Any number of points may be given and a True/False tuple will be generated that corresponds
    to significant differences (True if different) between different pairs of groups.
    
    Output is a tuple that corresponds to significant diffs between groups.
    """
    codes = {
        (False, False, False): ("a", "a", "a"),
        (False, False, True): ("ab", "a", "b"),
        (False, True, False): ("a", "ab", "b"),
        (True, False, False): ("a", "b", "ab"),
        (True, True, False): ("a", "bc", "c"), #################################### Shouldn't this be (a,b,b)?????????????????????????????????
        (True, False, True): ("a", "b", "a"),
        (False, True, True): ("a", "a", "b"),
        (True, True, True): ("a", "b", "c")}
    
    
    
    
    list_count = len(data_list_of_lists)
    if list_count == 3:
        f_statistic, p_value = f_oneway(data_list_of_lists[0], data_list_of_lists[1], data_list_of_lists[2])
        
        if p_value < 0.05:
            print("There are significant differences among the groups.")
            
            # All data will be combined and labels will be assigned for post-hoc test
            all_data = np.concatenate((data_list_of_lists[0], data_list_of_lists[1], data_list_of_lists[2]))
            labels = ["Group 1"]*len(data_list_of_lists[0]) + ["Group 2"]*len(data_list_of_lists[1]) + ["Group 3"]*len(data_list_of_lists[2])
            
            # Tukey results
            tukey_results = pairwise_tukeyhsd(all_data, labels, alpha = 0.05)
            print(tukey_results)
        
            reject = tuple(tukey_results.reject)
            output = codes[reject]
            return output
        else:
            return ("a", "a", "a")
        
    else:
        print("FUNCTIONALITY OF MORE OR FEWER THAN 3 GROUPS NOT YET SUPPORTED FOR TUKEY'S HSD")
    
    
def try_cmap(cmap_name):
    try:
        cmap = plt.get_cmap(cmap_name)
    except ValueError:
        print(f"Error: '{cmap_name}' is not a valid colormap name.")
        print("Valid colormap names:")
        print(plt.colormaps())
    else:
        print(f"Successfully used colormap: {cmap_name}")



def error(err_code):
    print("ERROR")
    print(err_code)

def truefalse(state, value_if_true, value_if_false):
    """
    Parameters
    ----------
    state : BOOL
        True or False value corresponding to the state of a checkbox (checked or not).
    value_if_true : str or int
        Value that will be returned if state is True.
    value_if_false : str or int
        Value that will be returned if state is False.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    
    if state == True:
        return value_if_true
    elif state == False:
        return value_if_false

def bar_plot_3_groups_extractions_combined(sample_class, polymer_class, dataframe_subplots, group_ranges:list):
    save_barplots = False
    
    combine_extractions_barplot = True
    
    display_full_probe_name = True
    display_additional_probe_code = True # This will display the antibody code underneath the target name on the x-axis

    if True: # Indicates which polymers to plot on bar plots based on polymer class supplied to this function
        if polymer_class.lower() == "pectin":
            label_x_pos = 0.07
            polymers = [ 
                "JIM5",
                "JIM7",
                "LM18",
                "LM19",
                "LM20",
                "LM7",
                "RU2",
                "RU1",
                "LM5",
                "LM6",
                "LM16"]
            
        elif polymer_class.lower() == "pectin_for_viz_figure":
            polymers = [ 
                
                "JIM7",
                "LM18",
                "LM20",
                
                
                "RU1",
                "LM5",
                "LM6",
                ]
            
        elif polymer_class.lower() == "agp":
            label_x_pos = 0.3
            polymers = [
                "JIM8",
                "JIM13",
                "JIM15",
                "JIM16",
                "LM2"]
        elif polymer_class.lower() == "hemicellulose":
            polymers = [
                "LM21",
                "LM15",
                "LM24",
                "LM25",
                "LM11",
                "LM23"]
        elif polymer_class.lower() == "cellulose and ß-glucan":
            polymers = [
                "BS-400-2",
                "BS-400-3",
                "CBM3a"]
        elif polymer_class.lower() == "cellulose and hemicellulose":
            label_x_pos = 0.11
            polymers = [
                "CBM3a",
                "BS-400-2",
                "BS-400-3",           
                "LM21",
                "LM15",
                "LM24",
                "LM25",
                "LM11",
                "LM23"]
        elif polymer_class.lower() == "hg & rg-i":
            polymers = ["JIM5", "LM16"]
        elif polymer_class.lower() == "lm25+lm15":
            polymers = ["LM15", "LM25"]
        else:
            print("INVALID POLYMER CLASS!")
            print(polymer_class, "not found")
            sys.exit()

    
    ##########################################################################
    #                   D A T A     S E T T I N G S                          #
    ##########################################################################


    # This will be used to plot a varable number of probes on a single plot without needing manual adjustment
    index = np.arange(len(polymers))
    
    group1_range = group_ranges[0]
    group2_range = group_ranges[1]
    group3_range = group_ranges[2]
    
    
    group1_data = []
    group1_error = []
    group1_negative_data = []
    group1_negative_error = []
    
    group2_data = []
    group2_error = []
    group2_negative_data = []
    group2_negative_error = []
    
    
    group3_data = []
    group3_error = []
    group3_negative_data = []
    group3_negative_error = []



    # TUKEY'S HSD TESTING
    CDTA_tukey_groups = []
    NaOH_tukey_groups = []
    
    
    
    for df_no, dataframe in enumerate(dataframe_subplots):
        
        if display_full_probe_name:
            probe_names = []
        else:
            probe_names = polymers[:]
        
        
        
        for column in polymers:  # Use the alternative version of this loop-conditional pair if anything goes wrong with sequence and labels on bar plots
            if column in dataframe_subplots[dataframe]:
                
        # for column in dataframe: ###### REUSE THIS if anything is broken with bar plots and order
        #     if column in polymers:
                
            
                data_1 = []
                data_2 = []
                data_3 = []
                error_1 = []
                error_2 = []
                error_3 = []
                
                ### Here we get out full epitope name out of the list of antibodies defined at the start of the program
                if display_full_probe_name:
                    name_check = column + ")"
                    for name in antibodies:
                        if name_check in name:
                            trimmed_name = name[:-len(name_check)-2]
                            probe_names.append(trimmed_name)

                
                
                
                #Group1
                group1 = dataframe_subplots[dataframe][column].iloc[group1_range[0]: group1_range[1]]
                group1_list = list(group1)
                if column == "CBM3a":
                    group1_list = [0,0,0,0]
                    
                # Here we account for outliers that are marked as "-1"
                if -1 in group1_list:
                    #input("STOP STOP STOP STOP")
                    while -1 in group1_list:
                        group1_list.remove(-1)
                    print("group 1:", group1_list)
                    #input("STOP STOP STOP STOP")
                
                standard_error_group1 = stats.sem(group1_list)
                group1_average = np.mean(group1_list)
                data_1.append(group1_average)
                error_1.append(standard_error_group1)
                #print(column, group1_name+" average:", group1_average)
                
                
                #Group2
                group2 = dataframe_subplots[dataframe][column].iloc[group2_range[0]: group2_range[1]]
                group2_list = list(group2)
                if column == "CBM3a":
                    group2_list = [0,0,0,0]
                
                # Here we account for outliers that are marked as "-1"
                if -1 in group2_list:
                    #input("STOP STOP STOP STOP")
                    while -1 in group2_list:
                        group2_list.remove(-1)
                    print("group 2:", group2_list)
                    #input("STOP STOP STOP STOP")
                    
                group2_average = np.mean(group2_list)
                standard_error_group2 = stats.sem(group2_list)
                data_2.append(group2_average)
                error_2.append(standard_error_group2)
                #print(column, group2_name+" average:", group2_average)
                
                
                
                #Group3
                group3 = dataframe_subplots[dataframe][column].iloc[group3_range[0]: group3_range[1]]
                group3_list = list(group3)
                if column == "CBM3a":
                    group3_list = [0,0,0,0]
                
                # Here we account for outliers that are marked as "-1"
                if -1 in group3_list:
                    #input("STOP STOP STOP STOP")
                    while -1 in group3_list:
                        group3_list.remove(-1)
                    print("group 3:", group3_list)
                    #input("STOP STOP STOP STOP")
                
                standard_error_group3 = stats.sem(group3_list)
                group3_average = np.mean(group3_list)
                data_3.append(group3_average)
                error_3.append(standard_error_group3)
                #print(column, group1_name+" average:", group3_average)
                

                
                # This removes really small values from significance indication on bar plots
                for i in range(len(group1_list)):
                    if group1_list[i] < 5:
                        group1_list[i] = 0
                for i in range(len(group2_list)):
                    if group2_list[i] < 5:
                        group2_list[i] = 0
                for i in range(len(group3_list)):
                    if group3_list[i] < 5:
                        group3_list[i] = 0
                
                
                
                if df_no == 0:
                    group1_data += data_1
                    group1_error += error_1
                    group2_data += data_2
                    group2_error += error_2
                    group3_data += data_3
                    group3_error += error_3
                    
                    tukey_codes = tukey_HSD([group1_list, group2_list, group3_list])
                    CDTA_tukey_groups.append(tukey_codes)
                    
                    
                elif df_no == 1:
                    group1_negative_data += data_1
                    group1_negative_error += error_1
                    group2_negative_data += data_2
                    group2_negative_error += error_2
                    group3_negative_data += data_3
                    group3_negative_error += error_3
                    
                    tukey_codes = tukey_HSD([group1_list, group2_list, group3_list])
                    NaOH_tukey_groups.append(tukey_codes)
            else:
                print(column, " ERROR (*o*)")
    
    
    

    # The following will take each group (the CDTA and NaOH data for root, stem, and leaf for a SINGLE CULTIVAR and normalise it to 100)
    if normalise:
        for i in range(len(group1_data)):
            multiplier = 100 / (max([group1_data[i], group2_data[i], group3_data[i], group1_negative_data[i], group2_negative_data[i], group3_negative_data[i]]))
            group1_data[i] *= multiplier
            group2_data[i] *= multiplier
            group3_data[i] *= multiplier
            group1_negative_data[i] *= multiplier 
            group2_negative_data[i] *= multiplier
            group3_negative_data[i] *= multiplier
    
            group1_error[i] *= multiplier
            group2_error[i] *= multiplier
            group3_error[i] *= multiplier
            group1_negative_error[i] *= multiplier 
            group2_negative_error[i] *= multiplier
            group3_negative_error[i] *= multiplier
    
    
    
    # Here all NaOH values are multiplied by (-1) so that we can plot it as negative values on our bar plots
    for data_list in [group1_negative_data, group2_negative_data, group3_negative_data]:
        for i in range(len(data_list)):
            data_list[i] *= 1
    
    # Here we find the upper and lower limits of the bar plot figure by determining the max height of the data and error bars combines (allowing for +10 for significance codes to be added)
    y_upper = 135
    y_lower = -135
    for i in range(len(group1_data)):
        if group1_data[i] + group1_error[i] + 10 > y_upper:
            y_upper = group1_data[i] + group1_error[i] + 10
        if group2_data[i] + group2_error[i] + 10 > y_upper:
            y_upper = group2_data[i] + group2_error[i] + 10
        if group3_data[i] + group3_error[i] + 10 > y_upper:
            y_upper = group3_data[i] + group3_error[i] + 10

    for i in range(len(group1_data)):
        if group1_negative_data[i] - group1_negative_error[i] - 10 < y_lower:
            y_upper = group1_negative_data[i] - group1_negative_error[i] - 10
        if group2_negative_data[i] - group2_negative_error[i] - 10 < y_lower:
            y_upper = group2_negative_data[i] - group2_negative_error[i] - 10
        if group3_negative_data[i] - group3_negative_error[i] - 10 < y_lower:
            y_upper = group3_negative_data[i] - group3_negative_error[i] - 10

    


            
    ##########################################################################
    #               F I G U R E     S E T T I N G S                          #
    ##########################################################################
            
    plt.figure()  
    fig, (ax1, ax2) = plt.subplots(2,1, sharex = True, figsize = (len(polymers)*1.5, 5)) # This will auto-adjust the figure to the right size. Play around with this if you want other figure dimentions

    #plt.style.use("default")
    #fig.patch.set_facecolor("white")
    
    #print("group 1 error", group1_error)
    #print("group 1 data", group1_data)


    name1 = "Leaf" 
    name2 = "Stem"
    name3 = "Root"
    
    bar_colours = ["darkseagreen",
                   "orange",
                   "cornflowerblue"]
    
    # cute pastel purple and blue
    bar_colours = ["#534B62",
                   "#A499B3",
                   "#D0BCD5"]
    
    bar_colours = ["#06A77D",
                   "#005377",
                   "#052F5F"]
    
       
    bar_colours = ["#5C80BC",
                   "#30323D",
                   "#4D5061"]
    
    # works well with #DBA159 (earth yellow - BEST SO FAR) as NaOH background
    # and UCLA blue #2274A5
    # and sepia #6D3D14
    # SECOND BEST - #8B9EB7 cadet grey
    bar_colours = ["#81C14B", 
                   "#2E933C",
                   "#297045"]

    bar_width = 0.22            
    ax1.bar(index - bar_width,  group1_data,  bar_width, yerr = group1_error, color = bar_colours[0], label = name1)
    ax1.bar(index,  group2_data,  bar_width, yerr = group2_error, color = bar_colours[1], label = name2)          
    ax1.bar(index + bar_width,  group3_data,  bar_width, yerr = group3_error, color = bar_colours[2], label = name3)
    #ax1.axhline(y=0, color="black", linestyle="-", linewidth = 0.4)
    


    
    # EXPERIMENTAL - adding negative equivalents
    if combine_extractions_barplot:
        bar_width = 0.22
        ax2.bar(index - bar_width,  group1_negative_data,  bar_width, yerr = group1_negative_error, color = bar_colours[0], label = name1 + " negative")
        ax2.bar(index,  group2_negative_data,  bar_width, yerr = group2_negative_error, color = bar_colours[1], label = name2 + " negative")
        ax2.bar(index + bar_width,  group3_negative_data,  bar_width, yerr = group3_negative_error, color = bar_colours[2], label = name3 + " negative")
    
    
    # Background colour
    #y_lower, y_upper = ax.get_ylim()
    # y_lower = -200
    # y_upper = 200
    #ax.axhspan(y_lower, 0, facecolor="#5EB1BF", alpha=0.5, zorder=0)
    #ax.axhspan(0, y_upper, facecolor="white", alpha=0.5, zorder=0)

    
    title = "Relative levels of " + polymer_class + " in " + sample_class + " organs"
    ax1.set_title(title, fontdict={'fontsize': 16, 'fontweight': 'bold', 'fontname': 'Cambria'})
    #ax1.set_ylabel("Absorbance at 450 nm", fontdict={'fontsize': 12, 'fontweight': 'bold', 'fontname': 'Cambria'})
    
    
    ax1.set_ylim(0, 136)
    ax2.set_ylim(0, 136)
    
    
    
    # Add text for CDTA/NaOH
    plt.text(len(polymers)-label_x_pos, 190, "CDTA", rotation=90)
    plt.text(len(polymers)-label_x_pos, 60, "NaOH", rotation=90)
    
    ### Here we add the significance codes on the bar plots
    for i in range(len(polymers)):
        text_CDTA = str(CDTA_tukey_groups[i][0]) + " " + str(CDTA_tukey_groups[i][1]) + " " + str(CDTA_tukey_groups[i][2])
        text_NaOH = str(NaOH_tukey_groups[i][0]) + " " + str(NaOH_tukey_groups[i][1]) + " " + str(NaOH_tukey_groups[i][2])
        #plt.text(index[i], 20, text_CDTA, ha = "center", va = "bottom", fontdict = {'fontsize': 20, 'fontweight': 'bold', 'fontname': 'Cambria'})
        #plt.text(index[i], -20, text_NaOH, ha = "center", va = "bottom", fontdict = {'fontsize': 20, 'fontweight': 'bold', 'fontname': 'Cambria'})
        
        

        ### Here we plot the relevant letter codes on the bar plots to indicate significance


        #group1
        code = CDTA_tukey_groups[i][0]
        #plt.text(index[i] - bar_width, group1_data[i]+group1_error[i], code, ha = "center", va = "bottom", fontdict = {'fontsize': 10, 'fontweight': 'bold', 'fontname': 'Cambria'})
        ax1.annotate(code, (index[i] - bar_width, group1_data[i]+5+group1_error[i]*1.5), ha="center",color='black', size=11)
        #group1_negative
        code = NaOH_tukey_groups[i][0]
        #plt.text(index[i] - bar_width, group1_negative_data[i]-10 - group1_negative_error[i]*1.5, code, ha = "center", va = "bottom", fontdict = {'fontsize': 10, 'fontweight': 'bold', 'fontname': 'Cambria'})
        ax2.annotate(code, (index[i] - bar_width, group1_negative_data[i]+5 + group1_negative_error[i]*1.5), ha="center",color='black', size=11)
        
        #group2
        code = CDTA_tukey_groups[i][1]
        #plt.text(index[i], group2_data[i]+group2_error[i], code, ha = "center", va = "bottom", fontdict = {'fontsize': 10, 'fontweight': 'bold', 'fontname': 'Cambria'})
        ax1.annotate(code, (index[i], group2_data[i]+5+group2_error[i]*1.5), ha="center",color='black', size=11)
        #group2_negative
        code = NaOH_tukey_groups[i][1]
        #plt.text(index[i], group2_negative_data[i]-10 - group2_negative_error[i]*1.5, code, ha = "center", va = "bottom", fontdict = {'fontsize': 10, 'fontweight': 'bold', 'fontname': 'Cambria'})
        ax2.annotate(code, (index[i], group2_negative_data[i]+5 + group2_negative_error[i]*1.5), ha="center",color='black', size=11)
        #group3
        code = CDTA_tukey_groups[i][2]
        #plt.text(index[i] + bar_width, group3_data[i]+group3_error[i], code, ha = "center", va = "bottom", fontdict = {'fontsize': 10, 'fontweight': 'bold', 'fontname': 'Cambria'})
        ax1.annotate(code, (index[i] + bar_width, group3_data[i]+5+group3_error[i]*1.5), ha="center",color='black', size=11)
        #group3_negative
        code = NaOH_tukey_groups[i][2]
        #plt.text(index[i] + bar_width, group3_negative_data[i]-10 - group3_negative_error[i]*1.5, code, ha = "center", va = "bottom", fontdict = {'fontsize': 10, 'fontweight': 'bold', 'fontname': 'Cambria'})
        ax2.annotate(code, (index[i] + bar_width, group3_negative_data[i]+5 + group3_negative_error[i]*1.5), ha="center",color='black', size=11)
    
    
    x = 0.025
    if polymer_class.lower() == "cellulose and hemicellulose":
        x = 0.07
    elif polymer_class.lower() == "pectin":
        x = 0.07
    fig.text(x, 0.5, 'Relative absorbance', va='center', rotation='vertical', fontdict={'fontsize': 15, 'fontweight': 'bold', 'fontname': 'Cambria'})
    
    
    # The code below adds the antibody abbreviation alongside the target
    if display_additional_probe_code:
        for i in range(len(probe_names)):
            spaces = " " * round(len(probe_names[i])/2)
            #print(spaces, probe_names[i])
            temp_name = probe_names[i] + "\n" + spaces + "(" + polymers[i] + ")"
            probe_names[i] = temp_name
    
    if not len(probe_names) < len(index):
        plt.xticks(index, probe_names)
        if display_full_probe_name:
            ax2.set_xticklabels(probe_names, rotation = 50, ha = "right",fontdict={'fontsize': 12, 'fontweight': 'light', 'fontname': 'Cambria'})
        else:
            ax2.set_xticklabels(probe_names)
    
    
    """
    ADD functionality of limiting y-axis to height of graph + error
    
    ADD Tukey significance codes here
    
    Old format:
    plt.ylim(0, max(max(group1_data), max(group2_data), max(group3_data), max(group4_data))*1.4)
    """
    
    
        
    ax1.legend(loc = "upper left", bbox_to_anchor = (1,1)) #USE bbox_to_anchor = (1,1) for legend outside of plot frame; USE loc = "upper/lower left/right" to change location
    
    if polymer_class == "AGP" or polymer_class == "cellulose and hemicellulose":
        ax1.set_ylim(0, 125)
        ax2.set_ylim(0, 125)
    else:
        ax1.set_ylim(0, y_upper)
        ax2.set_ylim(0, y_upper)
    
    if save_barplots:
        settings = "__normalised__"

        fig_name = sample_class + " " + polymer_class + " - " + settings + ".png"    
        
        create_folder_if_not_exists(sample_class)
        folder = os.getcwd() + "\\" + sample_class + "\\"
        name = folder + "____FINAL" + fig_name

        plt.savefig(name, dpi = 300, bbox_inches = "tight")

def bar_plot_2_groups_extractions_combined(sample_class, polymer_class, dataframe_subplots, group_ranges:list):
    save_barplots = True
    
    combine_extractions_barplot = True
    
    display_full_probe_name = True
    display_additional_probe_code = True # This will display the antibody code underneath the target name on the x-axis
    
    label_x_pos = 0.05
    
    if True: # Indicates which polymers to plot on bar plots based on polymer class supplied to this function
        if polymer_class.lower() == "pectin": # RU1 removed
            label_x_pos = 0.5
            polymers = [ 
                "JIM5",
                "JIM7",
                "LM18",
                "LM19",
                "LM20",
                "LM7",
                "RU2",
                "LM5",
                "LM6",
                "LM16"]
            
            
        elif polymer_class.lower() == "agp":
            label_x_pos = 0.7
            polymers = [
                "JIM8",
                "JIM13",
                "JIM15",
                "JIM16",
                "LM2"]

        elif polymer_class.lower() == "cellulose and hemicellulose": #BS-400-3, LM24, and CBM3a removed
            label_x_pos = 0.65    
            polymers = [
                "CBM3a",
                "BS-400-2",
         
                "LM21",
                "LM15",

                "LM25",
                "LM11",
                "LM23"]
        else:
            print("INVALID POLYMER CLASS!")
            print(polymer_class, "not found")
            sys.exit()

    
    ##########################################################################
    #                   D A T A     S E T T I N G S                          #
    ##########################################################################


    # This will be used to plot a varable number of probes on a single plot without needing manual adjustment
    index = np.arange(len(polymers))
    
    group1_range = group_ranges[0]
    group2_range = group_ranges[1]

    
    
    group1_data = []
    group1_error = []
    group1_negative_data = []
    group1_negative_error = []
    
    group2_data = []
    group2_error = []
    group2_negative_data = []
    group2_negative_error = []


    # TUKEY'S HSD TESTING
    CDTA_tukey_groups = []
    NaOH_tukey_groups = []
    significance_pos = []
    significance_neg = []
    
    
    
    
    for df_no, dataframe in enumerate(dataframe_subplots):
        
        if display_full_probe_name:
            probe_names = []
        else:
            probe_names = polymers[:]
        
        
        for column in polymers:  # Use the alternative version of this loop-conditional pair if anything goes wrong with sequence and labels on bar plots
            if column in dataframe_subplots[dataframe]:
                
        # for column in dataframe: ###### REUSE THIS if anything is broken with bar plots and order
        #     if column in polymers:
                
            
                data_1 = []
                data_2 = []

                error_1 = []
                error_2 = []

                
                ### Here we get out full epitope name out of the list of antibodies defined at the start of the program
                if display_full_probe_name:
                    name_check = column + ")"
                    for name in antibodies:
                        if name_check in name:
                            trimmed_name = name[:-len(name_check)-2]
                            probe_names.append(trimmed_name)

                
                
                
                #Group1
                group1 = dataframe_subplots[dataframe][column].iloc[group1_range[0]: group1_range[1]]
                group1_list = list(group1)
                
                # Here we account for outliers that are marked as "-1"
                if -1 in group1_list:
                    #input("STOP STOP STOP STOP")
                    while -1 in group1_list:
                        group1_list.remove(-1)
                    print("group 1:", group1_list)
                    #input("STOP STOP STOP STOP")
                
                standard_error_group1 = stats.sem(group1_list)
                group1_average = np.mean(group1_list)
                data_1.append(group1_average)
                error_1.append(standard_error_group1)

                
                
                #Group2
                group2 = dataframe_subplots[dataframe][column].iloc[group2_range[0]: group2_range[1]]
                group2_list = list(group2)
                
                # Here we account for outliers that are marked as "-1"
                if -1 in group2_list:
                    #input("STOP STOP STOP STOP")
                    while -1 in group2_list:
                        group2_list.remove(-1)
                    print("group 2:", group2_list)
                    #input("STOP STOP STOP STOP")
                    
                group2_average = np.mean(group2_list)
                standard_error_group2 = stats.sem(group2_list)
                data_2.append(group2_average)
                error_2.append(standard_error_group2)
                #print(column, group2_name+" average:", group2_average)
                
                
                
                
                
                if True: # P-VALUE CALCULATIONS
                    #### Below we determine SIGNIFICANT DIFFERENCES
                    t_stat, p_value = stats.ttest_ind(group1_list, group2_list, equal_var=False)
                    

                    if df_no == 0:
                        if p_value < 0.01:
                            asterisk = "**"
                            significance_pos += ["**"]
                        elif p_value < 0.05:
                            asterisk = "*"
                            significance_pos += ["*"]
                        elif p_value > 0.05 and p_value < 0.1: # This condition is useful for finding values that almost made the cut-off. Perhaps evaluate the raw data and remove outliers if applicable.
                            significance_pos += [""]
                        else:
                            asterisk = ""
                            significance_pos += [""]
                    
                    elif df_no == 1:
                        # Negative bars 
                        t_stat, p_value = stats.ttest_ind(group1_list, group2_list, equal_var=False)
                        
                        # DEBUGGING
                        """
                        significance += [str(p_value)[:4]]
                        """
                        if p_value < 0.01:
                            asterisk = "**"
                            significance_neg += ["**"]
                        elif p_value < 0.05:
                            asterisk = "*"
                            significance_neg += ["*"]
                        elif p_value > 0.05 and p_value < 0.1: # This condition is useful for finding values that almost made the cut-off. Perhaps evaluate the raw data and remove outliers if applicable.
                            significance_neg += [""]
                        else:
                            asterisk = ""
                            significance_neg += [""]

                

                
                # This removes really small values from significance indication on bar plots
                for i in range(len(group1_list)):
                    if group1_list[i] < 5:
                        group1_list[i] = 0
                for i in range(len(group2_list)):
                    if group2_list[i] < 5:
                        group2_list[i] = 0

                              
                if df_no == 0:
                    group1_data += data_1
                    group1_error += error_1
                    group2_data += data_2
                    group2_error += error_2

                    tukey_codes = tukey_HSD([group1_list, group2_list])
                    CDTA_tukey_groups.append(tukey_codes)
                    
                    
                elif df_no == 1:
                    group1_negative_data += data_1
                    group1_negative_error += error_1
                    group2_negative_data += data_2
                    group2_negative_error += error_2

                    
                    tukey_codes = tukey_HSD([group1_list, group2_list])
                    NaOH_tukey_groups.append(tukey_codes)
            else:
                print(column, " ERROR (*o*)")
           
    

    # The following will take each group (the CDTA and NaOH data for root, stem, and leaf for a SINGLE CULTIVAR and normalise it to 100)
    if normalise:
        for i in range(len(group1_data)):
            multiplier = 100 / max([group1_data[i], group2_data[i], group1_negative_data[i], group2_negative_data[i]])
            group1_data[i] *= multiplier
            group2_data[i] *= multiplier

            group1_negative_data[i] *= multiplier 
            group2_negative_data[i] *= multiplier

    
            group1_error[i] *= multiplier
            group2_error[i] *= multiplier

            group1_negative_error[i] *= multiplier 
            group2_negative_error[i] *= multiplier

    

    
    # Here all NaOH values are multiplied by (-1) so that we can plot it as negative values on our bar plots
    for data_list in [group1_negative_data, group2_negative_data]:
        for i in range(len(data_list)):
            data_list[i] *= 1
    
    # Here we find the upper and lower limits of the bar plot figure by determining the max height of the data and error bars combines (allowing for +10 for significance codes to be added)
    

    y_upper = 120
    y_lower = -120

    for i in range(len(group1_data)):
        if group1_data[i] + group1_error[i] + 10 > y_upper:
            y_upper = group1_data[i] + group1_error[i] + 10
        if group2_data[i] + group2_error[i] + 10 > y_upper:
            y_upper = group2_data[i] + group2_error[i] + 10


    for i in range(len(group1_data)):
        if group1_negative_data[i] - group1_negative_error[i] - 10 < y_lower:
            y_upper = group1_negative_data[i] - group1_negative_error[i] - 10
        if group2_negative_data[i] - group2_negative_error[i] - 10 < y_lower:
            y_upper = group2_negative_data[i] - group2_negative_error[i] - 10



        
            
    ##########################################################################
    #               F I G U R E     S E T T I N G S                          #
    ##########################################################################
            
    plt.figure()  
    fig, (ax1,ax2) = plt.subplots(2,1, sharex = True, figsize = (len(polymers)*1.5, 5)) # This will auto-adjust the figure to the right size. Play around with this if you want other figure dimentions

    #plt.style.use("default")
    #fig.patch.set_facecolor("white")
    
    #print("group 1 error", group1_error)
    #print("group 1 data", group1_data)


    name1 = "Skin" 
    name2 = "Flesh"
    name3 = "Root"
    
    bar_colours = ["darkseagreen",
                   "orange",
                   "cornflowerblue"]
    
    # cute pastel purple and blue
    bar_colours = ["#534B62",
                   "#A499B3",
                   "#D0BCD5"]
    
    bar_colours = ["#06A77D",
                   "#005377",
                   "#052F5F"]
    
       
    bar_colours = ["#5C80BC",
                   "#30323D",
                   "#4D5061"]
    
    # works well with #DBA159 (earth yellow - BEST SO FAR) as NaOH background
    # and UCLA blue #2274A5
    # and sepia #6D3D14
    # SECOND BEST - #8B9EB7 cadet grey
    bar_colours = ["purple", 
                   "white",
                   "#297045"]


    bar_width = 0.22            
    ax1.bar(index - 0.5*bar_width,  group1_data,  bar_width, yerr = group1_error, color = bar_colours[0], label = name1)
    ax1.bar(index + 0.5*bar_width,  group2_data,  bar_width, yerr = group2_error, color = bar_colours[1], label = name2, edgecolor="black")          
    #ax.axhline(y=0, color="black", linestyle="-", linewidth = 0.4)
    
    # Add text for CDTA/NaOH
    plt.text(len(polymers)-label_x_pos, 190, "CDTA", rotation=90)
    plt.text(len(polymers)-label_x_pos, 60, "NaOH", rotation=90)
    
    
    #Remove axis spines
    # ax.spines["bottom"].set_visible(False)
    # ax.spines["right"].set_visible(False)
    
    
    # EXPERIMENTAL - adding negative equivalents
    if combine_extractions_barplot:
        bar_width = 0.22
        ax2.bar(index - 0.5*bar_width,  group1_negative_data,  bar_width, yerr = group1_negative_error, color = bar_colours[0])
        ax2.bar(index + 0.5*bar_width,  group2_negative_data,  bar_width, yerr = group2_negative_error, color = bar_colours[1], edgecolor="black")
    
    # ADD SIGNIFICANCE INDICATORS
    # for i in range(len(polymers)): # This add asterisks to the significant pos groups
    #     plt.text(index[i], max(group2_data[i], group1_data[i]) + max(group1_error[i], group2_error[i])/2, significance_pos[i], ha = "center", va = "bottom", fontdict = {'fontsize': 20, 'fontweight': 'bold', 'fontname': 'Cambria'})
    # for i in range(len(polymers)): # This add asterisks to the significant neg groups
    #     plt.text(index[i], min(group2_negative_data[i], group1_negative_data[i]) - 40, significance_neg[i], ha = "center", va = "bottom", fontdict = {'fontsize': 20, 'fontweight': 'bold', 'fontname': 'Cambria', 'color': 'black'})
    
    
    for i in range(len(polymers)): # This add asterisks to the significant pos groups
        #plt.text(index[i], max(group2_data[i], group1_data[i]) + max(group1_error[i], group2_error[i])/2, significance_pos[i], ha = "center", va = "bottom", fontdict = {'fontsize': 20, 'fontweight': 'bold', 'fontname': 'Cambria'})
        ax1.annotate(significance_pos[i], (index[i], max(group2_data[i], group1_data[i]) + max(group1_error[i], group2_error[i])/2), ha="center",color='black', size=16)
    for i in range(len(polymers)): # This add asterisks to the significant neg groups
        #plt.text(index[i], min(group2_negative_data[i], group1_negative_data[i]) - 40, significance_neg[i], ha = "center", va = "bottom", fontdict = {'fontsize': 20, 'fontweight': 'bold', 'fontname': 'Cambria', 'color': 'black'})
        ax2.annotate(significance_neg[i], (index[i],   max(group2_negative_data[i], group1_negative_data[i]) + max(group1_negative_error[i], group2_negative_error[i])/2), ha="center",color='black', size=16)
        
    
    # Background colour
    #y_lower, y_upper = ax.get_ylim()
    #y_lower = -145
    y_upper = 145
    #ax.axhspan(y_lower, 0, facecolor="#5EB1BF", alpha=0.5, zorder=0)
    #ax.axhspan(0, y_upper, facecolor="white", alpha=0.5, zorder=0)

    
    title = "Relative levels of " + polymer_class + " in eggplant fruit (BB)"
    ax1.set_title(title, fontdict={'fontsize': 17, 'fontweight': 'bold', 'fontname': 'Cambria'})
    #ax.set_ylabel("Relative absorbance at 450 nm", fontdict={'fontsize': 12, 'fontweight': 'bold', 'fontname': 'Cambria'})
    ax1.tick_params(labeltop=False, top=False)
    
    
    xxx = 0.05
    if polymer_class.lower() == "pectin":
        xxx = 0.075
    elif polymer_class.lower() == "agp":
        xxx = 0.04
    fig.text(xxx, 0.5, 'Relative absorbance', va='center', rotation='vertical', fontdict={'fontsize': 15, 'fontweight': 'bold', 'fontname': 'Cambria'})
    
    # The code below adds the antibody abbreviation alongside the target
    if display_additional_probe_code:
        for i in range(len(probe_names)):
            spaces = " " * round(len(probe_names[i])/2)
            #print(spaces, probe_names[i])
            temp_name = probe_names[i] + "\n" + spaces + "(" + polymers[i] + ")"
            probe_names[i] = temp_name
    
    if not len(probe_names) < len(index):
        plt.xticks(index, probe_names)
        if display_full_probe_name:
            ax2.set_xticklabels(probe_names, rotation = 50, ha = "right", fontdict={'fontsize': 11, 'fontweight': 'light', 'fontname': 'Cambria'})
        else:
            ax2.set_xticklabels(probe_names)
    
    
    """
    ADD functionality of limiting y-axis to height of graph + error
    
    ADD Tukey significance codes here
    
    Old format:
    plt.ylim(0, max(max(group1_data), max(group2_data), max(group3_data), max(group4_data))*1.4)
    """
    
    
        
    ax1.legend(loc = "upper left", bbox_to_anchor = (1,1)) #USE bbox_to_anchor = (1,1) for legend outside of plot frame; USE loc = "upper/lower left/right" to change location
    ax1.set_ylim(0, 130)
    ax2.set_ylim(0, 130)
    
    if save_barplots:
        settings = "__normalised__"

        fig_name = sample_class + " " + polymer_class + " - " + settings + ".png"    
        
        create_folder_if_not_exists(sample_class)
        folder = os.getcwd() + "\\" + sample_class + "\\"
        name = folder + "__FINAL__" + fig_name

        plt.savefig(name, dpi = 300, bbox_inches = "tight")

def initiate_barplots():
    #TODO - use this function to call barplot_paired_3_groups, 1, 2 depending on barplot_bar_count
    #TODO - add  .xlsx file on which to indicate groups (different sheets)
    global barplot_df_dict
    
    xl = pd.ExcelFile("barplot_groups.xlsx")

    sheet_list = xl.sheet_names # see all sheet names
    barplot_df_dict = {}
    for i, name in enumerate(sheet_list):
        barplot_df_dict[name] = pd.read_excel("barplot_groups.xlsx", sheet_name = i, engine="openpyxl")
    
    global polymer_class_names
    polymer_class_names = []
    for key in barplot_df_dict:
        print(barplot_df_dict[key].iloc[:, 0].tolist())
        polymer_class_names += [key]
    
    
    
    # Here we decide how many pairs of plots to make
    if barplot_bar_count == 2:
        # indexes correspond to [(bar1), (bar2)]
        sample_indexes_leaf_stem_root = {
            "Black Beauty": [(0,5), (5,10)]
            }
    
        plt.style.use("default") # Make sure to call this function once before calling functions using matplotlib
        
        ### 3 groups compared per antibody in each antibody class
        for key in sample_indexes_leaf_stem_root:
            sample_class = key
            group_ranges = sample_indexes_leaf_stem_root[key]
            for polymer in ["pectin", "cellulose and hemicellulose", "AGP"]:#polymer_class_names:#polymer_classes:
                bar_plot_2_groups_extractions_combined(sample_class, polymer, subplots, group_ranges)
    
    elif barplot_bar_count == 3:
    
        # indexes correspond to [(Leaf), (Stem), (Root)]
        sample_indexes_leaf_stem_root = {
            "group1": [(0,4), (8,12), (16, 20)], 
            "group2":   [(4,8), (12,16), (20,24)]
            }        
        
        plt.style.use("default") # Make sure to call this function once before calling functions using matplotlib
        ### 3 groups compared per antibody in each antibody class
        for key in sample_indexes_leaf_stem_root:
            sample_class = key
            ranges = sample_indexes_leaf_stem_root[key]
            for polymer in ["pectin", "AGP", "cellulose and hemicellulose"]:
                bar_plot_3_groups_extractions_combined(sample_class, polymer, subplots, ranges)

if __name__ == "__main__":
    app = MyApp()
    print(logo)
    app.protocol("WM_DELETE_WINDOW", app.on_closing)
    app.mainloop()
