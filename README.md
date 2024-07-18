Welcome to VizELISA 1.2!
A tool for automating the file handling, data processing, statistics and visualisation of ELISA data


.::         .::            .::::::::.::      .::  .:: ::        .:       
 .::       .::  .:         .::      .::      .::.::    .::     .: ::     
  .::     .::     .:::: .::.::      .::      .:: .::          .:  .::    
   .::   .::   .::     .:: .::::::  .::      .::   .::       .::   .::   
    .:: .::    .::   .::   .::      .::      .::      .::   .:::::: .::  
     .::::     .::  .::    .::      .::      .::.::    .:: .::       .:: 
      .::      .::.::::::::.::::::::.::::::::.::  .:: ::  .::         .::   

For the best reading experience of this text file, turn off text wrapping. 





This program is designed to process and visualise ELISA glycan array data. 
The program is split into 4 sections:
	* Data Input
			- this is where instructions and links to important files can be found. You may also validate their data to make sure the formatting is correct.
	* Quick Parameters
			- this section allows for the quick use of checkboxes for you to indicate preferences such as outlier removal or bolding of figure titles.
	* Graph Specifications
			- this section allows the user to manually enter/change values relating to the physical characteristics of the output figures.
	* Standard Output
			- this section houses Python's standard output and updates/errors will be given here along every step of the analysis.





##############################
		    USER GUIDE
##############################

Note: Testing data has been provided under the data folder, mapping file, bar plot groups and excluded antibodies for testing purposes.


Before use, you will need to install the following: 

Python 3.11 or higher and the following libraries:
	* customtkinter
	* matplotlib
	* numpy
	* pandas
	* scipy
	* seaborn
	* statsmodels
	* tkinter


Standard Operating Procedure:

	1) Click on the Open Data Mapping File (under Data Input) and indicate the locations of your samples (and blanks) on the 96-well plate map. Note that different sheets are used for different extractions (e.g., CDTA and NaOH)
		* Note that numbered suffixes are used to indicate extractions. E.g., "LM18_1" and "LM18_2" could refer to cell walls extracted using CDTA and NaOH, respectively, before being tested with LM18.
	2) Change all the parameters necessary for the final heatmap output. This will be under Data Input and Quick Parameters. You may validate your data here (optional).
	3) Click the Load Files button under Quick Parameters. This will extract the numerical information from all files in the "data" folder. Common errors will be met with advice on how to solve it.
		* Note: a new folder will be generated each time Load Data is clicked.
	4) Adjust the parameters as required under Graph Specifications.
	5) Click on the Confirm Parameters button. This will generate a low-resolution figure (as a pop-up) for you to evaluate your chosen parameters. Note: close the pop-up before returning to VizELISA or the program will crash.
	6) Click the Generate Final Figure button under Standard Output to generate your heatmap.




##############################
	  EXPERIMENTAL FEATURES
##############################

Paired bar plots (2 or 3 grouped bars per probe) can be generated using the same parameter input fields as heatmaps. Click on Generate Barplots (under Graph Specifications) after loading files under Quick Parameters.
Note: It is not necessary to click on Generate Final Figure for barplots. However, after clicking on Generate Barplots, VizELISA will freeze temporarily, however the figures will still be output in the output folder marked with the time of analysis.






##############################
	 User-changeable parameters
##############################

DATA INPUT

Instructions                |   A button that opens this help text file
Open data mapping file      |   A button that opens the template file on which sample names should be indicated
Check data validity         |   (OPTIONAL) This checks whether data provided is in the right formats and is ready to be processed
Antibody list               |   This opens a text file containing common epitopes and their corresponding probes. Use the same format to add any antibodies that are missing, i.e., epitopeName(probeCode)
Polymer class groupings     |   This opens a speadsheet where different sheets represent polymer classes (e.g., pectin, AGP, etc.). Feel free to change the sheet name or lists of antibodies.
Exclude antibodies          |   This will temporarily remove the indicated probes from the generation of heatmaps and bar plots
Custom antibody order       |   This allows the user to indicate a custom order in which to display probes on the generated heatmap
Polymer class colours       |   This opens a spreadsheet where the sheet name represents the colour in which the probe name will be written on the heatmap if the user checks the "Colour code antibodies by class" checkbox. 
Bulk file renaming settings |   This text file can be used to set parameters for renaming files, i.e., adding or removing prefixes or suffixes to a particular file type.
FINAL file renaming         |   This commits the above renaming settings and applies them to the relevant files in the subfolder "to_rename"
List of valid colours       |   This is a list of colour schemes that may be used on the generated heatmap. 
Global normalisation        |   This will apply CoMPP-style normalisation to heatmap data, i.e., the max in the whole table will be assigned a value of 100 and all other values will be adjusted relatively thereto.
Local normalisation         |   This will normalise each column (for all extractions) individually for better comparisons between samples without high signals drowning out lower signals.
Remove extraction labels    |   This will remove the indicated extraction labels from sample names, e.g., CDTA_leaf1 -> leaf1
A1 position on 96-well plate|   This is the value on each data file that corresponds to the A1 well on a 96-well plate. E.g., B11, C13
Plate outlier label         |   This is a text label that may be used to remove outlier technical or biological outliers from raw data files. Replace an outlier with this label and VizELISA will handle it accordingly.
Number of extractions       |   The number of extractions being compared, e.g., this value would be 2 if CDTA and NaOH extractions were used
Number of template sheets   |   The number of sheets correspond to the number of plates of which all samples are distributed. Create more sheets if necessary by copying sheet 1, 2, or 3
Extraction label 1          |   Label used to label samples, e.g., this would be "NaOH" if samples are labeled "NaOH_leaf1", "NaOH_leaf2"
Extraction label 2          |   Label used to label samples, e.g., this would be "CDTA" if samples are labeled "CDTA_leaf1", "CDTA_leaf2"
Extraction label 3          |   Label used to label samples, e.g., this would be "NaOH" if samples are labeled "NaOH_leaf1", "NaOH_leaf2"
% of blank cut-off          |   This setting will remove all data that is not at least x% higheer than the background reading, E.g., blank=0.1, 20% cut-off means that values of 0.12 and lower will be excluded from analysis.
Remove values below         |   This allows you to remove values on the heatmap below a particular value, typically 5 as per CoMPP
P<0.05                      |   A symbol that indicates a significant difference with 0.01 < P < 0.05, e.g., *
P<0.01                      |   A symbol that indicates a significant difference with P < 0.0, e.g., **
Sample indices              |   Barplot groups will be assigned using these values. Each bar is represented by a set of () where the two values represents the start and end position of the samples on the indicated sample order list (on the mapping file). E.g., for samples = leaf, leaf, leaf, root, root, root -> bars = [(1,3), (4,6)]. Only use the number of brackets corresponding the number of bars you want to plot.
Correlation range           |   This will select a range of samples on which to do correlation analysis followed by visualisation per correlation heatmap



QUICK PARAMETERS

Remove numbers from heatmap             |   Hides numerical values from heatmaps, leaving only colours to indicate signal intensity
Use median instead of mean              |   The median between technical replicates will  be used 
Pearson correlation (default: Spearman) |   Indicat which method should be used to determine correlation
Full epitope name on label              |   This will include the target epitope instead of just the probe code
Colour code antibodies by class         |   Colours will be used to distinguish polymer classes. Colours may be changed by clicking "Polymer class colours"
Bold graph title                        |   Bolds the title for heatmaps and barplots
Bold axis titles                        |   Bolds the axis labels for heatmaps and barplots
Load Files                              |   This initiates the data processing and creates a data log and output folder.



GRAPH SPECIFICATIONS

Graph title             |   Title to be used for heatmaps and barplots
Fig size                |   Dimensions of the final heatmap in the format (x, y)
X-axis title            |   Label to be used for the x-axis
Y-axis title            |   Label to be used for the y-axis
Title font size         |   Font size to be used for figure titles
Axis font size          |   Font size to be used for figure axis labels
X-label rotation        |   Degrees by which the x-axis is rotated for readability 
Y-label rotation        |   Degrees by which the y-axis is rotated for readability
Colour scheme           |   Colour scheme of the generated heatmap. This can be changed under Data Input
BARPLOT #pairs          |   This corresponds the the number of barplot sample indices indicated under Data Input. E.g., this value would be 2 if only 2 of the 3 brackets were used -> [(1,5), (6,10), (x,x)]
Save image as           |   Here a custom name can be chosen for generated heatmaps
Image resolution        |   The DPI of the final high-quality rendered figure. Values above 200 are recommended
Confirm parameters      |   Button that generates a TEMPORARY low-resolution depiction of the final heatmap for the user to adjust parameters before generating the final figure
Generate barplots       |   (EXPERIMENTAL) Barplots can be generated for 2 or 3 groups for be compared once files have been loaded
Correlation table       |   (EXPERIMENTAL) A correlation heatmap can be generated for a custom set of samples using either Spearman or Pearson correlation, depending on the user's preference.


STANDARD OUTPUT

Generate Final Figure   |   This will generate a final, high-quality heatmap and save it in the output folder.






Author: Eugene Graham Badenhorst

Version: 1.1.3

Date: 31/12/2024

License: MIT

Contact: ebaden@sun.ac.za
