# Single-cell analysis of the cellular heterogeneity and interactions in the injured mouse spinal cord

This repository accompanies the paper **Single-cell analysis of the cellular heterogeneity and interactions in the injured mouse spinal cord** by Milich, Choi, et al. Paper can be found here: [link](https://doi.org/10.1084/jem.20210040). 

Code used for analysis in producing results presented in the paper can be found in the *scripts* folder. Code used to build the shinyApp site can be found in the *mouseSCI_2021* folder. (Apologies for the disorganized code, it will eventually be reorganized.)

If you have any questions regarding the study, data, or code, feel free to contact James Choi or Dr Jae Lee ([lab website](https://www.jaeleelab.com/)).


## Data Availability

scRNAseq data for this study is available at NCBI GEO accession no. GSE162610.

Full gene x cell count matrix for this study can be found here: https://drive.google.com/drive/u/0/folders/1AUbh2A3ogbDBaWlhf-x3UGjLno-__M-z

* **x_sci.rds** is a dgCMatrix class object that can be read into R. The matrix contains the full count matrix (un-normalized, raw counts). Column names contain unique cell barcodes and row names contain gene names.
* **obs_sci.csv** contains cell-level metadata such as sample ID, injury time-point, cell-type classifications, etc. All cluster analysis results from the paper can be found in this spreadsheet.
* **vars_sci.csv** contains gene-level metadata (empty).


## Using the MouseSCI_2021 portal

The mouseSCI scRNAseq data can been viewed in this [browser](https://jaeleelab.shinyapps.io/mouseSCI_2021/). Please be patient since data sizes are large (> 1Gb). Chrome users may have receive a "Page unresponsive" message - click "Wait" and site will eventually load.

NOTE: We are using a free version of shinyapp.io so there are restrictions on data upload size. **25% of cells are sampled and displayed**. 
