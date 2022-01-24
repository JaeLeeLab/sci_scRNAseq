# Single-cell analysis of the cellular heterogeneity and interactions in the injured mouse spinal cord

This repository accompanies the paper **Single-cell analysis of the cellular heterogeneity and interactions in the injured mouse spinal cord** by Milich and Choi, et al. Paper can be found here: [link](https://doi.org/10.1084/jem.20210040). 

Code used for analysis in producing results presented in the paper can be found in the *scripts* folder. Code used to build the shinyApp site can be found in the *mouseSCI_2021* folder. 

If you have any questions regarding this repository, study, data, or code, please contact James Choi or Dr Jae Lee (corresponding author, [lab website](https://www.jaeleelab.com/)).


## Data Availability

Raw data are available from the SRA (Sequence Read Archive) database under study accession [SRP295673](https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP295673). Processed scRNAseq data are available at NCBI GEO accession [GSE162610](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162610).

Direct downloads for Seurat object and gene x cell count matrix can be found here: https://drive.google.com/drive/u/0/folders/1AUbh2A3ogbDBaWlhf-x3UGjLno-__M-z

* **x_sci.rds** is a dgCMatrix class object that can be read into R. The matrix contains the full count matrix (un-normalized, raw counts). Column names contain unique cell barcodes and row names contain gene names.
* **obs_sci.csv** contains cell-level metadata such as sample ID, injury time-point, cell-type classifications, etc. All cluster analysis results from the paper can be found in this spreadsheet.
* **vars_sci.csv** contains gene-level metadata (empty).
* **sci.rds** contains the SeuratObject for the full SCI dataset.


## SCI_singlecell web portal

A web portal is available to more easily browse through the gene expression data. Portal is available at https://jaeleelab.shinyapps.io/SCI_singlecell
 
Code repository to build the ShinyApp portal is avaiable here: https://github.com/JaeLeeLab/SCI_SingleCell_portal


## Acknowledgements

This study was funded by National Institute of Neurological Disorders and Stroke grant R01NS081040 (to J.K. Lee), University of Miami SAC Award UM SIP 2019-2 (to J.K. Lee), National Institute of Neurological Disorders and Stroke grant F31NS115225 (to C. Ryan), the Miami Project to Cure Paralysis, the Buoniconti Fund to Cure Paralysis, a generous gift from the Craig H. Nielsen Foundation (J.K. Lee), and National Institutes of Health grant 1S10OD023579-01 for the VS120 Slide Scanner housed at the University of Miami Miller School of Medicine Analytical Imaging Core Facility.
