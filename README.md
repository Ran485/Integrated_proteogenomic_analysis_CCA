# Integrated Proteogenomic Characterization of Cholangiocarcinoma

![Schematic representation of the experimental design](https://github.com/Ran485/Integrated_proteogenomic_analysis_CCA/blob/main/Schematic_workflow.png)

This github repository contains the data files and analysis code used to generate the figures for the manuscript **"Proteogenomic Characterization of Cholangiocarcinoma"** published in _Hepatology_. The Schematic representation of the experimental design and data analysis is shown above.

## Dependencies

Data analysis was performed using Mac platform in Python (version 3.8.6) and R (version 4.1.0) with open-source libraries. The detail library dependencies are not listed here, they should be resolved automatically when you install the dependencies we mentioned in specific analysis scripts.

## Folders Structure

The files are organised into five folders:

- [***data***](https://github.com/Ran485/Integrated_proteogenomic_analysis_CCA/tree/main/data): which contains all the genomic, transcriptomic, proteomics, phosphoproteomic and clinical patient informations required to perform the analyses described in the paper. The data files is currently deposited to the zenodo repository and can be available on [*Supplementary_Tables*](https://zenodo.org/record/6536180#.YnsscxMza1s) link.
- [***R***](https://github.com/Ran485/Integrated_proteogenomic_analysis_CCA/tree/main/R): which contains the R code to reproduce all analyses and generate the the figures in our manuacript.
- [***python***](https://github.com/Ran485/Integrated_proteogenomic_analysis_CCA/tree/main/python): there are python modules that contain the bulk of the code mainly used for data integration and preprocessing before performing further data analysis, and also including a few figures generation. 
- [***shell***](https://github.com/Ran485/Integrated_proteogenomic_analysis_CCA/tree/main/shell): which contains the shell scripts to calling the gene fusions.
- [***resources***](https://github.com/Ran485/Integrated_proteogenomic_analysis_CCA/tree/main/resources): which contains additional public databases (eg: CancerDriveGene lists, DrugTarget, Kinase_Substrate_Dataset and etc.) which are required by the R and python scripts.

## Reference

1. Deng, M., Ran, P., Chen, L., Wang, Y., Yu, Z., Cai, K., Feng, J., Qin, Z., Yin, Y., Tan, S., et al. (2022). Proteogenomic characterization of cholangiocarcinoma. Hepatology, hep.32624. 10.1002/hep.32624.
