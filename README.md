# Spatial-measures-of-TILs-maps-heterogeneity-to-predict-cancer-survival

Following code was made as a part of master thesis reaserch. TIL (Tumor-Infiltrating Lymphocytes) are belived to have an important impact on patients survival in
cancer prognosis.
However simple densities (of TILs) are not a sufficient factor because other factors such as TIL localization or it they form any types of structures (e.g. they surround the cancer cells) have influence on the results. Having that in mind, this work tries to "rank" the effects of TIL on patients results (if any)
by utilizing 3 different methods for structure analysis. 
The used methods are: Grey Level Co-occurance Matrix, Spatial Chaos and Affinity propagarion. The research was conducted on 7 cancer types: BLCA, BRCA, CESC, LUAD, LUSC, PAAD and PRAD


Materials and other informations needed by the algorithm:

- TIL maps (are available under the adress: https://wiki.cancerimagingarchive.net/pages/viewpage.action?pageId=33948919)

- The excel file that contiants information about patiets e.g. vital status, age etc. (available under the adress: https://gdc.cancer.gov/about-data/publications/PanCan-Clinical-2018). 

- The naming needs to be unchanged (both names of TIL images as well as the collumns in the excel file), otherwise additonal changes in the code will be required. 

- The matlab code requires the Statistical Toolbox to be installed in order to run it and the R code needs clusterCrit and apcluster.

Both R and Matlab codes requires to provide a path to the files (i.e. TIL maps location or path to excel files),
the correct path should be put in appropriate fragments of code. For the R code its specifing where to save the images (coordinates and the ap results)
and the indices results.
The naming scheme of TIL maps should be unchanged otherwise the substr will need to be altered


Additional informations

Chunks of code in the GLCM file are taken from EXIMS research conducted by As Chalini D. et al. (https://sourceforge.net/projects/exims/files/)
Code from the folder called "Spatial chaos original source code" is taken directly from authors (Andreas Bartels, Theodore Alexandrov, University of Bremen) github pageand its available under the adress (https://github.com/alexandrovteam/IMS_quality/tree/master/codebase), 
the reason why those files are included in this repository is to make sure all the necessery files are stored in one place. 
However if the authors are not happy about me posting their code here, please contanct me and i will make appropiate changes to the repository.

Results are presented in the folder called "Results", each folder hold 3 excel files - one for spatial Chaos method, 
one for GLCM method and one for Affinity Propagation. More details about the program as well as the results disscusion can be found in the thesis 
pdf file located under the "Thesis" folder
