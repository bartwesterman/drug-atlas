# drug-atlas
Synergy prediction and the drug atlas

## Motivation and model
The objective is to predict drug synergy by using not only sensitivity data of different cell lines to individual drugs, but also the drug atlas distance. It may be also of interest to compare predictions with and without using the drug atlas distance, so as to demonstrate its added value.

## The data
Each row in the data corresponds to values for one cell line and a pair of drugs. As such, there is a lot of structure in the data, via the drug pairs. In addition to the distance between the two drugs from the drug atlas, the data includes individual drug sensitivity value as well as target information relating to each drug of the pair, where target information indicates whether or not the drug targets a gene known to be affected/mutated in the cell line at hand. Finally, it includes an indicator variable if the two drugs are known to display synergy.

## Application
The idea is to fit the model to the cell lines for which synergy is known. Note that there is no information if no synergy is also known. We will only use the cell lines for which sensitivity is available for both drugs. At this point, the other ones are not informative.

~~~ mydir.base <- "/media/renee/Seagate Expansion Drive/Projects/" 
mywd <- paste(mydir.base, "westerman/NatComm model", sep="")
source(paste(mydir.base,"functions_wilcoxon_t_logistic_fdr.R",sep="/")) # loads functions for FDR and graphs
source(paste(mydir.base,"mysplit.R",sep="/")) # loads functions for mt cor and graphs
source(paste(mydir.base,"plot_density_percolumn.R",sep="/")) # function to make graphs of densities per column of a data matrix, with proper limits  
source(paste(mydir.base,"var_in_colour.R",sep="/")) # loads functions for mt cor and graphs
source(paste(mydir.base,"functions_sens.R",sep="/")) # loads functions for mt cor and graphs
source(paste(mydir.base,"f_roc.R",sep="/")) # loads functions for mt cor and graphs

library(gplots)
library(DescTools)

cl.data <- read.delim(paste(mywd, "/distance_sensitivity-model_input_data_literaturevsbreast.csv", sep=""), sep=",")
cl.data$Ward_Drug_Distance <- NULL
~~~ 
We start by extracting the tissue, cell line ID and drug names.
~~~ 
id <- as.character(cl.data$Identifyer)
f.tissue <- sapply(1:nrow(cl.data), mysplit, char.vector = id, split.by = "_", result.sel = 1)
f.cl  <- sapply(1:nrow(cl.data), mysplit, char.vector = id, split.by = "_", result.sel = 2)
f.d1 <-  sapply(1:nrow(cl.data), mysplit, char.vector = id, split.by = "_", result.sel = 4)
f.d2 <-  sapply(1:nrow(cl.data), mysplit, char.vector = id, split.by = "_", result.sel = 5)
~~~ 

First we notice that 3 rows yield a warning (not shown), for all variables. All results have the correct length, so the rows with problems likely generated NAs. Indeed, there are 3 NAs in the variables, and the corresponding data rows are given below:

~~~ 
cl.data[is.na(f.tissue), ]
~~~ 
 
~~~
##                                              Identifyer
## 189      breast_BT-474_breast_Trastuzumab_pertuzumab\xa0
## 255         breast_CAL-51_breast_ABT737_\xa02'-Hydroxy-4
## 400 lung_A549_lung: NSCLC: NOS_\xdf-Lapachone_paclitaxel
##     Ward_Target_Distance Delta1   Delta2 Targeted Label Coverage
## 189                   NA     NA       NA     1.00     1        0
## 255                   NA     NA       NA     0.01     1        0
## 400                   NA     NA 1.270298     0.01     1        1
~~~
Since these rows do not contain sensitivity info anyway, they are left out of the analysis.
~~~
cl.data$Tissue <- f.tissue
cl.data$cl <- f.cl
cl.data$drug1 <- tolower(f.d1)
cl.data$drug2 <- tolower(f.d2)
cl.data <- cl.data[!is.na(f.tissue), ]
~~~
A total of 128 different cell lines corresponding to 16 tissues are included. The number of observations per tissue is:
~~~
table(cl.data$Tissue)
~~~
 
~~~
## 
##             bladder               blood                bone 
##                  45                  19                   1 
##              breast                 CNS            GI tract 
##              129489                  48                  30 
##              kidney                lung               other 
##                  20                  34                  16 
##               ovary            pancreas            prostate 
##                  29                   4                   1 
##                skin         soft tissue             thyroid 
##                   5                   2                   3 
## upper aerodigestive 
##                   6
~~~
