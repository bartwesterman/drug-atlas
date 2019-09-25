# drug-atlas
Synergy prediction and the drug atlas R script

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
## 400 Lung_A549_lung: NSCLC: NOS_\xdf-Lapachone_paclitaxel
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
## Logistic regression
The data contains an indicator variable Label that is 1 for cases where the two drugs display synergy, and 0 otherwise.

Per tissue, we will fit a logistic regression using all rows for which Label is 1, and then choosing at random the same number of cell lines (observations) for that tissue for which Label is 0.

I want to compare model fits. <Model IA> relates only the target information (at least one drug targets the cell line modification) is:

logit(θij)=α+βtijk+eijk.

Model IB relates both target information and the drug atlas distance:

logit(θij)=α+βtijk+δdij+eijk.
This model can be used to assess the added value of the distance by comparing its results to those using model IA.

The second pair of models involves the target information, as well as the individual sensitivities of the cell line to the drugs. Model IIA is:

logit(θij)=α+βtijk+γikSik+γjkSjk+eijk.

Finally, model IIB is the same as model IIA, but also includes the drug atlas distance:

logit(θij)=α+βtijk+δdij+γikSik+γjkSjk+eijk.

This model can be used to assess the added value of the distance by comparing its results to those using model IIA.

~~~
# We define selection variables to enable running the regressions
# First for breast cancer as it is the largest 
myt <- "panCancer"
#tsel <- cl.data$Tissue == myt
mydata <- cl.data
~~~~
Note that a large number (63609) of cell lines have no drug atlas distance available. In addition, 82834 cell lines have no information for sensitivity to drug 1, and 77768 have no information for sensitivity to drug 2. The number of cell lines with information available on all three variables (drug atlas distance, sensitivity to drug 1 and to drug 2) is 11519. For convenience, and to guarantee that results are not biased (which can happen if one group, say with synergy information, is larger than the other group), we leave those out of the analyses.
~~~
sel.nonas.dist <- !( is.na(mydata$Ward_Target_Distance) )
sel.nonas.sens <- !( is.na(mydata$Delta1)  | is.na(mydata$Delta2) )
sel.nonas.sens.dist <- !( is.na(mydata$Ward_Target_Distance) | is.na(mydata$Delta1)  | is.na(mydata$Delta2) )
mydata1a <- mydata
mydata1b <- mydata[sel.nonas.dist, ]
#mydata2 <- mydata[sel.nonas.sens, ]
mydata2a <- mydata[sel.nonas.sens, ]
mydata2b <- mydata[sel.nonas.sens.dist, ]
# Now make selection of cell lines to be included
set.seed <- 45621
sel1a <- sel.samples(mydata1a)
sel1b <- sel.samples(mydata1b)
sel2a <- sel.samples(mydata2a)
sel2b <- sel.samples(mydata2b)
# Models IA, IB
# Use the same basic data for both in order to be able to make the graphs
dat1a <- mydata1a[sel1a, ]
m1a <- glm(Label ~ Targeted, data = dat1a, na.action = na.exclude)
dat1b <- dat1a # dat1b <- mydata1b[sel1b, ]
m1b <- glm(Label ~ Targeted + Ward_Target_Distance, data = dat1a, 
           na.action = na.exclude)
# Models IIA, IIB
# Use the same basic data for both in order to be able to make the graphs
dat2a <- mydata2a[sel2a, ]
m2a <- glm(Label ~ Targeted + Delta1 + Delta2, data = dat2a, na.action = na.exclude)
dat2b <- dat2a # dat2b <- mydata2a[sel2a, ]
m2b <- glm(Label ~ Targeted + Delta1 + Delta2 + Ward_Target_Distance, 
           data = dat2a, na.action = na.exclude)
# Predictions
dat1a$f1a <- predict(m1a, type = "response")
dat1a$f1b <- predict(m1b, type = "response")
dat2a$f2a <- predict(m2a, type = "response")
dat2a$f2b <- predict(m2b, type = "response")
# Random response
dat1ar <- dat1a
dat1ar$Label <- sample(dat1ar$Label)
m1ar <- glm(Label ~ Targeted, data = dat1ar, na.action = na.exclude)
dat1br <- dat1ar # dat1b
#dat1br$Label <- sample(dat1br$Label)
m1br <- glm(Label ~ Targeted + Ward_Target_Distance, data = dat1br, na.action = na.exclude)
dat2ar <- dat2a
dat2ar$Label <- sample(dat2ar$Label)
m2ar <- glm(Label ~ Targeted + Delta1 + Delta2, data = dat2ar, na.action = na.exclude)
dat2br <- dat2ar # dat2b
#dat2br$Label <- sample(dat2br$Label)
m2br <- glm(Label ~ Targeted + Delta1 + Delta2 + Ward_Target_Distance, 
           data = dat2br, na.action = na.exclude)
dat1ar$f1a <- predict(m1ar, type = "response")
dat1ar$f1b <- predict(m1br, type = "response")
dat2ar$f2a <- predict(m2ar, type = "response")
dat2ar$f2b <- predict(m2br, type = "response")
~~~
For completeness, the model fit coefficients tables are given below:
~~~
round(summary(m1a)$coefficients, 3)
##             Estimate Std. Error t value Pr(>|t|)
## (Intercept)    0.512      0.017  29.739    0.000
## Targeted      -0.108      0.054  -2.019    0.044
~~~
round(summary(m1b)$coefficients, 3)
~~~
##                      Estimate Std. Error t value Pr(>|t|)
## (Intercept)             0.159      0.039   4.095    0.000
## Targeted                0.227      0.076   2.972    0.003
## Ward_Target_Distance    0.006      0.001   3.953    0.000
~~~
round(summary(m2a)$coefficients, 3)
~~~
##             Estimate Std. Error t value Pr(>|t|)
## (Intercept)    0.460      0.045  10.338    0.000
## Targeted      -0.248      0.150  -1.657    0.100
## Delta1        -0.061      0.015  -3.999    0.000
## Delta2        -0.015      0.015  -0.959    0.339
~~~
round(summary(m2b)$coefficients, 3)
~~~~
##                      Estimate Std. Error t value Pr(>|t|)
## (Intercept)             0.172      0.088   1.959    0.055
## Targeted                0.055      0.202   0.273    0.786
## Delta1                 -0.074      0.023  -3.260    0.002
## Delta2                  0.083      0.024   3.381    0.001
## Ward_Target_Distance    0.006      0.004   1.538    0.129
~~~
~~~
# Plotting results
dat1s <- dat1a[order(dat1a$Label, dat1a$f1a), ]
lcol <- rep("blue", nrow(dat1a))
lcol[dat1s$Label == 0] <- "red"
#par(mfrow = c(1, 2))
plot(dat1s$Label, pch = 21, col = lcol, cex=2, 
     main=paste(myt, ", sorted by pred model IA"), xlab = "cases", 
     ylab = "observed and predicted synergy")
points(dat1s$f1a, pch = 2, col = lcol)
points(dat1s$f1b, pch = 3, col = lcol)
legend("topleft", legend = c("observed", "model IA", "model IB"), pch = c(21, 2, 3) )
legend("bottomright", legend = c("synergy obs", "no synergy obs"), fill = c("blue", "red"))
segments(0,0.5, nrow(dat1a), 0.5, lty = "dashed", col = "grey", lwd = 2)
~~~
~~~
#lcol1 <- lcol

dat2s <- dat2a[order(dat2a$Label, dat2a$f2a), ]
dat2s <- dat2s[!(is.na(dat2s$f2a)), ]
lcol <- rep("blue", nrow(dat2a))
lcol[dat2s$Label == 0] <- "red"
plot(dat2s$Label, pch = 21, col = lcol, cex=2,
     main=paste(myt, ", sorted by pred model IIA"), xlab = "cases", 
     ylab = "observed and predicted synergy")

points(dat2s$f2a, pch = 2, col = lcol)
points(dat2s$f2b, pch = 3, col = lcol)
legend("topleft", legend = c("observed", "model IIA", "model IIB"), pch = c(21, 2, 3) )
segments(0,0.5, nrow(dat2a), 0.5, lty = "dashed", col = "grey", lwd = 2)
legend("bottomright", legend = c("synergy obs", "no synergy obs"), fill = c("blue", "red"))
~~~
~~~
#lcol2 <- lcol

dat2s <- dat2a[order(dat2a$Label, dat2a$f2b), ]
dat2s <- dat2s[!(is.na(dat2s$f2a)|is.na(dat2s$f2b)), ]
lcol <- rep("blue", nrow(dat2s))
lcol[dat2s$Label == 0] <- "red"
plot(dat2s$Label, pch = 21, col = lcol, cex=2,
     main=paste(myt, ", sorted by pred model IIB"), xlab = "cases", 
     ylab = "observed and predicted synergy")

points(dat2s$f2a, pch = 2, col = lcol)
points(dat2s$f2b, pch = 3, col = lcol)
legend("topleft", legend = c("observed", "model IIA", "model IIB"), pch = c(21, 2, 3) )
segments(0,0.5, nrow(dat2a), 0.5, lty = "dashed", col = "grey", lwd = 2)
legend("bottomright", legend = c("synergy obs", "no synergy obs"), fill = c("blue", "red"))
~~~~~
The graphs above display predicted probabilities for the same basic model per graph, without the drug distance (triangle pointing upwards) and with the drug distance (triangle pointing downwards). In both cases, we can see an improvement (decrease in false negatives as well as false positives) when drug distance is included, compared to the model without it. Specifically, for models IIA and IIB, the number of false predictions are:
~~~
t1a <- table(factor(dat1a$Label, labels = c("no synergy obs", "synergy obs")),
      factor(dat1a$f1a>0.5, labels = c("no synergy pred IA", "synergy pred IA")) )
t1a
##                 
##                  no synergy pred IA synergy pred IA
##   no synergy obs                 59             420
##   synergy obs                    40             439
round(prop.table(t1a, margin = 1), 2)
##                 
##                  no synergy pred IA synergy pred IA
##   no synergy obs               0.12            0.88
##   synergy obs                  0.08            0.92
t1b <- table(factor(dat1a$Label, labels = c("no synergy obs", "synergy obs")),
      factor(dat1a$f1b>0.5, labels = c("no synergy pred IB", "synergy pred IB")) )
t1b
##                 
##                  no synergy pred IB synergy pred IB
##   no synergy obs                272               5
##   synergy obs                   112               8
round(prop.table(t1b, margin = 1), 2)
##                 
##                  no synergy pred IB synergy pred IB
##   no synergy obs               0.98            0.02
##   synergy obs                  0.93            0.07
t2a <- table(factor(dat2a$Label, labels = c("no synergy obs", "synergy obs")),
      factor(dat2a$f2a>0.5, labels = c("no synergy pred IIA", "synergy pred IIA")) )
t2a
##                 
##                  no synergy pred IIA synergy pred IIA
##   no synergy obs                  48               22
##   synergy obs                     27               43
round(prop.table(t2a, margin = 1), 2)
##                 
##                  no synergy pred IIA synergy pred IIA
##   no synergy obs                0.69             0.31
##   synergy obs                   0.39             0.61
t2b <- table(factor(dat2a$Label, labels = c("no synergy obs", "synergy obs")),
      factor(dat2a$f2b>0.5, labels = c("no synergy pred IIB", "synergy pred IIB")) )
t2b
##                 
##                  no synergy pred IIB synergy pred IIB
##   no synergy obs                  34                4
##   synergy obs                     10               15
round(prop.table(t2b, margin = 1), 2)
##                 
##                  no synergy pred IIB synergy pred IIB
##   no synergy obs                0.89             0.11
##   synergy obs                   0.40             0.60
~~~
So, when we compare predictions with model IB with those from model IA, we see that false negatives and positive ratios decrease from 0.88 and 0.08 with IA to 0.02 and 0.93 with IB. Similarly, comparing predicted probabilities with model IIB with those from model IIA, false negatives and positive ratios decrease from 0.31 and 0.39 with IIA to 0.11 and 0.4 with IIB. Thus, predicted synergy improves when drug atlas distance is included in these models.

ROC curves breast cell lines
~~~
mgrid <- seq(from = 0.01, to = 0.99, by = 0.01)
roc.1a <- get.roc(status = dat1a$Label, predicted = dat1a$f1a, plot = FALSE, grid = mgrid)
roc.1b <- get.roc(status = dat1a$Label, predicted = dat1a$f1b, plot = FALSE, grid = mgrid)
roc.2a <- get.roc(status = dat2a$Label, predicted = dat2a$f2a, plot = FALSE, grid = mgrid)
roc.2b <- get.roc(status = dat2a$Label, predicted = dat2a$f2b, plot = FALSE, grid = mgrid)
# And for the random response
roc.1ar <- get.roc(status = dat1ar$Label, predicted = dat1ar$f1a, plot = FALSE, grid = mgrid)
roc.1br <- get.roc(status = dat1ar$Label, predicted = dat1ar$f1b, plot = FALSE, grid = mgrid)
roc.2ar <- get.roc(status = dat2ar$Label, predicted = dat2ar$f2a, plot = FALSE, grid = mgrid)
roc.2br <- get.roc(status = dat2ar$Label, predicted = dat2ar$f2b, plot = FALSE, grid = mgrid)
To compare the predictions from model IB to those obtained with model IA for various cut-offs on the predicted synergy, we make ROC curves for both cases (left-hand side graph below).

par(mfrow=c(1, 2))
mcols <- rep(c("blue", "purple"), 2)
plot(1-roc.1a[, "spec"], roc.1a[, "sens"], type = "l", xlim = c(0, 1), ylim = c(0, 1),
     col = mcols[1], main = "ROC curves models IA and IB", xlab = "1-specificity",
     ylab = "sensitivity", lwd = 2)
lines(1-roc.1b[, "spec"],  roc.1b[, "sens"],  col = mcols[2], lwd = 2)
lines(1-roc.1ar[, "spec"], roc.1ar[, "sens"], col = mcols[3], lwd = 2, lty = "dotted")
lines(1-roc.1br[, "spec"], roc.1br[, "sens"], col = mcols[4], lwd = 2, lty = "dotted")
segments(0, 0, 1, 1, lty = "dashed", col = "gray", lwd = 2)
legend("bottomright", 
       legend = c("IA - no drug distance", "IB - with drug distance", "IA - random", "IB - random"),
       lty = c(rep("solid", 2), rep("dotted", 2)), 
       col = mcols, lwd = 2)

plot(1-roc.2a[, "spec"], roc.2a[, "sens"], type = "l", xlim = c(0, 1), ylim = c(0, 1),
     col = mcols[1], main = "ROC curves models IIA and IIB", xlab = "1-specificity",
     ylab = "sensitivity", lwd = 2)
lines(1-roc.2b[, "spec"],  roc.2b[, "sens"],  col = mcols[2], lwd = 2)
lines(1-roc.2ar[, "spec"], roc.2ar[, "sens"], col = mcols[3], lwd = 2, lty = "dotted")
lines(1-roc.2br[, "spec"], roc.2br[, "sens"], col = mcols[4], lwd = 2, lty = "dotted")
segments(0, 0, 1, 1, lty = "dashed", col = "gray", lwd = 2)
legend("bottomright", 
       legend = c("IIA - no drug distance", "IIB - with drug distance", 
                  "IIA - random", "IIB - random"),
       lty = c(rep("solid", 2), rep("dotted", 2)), 
       col = mcols, lwd = 2)
~~~
