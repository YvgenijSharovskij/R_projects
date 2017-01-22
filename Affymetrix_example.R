###############################################################################################################
# Affymetrix chip data processing 
# with source("https://bioconductor.org/biocLite.R")
#
# based on the common workflow:
# 1. read the probe level CEL-data files
# 2. background correction
# 3. normalisation
# 4. probe specific background correction
# 5. ouput summary into an expression matrix
#
# summary:
# 1. reading the probe level CEL-data files
# 2. normalisation of the probesets
# 3. individual gene analysis of differentially expressed genes
#    as change of fold
#
###############################################################################################################

# load the Affymetrix library with Biobase
library("affy") 
require("Biobase")

###############################################################################################################
# 1. reading the probe level CEL-data files:
#
# read the Affymetrix CEL files, which contain the raw probe data with 
# hybridisation intensities and locations, into affy.data
#
###############################################################################################################
#
# alternatively, use GUI:
# affy.data <- ReadAffy(widget=TRUE) ## select the CEL files with the tkWidgets GUI widget


# set the working dir
setwd("C:/Users/admin/Documents/R/examples/affy_su")

# list the CEL-files in the set dir
files <- list.files(pattern = "\\.CEL$")
files

# read all of the listed CEL-files
affy.data <- ReadAffy(filenames = files)


###############################################################################################################
# 2. normalisation of the probesets
#
# note: 
# alternatively, can use:
# eset.rma = justRMA() ## includes log2-transformation
#
# or, the Biobase expresso class:
# eset <- expresso(affy.data, normalize.method="qspline",
#         bgcorrect.method="rma",pmcorrect.method="pmonly",
#         summary.method="liwong")
###############################################################################################################

# normalise with MAS5 and write into a  tab-delimited file
eset.mas5 = mas5(affy.data)
write.exprs(eset.mas5, file="affy_mas5_normalised_data.txt")

# set expression matrix:
# [ row     :    gene
#   column  :    sample ]
exprSet.nologs = exprs(eset.mas5)
# list the column names
colnames(exprSet.nologs)
# set new column names
colnames(exprSet.nologs) = c("SAMPLE_A", "SAMPLE_B", 
                             "SAMPLE_C", "SAMPLE_D",
                             "SAMPLE_E", "SAMPLE_F", 
                             "SAMPLE_G", "SAMPLE_H")
# list the new column names
colnames(exprSet.nologs)

# log2-transformed expression values (+/- 1)
exprSet = log(exprSet.nologs, 2)

# print to file
write.table(exprSet, file="affy_mas5_matrix.txt", quote=F, sep="\t")


###############################################################################################################
# 3. individual gene analysis of differentially expressed genes
# as change of fold
#
###############################################################################################################

# get the means of each set of chips, according to columns 1-8
sampleAB.mean = apply(exprSet[, c(1,2)], 1, mean)
sampleCD.mean = apply(exprSet[, c(3,4)], 1, mean)
sampleEF.mean = apply(exprSet[, c(5,6)], 1, mean)
sampleGH.mean = apply(exprSet[, c(7,8)], 1, mean)

# get the mean ratios;
# since the data is log2-transformed, use the general formula:
# log(A / B) = log(A) - log(B)
sampleCD.sampleAB.ratio = sampleCD.mean - sampleAB.mean
sampleGH.sampleEF.ratio = sampleGH.mean - sampleEF.mean

# column bind for all data
all.data = cbind(exprSet, sampleAB.mean, sampleCD.mean, sampleEF.mean, sampleGH.mean,
                 sampleCD.sampleAB.ratio, sampleGH.sampleEF.ratio)
# list column names
colnames(all.data)

# print to a single file
write.table(all.data, file="Affy_Microarray_Analaysis_Data_JAN-2017.txt", quote=F, sep="\t")


# scatter plot
colnames(all.data)
x.data = all.data[, "sampleAB.mean"]
y.data = all.data[, "sampleCD.mean"]
plot(x.data, y.data)
plot(x.data, y.data, main = "Log2 expression sampleCD vs. sampleAB",
     xlab="sampleAB", ylab="sampleCD", col="purple", cex=0.5)

