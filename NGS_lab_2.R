# ############################################################################################################
# NGS lab 2: Sekvensanalyse og rapport
# ############################################################################################################

# 1 remove previous files in R
rm(list=ls())

# 2
# set wd C:\Users\admin\Documents\NGS_1
setwd("C:/Users/admin/Documents/NGS_1")

library(reshape2)

# 3 load csv table
ct<-read.table(file="CountTable.csv", header=T, sep=",")

# 4 string output table
str(ct)

# 5  
rownames(ct) <- ct$X
ct$X <- NULL

str(ct)

# 6 nr of sequence pairs
rowSums(ct)

# 7 nr of sequence pairs per sample
colSums(ct)

# 8 transform the table
t(ct)

ct_l <- melt(t(ct))

colnames(ct_l) <- c("SampleID", "HPVtype", "count")

# Hvilke prøver og har flere enn 100 HPV sekvenspar og for hvor mange ulike HPV typer?
ct_l_cut10 <- ct_l[ which(ct_l$count>99), ]


# 9 plot 
library(ggplot2)

c <- ggplot(ct_l_cut10, aes(factor(HPVtype)))

# c + geom_bar(stat="identity", position=position_dodge(), colour="seashell")
c + geom_bar()
c + geom_bar(position=position_dodge(), fill="green") + ggtitle("HPV dist")




