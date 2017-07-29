# ############################################################################################################
# NGS lab 1:
# ############################################################################################################

# from bioconductor.R
library(ShortRead)

# 1
# Opprydding før du starter, fjerne tidligere filer i R (om det skulle være noen)

rm(list=ls())

# 2
# set wd C:\Users\admin\Documents\NGS_1
setwd("C:/Users/admin/Documents/NGS_1")

# 3
# load the fastq file
fq <- readFastq("MS751-20161102-method-validation-MGP-1H9-SA709-SA508_S72_L001_R1_001.fastq.gz")

# 4 check that the fastq file is loaded
fq

# 5 read the header of fq
head(sread(fq))

# 6 look at the quality of fq header
head(quality(fq))

# 7 summarise the fq sequences in a table
tbls <- tables(fq)

# 8 how many different sequences are there in fq?
head(tbls$distribution)

# 9 plot the tbls distribution
plot(tbls$distribution)

# 10 find the most frequent frequencies in tbls
tbls$top[1:5]

# 11 BLAST of the 1st sequence gives:
# https://blast.ncbi.nlm.nih.gov/Blast.cgi
# A: the sequence comes from Humans and C. carpio
# with a >90% sequence homology


# 12 load a QA summary
Summary <- qa("C:/Users/admin/Documents/NGS_1", "*.fastq.gz", type = "fastq")

# 13 how many sequences are ther in each loaded summary file?
head(Summary[["readCounts"]])




