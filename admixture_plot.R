tbl <- read.table("/Users/robbi/Documents/DataAnalysis/Admixture_small/for_smallAdmixture_4.5.Q")
#indTable <- read.table("/Users/robbi/Documents/DataAnalysis/Admixture_small/YorSard_poplabels.txt",
  #                    col.names = c("Sample", "Pop"))
labTable <- read.table("/Users/robbi/Documents/DataAnalysis/Admixture_small/smallAdmixture_4labels.txt",
                       col.names=c("Num", "XLabel"))
#you don't need below command, this is if you ran a big data set and wanted to extract specific groups
#popGroups <- read.table("", col.names=c("Pop", "PopGroup"))

#this binds (column-wise) your admixture proportions and your population names taken from the .ind file
mergedAdmixtureTable = cbind(tbl, labTable)

#don't need below command
#mergedAdmWithPopGroups = merge(mergedAdmixtureTable, popGroups, by="Pop")

#this orders your table by your populations
#ordered <- mergedAdmixtureTable[order(mergedAdmixtureTable$Pop),]

#time to hatch the plot mwahahaha
barplot(t(as.matrix(subset(mergedAdmixtureTable, select=V1:V5))), 
        col=rainbow(11), space =0.05, border=NA, ylab="Ancestry", xlab="K=5", names.arg=mergedAdmixtureTable$XLabel, 
        las = 3)



