#read in admixture *Q output
tbl <- read.table("<path>/<name>.5.Q")

labTable <- read.table("<path>/<name>labels.txt", col.names=c("Num", "XLabel"))

#this binds (column-wise) your admixture proportions and your population names taken from the .ind file
mergedAdmixtureTable = cbind(tbl, labTable)

#this orders your table by your populations
ordered <- mergedAdmixtureTable[order(mergedAdmixtureTable$XLabel),]

#plot
barplot(t(as.matrix(subset(mergedAdmixtureTable, select=V1:V5))), 
        col=rainbow(5), space =0.05, border=NA, ylab="Ancestry", xlab="K=5", 
        names.arg=mergedAdmixtureTable$XLabel, las = 3)
