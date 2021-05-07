#read in data
fn = "/Users/robbi/Documents/DataAnalysis/Inca_PCA/Inca_PCA1/IncaMod+HO_NoAdm+TIW.pca.evec.txt"
evecDat2 = read.table(fn, col.names=c("Sample", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", 
                                     "PC7", "PC8", "PC9", "PC10", "Pop"))

#read in population groups to table for labelling purposes
popGroups2=read.table("/Users/robbi/Documents/DataAnalysis/Inca_PCA/Inca_PCA1/IncaMod+HO_NoAdm+TIW.ind", 
                     col.names=c("Sample", "Sex", "PopGroup.a"))

#file for population grouping in legend
populations2 = read.table("/Users/robbi/Documents/DataAnalysis/Inca_PCA/Inca_PCA1/IncaMod+HO_NoAdm+TIW.popGroups.txt", 
                         col.names=c("PopGroup.b", "Population") )
mergedPopDat2 = cbind(popGroups2, populations2)

#add population column 
mergedEvecDat2 = merge(mergedPopDat2, evecDat2, by="Sample")

#plot
plot(mergedEvecDat2$PC1, mergedEvecDat2$PC2, col=mergedEvecDat2$Population, xlab="PC1", ylab="PC2", pch=1)
legend("left", legend=levels(mergedEvecDat2$Population), col=1:length(levels(mergedEvecDat2$Population)), pch=20)


#to zoom to a part of the plot add this into the plot line
#xlim=c(-0.01,0.01), ylim=c(-0.03,-0.01)