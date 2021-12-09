#read in data
fn = "<path>/<name>.pca.evec.txt"
evecDat = read.table(fn, col.names=c("Sample", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", 
			"PC7", "PC8", "PC9", "PC10", "Population"))

#read in population groups to table for labelling purposes
populations = read.table("<path>/<name>.popLabels.txt",col.names=c("Sample", "Label"))
ind = read.table("<path>/<name>.ind", 
                col.names=c("Sample", "Sex", "Population"))

#bind population colums (horizontally) for labelling first
mergedPopDat = cbind(ind, populations)
#merge with evec data
mergedEvecDat3 = merge(mergedPopDat3, evecDat3, by="Sample")

#plot
plot(mergedEvecDat$PC1, mergedEvecDat$PC2, col=mergedEvecDat$Label, 
			pch=as.integer(mergedEvecDat$Label) %% 24, xlab="PC1", ylab="PC2")
legend("topright", xpd=TRUE, legend=levels(mergedEvecDat$Label), 
			col=1:length(levels(mergedEvecDat$Label)), pch=1:length(levels(mergedEvecDat$Label)))
