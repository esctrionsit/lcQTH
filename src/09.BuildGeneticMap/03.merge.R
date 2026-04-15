args <- commandArgs(trailingOnly = TRUE)

chrcount <- as.numeric(args[1])

library(qtl)

fullmap <- read.cross("csv", "GM/", "data.merged.csv", genotypes=c("GG", "GT", "TT"), alleles=c("G", "T"))
chrs <- chrnames(fullmap)

for (chrid in chrs){
	ssug <- read.cross("csv", "GM/", paste("GM.", chrid, ".csv", sep=""), genotypes=c("GG", "GT", "TT"), alleles=c("G", "T"))
	ssug2 <- subset(ssug, chr=chrid)
	write.cross(ssug2, format="csv", filestem=paste("GM/GM.s.", chrid, sep=""))
}


data <- read.csv(paste("GM/GM.s.", chrs[1],".csv", sep=""), header=F, stringsAsFactors = F)
if (length(chrs) > 1){
	for (chridx in 2:length(chrs)){
		data2 <- read.csv(paste("GM/GM.s.", chrs[chridx], ".csv", sep=""),header=F, stringsAsFactors = F)
		data <- cbind(data, data2[-c(1)])
	}
}
write.table(x=data, file="GM/GM.full.csv", append = FALSE, quote = F, sep = ",",
	eol = "\n", na = "NA", dec = ".", row.names = F,
	col.names = F)
