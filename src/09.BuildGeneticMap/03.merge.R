args <- commandArgs(trailingOnly = TRUE)

chrcount <- as.numeric(args[1])

library(qtl)
for (chrid in 1:chrcount){
	ssug <- read.cross("csv", "GM/", paste("GM.", as.character(chrid), ".csv", sep=""), genotypes=c("GG", "GT", "TT"), alleles=c("G", "T"))
	ssug2 <- subset(ssug, chr=chrid)
	write.cross(ssug2, format="csv", filestem=paste("GM/GM.s.", as.character(chrid), sep=""))
}


data <- read.csv("GM/GM.s.1.csv", header=F, stringsAsFactors = F)
if (chrcount > 1){
	for (chrid in 2:chrcount){
		data2 <- read.csv(paste("GM/GM.s.", as.character(chrid), ".csv", sep=""),header=F, stringsAsFactors = F)
		data <- cbind(data, data2[-c(1)])
	}
}
write.table(x=data, file="GM/GM.full.csv", append = FALSE, quote = F, sep = ",",
	eol = "\n", na = "NA", dec = ".", row.names = F,
	col.names = F)