library(qtl)
sug <- read.cross("csv", "GM/", "GM.full.csv",
                  genotypes=c("GG", "GT", "TT"), alleles=c("G", "T"))

png("GM/plot.png", width = 500, height = 300, units = "px")
plot.map(sug)
dev.off()
