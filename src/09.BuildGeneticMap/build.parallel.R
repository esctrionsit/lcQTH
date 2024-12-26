args <- commandArgs(trailingOnly = TRUE)
WP <- args[1]

library(qtl)
sug <- read.cross("csv", WP, "data.merged.csv",
                  genotypes=c("GG", "GT", "TT"), alleles=c("G", "T"))

# png("LD.heatmap.png", width = 1200, height = 1200, units = "px")
# sug2 <- est.rf(sug)
# plotRF(sug2, alternate.chrid=TRUE)
# nm <- est.map(sug, error.prob=0.001)
# plot(nm)
# dev.off()

nm <- est.map(sug, chr=args[2], error.prob=0.001, map.function = "kosambi", verbose=FALSE)
sug <- replace.map(sug, nm)

# png("Genetic.map.png", width = 1200, height = 1200, units = "px")
# plot.map(sug)
# dev.off()

write.cross(sug, format="csv", filestem=paste(WP, "/GM.", args[2], sep=""))