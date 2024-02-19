args <- commandArgs(trailingOnly = TRUE)
WP <- args[1]

library(qtl)
# sug <- read.cross("csv", "../11.QTL定位/", "data.merged.csv",
#                   genotypes=c("GG", "GT", "TT"), alleles=c("G", "T"))
sug <- read.cross("csv", WP, "data.merged.csv",
                  genotypes=c("GG", "GT", "TT"), alleles=c("G", "T"))
# out.SL <- scanone(sug, pheno.col=1, method="hk")
# out.SL.em <- scanone(sug, pheno.col=1, method="em")
# plot(out.SL.em, col="blue")
# plot(out.SL, col="red", add=TRUE)

# png("LD.heatmap.png", width = 1200, height = 1200, units = "px")
# sug2 <- est.rf(sug)
# plotRF(sug2, alternate.chrid=TRUE)
# nm <- est.map(sug, error.prob=0.001)
# plot(nm)
# dev.off()

# for (i in 1:21){
#     print(i)
#     nm <- est.map(sug, chr=as.character(i), error.prob=0.05, map.function = "kosambi", verbose=FALSE)
#     sug <- replace.map(sug, nm)
# }
nm <- est.map(sug, chr=args[2], error.prob=0.001, map.function = "kosambi", verbose=FALSE)
sug <- replace.map(sug, nm)
# png("Genetic.map.png", width = 1200, height = 1200, units = "px")
# plot.map(sug)
# dev.off()
write.cross(sug, format="csv", filestem=paste(WP, "/GM.", args[2], sep=""))