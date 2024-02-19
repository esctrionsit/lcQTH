
library(qtl)
sug <- read.cross("csv", "./", "data.merged.csv",
                  genotypes=c("GG", "GT", "TT"), alleles=c("G", "T"))
# out.hk <- scanone(sug, pheno.col=1, method="hk")
out.em <- scanone(sug, pheno.col=1, method="em")

png("plot.merged.png", width = 900, height = 300, units = "px")
plot(out.em, chr="11", col="black")
# plot(out.hk, col="red", add=TRUE)
# abline(v=267.6114686,col="red")
dev.off()

write.csv(out.em,file='plot.merged.csv')