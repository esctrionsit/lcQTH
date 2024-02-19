data <- read.table("data.txt", header=F, stringsAsFactors = F)
pos <- data[1,]
pos <- as.character(pos[-c(1, length(pos))])
data <- data[-c(1),]
SMlst <- data[,1]
data <- data[,-c(1)]
ncol(data)
pheno <- data[,ncol(data)]
data <- data[,-c(ncol(data))]
data <- as.matrix(data)
rownames(data) <- SMlst
colnames(data) <- pos
library(ComplexHeatmap)
cm = c("A" = "#c55a11", "B" = "#2e75b6", "C" = "#dddddd", "D" = "#ffffff")

colanno <- read.table("colAnno.txt", header=F, stringsAsFactors = F)

pdf("plot.pdf", width = 10, height = nrow(data)*0.15+3)
Heatmap(data, col = cm, cluster_rows = F, cluster_columns = F, right_annotation = rowAnnotation(Pheno = anno_lines(as.numeric(pheno))), top_annotation = HeatmapAnnotation(Highlight = anno_barplot(as.numeric(colanno[,1]))))
dev.off()

svg("plot.svg", width = 10, height = nrow(data)*0.15+3)
Heatmap(data, col = cm, cluster_rows = F, cluster_columns = F, right_annotation = rowAnnotation(Pheno = anno_lines(as.numeric(pheno))), top_annotation = HeatmapAnnotation(Highlight = anno_barplot(as.numeric(colanno[,1]))))
dev.off()