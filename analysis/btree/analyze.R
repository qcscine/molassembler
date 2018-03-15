filedata <- read.csv("btree.csv", header=FALSE, sep=",")

pdf(file="btree.pdf")

matrixData <- matrix(as.numeric(unlist(filedata)),nrow=nrow(filedata))
#matrixData <- t(matrixData[ncol(matrixData):1,])

imageplot <- image(
  matrixData,
  axes=FALSE,
  ylab="Requested size",
  xlab="Minimum degree",
  main="Space efficiency of constexpr BTree",
  col=grey.colors(ncol(matrixData), start = 0, end = 1, alpha = NULL),
  useRaster=TRUE
)

axis(1, at=c(0:18)/18, labels=c(2:20))
axis(2, at=c(0:4)/4, labels=c(10, 100, 1e3, 1e4, 1e5))

box()

dev.off()
