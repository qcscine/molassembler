# TODO modularize for any number of csv files
prefixes <- c("pre-", "post-")

symmetryNames <- c(
  "linear", # 2
  "bent",
  "T-shaped", # 3
  "trigonal-planar", 
  "trigonal-pyramidal", # 4
  "tetrahedral",
  "seesaw",
  "square-planar",
  "trigonal-bipyramidal", # 5
  "pentagonal-planar",
  "square-pyramidal",
  "trigonal-prismatic", # 6
  "pentagonal-pyramidal",
  "octahedral",
  "pentagonal-bipyramidal", # 7
  "square-antiprismatic" # 8
)

nFiles <- length(symmetryNames)

pdf(
  file="matrices.pdf",
  width=4*(length(prefixes) + 1),
  height=4*nFiles
)

par(
  mfrow=c(nFiles, length(prefixes) + 1), 
  mar=c(2, 2, 2, 2)
)

readMatrixData <- function(csvFilename) {
  rawMatrixData <- read.csv(
    file=csvFilename, 
    header=FALSE
  )

  matrixData <- matrix(as.numeric(unlist(rawMatrixData)),nrow=nrow(rawMatrixData))
  matrixData <- t(matrixData[ncol(matrixData):1,])

  return(matrixData)
}

makeImagePlot <- function(matrixData, titleString) {
  image(
    matrixData,
    axes=FALSE,
    xlab="",
    ylab="",
    main=titleString,
    col=rev(heat.colors(100))
  )

  box()
}

for(symmetryName in symmetryNames) {
  for(prefix in prefixes) {
    makeImagePlot(
      readMatrixData(paste(prefix, symmetryName, ".csv", sep="")),
      paste(prefix, symmetryName, sep="")
    )
  }

  makeImagePlot(
    readMatrixData(paste(prefixes[1], symmetryName, ".csv", sep=""))
    - readMatrixData(paste(prefixes[2], symmetryName, ".csv", sep="")),
    paste(symmetryName, "-diff", sep="")
  )
}

dev.off()
