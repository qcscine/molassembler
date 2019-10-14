shapeNames <- c(
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

coordinationNumber <- c(2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 8)

# file format is now:
# V1 - symmetry index (zero-based)
# V2 - error value

# read data
filedata <- read.csv(
  file="DGRefinementProblem-symmetric-ensemble-errors.csv",
  header=FALSE,
  sep=","
)

symmetryCol <- as.numeric(filedata$V1)
errorCol <- as.numeric(filedata$V2)

symmetryErrors <- list()
filedata <- list()
for(i in seq(1:length(shapeNames) - 1)) {
  # add into list
  symmetryErrors[[i]] <- errorCol[
    which(
      symmetryCol + 1 == i
    )
  ]
}

# 1 Histograms 4x4 plot
pdf("histograms.pdf", width=14, height=14)
par(mfrow=c(4,4))

for(i in seq(1:length(shapeNames) - 1)) {
  hist(
    symmetryErrors[[i]],
    main = shapeNames[i],
    xlab="errf"
    #breaks=20,
  )
}

dev.off()

# 2 Coordination Number vs average errf plot
averages <- c()
for(i in seq(1:length(shapeNames) - 1)) {
  averages <- c(averages, mean(symmetryErrors[[i]]))
}

pdf("CNvsErrf.pdf", width=7, height=7)
par(mfrow=c(1, 1))

bp <- barplot(
  averages,
  xlab="Average error function value",
  ylab="Coordination number",
  horiz=TRUE,
  col="black"
)

axis(2, at=bp, labels=coordinationNumber)

for(i in 1:length(shapeNames)) {
  if(averages[i] < 0.5 * max(averages)) { # this boundary must be varied if values change
    text(
      averages[i],
      bp[i],
      labels=shapeNames[i],
      pos=4
    )
  } else {
    text(
      averages[i],
      bp[i],
      labels=shapeNames[i],
      pos=2,
      col="white"
    )
  }
}

dev.off()
