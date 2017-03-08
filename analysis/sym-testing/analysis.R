filenames <- c(
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

# read data
filedata <- list()
upperlimit <- 0
for(i in seq(1:length(filenames) - 1)) {
  # add into list
  filedata[[i]] <- as.numeric(read.csv(
    file=paste(filenames[i], ".csv", sep=""),
    header=FALSE,
    sep=","
  )$V1)

  # update maximum
  upperlimit <- max(upperlimit, max(filedata[[i]]))
}

# 1 Histograms 4x4 plot
pdf("histograms.pdf", width=14, height=14)
par(mfrow=c(4,4))

for(i in seq(1:length(filenames) - 1)) {
  hist(
    filedata[[i]],
    main = filenames[i],
    xlab="errf"
    #breaks=20,
  )
}

dev.off()

# 2 Coordination Number vs average errf plot
averages <- c()
for(i in seq(1:length(filenames) - 1)) {
  averages <- c(averages, mean(filedata[[i]]))
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

for(i in 1:length(filenames)) {
  if(averages[i] < 0.04) { # this boundary must be varied if values change
    text(
      averages[i],
      bp[i],
      labels=filenames[i],
      pos=4
    )
  } else {
    text(
      averages[i],
      bp[i],
      labels=filenames[i],
      pos=2,
      col="white"
    )
  }
}

dev.off()
