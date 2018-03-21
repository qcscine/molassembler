#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
filename <- toString(args[1])

filedata <- read.csv(
  filename,
  header=FALSE,
  sep=",",
  colClasses=c(rep("numeric", 2))
)

# format is
# V1 - distanceError
# V2 - chiralError
# V3 - fourthDimError
# V4 - gradient norm
# V5 - compress
# V6 - fraction correct

distanceError <- filedata$V1
chiralError <- filedata$V2
fourthDimError <- filedata$V3
gradientNorm <- filedata$V4
compress <- filedata$V5
proportionCorrect <- filedata$V6

error <- distanceError + chiralError + fourthDimError


xSeq <- seq(1, length(error))

if(substr(filename, 1, 2) == "./") {
  substr(filename, 1,2) <- ""
  filename <- substr(filename, 3, nchar(filename))
  print(filename)
}

baseName <- strsplit(filename, "\\.")
pdfFilename <- paste(baseName[[1]][length(baseName[[1]]) - 1], ".pdf", sep="")

pdf(
  pdfFilename,
  width=14,
  height=8
)
par(
  mar=c(0, 4, 0, 1),
  oma=c(4, 0, 1, 0)
)
layout(matrix(c(1,1,2,2,3,3,4), nrow=7, ncol=1))

# error plot
plot(
  xSeq,
  error,
  type="n",
  xlab="",
  ylim=c(
    min(
      distanceError[which(distanceError > 0)],
      chiralError[which(chiralError > 0)],
      fourthDimError[which(fourthDimError > 0)]
    ),
    max(
      distanceError[which(distanceError > 0)],
      chiralError[which(chiralError > 0)],
      fourthDimError[which(fourthDimError > 0)]
    )
  ),
  ylab="Error (log scale)",
  log="y",
  xaxt="n"
)

lines(
  xSeq,
  error,
  lwd=2
)

lines(
  xSeq,
  distanceError,
  col="steelblue",
  lwd=2
)

lines(
  xSeq,
  chiralError,
  col="tomato",
  lwd=2
)

lines(
  xSeq,
  fourthDimError,
  col="olivedrab",
  lwd=2
)

# proportion of error development
plot(
  xSeq,
  rep(1, length(error)),
  ylim=c(0, 1),
  type="n",
  xlab="",
  ylab="Proportion of error contrib.",
  xaxt="n"
)

lines(
  xSeq,
  distanceError / error,
  col="steelblue",
  lwd=2
)

lines(
  xSeq,
  chiralError / error,
  col="tomato",
  lwd=2
)

lines(
  xSeq,
  fourthDimError / error,
  col="olivedrab",
  lwd=2
)

# gradient norm plot
plot(
  xSeq,
  gradientNorm,
  type="n",
  xlab="",
  ylab="Gradient norm (log scale)",
  log="y",
  xaxt="n"
)

lines(
  xSeq,
  gradientNorm,
  lwd=2
)

# proportion correct chirality constraints
plot(
  xSeq,
  proportionCorrect,
  type="n",
  xlab="Step number",
  ylab="Frac. correct chir."
)

lines(
  xSeq,
  proportionCorrect,
  lwd=2
)

dev.off()
