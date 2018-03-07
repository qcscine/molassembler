args <- commandArgs(TRUE)
filename <- toString(args[1])
splitFilename <- strsplit(filename, "\\.")
baseName <- splitFilename[[1]][length(splitFilename[[1]]) - 1]
pdfFilename <- paste(baseName, ".pdf", sep="")

values <- read.csv(file=filename, header=FALSE)
meta <- read.csv(file=paste(baseName, "-meta.csv", sep=""), header=FALSE)

edge_lengths <- as.numeric(meta[1,])
guessC <- as.numeric(meta[2,])[1]
determinedC <- as.numeric(meta[2,])[2]

x_seq <- values$V1
y_vals <- values$V2

pdf(pdfFilename, width=14, height=7)
par(
  mar=c(4, 4.5, 4, 1),
  mfrow=c(1, 2)
)

plot(
  x_seq,
  y_vals,
  type="n",
  xlab="c",
  ylab=expression(
    paste(
      "central angle sum deviation from ",
      2 * pi
    )
  )
)

edgesString <- paste(
  "{",
  paste(round(edge_lengths, 2), collapse=", "), 
  "}",
  sep=""
)

roundedCircumradius <- round(determinedC, 2)
title(edgesString)

abline(h=0)

lines(x_seq, y_vals)

abline(v=guessC, col="red")
abline(v=determinedC, col="blue")

# try constructing the circle
circumradius <- determinedC
plot(
  1,
  asp=1,
  type="n",
  xlab="x",
  ylab="y",
  xlim=c(-circumradius, circumradius), 
  ylim=c(-circumradius, circumradius)
)

roundedCircumradius <- round(circumradius, 2)
title(
  bquote(r == .(roundedCircumradius))
)

circlePoints <- 1000
circleX <- c(circumradius)
circleY <- c(0)
for(i in 1:circlePoints) {
  angle <- i * 2 * pi / circlePoints
  circleX <- c(circleX, cos(angle) * circumradius)
  circleY <- c(circleY, sin(angle) * circumradius)
}

lines(circleX, circleY)

currentAngle <- 0
currentX <- circumradius
currentY <- 0
for(i in 1:length(edge_lengths)) {
  centerpointAngle <- acos(1 - edge_lengths[i]^2 / (2 * circumradius^2))

  # compute next x,y coordinates
  nextX <- cos(currentAngle + centerpointAngle) * circumradius
  nextY <- sin(currentAngle + centerpointAngle) * circumradius

  # draw line segment
  segments(currentX, currentY, nextX, nextY)

  # overwrite for next iteration
  currentX <- nextX
  currentY <- nextY
  currentAngle <- currentAngle + centerpointAngle
}

dev.off()
