args <- commandArgs(TRUE)
filename <- toString(args[1])
splitFilename <- strsplit(filename, "\\.")
baseName <- splitFilename[[1]][length(splitFilename[[1]]) - 1]
pdfFilename <- paste(baseName, ".pdf", sep="")

transformY <- TRUE

values <- read.csv(file=filename, header=FALSE)
meta <- read.csv(file=paste(baseName, "-meta.csv", sep=""), header=FALSE)

pdf(pdfFilename, width=14, height=7)
par(
  mar=c(4, 4.5, 4, 1),
  mfrow=c(1, 2)
)

edge_lengths <- as.numeric(meta[1,])
strategies <- as.numeric(meta[2,])[1:2]
rhoDistribution <- as.numeric(meta[3,])[1:2]
rhoRoot <- as.numeric(meta[4,])[1]
validRoot <- as.numeric(meta[4,])[2]

x_seq <- values$V1
y_vals <- values$V2

#if(rhoRoot != 0 && rhoRoot < 0.1 * max(x_seq)) {
#  which_indices <- which(x_seq <= 2 * rhoRoot)
#  x_seq <- x_seq[which_indices]
#  y_vals <- y_vals[which_indices]
#}

if(transformY) { # transform to logarithmic values
  y_vals <- log(1 + abs(y_vals))
  y_lims <- c(min(y_vals), max(y_vals))

  y_label <- expression(paste("log(1 + |", A[5]^2 - B[5]^2*Delta[5], "|)"))
} else {
  y_lims <- c(
    max(
      min(y_vals),
      -0.01
   ),
    min(
      max(y_vals),
      0.01
    )
  )
  y_label <- expression(A[5]^2 - B[5]^2*Delta[5])
}

plot(
  x_seq,
  y_vals,
  xlim=c(min(x_seq), max(x_seq)),
  ylim=y_lims,
  log="x",
  type="n",
  xlab=expression(rho),
  ylab=y_label
)

titleCol <- "red"
if(validRoot) {
  titleCol <- "olivedrab"
}

title(
  paste(
    paste(round(edge_lengths, 2), collapse=" - "), 
    " -> ", 
    rhoRoot, 
    sep=""
  ),
  col.main=titleCol
)

abline(h=0)

lines(x_seq, y_vals)

if(strategies[1] == 1) {
  abline(v = rhoDistribution[0], col="blue", lty=2)
  if(strategies[2] == 0) { # found with strategy 1
    abline(v = rhoRoot, col="blue")
  }
}

if(strategies[2] == 1) {
  abline(v = rhoRoot, col="red")
}

if(rhoRoot != 0) {
  # try constructing the circle
  circumradius <- sqrt(1 / rhoRoot)
  plot(
    1,
    asp=1,
    type="n",
    xlab="",
    ylab="",
    xlim=c(-circumradius, circumradius), 
    ylim=c(-circumradius, circumradius)
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
  for(i in 1:5) {
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

}

dev.off()
