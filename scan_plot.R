args <- commandArgs(TRUE)
filename <- toString(args[1])
splitFilename <- strsplit(filename, "\\.")
baseName <- splitFilename[[1]][length(splitFilename[[1]]) - 1]
pdfFilename <- paste(baseName, ".pdf", sep="")

values <- read.csv(file=filename, header=FALSE)
meta <- read.csv(file=paste(baseName, "-meta.csv", sep=""), header=FALSE)

pdf(pdfFilename, width=14, height=7)
par(
  mar=c(4, 4.5, 4, 1)
)

edge_lengths <- as.numeric(meta[1,])
strategies <- as.numeric(meta[2,])[1:2]
rhoDistribution <- as.numeric(meta[3,])[1:2]
rhoRoot <- as.numeric(meta[4,])[1]

x_seq <- values$V1
y_vals <- values$V2

#if(rhoRoot != 0 && rhoRoot < 0.1 * max(x_seq)) {
#  which_indices <- which(x_seq <= 2 * rhoRoot)
#  x_seq <- x_seq[which_indices]
#  y_vals <- y_vals[which_indices]
#}

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

plot(
  x_seq,
  y_vals,
  xlim=c(min(x_seq), max(x_seq)),
  ylim=y_lims,
  type="n",
  xlab=expression(rho),
  ylab=expression(A[5]^2 - B[5]^2*Delta[5])
)

title(paste(round(edge_lengths, 2), collapse=" - "))

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

dev.off()
