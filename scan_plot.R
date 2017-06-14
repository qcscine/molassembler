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
searchBounds <- as.numeric(meta[2,])[1:4]

x_seq <- values$V1
y_vals <- values$V2

plot(
  x_seq,
  y_vals,
  xlim=c(min(x_seq), max(x_seq)),
  ylim=c(-300, 300),
  type="n",
  xlab=expression(rho),
  ylab=expression(A[5]^2 - B[5]^2*Delta[5])
)

title(paste(round(edge_lengths, 2), collapse=" - "))

lines(x_seq, y_vals)

abline(h=0)

abline(v=searchBounds[2])
abline(v=searchBounds[4], lty=3)
abline(v=c(searchBounds[1], searchBounds[3]), lty=2)

dev.off()
