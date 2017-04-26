args <- commandArgs(TRUE)
filename <- toString(args[1])

filedata <- read.csv(
  filename,
  header=FALSE,
  sep=",",
  colClasses=c(rep("numeric", 2))
)

# format is
# V1 - error
# V2 - gradient norm
# V3 - compress
# V4 - sum abs 4D

error <- filedata$V1
gradientNorm <- filedata$V2
compress <- filedata$V3
totalAbs4D <- filedata$V4

lastUncompressed <- tail(which(compress == 0), n=1)
firstCompressed <- head(which(compress == 1), n=1)

xSeq <- seq(1, length(error))

baseName <- strsplit(filename, "\\.")
pdf(
  paste(baseName[[1]][1], ".pdf", sep=""), 
  width=14,
  height=7
)
par(
  mfrow=c(3, 1),
  mar=c(0, 4, 0, 1),
  oma=c(4, 0, 1, 0)
)

# error plot
plot(
  xSeq,
  error,
  type="n",
  xlab="",
  ylab="Error (log scale)",
  log="y",
  xaxt="n"
)

lines(
  xSeq,
  error
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
  gradientNorm
)

# sum abs 4D
plot(
  xSeq,
  totalAbs4D,
  type="n",
  xlab="Step number",
  ylab="Total absolute fourth dimension component"
)

lines(
  xSeq,
  totalAbs4D
)

abline(v=mean(lastUncompressed, firstCompressed), col="blue")

dev.off()
