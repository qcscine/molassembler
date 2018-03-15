filedata <- read.csv("asin.csv", header=FALSE, sep=",", colClasses=rep("numeric", 3))

pdf(file="asin.pdf")

# value, asin, std::asin
value <- filedata$V1
asinValues <- filedata$V2
stdAsinValues <- filedata$V3

diff <- stdAsinValues - asinValues

plot(
  value,
  diff,
  type="n",
  xlab="x",
  ylab="Deviation of asin from std::asin"
)

lines(
  value,
  diff
)

abline(h=0, col="grey60")

dev.off()
