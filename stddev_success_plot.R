values <- read.csv(file="edge-lengths-variance-success.csv", header=FALSE)

variance <- as.numeric(values[,1])
success <- as.numeric(values[,2])

pdf("variance-success.pdf", width=7, height=7)
par(mar=c(4, 4.5, 1, 1))

xLabelString <- "of sampled trunc. norm. distr."

plot(
  variance,
  success,
  type="n",
  ylim=c(0, 100),
  xlab=bquote(sigma / mu ~ .(xLabelString)),
  ylab="Root search success rate in %"
)

lines(
  variance,
  success
)

dev.off()
