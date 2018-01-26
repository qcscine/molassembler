#!/usr/bin/env Rscript

require(plyr)
library(boot)

x.errorbars <- function(x,y,xbar,ebl=0.02,...){
  arrows(x-xbar,y,x+xbar,y,code=3,angle=90,length=ebl,...)
}
y.errorbars <- function(x,y,ybar,ebl=0.02,...){
  arrows(x,y-ybar,x,y+ybar,code=3,angle=90,length=ebl,...)
}
xy.errorbars <- function(x,y,xbar,ybar,ebl=0.02,...){
  arrows(x-xbar,y,x+xbar,y,code=3,angle=90,length=ebl,...)
  arrows(x,y-ybar,x,y+ybar,code=3,angle=90,length=ebl,...)
}

filename <- "graph_timings.csv"

filedata <- read.csv(
  filename,
  header=TRUE,
  sep=",",
  quote="\"",
  colClasses=c(rep("numeric", 2)),
  check.names=FALSE,
  row.names=NULL
)

N <- filedata[,"N"]
B <- filedata[,"E"]

pdf(
  "graph_timings.pdf",
  width=10,
  height=10
)
par(
  mar=c(5, 4, 2, 1)
)

layout(
  matrix(
    c(
      1, 2,
      3, 4
    ),
    2, 2, byrow=TRUE
  )
)

color_list <- c("black", "tomato", "steelblue")

nPoints <- 100
low <- min(filedata[, "N"])
high <- max(filedata[, "N"])

lseq <- function(from=1, to=100000, length.out=6) {
  # logarithmic spaced sequence
  exp(seq(log(from), log(to), length.out = length.out))
}

linearModel <- function(D, d, which) {
  E = D[d,]

  t_averages <- E[, which]
  N <- E[, "N"]

  model <- lm(log(t_averages) ~ log(N))

  prefactor <- exp(coef(summary(model))["(Intercept)", "Estimate"])
  exponent <- coef(summary(model))["log(N)", "Estimate"]

  values <- prefactor * (lseq(low, high, length.out=nPoints) ** exponent)

  c(
    prefactor,
    exponent,
    values
  )
}

make_regression_plot <- function(which, col) {

  b = boot(filedata, linearModel, R = 500, which=which)

  cat("Average prefactor ", mean(b$t[,1]), "\n")
  cat("Average exponent ", mean(b$t[,2]), "\n")

  prefactor.ci <- boot.ci(b, type="basic", index=c(1))
  exponent.ci <- boot.ci(b, type="basic", index=c(2))

  cat(
    "Prefactor ci ",
    prefactor.ci$basic[1, 4], 
    " - ",
    prefactor.ci$basic[1, 5],
    "\n"
  )

  cat(
    "Exponent ci ",
    exponent.ci$basic[1, 4], 
    " - ",
    exponent.ci$basic[1, 5],
    "\n"
  )

  get_top_ci <- function(index) {
    boot.ci(b, type="basic", index=c(index))$basic[1, 5]
  }

  get_bottom_ci <- function(index) {
    boot.ci(b, type="basic", index=c(index))$basic[1, 4]
  }

  # top confidence interval line
  x.seq <- lseq(low, high, length.out=nPoints)
  top.ci <- sapply(seq(3, nPoints + 2), get_top_ci)
  bottom.ci <- sapply(seq(3, nPoints + 2), get_bottom_ci)

  y.max <- max(top.ci)
  y.min <- min(bottom.ci[which(bottom.ci > 1)])

  plot(
    filedata[,"N"],
    filedata[,which],
    type="n",
    xlab="N (log scale)",
    ylab="t / ns (log scale)",
    xlim=c(low, high),
    ylim=c(y.min, y.max),
    log="xy",
    main=which
  )

  prefactor.mean <- mean(b$t[,1])
  exponent.mean <- mean(b$t[,2])

  lines(
    x.seq,
    top.ci,
    col="gray60",
    lty=2
  )

  lines(
    x.seq,
    bottom.ci,
    col="gray60",
    lty=2
  )

  # regression line
  segments(
    low,
    prefactor.mean * (low ** exponent.mean),
    high,
    prefactor.mean * (high ** exponent.mean),
    col=col,
    lwd=2
  )

  points(
    filedata[, "N"],
    filedata[, which],
    pch=21,
    bg="white"
  )

  legend(
    "topleft",
    legend = c(
      paste(
        "t = ",
        round(prefactor.mean, 2), 
        " * (N ^ ",
        round(exponent.mean, 2),
        ") ns",
        sep = ""
      ),
      "Bootstrap 95% confidence intervals"
    ),
    lty=c(1, 2),
    col=c(col, "gray60")
  )

  c(prefactor.mean, exponent.mean)
}

dbm <- make_regression_plot("Floyd-Warshall & DBM", color_list[1])
lg <- make_regression_plot("Gor & LG", color_list[2])
spg <- make_regression_plot("Gor & SPG", color_list[3])

timeLines <- c(1, 1e3, 1e6, 1e9, 3.6e12, 8.64e13)
timeLabels <- c("1 ns", expression(paste("1 ", mu, "s", sep="")), "1 ms", "1 s", "1h", "1d")

comparison.xlims <- c(1, 1e4)
dbm.y <- dbm[1] * (comparison.xlims ** dbm[2])
lg.y <- lg[1] * (comparison.xlims ** lg[2])
spg.y <- spg[1] * (comparison.xlims ** spg[2])

min.y <- min(dbm.y, lg.y, spg.y)
max.y <- max(dbm.y, lg.y, spg.y)

plot(
  1,
  type="n",
  xlab="N (log scale)",
  ylab="t (log scale)",
  xlim=comparison.xlims,
  ylim=c(min.y, min(max.y, 9e13)),
  yaxt="n",
  log="xy",
  main="Comparison"
)

axis(2, at=timeLines, labels=timeLabels)

abline(h = timeLines, col="gray60")

lines(comparison.xlims, dbm.y, col=color_list[1], lwd=2)
lines(comparison.xlims, lg.y, col=color_list[2], lwd=2)
lines(comparison.xlims, spg.y, col=color_list[3], lwd=2)

legend(
  "topleft",
  legend = c(
    "Floyd-Warshall & DBM",
    "Gor & LG",
    "Gor & SPG"
  ),
  lty = c(1, 1, 1),
  col = color_list[1:3],
  bg = "white"
)

dev.off()
