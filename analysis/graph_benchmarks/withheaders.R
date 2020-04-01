#!/usr/bin/env Rscript

require(plyr)

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
  width=8,
  height=12
)
par(
  mar=c(5, 4, 2, 1)
  #oma=c(4, 0, 1, 0),
)

layout(
  matrix(
    c(
      1, 4,
      2, 4,
      3, 4
    ), 
    3, 2, byrow=TRUE)
)

color_list <- c("black", "tomato", "steelblue", "forestgreen", "purple")

make_regression_plot <- function(t_averages, t_stddevs, name, color, linetype) {
  if(any(t_stddevs > t_averages)) {
    print("Sigma is larger than average value! This will not work on a log-log plot.")
  }

  plot(
    N,
    t_averages,
    type="n",
    xlab="N (log scale)",
    ylab="t / ns (log scale)",
    log="xy",
    main=name
  )

  model <- lm(log(t_averages) ~ log(N))

  prefactor <- exp(coef(summary(model))["(Intercept)", "Estimate"])
  exponent <- coef(summary(model))["log(N)", "Estimate"]

  xSeq <- seq(1, max(N))
  predicted <- prefactor * (xSeq ^ exponent)

  lines(
    xSeq,
    predicted,
    lwd=1.5,
    col=color,
    lty=linetype
  )

  y.errorbars(
    N,
    t_averages,
    t_stddevs
  )

  points(
    N,
    t_averages,
    pch=21,
    col="black",
    bg="white"
  )

  legend(
    "topleft",
    legend = c(paste("t = ", round(prefactor, 2), " * (N ^ ", round(exponent, 2), ") ns", sep="")),
    lty = linetype,
    col=color
  )

  list(prefactor, exponent, name, color, linetype)
}

print(paste(ncol(filedata), ncol(filedata), length(colnames(filedata)), length(color_list)), 2)

# DBM
regressions <- rbind(
  make_regression_plot(filedata[,3], filedata[,4], colnames(filedata)[3], color_list[1], 1)
)

for(i in seq(1:(ncol(filedata)-4))) {
  if(i %% 2 == 1) {
    print(paste(i + 4, i + 5, i + 4, 1 + ceiling(i / 2) %% 6, 1))

    regressions <- rbind(
      regressions,
      make_regression_plot(
        filedata[,i + 4], 
        filedata[,i + 5], 
        colnames(filedata)[i + 4],
        color_list[1 + ceiling(i / 2) %% 6],
        1
      )
    )
  }
}

# Summary graph
comparison.xlims <- c(1, 1e4)

make_y_values <- function(prefactor, exponent) {
  prefactor * (comparison.xlims ^ exponent)
}

# plot regressions
yValues <- adply(
  regressions,
  1,
  function(row) make_y_values(row[[1]], row[[2]])
)

lowestNColumn <- as.numeric(yValues$V1)
minY <- min(lowestNColumn)

highestNColumn <- as.numeric(yValues$V2)
maxY <- max(highestNColumn)

plot(
  1,
  type="n",
  xlab="N (log scale)",
  ylab="t / ns (log scale)",
  xlim=comparison.xlims,
  ylim=c(minY, min(maxY, 9e13)),
  yaxt="n",
  log="xy",
  main="Comparison"
)

timeLines <- c(1, 1e3, 1e6, 1e9, 3.6e12, 8.64e13)

axis(2, at=timeLines, labels=c("1 ns", expression(paste("1 ", mu, "s", sep="")), "1 ms", "1 s", "1h", "1d"))

abline(h = timeLines, col="gray60")

add_lines <- function(prefactor, exponent, color, linetype) {
  y_seq <- prefactor * (comparison.xlims ^ exponent)
  lines(
    comparison.xlims,
    y_seq,
    lwd = 1.5,
    col=color,
    lty=linetype
  )
}

adply(
  regressions,
  1,
  function(row) add_lines(row[[1]], row[[2]], row[[4]], row[[5]])
)

print(regressions)

legend(
  "topleft",
  legend = regressions[,3],
  col = as.character(regressions[,4]),
  lty = as.numeric(regressions[,5]),
  bg="white"
)

dev.off()
