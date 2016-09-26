pdf("both.pdf", width=8, height=5)
par(
    mfrow = c(2,1),
# margin order is bottom, left, top, right
    oma = c(5,0,0,0) + 0.1,
    mar = c(-0.1,4,1,1) + 0.1
)

data <- read.table("tanimoto_times.csv", header=FALSE, sep=",", quote="\"")
for(i in 2:4) {
    data[,i] <- data[,i] * 1000
}
colors <- c("tomato", "steelblue", "black")
linewidth <- 3

plot(
     data[,1], 
     data[,2], 
     type="n", 
     ylim=c(0, max(data[,2], data[,3], data[,4])), 
     ylab=expression(paste("Time [", mu, "s]", sep="")),
     xaxt="n"
)

for(i in 1:3) {
    lines(data[,1], data[,i+1], col=colors[i], lwd=linewidth)
    points(data[,1], data[,i+1], pch=21, bg="white")
}

legend("topleft", c("string", "bool array", "int array"), lty=1, lwd=linewidth, col=c("tomato", "steelblue", "black"))

plot(
    data[,1], 
    data[,2], 
    type="n", 
   # ylim=c(0.0004, upper), 
    ylim=c(0.4, 500),
    xlab="Fingerprint length in bits", 
    ylab=expression(paste("Time [", mu, "s] (log scale)")), 
    log="y", 
    xaxt="n"
)

for(i in 1:3) {
    lines(data[,1], data[,i+1], col=colors[i], lwd=linewidth) # bv_str
    points(data[,1], data[,i+1], pch=21, bg="white")
}

axis(1, c(64, 128, 256, 512, 1024, 2048), outer=TRUE)
#legend("topleft", c("string", "bool array", "int array"), lty=1, lwd=2, col=colors, bg="white")

title(xlab="Fingerprint length in bits", ylab="Time[ms]", outer=TRUE)

dev.off()
