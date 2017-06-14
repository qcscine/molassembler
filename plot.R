values <- read.csv(file="build/scan.csv", header=FALSE)

pdf("rho_pentagon.pdf", width=14, height=7)
par(
  mar=c(4, 4.5, 1, 1)
)

roots <- c(26.385, 16.512, 17.026, 17.595, 17.991, 18.335, 18.651)
rho_roots <- 1 /( roots ^ 2)

x_seq <- values$V1
y_vals <- values$V2

plot(
  x_seq,
  y_vals,
  xlim=c(min(x_seq), max(x_seq)),
  type="n",
  xlab=expression(rho),
  ylab=expression(A[5]^2 - B[5]^2*Delta[5])
)

lines(x_seq, y_vals)

abline(v=rho_roots, lty=2)
abline(h=0)

abline(v=c(0.00117, 0.00170), lty=3, col="blue")

dev.off()
