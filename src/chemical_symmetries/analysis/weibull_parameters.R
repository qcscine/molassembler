library(fitdistrplus)
source("elements.R")

# arrays are in elementArrays
# array's names are in elementArrayNames

puts <- function(x) {
  write(x, "")
}

remove_zeroes <- function(x) {
  x[which(x != 0.0)]
}

elementNames <- names(elementArrays)

for(i in 1:length(elementNames)) {
  elementArrayName <- elementNames[i]
  puts(paste("static constexpr std::array<Weibull, ", nrow(elementArray), "> ", elementArrayName, " {{", sep=""))
  elementArray <- elementArrays[[elementArrayName]]

  for(j in 1:nrow(elementArray)) {
    mu <- mean(elementArray[j,])
    sigma <- sd(elementArray[j,])
    if(sigma < 1e-4) {
      puts(paste("  --- ", j, " is ", mu, " +- ", sigma, " ---", sep=""))
    } else {
      fit.weibull <- fitdist(remove_zeroes(elementArray[j,]), "weibull")
      puts(
        paste(
          "  {",
          fit.weibull$estimate["shape"],
          ", ",
          fit.weibull$estimate["scale"],
          "},"
        )
      )
    }
  }

  puts("}};")
}
