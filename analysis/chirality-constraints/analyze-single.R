# symmetries in same order as Symmetries::allNames from symmetry_information
symmetries <- c(
  "linear", # 2
  "bent",
  "trigonal-planar", # 3
  "trigonal-pyramidal",
  "T-shaped", 
  "tetrahedral", # 4
  "square-planar",
  "seesaw",
  "square-pyramidal", # 5
  "trigonal-bipyramidal",
  "pentagonal-planar",
  "octahedral", #6
  "trigonal-prismatic",
  "pentagonal-pyramidal",
  "pentagonal-bipyramidal", # 7
  "square-antiprismatic" # 8
)

coordinationNumber <- c(2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 8)

# Columns are now
# V1 -> index of symmetry in vector of symmetries (0-based)
# V2 -> assignment
# V3 -> angle deviation
# V4 -> oneThreeDeviation
# V5 -> chiralityDeviation

makePerSymmetryPlots <- FALSE
makeAngleDeviationVsOneThreePlots <- FALSE
makeMultiplicityHistograms <- TRUE
makeMultiplicityBarPlot <- TRUE

if(makePerSymmetryPlots) {
  colors <- c("steelblue", "tomato", "forestgreen", "darkorchid", "darkorange")

  makePlotOnDiagonal <- FALSE

  for(i in 1:length(symmetries)) {

    pdf(
      paste(symmetries[i], ".pdf", sep=""),
      width=10,
      height=10
    )
    par(mfrow=c(1, 1))

    filename <- paste(symmetries[i], ".csv", sep="")

    filesize <- length(readLines(filename))

    if(filesize > 0) {
      filedata <- read.csv(
        file=filename,
        header=FALSE,
        sep=","
      )

      angleDeviations <- as.numeric(
        filedata$V3
      )

      oneThreeDeviations <- as.numeric(
        filedata$V4
      )

      targetSymmetry <- as.numeric(
        filedata$V1
      ) + 1

      targetSymmetryUniques <- unique(targetSymmetry)

      if(makePlotOnDiagonal) {
        minCoords <- c(min(angleDeviations), min(oneThreeDeviations))
        maxCoords <- c(max(angleDeviations), max(oneThreeDeviations))

        plotLimits <- c(min(minCoords), max(maxCoords))

        plot(
          angleDeviations,
          oneThreeDeviations,
          xlab="Angle deviations in radians (log scale)",
          ylab="1-3 deviations in Angstrom (log scale)",
          xlim=plotLimits,
          ylim=plotLimits,
          type="n",
          log="xy"
        )

      } else {

        plot(
          angleDeviations,
          oneThreeDeviations,
          xlab="Angle deviations in radians (log scale)",
          ylab="1-3 deviations in Angstrom (log scale)",
          type="n",
          log="xy"
        )

      }

      title(
        main=paste("Fit onto generated", symmetries[i])
      )

      abline(0, 1, col="gray50")

      if(length(angleDeviations) > 1000) { # reduce set to sample
        indicesSubset <- sample(
          seq(1, length(angleDeviations)),
          1000
        )

        angleDeviations <- angleDeviations[indicesSubset]
        oneThreeDeviations <- oneThreeDeviations[indicesSubset]
        targetSymmetry <- targetSymmetry[indicesSubset]
      }

      for(i in 1:length(targetSymmetryUniques)) { # make colored points groups
        subsetIndices <- which(targetSymmetry == targetSymmetryUniques[i])

        points(
          angleDeviations[subsetIndices],
          oneThreeDeviations[subsetIndices],
          pch=21,
          col="black",
          bg=colors[i]
        )
      }

      legend(
        "topleft",
        legend = symmetries[targetSymmetryUniques],
        col = colors[seq(1, length(targetSymmetryUniques))],
        lty = 1,
        lwd = 3
      )
    }
  }

  dev.off()
}

if(makeAngleDeviationVsOneThreePlots) {

  makePlotOnDiagonal <- FALSE

  nTotalPlots <- 0
  for(i in 1:length(symmetries)) {
    nTotalPlots <- nTotalPlots + length(
      which(coordinationNumber == coordinationNumber[i])
    )
  } # 42 currently

  print(nTotalPlots)

  nRows <- 6
  nCols <- 7

  pdf("xyPlot.pdf", width=nCols*4, height=nRows*4)
  par(mfrow=c(nRows, nCols))

  for(i in 1:length(symmetries)) {
    # print(paste("i = ", i, ", symmetries[i]=", symmetries[i], ", coordinationNumber[i]=", coordinationNumber[i], sep=""))
    filename <- paste(symmetries[i], ".csv", sep="")

    filesize <- length(readLines(filename))

    if(filesize > 0) {
      filedata <- read.csv(
        file=filename,
        header=FALSE,
        sep=","
      )


      for(j in which(coordinationNumber == coordinationNumber[i])) {
        # print(paste("From ", symmetries[i], " (size ", coordinationNumber[i], ") to", symmetries[j], "(size ", coordinationNumber[j], ")", sep=""))
        # We're in file symmetries[i].csv -> whose V1 references the 0 based symmetry name
        toJData <- filedata[
          which(as.numeric(filedata$V1) == j - 1),
        ]

        print(length(toJData$V3))

        angleDeviations <- as.numeric(
          toJData$V3
        )

        oneThreeDeviations <- as.numeric(
          toJData$V4
        )


        if(makePlotOnDiagonal) {
          minCoords <- c(min(angleDeviations), min(oneThreeDeviations))
          maxCoords <- c(max(angleDeviations), max(oneThreeDeviations))

          plotLimits <- c(min(minCoords), max(maxCoords))

          plot(
            angleDeviations,
            oneThreeDeviations,
            xlab="Angle deviations in radians",
            ylab="1-3 deviations in Angstrom",
            xlim=plotLimits,
            ylim=plotLimits,
            type="n"
          )

        } else {

          plot(
            angleDeviations,
            oneThreeDeviations,
            xlab="Angle deviations in radians",
            ylab="1-3 deviations in Angstrom",
            type="n"
          )

        }

        title(
          main=paste("Fit", symmetries[j], "onto generated", symmetries[i])
        )

        abline(0, 1, col="gray50")

        pointscolor <- rgb(0, 0, 0, 0.5)

        if(length(angleDeviations) > 1000) {
          indicesSubset <- sample(
            seq(1, length(angleDeviations)),
            1000
          )
          points(
            angleDeviations[indicesSubset],
            oneThreeDeviations[indicesSubset],
            pch=21,
            col=pointscolor,
            bg=pointscolor
          )
        } else {
          points(
            angleDeviations,
            oneThreeDeviations,
            pch=21,
            col=pointscolor,
            bg=pointscolor
          )
        }
      }
    } 
  }

  dev.off()
}

# histogram of multiplicity of lowest total deviation for symmetric deviations

if(makeMultiplicityHistograms) {
  pdf("multiplicityHistograms.pdf", width=7, height=7)
  par(mfrow=c(1, 1))

  for(i in 1:length(symmetries)) {
    filename <- paste(symmetries[i], ".csv", sep="")
    filesize <- length(readLines(filename))
    
    print(paste(filename, "has length", filesize))

    if(filesize > 0) {
      filedata <- read.csv(
        file=filename,
        header=FALSE,
        sep=","
      )

      multiplicityOfLowestDeviation <- c()

      maxLines <- 10000

      # collect the data for same-symmetry fits
      assignmentsDeviation <- c()

      for(lineNum in 1:min(filesize, maxLines)) {
        symmetryNum <- as.numeric(filedata[lineNum,]$V1) + 1

        if(symmetryNum == i) {
          assignmentNum <- as.numeric(filedata[lineNum,]$V2)

          totalDeviation <- as.numeric(filedata[lineNum,]$V3) +
            as.numeric(filedata[lineNum,]$V4) +
            as.numeric(filedata[lineNum,]$V5)

          if( 
            ( # if we have just finished reading a full set
             assignmentNum == 0 && length(assignmentsDeviation) > 0
            ) || ( # or this is the last line
              lineNum == min(filesize, maxLines) 
            )
          ) {
            # add the current one if we're at the last line only
            if(lineNum == min(filesize, maxLines)) {
              assignmentsDeviation <- c(
                assignmentsDeviation,
                totalDeviation
              )
            }

            multiplicityOfLowestDeviation <- c(
              multiplicityOfLowestDeviation,
              length(
                which(
                  assignmentsDeviation == min(assignmentsDeviation)
                )
              )
            )

            # empty the list
            assignmentsDeviation <- c()
          } 

          # add to vector
          assignmentsDeviation <- c(
            assignmentsDeviation,
            totalDeviation
          )
        }
      }

      hist(
        multiplicityOfLowestDeviation,
        main=paste(
          "Histogram of multiplicity of lowest deviation of ", 
          symmetries[i],
          " fit to generated ",
          symmetries[i],
          " structure.",
          sep=""
        )
      )
    }
  }

  dev.off()
}

if(makeMultiplicityBarPlot) {
  pdf("multiplicityBarPlot.pdf", width=7, height=7)
  par(mfrow=c(1, 1))

  symmetryMultiplicities <- c()
  symmetryMultiplicityAllSame <- c()


  for(i in 1:length(symmetries)) {
    filename <- paste(symmetries[i], ".csv", sep="")
    filesize <- length(readLines(filename))
    
    print(paste(filename, "has length", filesize))

    if(filesize > 0) {
      filedata <- read.csv(
        file=filename,
        header=FALSE,
        sep=","
      )

      multiplicityOfLowestDeviation <- c()

      # collect the data for same-symmetry fits
      assignmentsDeviation <- c()

      for(lineNum in 1:min(filesize, maxLines)) {
        symmetryNum <- as.numeric(filedata[lineNum,]$V1) + 1

        if(symmetryNum == i) {
          assignmentNum <- as.numeric(filedata[lineNum,]$V2)

          totalDeviation <- as.numeric(filedata[lineNum,]$V3) +
            as.numeric(filedata[lineNum,]$V4) +
            as.numeric(filedata[lineNum,]$V5)

          if( 
            ( # if we have just finished reading a full set
             assignmentNum == 0 && length(assignmentsDeviation) > 0
            ) || ( # or this is the last line
              lineNum == min(filesize, maxLines) 
            )
          ) {
            # add the current one if we're at the last line only
            if(lineNum == min(filesize, maxLines)) {
              assignmentsDeviation <- c(
                assignmentsDeviation,
                totalDeviation
              )
            }

            multiplicityOfLowestDeviation <- c(
              multiplicityOfLowestDeviation,
              length(
                which(
                  assignmentsDeviation == min(assignmentsDeviation)
                )
              )
            )

            # empty the list
            assignmentsDeviation <- c()
          } 

          # add to vector
          assignmentsDeviation <- c(
            assignmentsDeviation,
            totalDeviation
          )
        }
      }

      symmetryMultiplicityAllSame <- c(
        symmetryMultiplicityAllSame, 
        length(
          which(
            multiplicityOfLowestDeviation == multiplicityOfLowestDeviation[1]
          )
        ) == length(multiplicityOfLowestDeviation)
      )

      # print(multiplicityOfLowestDeviation[1:min(100, length(multiplicityOfLowestDeviation))])

      symmetryMultiplicities <- c(
        symmetryMultiplicities,
        mean(multiplicityOfLowestDeviation)
      )
    } else {
      symmetryMultiplicities <- c(
       symmetryMultiplicities,
       0
      )
      symmetryMultiplicityAllSame <- c(
        symmetryMultiplicityAllSame,
        TRUE
      )
    }
  }

  if(!all(symmetryMultiplicityAllSame)) {
    print("NOT ALL MULTIPLICITIES ARE ALWAYS IDENTICAL!")
  } 

  bp <- barplot(
    symmetryMultiplicities,
    xlab="Average multiplicity of lowest deviation of symmetric assignment fit",
    ylab="Symmetry name",
    horiz=TRUE,
    col="black"
  )

  axis(2, at=bp, labels=coordinationNumber)

  for(i in 1:length(symmetries)) {
    if(symmetryMultiplicities[i] >= 1) {
      text(
        symmetryMultiplicities[i],
        bp[i],
        labels=symmetries[i],
        pos=2,
        col="white"
      )
    } else {
      text(
        symmetryMultiplicities[i],
        bp[i],
        labels=symmetries[i],
        pos=4,
        col="black"
      )
    }
  }

  dev.off()
}


# generated symmetry
#  -> list of attempted symmetries of same size
#  -> list of attempted assignments
#  -> tuple of angle deviations, 1-3 deviations
