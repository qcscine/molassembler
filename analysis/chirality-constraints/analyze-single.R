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

# Columns are
# V1 -> centerAtom
# V2 -> index of symmetry in vector of symmetries (0-based)
# V3 -> assignment
# V4 -> total number of assignments
# V5 -> angle deviation
# V6 -> oneThreeDeviation

makeExampleOfAngleDeviationVsOneThreeDeviation <- TRUE
makeMultiplicityHistograms <- TRUE
makeMultiplicityBarPlot <- TRUE

if(makeExampleOfAngleDeviationVsOneThreeDeviation) {

  makePlotOnDiagonal <- FALSE

  nTotalPlots <- 0
  for(i in 1:length(symmetries)) {
    nTotalPlots <- nTotalPlots + length(
      which(coordinationNumber == coordinationNumber[i])
    )
  } # 42 currently

  print(nTotalPlots)

  nRows <- 10
  nCols <- 3

  pdf("xyPlot.pdf", width=nCols*7, height=nRows*7)
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
        # We're in file symmetries[i].csv -> whose V2 references the 0 based
        toJData <- filedata[
          # rows with $V2 = j - 1 ( -1 because 
          which(as.numeric(filedata$V2) == j - 1),
        ]

        print(length(toJData$V5))

        angleDeviations <- as.numeric(
          toJData$V5
        )

        oneThreeDeviations <- as.numeric(
          toJData$V6
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

      # collect the data for same-symmetry fits
      assignmentsDeviation <- c()
      for(lineNum in 1:min(filesize, 10000)) {
        symmetryNum <- as.numeric(filedata[lineNum,]$V2) + 1

        if(symmetryNum == i) {
          assignmentNum <- as.numeric(filedata[lineNum,]$V3) + 1
          totalAssignments <- as.numeric(filedata[lineNum,]$V4)
          totalDeviation <- as.numeric(filedata[lineNum,]$V5) +
            as.numeric(filedata[lineNum,]$V6)

          assignmentsDeviation <- c(
            assignmentsDeviation,
            totalDeviation
          )

          if(assignmentNum == totalAssignments) { 
            # new fit set -> process the existing ones

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

      for(lineNum in 1:min(filesize, 10000)) {
        symmetryNum <- as.numeric(filedata[lineNum,]$V2) + 1

        if(symmetryNum == i) {
          assignmentNum <- as.numeric(filedata[lineNum,]$V3) + 1
          totalAssignments <- as.numeric(filedata[lineNum,]$V4)
          totalDeviation <- as.numeric(filedata[lineNum,]$V5) +
            as.numeric(filedata[lineNum,]$V6)


          assignmentsDeviation <- c(
            assignmentsDeviation,
            totalDeviation
          )


          if(assignmentNum == totalAssignments) { 
            # new fit set -> process the existing ones

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
