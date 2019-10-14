source("recognition_data.R")

# count elements for particular symmetry and recognizer
barplot_matrix <- function(symmetryIndex, recognizer) {
  # Count element occurrences function
  table_f <- function(array) { as.data.frame(table(array), stringsAsFactors = FALSE) }

  # Apply the counting to each column of the selected matrix. Naturally, counts
  # are reduced to the elements actually present, so if the list of results for
  # a particular distortion value only contains a single symmetry name index,
  # then none of the others of equal size are present in the table. This
  # complicates things slightly.
  tables <- apply(results[,,recognizer,symmetryIndex], c(2), table_f)
  # Figure out the expected symmetry indices in the data (0-based)
  expected_symmetry_indices <- which(symmetrySizes == symmetrySizes[symmetryIndex]) - 1

  # Create a matrix (symmetries, nDistortionValues) containing the counts of
  # each symmetry for each distortion value. This can be used for barplot()
  m <- matrix(0, length(expected_symmetry_indices), length(tables))
  for(i in 1:length(tables)) {
    for(j in 1:nrow(tables[[i]])) {
      rowSymmetryIndex <- as.numeric(tables[[i]][j,]$array)
      # a column for each distortion 
      m[which(expected_symmetry_indices == rowSymmetryIndex), i] <- tables[[i]][j,]$Freq
    }
  }

  # Collect the label positions (j-i sequenced since the offset has to be kept)
  xlabel_positions <- m
  for(j in 1:ncol(m)) {
    offset <- 0
    for(i in 1:nrow(m)) {
      xlabel_positions[i, j] <- offset + m[i, j] / 2
      offset <- offset + m[i, j]
    }
  }

  # Create the plot
  bp <- barplot(
    m,
    horiz=TRUE,
    main=paste(c(recognizers[recognizer], shapeNames[symmetryIndex])),
    xlab="Identified symmetry frequencies",
    ylab="Distortion vector norms",
    names.arg=x.values
  )
  # Add labels to each bar with the symmetry name
  xcolors <- c("white", rep("black", each=length(expected_symmetry_indices) - 1))
  # We have to i-j sequence label collection
  xlabels <- c()
  for(i in 1:nrow(m)) {
    for(j in 1:ncol(m)) {
      # Populate the current label only if the name can fit into the bar
      if(m[i, j] > strwidth(shapeNames[expected_symmetry_indices[i] + 1]) + 1) {
        xlabels <- c(xlabels, shapeNames[expected_symmetry_indices[i] + 1])
      } else {
        xlabels <- c(xlabels, "")
      }
    }
  }
  text(
    t(xlabel_positions),
    bp,
    labels=xlabels,
    col=rep(xcolors, each=ncol(m))
  )
}

multibarplot <- function(symmetryIndex) {
  recognizerCount <- length(recognizers)
  plotWidth <- 4
  plotHeight <- 4
  pdf(paste("multibarplot_", symmetryIndex, ".pdf", sep=""), width=recognizerCount * plotWidth, height=plotHeight)
  par(
    fin=c(plotWidth, plotHeight),
    mfrow=c(1, recognizerCount)
  )
  for(i in 1:recognizerCount) {
    barplot_matrix(symmetryIndex, i)
  }
  dev.off()
}

for(i in 1:length(shapeNames)) {
  multibarplot(i)
}

compare_algorithms_plot <- function(symmetryIndex) {
  gather_recognizer_data <- function(recognizerIndex) {
    table_f <- function(array) { as.data.frame(table(array), stringsAsFactors=FALSE) }
    tables <- apply(results[,,recognizerIndex,symmetryIndex], c(2), table_f)
    count_f <- function(table) {
      row_index <- which(table$array == symmetryIndex - 1)
      if(length(row_index) == 0) {
        return(0)
      }

      table[row_index,]$Freq
    }
    counts <- sapply(tables, count_f)
  }

  plot(
    x.values,
    rep(c(1), length(x.values)), 
    type="n",
    ylim=c(0, 100),
    xlab="Distortion vector norms",
    ylab="Frequency correctly identified symmetry",
    main=paste("Algorithm comparison", shapeNames[symmetryIndex])
  )
  colors <- c("black", "tomato", "steelblue", "slateblue", "yellowgreen")
  for(i in 1:length(recognizers)) {
    y.values <- gather_recognizer_data(i)
    lines(x.values, y.values, col=colors[i], lwd=2)
  }
  legend("topright", recognizers, lwd=2, col=colors)
}
