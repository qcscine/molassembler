pdf("speeds.pdf", width=8, height=10)

errorbars <- function(x_pos, y_pos, sdev, ...) {
    arrows(
        x_pos, # x0
        y_pos - sdev, # y0
        x_pos, # x1
        y_pos + sdev, # y1
        length = 0.05,
        angle = 90,
        code = 3,
        ...
    )
}

minimal <- read.csv(
    "minimal.csv",
    header=FALSE
)

fast <- read.csv(
    "fast.csv",
    header=FALSE
)

# fetch plot names
filter_out_sd_values <- function(string) {
    return(substr(string, 1, 3) != "sd_")
}
plot_names <- Filter(filter_out_sd_values, minimal[,1])

par(
    mfrow = c(5,2),
# margin order is bottom, left, top, right
    oma = c(5,0,0,0) + 0.1,
    mar = c(2,4,1,1) + 0.1
)

for(name in plot_names) {
    # which column has name?
    data_index <- which(minimal[,1] == name)
    # which column has sd_name
    minimal_stdev_index <- which(
        minimal[,1] == paste0(
            "sd_",
            substr(
                name,
                4,
                nchar(name)
            )
        )
    )
    # which column has sd_name
    fast_stdev_index <- which(
        fast[,1] == paste0(
            "sd_",
            substr(
                name,
                4,
                nchar(name)
            )
        )
    )

    # x sequence
    x_seq <- 10 + 55*seq(1, 19)

    # plot
    plot(
        x_seq,
        minimal[data_index, 2:length(minimal[data_index,])],
        type = "n",
        ylim = c(
            min(
                0,
                as.vector(unlist(
                    minimal[data_index, 2:length(minimal[data_index,])]
                    - minimal[minimal_stdev_index, 2:length(minimal[minimal_stdev_index,])]
                )),
                as.vector(unlist(
                    fast[data_index, 2:length(fast[data_index,])]
                    - fast[fast_stdev_index, 2:length(fast[fast_stdev_index,])]
                ))
            ),
            max(
                as.vector(unlist(
                    minimal[data_index, 2:length(minimal[data_index,])]
                    + minimal[minimal_stdev_index, 2:length(minimal[minimal_stdev_index,])]
                )),
                as.vector(unlist(
                    fast[data_index, 2:length(fast[data_index,])]
                    + fast[fast_stdev_index, 2:length(fast[fast_stdev_index,])]
                ))
            )
        ),
        xlab = "N",
        ylab = paste(substr(name, 4, nchar(name)), " [ns]")
    )

    lines(
        x_seq,
        minimal[data_index, 2:length(minimal[data_index,])], 
        col = "black"
    )
    polygon(
        c(x_seq, rev(x_seq)),
        c(
            as.vector(unlist(
                minimal[data_index, 2:length(minimal[data_index,])]
                + minimal[minimal_stdev_index, 2:length(minimal[minimal_stdev_index,])]
            )),
            rev(as.vector(unlist(
                minimal[data_index, 2:length(minimal[data_index,])]
                - minimal[minimal_stdev_index, 2:length(minimal[minimal_stdev_index,])]
            )))
        ),
        col=rgb(0,0,0,0.2),
        border=NA
    )

    lines(
        x_seq,
        fast[data_index, 2:length(minimal[data_index,])],
        col = "red"
    )
    polygon(
        c(x_seq, rev(x_seq)),
        c(
            as.vector(unlist(
                fast[data_index, 2:length(fast[data_index,])]
                + fast[fast_stdev_index, 2:length(fast[fast_stdev_index,])]
            )),
            rev(as.vector(unlist(
                fast[data_index, 2:length(fast[data_index,])]
                - fast[fast_stdev_index, 2:length(fast[fast_stdev_index,])]
            )))
        ),
        col=rgb(1,0,0,0.2),
        border=NA
    )
}
