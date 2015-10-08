#   Give us some nice colors
library("RColorBrewer")
#   Let us use mathematical notation in plot labels
library("grDevices")
require(graphics)
#   Read in the mutation rate file
mutation_rates <- read.table("Mutation_Rates.txt", header=T)
#   Get columns that actually contain rate information
r <- grep("Rate", names(mutation_rates), value=T)
#   We only want to plot the overall, the transition, and the transversion
r <- r[c(1:3)]
#   Then save those columns, and remove the first row, since that is M92-220
mut <- as.matrix(mutation_rates[,r])
#   Then, we take the means and standard deviations of the different classes
#   The groups are:
#       Rows 6-10, 12, 13: Gy16
#       Rows 11, 14, 15: Gy32
#       Rows 3, 4: GTTC
gy16_rows <- c(6, 7, 8, 9, 10, 12, 13)
gy32_rows <- c(11, 14, 15)
gttc_rows <- c(3, 4)
gy16_means <- apply(mut[gy16_rows,], 2, mean)
gy16_sds <- apply(mut[gy16_rows,], 2, sd)
gy32_means <- apply(mut[gy32_rows,], 2, mean)
gy32_sds <- apply(mut[gy32_rows,], 2, sd)
# wpt_means <- apply(mut[c(9:10),], 2, mean)
# wpt_sds <- apply(mut[c(9:10),], 2, sd)
#   put the means all together into a matrix for barplots
plot_data <- rbind(
    gy16_means,
    gy32_means,
    mut[3,],
    mut[4,]
    )
#   Convert negative values to juuuuuust below 0
plot_data[plot_data <= 0] <- -2e-10
pdf(file="Rate_Plot_Groups.pdf", 6, 8)
original_params <- par()
par(
    xaxt="n"
    )
plt <- barplot(
    plot_data,
    xlab="Mutational Class",
    ylab=expression(paste("Mutation Rate (", NULL %*% 10^{-8}, " bp", NULL^{-1}, ")")),
    beside=TRUE,
    col=brewer.pal(4, "Set1"),
    space=c(0, 1.5),
    main="Estimated Induced Mutational Spectrum",
    ylim=c(0, 1e-7),
    axes=F)
# abline(
#     h=5.9e-9,
#     lty=1,
#     lwd=1,
#     col="grey40"
#     )
#   Add the standard deviation bars
#   This adds the vertical lines for the bars: +/- 2.16*sd
segments(
    plt[1,],
    plot_data[1,] - gy16_sds,
    plt[1,],
    plot_data[1,] + gy16_sds,
    col="black",
    lwd=1
    )
#   And then this adds the horizontal marks, just for appearances
segments(
    plt[1,] - 0.2,
    plot_data[1,] - gy16_sds,
    plt[1,] + 0.2,
    plot_data[1,] - gy16_sds,
    col="black",
    lwd=1
    )
segments(
    plt[1,] - 0.2,
    plot_data[1,] + gy16_sds,
    plt[1,] + 0.2,
    plot_data[1,] + gy16_sds,
    col="black",
    lwd=1
    )
#   Repeat for the other bars!
#   Vertical segment
segments(
    plt[2,],
    plot_data[2,] - gy32_sds,
    plt[2,],
    plot_data[2,] + gy32_sds,
    col="black",
    lwd=1
    )
#   Lower ticks
segments(
    plt[2,] - 0.2,
    plot_data[2,] - gy32_sds,
    plt[2,] + 0.2,
    plot_data[2,] - gy32_sds,
    col="black",
    lwd=1
    )
#   Upper ticks
segments(
    plt[2,] - 0.2,
    plot_data[2,] + gy32_sds,
    plt[2,] + 0.2,
    plot_data[2,] + gy32_sds,
    col="black",
    lwd=1
    )
#   Vertical segment
# segments(
#     plt[3,],
#     plot_data[3,] - 2.16*wpt_sds,
#     plt[3,],
#     plot_data[3,] + 2.16*wpt_sds,
#     col="black",
#     lwd=1
#     )
# #   Lower ticks
# segments(
#     plt[3,] - 0.2,
#     plot_data[3,] - 2.16*wpt_sds,
#     plt[3,] + 0.2,
#     plot_data[3,] - 2.16*wpt_sds,
#     col="black",
#     lwd=1
#     )
# #   Upper ticks
# segments(
#     plt[3,] - 0.2,
#     plot_data[3,] + 2.16*wpt_sds,
#     plt[3,] + 0.2,
#     plot_data[3,] + 2.16*wpt_sds,
#     col="black",
#     lwd=1
#     )

axis(
    side=1,
    at=apply(plt, 2, mean),
    labels=FALSE,
    tick=FALSE)
text(
    apply(plt, 2, mean),
    -5e-9,
    adj=c(0.5, 0.5),
    srt=0,
    labels=c(
        "Overall",
        "Transition",
        "Transversion"
        # expression(paste("A:T", NULL%->%NULL, "G:C")),
        # expression(paste("G:C", NULL%->%NULL, "A:T")),
        # expression(paste("A:T", NULL%->%NULL, "C:G")),
        # expression(paste("A:T", NULL%->%NULL, "T:A")),
        # expression(paste("G:C", NULL%->%NULL, "T:A")),
        # expression(paste("C:G", NULL%->%NULL, "G:C"))
        ),
    xpd=TRUE,
    cex=0.8)
axis(
    side=2,
    at=seq(0, 1e-7, by=1e-8),
    labels=FALSE,
    tick=TRUE,
    lty=1,
    lwd=1
    )
text(
    0,
    y=seq(0, 1e-7, by=1e-8),
    labels=c("0.0", "1.0", "2.0", "3.0", "4.0", "5.0", "6.0", "7.0", "8.0", "9.0", "10.0"),
    xpd=TRUE,
    cex=0.75)
legend(
    "topright",
    as.vector(c("16 gray", "32 gray", "WPT389", "WPT391")),
    fill=brewer.pal(4, "Set1"),
    cex=1
    )
dev.off()
