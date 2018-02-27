#' @export
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

#' @export
#' @import ggplot2
barplot_with_ci <- function(names, data, ci1, ci2) {
    df <- data.frame(names=names, data=data, ci1=ci1, ci2=ci2)
    dodge <- position_dodge(width = 0.9)
    plot <- ggplot(data = df,
                            aes(x = names, y = data, fill = names))
    plot <- plot + geom_bar(stat = "identity", position = dodge)
    plot <- plot + geom_errorbar(position=dodge, width=0.25,
                                 aes(ymin=ci1, ymax=ci2))
    plot
}

#' @export
#' @import ggplot2
point_with_ci <- function(df) {
    dodge <- position_dodge(width = 0.9)
    plot <- ggplot(data = df,
                            aes(x = sample, y = data))
    plot <- plot + geom_point(stat = "identity", position = dodge)
    plot <- plot + geom_errorbar(position=dodge, width=0.25,
                                          aes(ymin=ci1, ymax=ci2))
    plot <- plot + theme_shared()
    plot
}


#' @export
#' @import ggplot2
#' @importFrom rworldmap getMap
extraction_map <- function() {
    newmap <- getMap(resolution = "low")
    plot(newmap, xlim = c(-130, 40), ylim = c(30, 65), asp = 2,
         lwd=1.55, col="white", bg="#96c1e1",border="#6C6C71")
    box(lwd=10, col="black")
    # germany_coords <- matrix(c(8.5, 10, 12,
    #                          52, 50, 52), ncol=2)
    # finland_coords <- matrix(c(26, 28, 28,
    #                          63, 62, 64), ncol=2)
    magna_coords <- matrix(c(10, 28, 34,
                             50, 64, 31), ncol=2)
    pulex_coords <- matrix(c(-80, -124, 43, 43), ncol=2)
    points(magna_coords, col="red", pch=4, cex=2, lwd=5)
    points(pulex_coords, col="blue", pch=5, cex=2, lwd=5)
}

#' @export
#' @import RColorBrewer
#' @import ggplot2
snv_barplot <- function(var_table) {
    require(RColorBrewer)
    snv_var_table <- var_table[which(var_table$class == "snv"),]

    nts <- length(which(var_table$ts_tv=="transition"))
    ntv <- length(which(var_table$ts_tv=="transversions"))
    a_t <- length(which((var_table$ref %in% c("A", "T")) &
                            (var_table$alt %in% c("A", "T"))))
    g_t <- length(which((var_table$ref %in% c("G", "T")) &
                            (var_table$alt %in% c("G", "T"))))
    g_c <- length(which((var_table$ref %in% c("G", "C")) &
                            (var_table$alt %in% c("G", "C"))))
    a_c <- length(which((var_table$ref %in% c("A", "C")) &
                            (var_table$alt %in% c("A", "C"))))
    a_g <- length(which((var_table$ref %in% c("A", "G")) &
                            (var_table$alt %in% c("A", "G"))))
    c_t <- length(which((var_table$ref %in% c("C", "T")) &
                            (var_table$alt %in% c("C", "T"))))
    ts_tv <- c(rep("Transversion", 4), "Transition", "Transition")
    snv_df <- data.frame(id=c("A-T", "G-T", "G-C", "A-C", "A-G",
                              "C-T"), ts_tv, n=c(a_t, g_t, g_c, a_c, a_g, c_t),
                         stringsAsFactors = FALSE)
    snv_df <- snv_df[with(snv_df, order(-n)), ]
    snv_df$id <- reorder(snv_df$id, rev(snv_df$n))
    ggplot(snv_df, aes(x=id, y=n)) +
        geom_bar(stat="Density",fill=gg_color_hue(3)[3]) +
        ylab("") +
        theme_shared() +
        labs(fill="")
}


# For subsampled analysis
# iso <- unique(line_info$isolate)
# bs <- unlist(lapply(1:length(base_surv), function(i) {
#     isonum <- length(which(line_info$isolate==iso[i]))
#     rep(base_surv[i], isonum)
# }))



#' @export
#' @import ggplot2
#' @importFrom dplyr filter
ins_del_stacked_bar <- function(var_table) {
    indel_var_table <- var_table[which(var_table$class %in%
                                           c("insertion", "deletion")),]
    ninsertions <- length(which(indel_var_table$class == "insertion"))
    ndeletions <- length(which(indel_var_table$class == "deletion"))
    ins_lengths <- nchar(filter(indel_var_table, class=="insertion")$alt) - 1
    del_lengths <- nchar(filter(indel_var_table, class=="deletion")$ref) - 1
    id <- c(rep("Deletion", length(del_lengths)),
            rep("Insertion", length(ins_lengths)))
    lengths <- c(sort(del_lengths, decreasing=TRUE),
                 sort(ins_lengths, decreasing=TRUE))
    df <- data.frame(id=factor(id), lengths=lengths)
    ggplot(df, ggplot2::aes(x=id, y=lengths, fill=id,
                                     label=lengths)) +
        geom_bar(stat="identity", color="white") +
        # geom_text(size = 7, color="white", position = position_stack(vjust = 0.5)) +
        # theme(axis.text.x=element_blank()) +
        guides(fill=FALSE) +
        ylab("Cumulative Length (bp)") +
        xlab("") +
        scale_fill_manual(values=c(gg_color_hue(3)[1],
                                            gg_color_hue(3)[2]))+
        scale_y_continuous(breaks=seq(0, max(c(sum(ins_lengths),
                                                        sum(del_lengths))), 5))+
        theme_shared() +
        theme(plot.margin=ggplot2::unit(c(0, 0, 0, 4), "cm")) +
        coord_flip()
}


#' @export
#' @import ggplot2
gg_pie <- function(d) {
    dcounts <- table(as.factor(d))
    df <- as.data.frame(dcounts)
    p <- ggplot(df, ggplot2::aes(x="", y=Freq, fill=Var1))+
        geom_bar(width=0.5, stat = "identity") +
    # p <- p + coord_polar("y", start=0)
        ylab("Count") +
        scale_y_continuous(breaks = seq(0, 150, by = 20)) +
        scale_x_discrete(breaks = NULL) +
        scale_fill_discrete(labels=c("Del", "Ins", "SNV")) +
        guides(fill=guide_legend(title=NULL)) +
        theme_shared() +
        theme(axis.text.x=element_blank(),
                       title = element_blank(),
                       panel.border = element_blank(),
                       legend.position="bottom")
}

#' @export
#' @import ggplot2
theme_shared <- function() {
    theme_bw(18) +
        theme(
            panel.grid.major.x = element_blank()
        )
}

#' @export
#' @import ggplot2
af_histogram <- function(line_info, var_table) {
    ggplot(var_table, ggplot2::aes(af)) +
        geom_histogram(breaks=seq(0, 1, by=0.02), position="dodge") +
        scale_x_continuous(breaks = seq(0, 1, by=0.1)) +
        xlab("Mutant Allele Frequency") +
        theme_shared() +
        ylab("Density")
}
