
#' @export
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

barplot_with_ci <- function(names, data, ci1, ci2) {
    df <- data.frame(names=names, data=data, ci1=ci1, ci2=ci2)
    dodge <- ggplot2::position_dodge(width = 0.9)
    plot <- ggplot2::ggplot(data = df,
                            ggplot2::aes(x = names, y = data, fill = names))
    plot <- plot + ggplot2::geom_bar(stat = "identity", position = dodge)
    plot <- plot + ggplot2::geom_errorbar(position=dodge, width=0.25,
                                 ggplot2::aes(ymin=ci1, ymax=ci2))
    plot
}

#' @export
point_with_ci <- function(df) {
    dodge <- ggplot2::position_dodge(width = 0.9)
    plot <- ggplot2::ggplot(data = df,
                            ggplot2::aes(x = sample, y = data))
    plot <- plot + ggplot2::geom_point(stat = "identity", position = dodge)
    plot <- plot + ggplot2::geom_errorbar(position=dodge, width=0.25,
                                          ggplot2::aes(ymin=ci1, ymax=ci2))
    plot <- plot + theme_shared()
    plot
}


#' @export
extraction_map <- function() {
    require(rworldmap)
    newmap <- rworldmap::getMap(resolution = "low")
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
    ggplot2::ggplot(snv_df, ggplot2::aes(x=id, y=n)) +
        ggplot2::geom_bar(stat="Density",fill=gg_color_hue(3)[3]) +
        ggplot2::ylab("") +
        theme_shared() +
        ggplot2::labs(fill="")
}


# For subsampled analysis
# iso <- unique(line_info$isolate)
# bs <- unlist(lapply(1:length(base_surv), function(i) {
#     isonum <- length(which(line_info$isolate==iso[i]))
#     rep(base_surv[i], isonum)
# }))

#' @export
custom_boxplot <- function(bootstrap_data, xlab, ylab, xval,
                           yscale="unscaled", fill=NULL) {
    # fill <- "#538cc6"
    p <- ggplot2::ggplot(bootstrap_data, ggplot2::aes(group))
    if (is.null(fill)) {
        p <- p + ggplot2::geom_boxplot(
            ggplot2::aes(x=group, ymin=ci1, lower=q25, middle=q50,
                         upper=q75, ymax=ci2), stat="identity", size=0.6)
    } else {
        p <- p + ggplot2::geom_boxplot(ggplot2::aes(x=group, ymin=ci1,
                                                    lower=q25, middle=q50,
                                                    upper=q75, ymax=ci2,
                                                    fill=fill),
                              stat="identity", size=0.6)
        p <- p + ggplot2::scale_fill_discrete("")
    }

    p <- p + ggplot2::scale_x_discrete(labels=xval)
    p <- p + ggplot2::xlab("")
    if (yscale=="log10") {
        p <- p + ggplot2::scale_y_log10()
    }

    p <- p + theme_shared()
    p <- p + ggplot2::labs(x = xlab, y = ylab)
    p
}


#' @export
ins_del_stacked_bar <- function(var_table) {
    indel_var_table <- var_table[which(var_table$class %in%
                                           c("insertion", "deletion")),]
    ninsertions <- length(which(indel_var_table$class == "insertion"))
    ndeletions <- length(which(indel_var_table$class == "deletion"))
    ins_lengths <- nchar(dplyr::filter(indel_var_table, class=="insertion")$alt) - 1
    del_lengths <- nchar(dplyr::filter(indel_var_table, class=="deletion")$ref) - 1
    id <- c(rep("Deletion", length(del_lengths)),
            rep("Insertion", length(ins_lengths)))
    lengths <- c(sort(del_lengths, decreasing=TRUE),
                 sort(ins_lengths, decreasing=TRUE))
    df <- data.frame(id=factor(id), lengths=lengths)
    ggplot2::ggplot(df, ggplot2::aes(x=id, y=lengths, fill=id,
                                     label=lengths)) +
        ggplot2::geom_bar(stat="identity", color="white") +
        # geom_text(size = 7, color="white", position = position_stack(vjust = 0.5)) +
        # theme(axis.text.x=element_blank()) +
        ggplot2::guides(fill=FALSE) +
        ggplot2::ylab("Cumulative Length (bp)") +
        ggplot2::xlab("") +
        ggplot2::scale_fill_manual(values=c(gg_color_hue(3)[1],
                                            gg_color_hue(3)[2]))+
        ggplot2::scale_y_continuous(breaks=seq(0, max(c(sum(ins_lengths),
                                                        sum(del_lengths))), 5))+
        theme_shared() +
        ggplot2::theme(plot.margin=ggplot2::unit(c(0, 0, 0, 4), "cm")) +
        ggplot2::coord_flip()
}



#' @export
gg_pie <- function(d) {
    dcounts <- table(as.factor(d))
    df <- as.data.frame(dcounts)
    p <- ggplot2::ggplot(df, ggplot2::aes(x="", y=Freq, fill=Var1))+
        ggplot2::geom_bar(width=0.5, stat = "identity") +
    # p <- p + coord_polar("y", start=0)
        ggplot2::ylab("Count") +
        ggplot2::scale_y_continuous(breaks = seq(0, 150, by = 20)) +
        ggplot2::scale_x_discrete(breaks = NULL) +
        ggplot2::scale_fill_discrete(labels=c("Del", "Ins", "SNV")) +
        ggplot2::guides(fill=ggplot2::guide_legend(title=NULL)) +
        theme_shared() +
        ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                       title = ggplot2::element_blank(),
                       panel.border = ggplot2::element_blank(),
                       legend.position="bottom")
}

#' @export
theme_shared <- function() {
    ggplot2::theme_bw(18) +
        ggplot2::theme(
            panel.grid.major.x = ggplot2::element_blank()
        )
}


#' @export
lollipop_plots <- function(var_table, mag_genbank, pul_genbank) {
    require(trackViewer)
    require(genoPlotR)
    for (gen in c(mag_genbank, pul_genbank)) {
        genes <- read_dna_seg_from_genbank(gen)
        if (gen==mag_genbank) {
            species <- "magna"
        } else {
            species <- "pulex"
        }
        vart <- filter(var_table, species==species)
        snp <- vart$pos
        gr <- GenomicRanges::GRanges(
            "mtDNA", IRanges::IRanges(snp, width=1), height=vart$af)
        feat <- GRanges::GRanges("mtDNA", devtools::IRanges(genes$start,
                                         width=genes$end-genes$start,
                                         names=genes$gene))


        p <- ggplot2::ggplot(vart) +
            ggplot2::geom_bar(data=vart, mapping=ggplot2::aes(x=pos, y=af),
                     stat="identity", width=20) +
            ggplot2::geom_rect(data=as.data.frame(genes),
                               mapping=ggplot2::aes(xmin=start,
                      xmax=end,
                      ymin=-0.05, ymax=0), inherit.aes=F, fill=gg_color_hue(1)[1])
        eval(p)
        print(p)

    }


}


#' @export
af_histogram <- function(line_info, var_table) {
    ggplot2::ggplot(var_table, ggplot2::aes(af)) +
        ggplot2::geom_histogram(breaks=seq(0, 1, by=0.02), position="dodge") +
        ggplot2::scale_x_continuous(breaks = seq(0, 1, by=0.1)) +
        ggplot2::xlab("Mutant Allele Frequency") +
        theme_shared() +
        ggplot2::ylab("Density")
}
