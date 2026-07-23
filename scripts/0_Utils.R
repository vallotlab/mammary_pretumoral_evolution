# --------------------------
### Libraries ----
# --------------------------

# General

library(dplyr)
library(tidyr)
library(cluster) # for silhouette metric
library(scCustomize)

library(Seurat)
library(SingleCellExperiment)

require(xlsx)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(wesanderson)
library(pheatmap)
library(MatrixGenerics)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)

# survival

library(dplyr)
library(survival)
library(survminer)
library(ggfortify)
library(ggplot2)

# GLM

library(purrr)
library(grr)
library(lmtest)

# epigenome

library(BSgenome.Hsapiens.UCSC.hg38)
library(Matrix)
library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(eulerr)
library(ChromSCape)
library(scCustomize)
library(Signac)
library(ggplot2)
library(patchwork)
library(GenomicRanges)
library(viridis)
library(stringr)
library(ArchR)
library(rGREAT)



# pathway enrichment
library(rGREAT)

# motif search
library(GenomicRanges)
library(magrittr)
library(universalmotif)
library(memes)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVARmotifs)
library(TFBSTools)


## bigwig

library(plotgardener)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


# --------------------------
### Directories ----
# --------------------------

# tools
tool_dir_macs2 <- ""

# general
output_dir_clinical_data <- "./data/clinical_data/"

# snH3K4me1
input_dir_raw <- "./data/snH3K4me1/raw/"
input_dir_metadata <- "./data/snH3K4me1/metadata/"
tmp_output_dir_sce <- "./data/snH3K4me1/tmp_filtered_sce_50k/"
output_dir_sce <- "./data/snH3K4me1/sce_50k/"
output_dir_gene_activity <- "./data/snH3K4me1/gene_activity/"
output_dir_frags <- "./data/snH3K4me1/filtered_fragments/"
output_dir_doublets <- "./results/snH3K4me1/doublets/"
output_dir_annot <- "./data/snH3K4me1/annotation_lineage/"
output_dir_bw <- "./data/snH3K4me1/bigWig_by_lineage/"
output_dir_peaks <- "./data/snH3K4me1/peaks_10k/"
output_dir_DA_10k <- "./data/snH3K4me1/DA_10k_peaks/"
output_dir_peaks_TME <- "./data/snH3K4me1/peaks_10k/TME/"
output_dir_peaks_tumor <- "./data/snH3K4me1/peaks_10k/Tumor/"
output_dir_TFs <- "./data/snH3K4me1/TFs_500bp/"
output_dir_peaks_500bp <- "./data/snH3K4me1/peaks_500bp/"
output_dir_snH3K4me1_TME <- "./data/snH3K4me1/TME/"
output_dir_snH3K4me1_TME_bed <- "./data/snH3K4me1/TME/bed/"
output_dir_snH3K4me1_TME_bed_NT_RD <- "./data/snH3K4me1/TME/bed_NT_RD/"
output_dir_snH3K4me1_TME_DA <- "./data/snH3K4me1/TME/DA/"

output_dir_snH3K4me1_bed_states <- "./data/snH3K4me1/bed_states/"
output_dir_snH3K4me1_states <- "./data/snH3K4me1/states/"
output_dir_snH3K4me1_peaks_states <- "./data/snH3K4me1/peaks_states/"
output_dir_snH3K4me1_TF_states <- "./data/snH3K4me1/TF_states/"
output_dir_snH3K4me1_DA_states <- "./data/snH3K4me1/DA_states/"


# snRNA
output_dir_snRNA <- "./data/snRNA/seurat_all_cells/"
output_dir_snRNA_states <- "./data/snRNA/states/"
output_dir_snRNA_umap <- "./data/snRNA/UMAP/"
output_dir_snRNA_TF <- "./data/snRNA/TF/"
output_dir_snRNA_TME <- "./data/snRNA/TME/"

# snRNA/snH3K4me1
output_dir_CCA <- "./data/integration/CCA_snRNA_H3K4me1/"

# scRNA
output_dir_scRNA <- "./data/scRNA/count_matrices/"
output_dir_scRNA_NMF <- "./data/scRNA/NMF/"
output_dir_scRNA_states <- "./data/scRNA/states/"
output_dir_scRNA_vis <- "./data/scRNA/visualization/"
output_dir_scRNA_TF <- "./data/scRNA/TF/"
output_dir_scRNA_CCC <- "./data/scRNA/CCC/"



# --------------------------
### Plot functions ----
# --------------------------


save_plot <- function(p, file, res = 300, width = 7, height = 7) {
    png(file, res = res, width = width, height = height, units = "cm")
    print(p)
    dev.off()
    paste0("Plot saved to: ", file)
}

generate_palette <- function(length) {
    set.seed(6)
    palette <- unique(c(
        wesanderson::wes_palette("Rushmore1")[3],
        wesanderson::wes_palette("Rushmore1")[5],
        wesanderson::wes_palette("Zissou1")[1],
        wesanderson::wes_palette("Zissou1")[3],
        wesanderson::wes_palette("Darjeeling1")[4],
        wesanderson::wes_palette("GrandBudapest2")[1],
        sample(c(wesanderson::wes_palette("Darjeeling1"), wesanderson::wes_palette("Darjeeling2")[1:4], wesanderson::wes_palette("Moonrise3")[1:3], wesanderson::wes_palette("GrandBudapest1"), wesanderson::wes_palette("GrandBudapest2")))
    ))
    if (length <= length(palette)) {
        palette <- palette[1:length]
    } else {
        palette <- rep(palette, times = round(length / length(palette)) + 1)
        palette <- palette[1:length]
    }
    return(palette)
}



plot_barplot_2_identities_df_cell_count <- function(df, annot_name_prop, annot_name_by, title = "Proportion plot", palette = NULL, order_annot_by = NULL, order_annot_prop = NULL) {
    prop <- lapply(unique(df[, annot_name_by]), function(x) prop.table(table(df[, annot_name_prop][which(df[, annot_name_by] == x)], useNA = "ifany")))
    names(prop) <- unique(df[, annot_name_by])
    Sample <- c()
    Cell_type <- c()
    Percentage <- c()
    for (sample in names(prop)) {
        Sample <- c(Sample, rep(sample, times = length(prop[[sample]])))
        Cell_type <- c(Cell_type, names(prop[[sample]]))
        Percentage <- c(Percentage, prop[[sample]])
    }
    data <- data.frame(
        Sample = Sample,
        Group = Cell_type,
        Percentage = Percentage * 100
    )

    # add number of cells by annotation
    df$annot_name_by <- df[, annot_name_by]
    tmp <- df %>% dplyr::count(annot_name_by)
    rownames(tmp) <- tmp$annot_name_by
    data$Nb_cells <- tmp[data$Sample, "n"]

    if (!is.null(order_annot_by)) {
        data$Sample <- factor(data$Sample, levels = rev(intersect(order_annot_by, unique(data$Sample))))
    }
    if (!is.null(order_annot_prop)) {
        data$Group <- factor(data$Group, levels = rev(intersect(order_annot_prop, unique(data$Group))))
    }

    p <- ggplot(data = data, aes(x = Sample, y = Percentage, fill = Group)) +
        geom_bar(stat = "identity", width = 0.7) +
        theme_minimal() +
        coord_flip() +
        xlab("") +
        ggtitle(title) +
        theme(plot.title = element_text(hjust = 0.5)) +
        geom_text(aes(y = 101, label = Nb_cells), hjust = -0.25) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        scale_y_continuous(breaks = c(0, 25, 50, 75, 100), limits = c(c(0, 109)))

    if (is.null(palette)) {
        p <- p + scale_fill_manual(breaks = unique(seurat@meta.data[, annot_name_prop]), values = generate_palette(length(unique(seurat@meta.data[, annot_name_prop]))), na.value = "lightgrey")
    } else {
        p <- p + scale_fill_manual(breaks = names(palette), values = palette, na.value = "lightgrey")
    }
    return(p)
}



plot_umap_annot <- function(df_umap, type = "Cell type", labels, title = NULL, colors = NULL, show_legend = T, pt.size = 2, pt.stroke = 0.5, pt.transp = 0.1) {
    #### Build plot ####
    df <- data.frame(
        x = df_umap[, 1],
        y = df_umap[, 2],
        Alias = labels
    )
    if (show_legend) {
        p <- ggplot(df[sample(1:nrow(df), size = nrow(df), replace = F), ], aes(x, y, colour = Alias)) +
            geom_point(shape = 21, colour = "white", aes(fill = Alias), size = pt.size, stroke = pt.stroke, alpha = pt.transp) +
            xlab("UMAP_1") +
            ylab("UMAP_2") +
            ggtitle(title) +
            theme_classic() +
            theme(plot.title = element_text(hjust = 0.5)) +
            scale_fill_manual(values = colors, breaks = unique(labels)) +
            labs(fill = type)
    } else {
        p <- ggplot(df[sample(1:nrow(df), size = nrow(df), replace = F), ], aes(x, y, colour = Alias)) +
            geom_point(shape = 21, colour = "white", aes(fill = Alias), size = pt.size, stroke = pt.stroke, alpha = pt.transp) +
            xlab("UMAP_1") +
            ylab("UMAP_2") +
            ggtitle(title) +
            theme_classic() +
            theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
            scale_fill_manual(values = colors, breaks = unique(labels)) +
            labs(fill = type)
    }

    print(p)
}



plot_umap_feature_continuous <- function(df_umap, type = "Cell type", labels, title = NULL, show_legend = T, pt.size = 2, pt.stroke = 0.5, pt.transp = 0.1, high_col = "#625793", low_col = "grey95") {
    #### Build plot ####
    df <- data.frame(
        x = df_umap[, 1],
        y = df_umap[, 2],
        Alias = labels
    )
    if (show_legend) {
        p <- ggplot(df[order(labels), ], aes(x, y, colour = Alias)) +
            geom_point(shape = 21, colour = "white", aes(fill = Alias), size = pt.size, stroke = pt.stroke, alpha = pt.transp) +
            xlab("UMAP_1") +
            ylab("UMAP_2") +
            ggtitle(title) +
            theme_classic() +
            theme(plot.title = element_text(hjust = 0.5)) +
            scale_fill_gradient(high = high_col, low = low_col) +
            labs(fill = type)
    } else {
        p <- ggplot(df[order(labels), ], aes(x, y, colour = Alias)) +
            geom_point(shape = 21, colour = "white", aes(fill = Alias), size = pt.size, stroke = pt.stroke, alpha = pt.transp) +
            xlab("UMAP_1") +
            ylab("UMAP_2") +
            ggtitle(title) +
            theme_classic() +
            theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
            scale_fill_gradient(high = high_col, low = low_col) +
            labs(fill = type)
    }

    print(p)
}


ggradar_custom <- function(plot.data, base.size = 15, font.radar = "sans", values.radar = c(
                               "0%",
                               "50%", "100%"
                           ), axis.labels = colnames(plot.data)[-1], grid.min = 0,
                           grid.mid = 0.5, grid.max = 1, centre.y = grid.min - ((1 / 9) *
                               (grid.max - grid.min)), plot.extent.x.sf = 1, plot.extent.y.sf = 1.2,
                           x.centre.range = 0.02 * (grid.max - centre.y), label.centre.y = FALSE,
                           grid.line.width = 0.5, gridline.min.linetype = "longdash",
                           gridline.mid.linetype = "longdash", gridline.max.linetype = "longdash",
                           gridline.min.colour = "white", gridline.mid.colour = "white",
                           gridline.max.colour = "white", grid.label.size = 6, gridline.label.offset = -0.1 *
                               (grid.max - centre.y), label.gridline.min = FALSE, label.gridline.mid = FALSE,
                           label.gridline.max = TRUE, axis.label.offset = 1.15, axis.label.size = 5,
                           axis.line.colour = "white", group.line.width = 1.5, group.point.size = 6,
                           group.colours = NULL, background.circle.colour = "#D7D6D1",
                           background.circle.transparency = 0.2, plot.legend = if (nrow(plot.data) >
                               1) {
                               TRUE
                           } else {
                               FALSE
                           }, legend.title = "", plot.title = "",
                           legend.text.size = 14, legend.position = "left", fill = FALSE,
                           fill.alpha = 0.5, draw.points = TRUE, point.alpha = 1, line.alpha = 1) {
    plot.data <- as.data.frame(plot.data)
    plot.data <- aggregate(x = plot.data[, -1], by = list(plot.data[
        ,
        1
    ]), FUN = "mean")
    if (!is.factor(plot.data[, 1])) {
        plot.data[, 1] <- as.factor(as.character(plot.data[
            ,
            1
        ]))
    }
    var.names <- colnames(plot.data)[-1]
    plot.extent.x <- (grid.max + abs(centre.y)) * plot.extent.x.sf
    plot.extent.y <- (grid.max + abs(centre.y)) * plot.extent.y.sf
    if (length(axis.labels) != ncol(plot.data) - 1) {
        stop("'axis.labels' contains the wrong number of axis labels",
            call. = FALSE
        )
    }
    if (min(plot.data[, -1]) < centre.y) {
        stop("plot.data' contains value(s) < centre.y", call. = FALSE)
    }
    if (max(plot.data[, -1]) > grid.max) {
        plot.data[, -1] <- (plot.data[, -1] / max(plot.data[
            ,
            -1
        ])) * grid.max
        warning("'plot.data' contains value(s) > grid.max, data scaled to grid.max",
            call. = FALSE
        )
    }
    plot.data.offset <- plot.data
    plot.data.offset[, 2:ncol(plot.data)] <- plot.data[, 2:ncol(plot.data)] +
        abs(centre.y)
    group <- NULL
    group$path <- CalculateGroupPath(plot.data.offset)
    axis <- NULL
    axis$path <- CalculateAxisPath(
        var.names, grid.min + abs(centre.y),
        grid.max + abs(centre.y)
    )
    axis$label <- data.frame(text = axis.labels, x = NA, y = NA)
    n.vars <- length(var.names)
    angles <- seq(from = 0, to = 2 * pi, by = (2 * pi) / n.vars)
    axis$label$x <- sapply(1:n.vars, function(i, x) {
        ((grid.max + abs(centre.y)) * axis.label.offset) * sin(angles[i])
    })
    axis$label$y <- sapply(1:n.vars, function(i, x) {
        ((grid.max + abs(centre.y)) * axis.label.offset) * cos(angles[i])
    })
    gridline <- NULL
    gridline$min$path <- funcCircleCoords(c(0, 0), grid.min +
        abs(centre.y), npoints = 360)
    gridline$mid$path <- funcCircleCoords(c(0, 0), grid.mid +
        abs(centre.y), npoints = 360)
    gridline$max$path <- funcCircleCoords(c(0, 0), grid.max +
        abs(centre.y), npoints = 360)
    gridline$min$label <- data.frame(
        x = gridline.label.offset,
        y = grid.min + abs(centre.y), text = as.character(grid.min)
    )
    gridline$max$label <- data.frame(
        x = gridline.label.offset,
        y = grid.max + abs(centre.y), text = as.character(grid.max)
    )
    gridline$mid$label <- data.frame(
        x = gridline.label.offset,
        y = grid.mid + abs(centre.y), text = as.character(grid.mid)
    )
    theme_clear <- theme_bw(base_size = base.size) + theme(
        axis.text.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), legend.key = element_rect(linetype = "blank")
    )
    if (plot.legend == FALSE) {
        legend.position <- "none"
    }
    base <- ggplot(axis$label) +
        xlab(NULL) +
        ylab(NULL) +
        coord_equal() +
        geom_text(
            data = subset(axis$label, axis$label$x < (-x.centre.range)),
            aes(x = x, y = y, label = text), size = axis.label.size,
            hjust = 1, family = font.radar
        ) +
        scale_x_continuous(limits = c(-1.5 *
            plot.extent.x, 1.5 * plot.extent.x)) +
        scale_y_continuous(limits = c(
            -plot.extent.y,
            plot.extent.y
        ))
    base <- base + geom_path(
        data = gridline$min$path, aes(
            x = x,
            y = y
        ), lty = gridline.min.linetype, colour = gridline.min.colour,
        linewidth = grid.line.width
    )
    # base <- base + geom_path(data = gridline$mid$path, aes(x = x,
    #                                                         y = y), lty = gridline.mid.linetype, colour = gridline.mid.colour,
    #                             linewidth = grid.line.width)
    base <- base + geom_path(
        data = gridline$max$path, aes(
            x = x,
            y = y
        ), lty = gridline.max.linetype, colour = gridline.max.colour,
        linewidth = grid.line.width
    )
    base <- base + geom_text(
        data = subset(axis$label, abs(axis$label$x) <=
            x.centre.range), aes(x = x, y = y, label = text), size = axis.label.size,
        hjust = 0.5, family = font.radar
    )
    base <- base + geom_text(
        data = subset(axis$label, axis$label$x >
            x.centre.range), aes(x = x, y = y, label = text), size = axis.label.size,
        hjust = 0, family = font.radar
    )
    base <- base + theme_clear
    base <- base + geom_polygon(data = gridline$max$path, aes(
        x,
        y
    ), fill = background.circle.colour, alpha = background.circle.transparency)
    base <- base + geom_path(data = axis$path, aes(
        x = x, y = y,
        group = axis.no
    ), colour = axis.line.colour)
    theGroupName <- names(group$path[1])
    if (length(line.alpha) == 1) {
        base <- base + geom_path(
            data = group$path, aes(
                x = .data[["x"]],
                y = .data[["y"]], group = .data[[theGroupName]],
                colour = .data[[theGroupName]]
            ), linewidth = group.line.width,
            alpha = line.alpha
        )
    } else {
        base <- base + geom_path(data = group$path, aes(
            x = .data[["x"]],
            y = .data[["y"]], group = .data[[theGroupName]],
            colour = .data[[theGroupName]]
        ), linewidth = group.line.width) +
            scale_alpha_manual(values = line.alpha)
    }
    if (draw.points) {
        if (length(point.alpha) == 1) {
            base <- base + geom_point(
                data = group$path, aes(
                    x = .data[["x"]],
                    y = .data[["y"]], group = .data[[theGroupName]],
                    colour = .data[[theGroupName]]
                ), size = group.point.size,
                alpha = point.alpha
            )
        } else {
            base <- base + geom_point(data = group$path, aes(
                x = .data[["x"]],
                y = .data[["y"]], group = .data[[theGroupName]],
                colour = .data[[theGroupName]]
            ), size = group.point.size) +
                scale_alpha_manual(values = point.alpha)
        }
    }
    if (fill == TRUE) {
        base <- base + geom_polygon(data = group$path, aes(
            x = .data[["x"]],
            y = .data[["y"]], group = .data[[theGroupName]],
            fill = .data[[theGroupName]]
        ), alpha = fill.alpha)
    }
    if (plot.legend == TRUE) {
        base <- base + labs(colour = legend.title, size = legend.text.size)
    }
    if (label.gridline.min == TRUE) {
        base <- base + geom_text(aes(x = x, y = y, label = values.radar[1]),
            data = gridline$min$label, size = grid.label.size *
                0.8, hjust = 1, family = font.radar
        )
    }
    if (label.gridline.mid == TRUE) {
        base <- base + geom_text(aes(x = x, y = y, label = values.radar[2]),
            data = gridline$mid$label, size = grid.label.size *
                0.8, hjust = 1, family = font.radar
        )
    }
    if (label.gridline.max == TRUE) {
        base <- base + geom_text(aes(x = x, y = y, label = values.radar[3]),
            data = gridline$max$label, size = grid.label.size *
                0.8, hjust = 1, family = font.radar
        )
    }
    if (label.centre.y == TRUE) {
        centre.y.label <- data.frame(x = 0, y = 0, text = as.character(centre.y))
        base <- base + geom_text(aes(x = x, y = y, label = text),
            data = centre.y.label, size = grid.label.size, hjust = 0.5,
            family = font.radar
        )
    }
    # if (!is.null(group.colours)) {
    #     colour_values <- rep(group.colours, length(unique(plot.data[,
    #                                                                 1]))/length(group.colours))
    # }
    # else {
    #     colour_values <- generate_color_values(length(unique(plot.data[,
    #                                                                    1])))
    # }
    base <- base + theme(
        legend.key.width = unit(3, "line"),
        text = element_text(size = 20, family = font.radar)
    ) +
        theme(
            legend.text = element_text(size = legend.text.size),
            legend.position = legend.position
        ) + theme(legend.key.height = unit(
            2,
            "line"
        )) + scale_colour_manual(values = group.colours, breaks = names(group.colours)) +
        theme(text = element_text(family = font.radar)) + theme(legend.title = element_blank())
    if (isTRUE(fill)) {
        base <- base + scale_fill_manual(
            values = group.colours, breaks = names(group.colours),
            guide = "none"
        )
    }
    if (legend.title != "") {
        base <- base + theme(legend.title = element_text())
    }
    if (plot.title != "") {
        base <- base + ggtitle(plot.title)
    }
    return(base)
}


CalculateAxisPath <- function(var.names, min, max) {
    # var.names <- c("v1","v2","v3","v4","v5")
    n.vars <- length(var.names) # number of vars (axes) required
    # Cacluate required number of angles (in radians)
    angles <- seq(from = 0, to = 2 * pi, by = (2 * pi) / n.vars)
    # calculate vectors of min and max x+y coords
    min.x <- min * sin(angles)
    min.y <- min * cos(angles)
    max.x <- max * sin(angles)
    max.y <- max * cos(angles)
    # Combine into a set of uniquely numbered paths (one per variable)
    axisData <- NULL
    for (i in 1:n.vars) {
        a <- c(i, min.x[i], min.y[i])
        b <- c(i, max.x[i], max.y[i])
        axisData <- rbind(axisData, a, b)
    }
    # Add column names + set row names = row no. to allow conversion into a data frame
    colnames(axisData) <- c("axis.no", "x", "y")
    rownames(axisData) <- seq(1:nrow(axisData))
    # Return calculated axis paths
    as.data.frame(axisData)
}

CalculateGroupPath <- function(df) {
    # Drop dead levels. This might happen if the data is filtered on the way
    # into ggradar.
    path <- forcats::fct_drop(df[, 1])
    # set the name of the variable that will be used for grouping
    theGroupName <- colnames(df)[1]

    ## find increment
    nPathPoints <- ncol(df) - 1
    angles <- seq(from = 0, to = 2 * pi, by = (2 * pi) / nPathPoints)
    ## create graph data frame
    nDataPoints <- ncol(df) * length(levels(path))
    graphData <- data.frame(
        seg = rep("", nDataPoints),
        x = rep(0, nDataPoints),
        y = rep(0, nDataPoints)
    )
    colnames(graphData)[1] <- theGroupName

    rowNum <- 1
    for (i in 1:length(levels(path))) {
        pathData <- subset(df, df[, 1] == levels(path)[i])
        for (j in c(2:ncol(df))) {
            graphData[rowNum, theGroupName] <- levels(path)[i]
            graphData$x[rowNum] <- pathData[, j] * sin(angles[j - 1])
            graphData$y[rowNum] <- pathData[, j] * cos(angles[j - 1])
            rowNum <- rowNum + 1
        }
        ## complete the path by repeating first pair of coords in the path
        graphData[rowNum, theGroupName] <- levels(path)[i]
        graphData$x[rowNum] <- pathData[, 2] * sin(angles[1])
        graphData$y[rowNum] <- pathData[, 2] * cos(angles[1])
        rowNum <- rowNum + 1
    }
    # Make sure that name of first column matches that of input data (in case !="group")
    graphData[, 1] <- factor(graphData[, 1], levels = levels(path)) # keep group order
    graphData # data frame returned by function
}

funcCircleCoords <- function(center = c(0, 0), r = 1, npoints = 100) {
    tt <- seq(0, 2 * pi, length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
}




plot_dotplot_cell_type_markers_ggplot <- function(expr_mat, annot, cell_type_markers, output_path = NULL, width = 10, height = 8, res = 300, annot_order = NULL, aggregate_markers_across_cell_types = F) {
    library(ggplot2)
    library(dplyr)
    library(tidyr)

    # Keep only genes present in the expression matrix, in the order of the list
    gene_order <- unique(unlist(cell_type_markers))
    gene_order <- gene_order[gene_order %in% rownames(expr_mat)]
    markers <- gene_order
    expr_sub <- expr_mat[rownames(expr_mat) %in% markers, , drop = FALSE]
    df <- as.data.frame(t(expr_sub))
    df$cell_type <- annot[rownames(df)]
    df_long <- df %>%
        pivot_longer(-cell_type, names_to = "gene", values_to = "expression")

    # Calculate percent expressing and average expression per cell type/gene
    plot_data <- df_long %>%
        dplyr::group_by(cell_type, gene) %>%
        dplyr::summarise(
            avg_expr = mean(expression, na.rm = TRUE),
            pct_expr = mean(expression > 0, na.rm = TRUE) * 100,
            .groups = "drop"
        )

    # Keep gene order and reverse cluster order
    plot_data$gene <- factor(plot_data$gene, levels = gene_order)
    if (is.null(annot_order)) {
        annot_order <- sort(unique(plot_data$cell_type))
    }
    plot_data$cell_type <- factor(plot_data$cell_type, levels = rev(annot_order))

    # Create a data frame mapping each gene to its associated cell type(s), only for genes present
    gene2celltype <- stack(cell_type_markers)
    colnames(gene2celltype) <- c("gene", "associated_cell_type")
    gene2celltype <- gene2celltype %>%
        dplyr::filter(gene %in% gene_order) %>%
        dplyr::group_by(gene) %>%
        dplyr::summarise(associated_cell_type = paste(unique(associated_cell_type), collapse = ", "))
    gene2celltype$gene <- factor(gene2celltype$gene, levels = gene_order)
    y_annot <- length(annot_order) + 0.5
    gene2celltype$y_annot <- y_annot

    # Plot
    p <- ggplot(plot_data, aes(x = gene, y = cell_type)) +
        geom_point(aes(size = pct_expr, color = avg_expr)) +
        # scale_color_gradientn(colors = c("white", "lightblue", "darkblue")) +
        scale_color_gradientn(colors = c("#ffffff", colorRampPalette(brewer.pal(9, "Purples"))(100)[15:100])) +
        scale_size(range = c(1, 8)) +
        theme_classic() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.margin = margin(60, 20, 40, 20)
        ) +
        labs(
            x = "Marker Genes",
            y = "",
            color = "Avg Expr",
            size = "% Expressing"
        ) +
        geom_text(
            data = gene2celltype,
            aes(x = gene, y = y_annot, label = associated_cell_type),
            inherit.aes = FALSE,
            angle = 45,
            vjust = 0,
            hjust = 0,
            size = 3
        ) +
        coord_cartesian(clip = "off")

    # Save or return
    if (!is.null(output_path)) {
        ggsave(filename = output_path, plot = p, width = width, height = height, dpi = res, units = "cm")
    } else {
        return(p)
    }
}


# --------------------------
### Scoring functions ----
# --------------------------

UCell_scores_per_cell <- function(raw_counts,
                                  sig) {
    # Keep genes in dataset
    sig_genes <- dplyr::intersect(unique(unlist(sig)), rownames(raw_counts)) # restrict to genes found in matrix
    sig <- lapply(sig, function(x) dplyr::intersect(x, sig_genes))
    # Compute Ucell scores
    scores_ucell <- as.data.frame(UCell::ScoreSignatures_UCell(raw_counts, features = sig))
    return(scores_ucell)
}


# --------------------------
### Survival functions ----
# --------------------------

survival_vs_celltype_proportions <- function(cell_annot,
                                             sample_annot,
                                             thresh = NULL,
                                             col_ID,
                                             col_annot,
                                             col_time = c("OS", "PFS"),
                                             col_impact = NULL,
                                             keep_annot) {
    ### Remove samples with too few cells
    nb_cells_by_sample <- table(cell_annot[, col_ID])
    keep_samples <- names(nb_cells_by_sample)[which(nb_cells_by_sample >= 20)]
    cell_annot <- cell_annot[which(cell_annot[, col_ID] %in% keep_samples), ]

    res <- list() # prepare list to contain results

    ### Compute cell annotation % by sample
    cell_annot[, col_annot] <- factor(cell_annot[, col_annot], levels = unique(cell_annot[, col_annot]))
    cell_annot_by_sample <- cell_annot %>%
        dplyr::group_by_at(c(col_ID, col_impact, col_annot), .drop = F) %>%
        dplyr::summarise(n = n()) %>%
        dplyr::mutate(prop_cells = n / sum(n)) %>%
        as.data.frame()
    cell_annot_by_sample <- cell_annot_by_sample[which(cell_annot_by_sample[, col_annot] == keep_annot), ]

    ### Classify % as high/low
    if (!is.null(thresh)) {
        cell_annot_by_sample$cat <- ifelse(cell_annot_by_sample$prop_cells >= thresh, "HIGH", "LOW")
    } else {
        cell_annot_by_sample$cat <- ifelse(cell_annot_by_sample$prop_cells >= median(cell_annot_by_sample$prop_cells), "HIGH", "LOW")
    }

    ### If all samples are in the same HIGH/LOW group: consider pval=1
    if (length(unique(cell_annot_by_sample$cat)) == 1) {
        pval <- 1
        p <- NULL
    } else {
        ### Add cell annotation % to sample information
        cell_annot_by_sample$time <- sample_annot[match(cell_annot_by_sample[, col_ID], sample_annot[, col_ID]), col_time]
        cell_annot_by_sample$status <- sample_annot[match(cell_annot_by_sample[, col_ID], sample_annot[, col_ID]), "status"]
        cell_annot_by_sample <- cell_annot_by_sample[which(!is.na(cell_annot_by_sample$time)), ] # remove unannotated samples

        ### Run survival analysis
        res_survival <- run_survival_logrank(
            data = cell_annot_by_sample,
            col_impact
        )
        pval <- res_survival$pvalue

        ### Plot curves
        surv_fit <- survival::survfit(survival::Surv(time, status) ~ cat, cell_annot_by_sample)
        require(ggfortify)
        p <- ggplot2::autoplot(surv_fit, size = 2, conf.int = F) +
            ggplot2::theme_classic() +
            ggplot2::ggtitle(paste0("Impact of ", keep_annot, " % on survival\npval=", format(pval, scientific = T, digits = 3))) +
            ggplot2::xlab(col_time) +
            ggplot2::ylab("Survival") +
            ggplot2::scale_color_manual(values = c("LOW" = "grey", "HIGH" = "#D357A1")) +
            ggplot2::theme(
                axis.text = element_text(size = 14), axis.title = element_text(size = 14), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                plot.title = element_text(size = 14),
                legend.title = element_text(size = 10), legend.text = element_text(size = 10), legend.key.size = unit(0.8, "cm")
            )
    }

    res <- list("pval" = pval, "plot" = p, table = cell_annot_by_sample)
    return(res)
}


survival_vs_celltype_proportions_ratio <- function(cell_annot,
                                                   sample_annot,
                                                   thresh = NULL,
                                                   col_ID,
                                                   col_annot_top,
                                                   col_annot_bottom,
                                                   col_time = c("OS", "PFS"),
                                                   col_impact = NULL,
                                                   keep_annot_top,
                                                   keep_annot_bottom) {
    ### Remove samples with too few cells
    nb_cells_by_sample <- table(cell_annot[, col_ID])
    keep_samples <- names(nb_cells_by_sample)[which(nb_cells_by_sample >= 20)]
    cell_annot <- cell_annot[which(cell_annot[, col_ID] %in% keep_samples), ]

    res <- list() # prepare list to contain results

    ### Compute cell annotation % by sample for top and bottom annotations
    cell_annot_tmp <- cell_annot
    cell_annot_tmp[, col_annot_top] <- factor(cell_annot_tmp[, col_annot_top], levels = unique(cell_annot_tmp[, col_annot_top]))
    cell_annot_by_sample <- cell_annot_tmp %>%
        dplyr::group_by_at(c(col_ID, col_impact, col_annot_top), .drop = F) %>%
        dplyr::summarise(n = n()) %>%
        dplyr::mutate(prop_cells = n / sum(n)) %>%
        as.data.frame()
    cell_annot_by_sample_top <- cell_annot_by_sample[which(cell_annot_by_sample[, col_annot_top] == keep_annot_top), ]

    cell_annot_tmp <- cell_annot
    cell_annot_tmp[, col_annot_bottom] <- factor(cell_annot_tmp[, col_annot_bottom], levels = unique(cell_annot_tmp[, col_annot_bottom]))
    cell_annot_by_sample <- cell_annot_tmp %>%
        dplyr::group_by_at(c(col_ID, col_impact, col_annot_bottom), .drop = F) %>%
        dplyr::summarise(n = n()) %>%
        dplyr::mutate(prop_cells = n / sum(n)) %>%
        as.data.frame()
    cell_annot_by_sample_bottom <- cell_annot_by_sample[which(cell_annot_by_sample[, col_annot_bottom] == keep_annot_bottom), ]

    ### Combine top and bottom annotations
    cell_annot_by_sample <- merge(cell_annot_by_sample_top, cell_annot_by_sample_bottom, by = col_ID, suffixes = c("_top", "_bottom"))

    ### compute ratio
    cell_annot_by_sample$logFC <- log2(cell_annot_by_sample$prop_cells_top / cell_annot_by_sample$prop_cells_bottom)

    ### Classify % as high/low
    if (!is.null(thresh)) {
        cell_annot_by_sample$cat <- ifelse(cell_annot_by_sample$logFC >= thresh, "HIGH", "LOW")
    } else {
        cell_annot_by_sample$cat <- ifelse(cell_annot_by_sample$logFC >= median(cell_annot_by_sample$logFC), "HIGH", "LOW")
    }

    ### If all samples are in the same HIGH/LOW group: consider pval=1
    if (length(unique(cell_annot_by_sample$cat)) == 1) {
        pval <- 1
        p <- NULL
    } else {
        ### Add cell annotation % to sample information
        cell_annot_by_sample$time <- sample_annot[match(cell_annot_by_sample[, col_ID], sample_annot[, col_ID]), col_time]
        cell_annot_by_sample$status <- sample_annot[match(cell_annot_by_sample[, col_ID], sample_annot[, col_ID]), "status"]
        cell_annot_by_sample <- cell_annot_by_sample[which(!is.na(cell_annot_by_sample$time)), ] # remove unannotated samples

        ### Run survival analysis
        res_survival <- run_survival_logrank(
            data = cell_annot_by_sample,
            col_impact
        )
        pval <- res_survival$pvalue

        ### Plot curves
        surv_fit <- survival::survfit(survival::Surv(time, status) ~ cat, cell_annot_by_sample)
        require(ggfortify)
        p <- ggplot2::autoplot(surv_fit, size = 2, conf.int = F) +
            ggplot2::theme_classic() +
            ggplot2::ggtitle(paste0("Impact of ", keep_annot_top, "/", keep_annot_bottom, " % on survival\npval=", format(pval, scientific = T, digits = 3))) +
            ggplot2::xlab(col_time) +
            ggplot2::ylab("Survival") +
            ggplot2::scale_color_manual(values = c("LOW" = "grey", "HIGH" = "#D357A1")) +
            ggplot2::theme(
                axis.text = element_text(size = 14), axis.title = element_text(size = 14), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                plot.title = element_text(size = 14),
                legend.title = element_text(size = 10), legend.text = element_text(size = 10), legend.key.size = unit(0.8, "cm")
            )
    }

    res <- list("pval" = pval, "plot" = p, table = cell_annot_by_sample)
    return(res)
}


run_survival_logrank <- function(data,
                                 col_impact) {
    ### Run survival analysis
    if (is.null(col_impact)) {
        surv_diff <- survival::survdiff(survival::Surv(time, status) ~ cat, data)
    } else if (length(col_impact) == 1) {
        data$col_impact <- data[, col_impact]
        surv_diff <- survival::survdiff(survival::Surv(time, status) ~ cat + col_impact, data)
    } else if (length(col_impact) == 2) {
        for (i in 1:length(col_impact)) {
            data[, paste0("col_impact", i)] <- data[, col_impact[i]]
        }
        surv_diff <- survival::survdiff(survival::Surv(time, status) ~ cat + col_impact1 + col_impact2, data)
    } else if (length(col_impact) == 3) {
        for (i in 1:length(col_impact)) {
            data[, paste0("col_impact", i)] <- data[, col_impact[i]]
        }
        surv_diff <- survival::survdiff(survival::Surv(time, status) ~ cat + col_impact1 + col_impact2 + col_impact3, data)
    }

    return(surv_diff)
}

# --------------------------
### Epigenome functions ----
# --------------------------


{
    imputeWeights_custom <- function(matSVD = NULL, reducedDims = "IterativeLSI", dimsToUse = NULL,
                                     scaleDims = NULL, corCutOff = 0.75, td = 3, ka = 4, sampleCells = 5000,
                                     nRep = 2, k = 15, epsilon = 1, useHdf5 = TRUE, randomSuffix = FALSE,
                                     threads = getArchRThreads(), seed = 1, verbose = TRUE, logFile = createLogFile("addImputeWeights")) {
        set.seed(seed)
        tstart <- Sys.time()

        matDR <- matSVD
        matDR <- matDR[, dimsToUse]

        N <- nrow(matDR)
        rn <- rownames(matDR)
        if (!is.null(sampleCells)) {
            if (sampleCells > nrow(matDR)) {
                sampleCells <- NULL
            }
        }
        if (is.null(sampleCells)) {
            binSize <- N
            nRep <- 1
        } else {
            cutoffs <- lapply(seq_len(1000), function(x) {
                N / x
            }) %>% unlist()
            binSize <- min(cutoffs[order(abs(cutoffs - sampleCells))[1]] +
                1, N)
        }

        weightList <- .safelapply(seq_len(nRep), function(y) {
            # .logDiffTime(sprintf("Computing Partial Diffusion Matrix with Magic (%s of %s)",
            #                      y, nRep), t1 = tstart, verbose = FALSE, logFile = logFile)
            if (!is.null(sampleCells)) {
                idx <- sample(seq_len(nrow(matDR)), nrow(matDR))
                blocks <- split(rownames(matDR)[idx], ceiling(seq_along(idx) / binSize))
            } else {
                blocks <- list(rownames(matDR))
            }
            blockList <- lapply(seq_along(blocks), function(x) {
                if (x %% 10 == 0) {
                    # .logDiffTime(sprintf("Computing Partial Diffusion Matrix with Magic (%s of %s, Iteration %s of %s)",
                    #                      y, nRep, x, length(blocks)), t1 = tstart,
                    #              verbose = FALSE, logFile = logFile)
                }
                ix <- blocks[[x]]
                Nx <- length(ix)
                knnObj <- nabor::knn(data = matDR[ix, ], query = matDR[ix, ], k = k)
                knnIdx <- knnObj$nn.idx
                knnDist <- knnObj$nn.dists
                rm(knnObj)
                if (ka > 0) {
                    knnDist <- knnDist / knnDist[, ka]
                }
                if (epsilon > 0) {
                    W <- Matrix::sparseMatrix(rep(seq_len(Nx), k),
                        c(knnIdx),
                        x = c(knnDist), dims = c(Nx, Nx)
                    )
                } else {
                    W <- Matrix::sparseMatrix(rep(seq_len(Nx), k),
                        c(knnIdx),
                        x = 1, dims = c(Nx, Nx)
                    )
                }
                W <- W + Matrix::t(W)
                if (epsilon > 0) {
                    W@x <- exp(-(W@x / epsilon^2))
                }
                W <- W / Matrix::rowSums(W)
                Wt <- W
                for (i in seq_len(td)) {
                    Wt <- Wt %*% W
                }
                rownames(Wt) <- ix
                colnames(Wt) <- ix
                rm(knnIdx)
                rm(knnDist)
                rm(W)
                gc()
                return(Wt)
            }) %>% SimpleList()
            names(blockList) <- paste0("b", seq_along(blockList))
            return(blockList)
        }, threads = threads) %>% SimpleList()
        names(weightList) <- paste0("w", seq_along(weightList))
        imputeWeights <- SimpleList(Weights = weightList, Params = list(reducedDims = reducedDims, td = td, k = k, ka = ka, epsilon = epsilon))

        return(imputeWeights)
    }



    .safelapply <- function(..., threads = 1, preschedule = FALSE) {
        if (tolower(.Platform$OS.type) == "windows") {
            threads <- 1
        }

        if (threads > 1) {
            .requirePackage("parallel", source = "cran")
            o <- mclapply(..., mc.cores = threads, mc.preschedule = preschedule)

            errorMsg <- list()

            for (i in seq_along(o)) { # Make Sure this doesnt explode!
                if (inherits(o[[i]], "try-error")) {
                    capOut <- utils::capture.output(o[[i]])
                    capOut <- capOut[!grepl("attr\\(\\,|try-error", capOut)]
                    capOut <- head(capOut, 10)
                    capOut <- unlist(lapply(capOut, function(x) substr(x, 1, 250)))
                    capOut <- paste0("\t", capOut)
                    errorMsg[[length(errorMsg) + 1]] <- paste0(c(paste0("Error Found Iteration ", i, " : "), capOut), "\n")
                }
            }

            if (length(errorMsg) != 0) {
                errorMsg <- unlist(errorMsg)
                errorMsg <- head(errorMsg, 50)
                errorMsg[1] <- paste0("\n", errorMsg[1])
                stop(errorMsg)
            }
        } else {
            o <- lapply(...)
        }

        o
    }


    .requirePackage <- function(x = NULL, load = TRUE, installInfo = NULL, source = NULL) {
        if (x %in% rownames(installed.packages())) {
            if (load) {
                suppressPackageStartupMessages(require(x, character.only = TRUE))
            } else {
                return(0)
            }
        } else {
            if (!is.null(source) & is.null(installInfo)) {
                if (tolower(source) == "cran") {
                    installInfo <- paste0('install.packages("', x, '")')
                } else if (tolower(source) == "bioc") {
                    installInfo <- paste0('BiocManager::install("', x, '")')
                } else {
                    stop("Unrecognized package source, available are cran/bioc!")
                }
            }
            if (!is.null(installInfo)) {
                stop(paste0("Required package : ", x, " is not installed/found!\n  Package Can Be Installed : ", installInfo))
            } else {
                stop(paste0("Required package : ", x, " is not installed/found!"))
            }
        }
    }

    .suppressAll <- function(expr = NULL) {
        suppressPackageStartupMessages(suppressMessages(suppressWarnings(expr)))
    }

    .getAssay <- function(se = NULL, assayName = NULL) {
        .assayNames <- function(se) {
            names(SummarizedExperiment::assays(se))
        }
        if (is.null(assayName)) {
            o <- SummarizedExperiment::assay(se)
        } else if (assayName %in% .assayNames(se)) {
            o <- SummarizedExperiment::assays(se)[[assayName]]
        } else {
            stop(sprintf("assayName '%s' is not in assayNames of se : %s", assayName, paste(.assayNames(se), collapse = ", ")))
        }
        return(o)
    }
}


differential_activation_custom <- function(scExp, by = c("cell_cluster", "sample_id")[1], groups = c("1"),
                                           verbose = TRUE, progress = NULL) {
    .par_chisq <- function(row) {
        return(chisq.test(
            matrix(c(row[1], row[2], row[3], row[4]), ncol = 2),
            simulate.p.value = FALSE
        )$p.value)
    }

    if (is.null(groups)) {
        groups <- unique(colData(scExp)[, by])
    }

    list_res <- list()

    mat <- SingleCellExperiment::counts(
        scExp[SingleCellExperiment::rowData(scExp)$top_feature, ]
    )
    bin_mat <- Matrix::Matrix(mat > 0 + 0, sparse = TRUE)
    feature <- as.data.frame(SummarizedExperiment::rowRanges(
        scExp[SingleCellExperiment::rowData(scExp)$top_feature, ]
    ))
    feature <- data.frame(
        ID = feature[, "ID"], chr = feature[, "seqnames"],
        start = feature[, "start"], end = feature[, "end"]
    )


    for (group in groups) {
        if (!is.null(progress)) {
            progress$inc(
                detail = paste0("Calculating differential activation - ", group, "..."), amount = 0.9 / length(groups)
            )
        }
        if (verbose) cat("ChromSCape::differential_activation - Calculating differential activation for", group, ".\n")
        cluster_bin_mat <- bin_mat[, which(SingleCellExperiment::colData(scExp)[, by] %in% group)]
        cluster_mat <- mat[, which(SingleCellExperiment::colData(scExp)[, by] %in% group)]
        reference_bin_mat <- bin_mat[, which(!SingleCellExperiment::colData(scExp)[, by] %in% group)]
        reference_mat <- mat[, which(!SingleCellExperiment::colData(scExp)[, by] %in% group)]

        rectifier <- mean(Matrix::colSums(cluster_mat)) / mean(Matrix::colSums(reference_mat))
        group_sum <- Matrix::rowSums(cluster_bin_mat)
        group_activation <- group_sum / ncol(cluster_bin_mat)
        group_corrected_activation <- group_activation / rectifier

        reference_sum <- Matrix::rowSums(reference_bin_mat)
        reference_activation <- reference_sum / ncol(reference_bin_mat)

        n_cell_cluster <- ncol(cluster_bin_mat)
        n_cell_reference <- ncol(reference_bin_mat)
        other_group <- n_cell_cluster - group_sum
        other_ref <- n_cell_reference - reference_sum

        chisq_mat <- cbind(group_sum, other_group, reference_sum, other_ref)

        suppressWarnings({
            pvalues <- apply(chisq_mat, 1, .par_chisq)
        })

        q.values <- p.adjust(pvalues, method = "BH")
        logFCs <- log2(group_corrected_activation / reference_activation)
        if (any(is.nan(logFCs))) logFCs[which(is.nan(logFCs))] <- 0
        if (any(logFCs == Inf)) {
            logFCs[which(logFCs == Inf)] <- max(
                logFCs[which(!is.infinite(logFCs))]
            )
        }
        if (any(logFCs == -Inf)) {
            logFCs[which(logFCs == -Inf)] <- min(
                logFCs[which(!is.infinite(logFCs))]
            )
        }

        res <- data.frame(
            logFC.gpsamp = logFCs,
            qval.gpsamp = q.values,
            group_activation.gpsamp = group_activation,
            reference_activation.gpsamp = reference_activation
        )
        colnames(res) <- gsub("gpsamp", group, colnames(res))
        list_res[[group]] <- res
    }
    gc()

    names(list_res) <- NULL
    res <- cbind(feature, do.call("cbind", list_res))

    return(res)
}


.computeEnrichment <- function(matches = NULL, compare = NULL, background = NULL) {
    matches <- .getAssay(matches, grep("matches", names(assays(matches)), value = TRUE, ignore.case = TRUE))

    # Compute Totals
    matchCompare <- matches[compare, , drop = FALSE]
    matchBackground <- matches[background, , drop = FALSE]
    matchCompareTotal <- Matrix::colSums(matchCompare)
    matchBackgroundTotal <- Matrix::colSums(matchBackground)

    # Create Summary DF
    pOut <- data.frame(
        feature = colnames(matches),
        CompareFrequency = matchCompareTotal,
        nCompare = nrow(matchCompare),
        CompareProportion = matchCompareTotal / nrow(matchCompare),
        BackgroundFrequency = matchBackgroundTotal,
        nBackground = nrow(matchBackground),
        BackgroundProporition = matchBackgroundTotal / nrow(matchBackground)
    )

    # Enrichment
    pOut$Enrichment <- pOut$CompareProportion / pOut$BackgroundProporition

    # Get P-Values with Hyper Geometric Test
    pOut$mlog10p <- lapply(seq_len(nrow(pOut)), function(x) {
        p <- -phyper(pOut$CompareFrequency[x] - 1, # Number of Successes the -1 is due to cdf integration
            pOut$BackgroundFrequency[x], # Number of all successes in background
            pOut$nBackground[x] - pOut$BackgroundFrequency[x], # Number of non successes in background
            pOut$nCompare[x], # Number that were drawn
            lower.tail = FALSE, log.p = TRUE
        ) # P[X > x] Returns LN must convert to log10
        return(p / log(10))
    }) %>%
        unlist() %>%
        round(4)

    # Minus Log10 Padj
    pOut$mlog10Padj <- pmax(pOut$mlog10p - log10(ncol(pOut)), 0)
    pOut <- pOut[order(pOut$mlog10p, decreasing = TRUE), , drop = FALSE]

    pOut
}

# --------------------------
### Compare proportions ----
# --------------------------


library(dplyr)
library(rstatix)
library(rlang)

wilcox_with_fold <- function(data, key, group_var = "group", value_var = "value", groups_order = NULL, paired = T) {
    # Accès dynamique aux colonnes
    group_sym <- sym(group_var)
    value_sym <- sym(value_var)

    # Extraire les 2 groupes
    if (is.null(groups_order)) {
        groups <- unique(data[[group_var]])
    } else {
        groups <- groups_order
    }

    if (length(groups) != 2) {
        return(NULL)
    }

    group1 <- groups[1]
    group2 <- groups[2]

    # Moyennes des deux groupes
    mean1 <- mean(data %>% filter(!!group_sym == group1) %>% pull(!!value_sym))
    mean2 <- mean(data %>% filter(!!group_sym == group2) %>% pull(!!value_sym))

    # Hypothèse directionnelle
    alternative <- if (mean2 < mean1) "greater" else "less"
    fold_change <- log2(mean2 / mean1)

    # Formule dynamique pour le test
    test_formula <- as.formula(paste(value_var, "~", group_var))

    # Appliquer le test
    test_result <- tryCatch(
        {
            data %>%
                rstatix::wilcox_test(formula = test_formula, alternative = alternative, paired = paired) %>%
                rstatix::add_significance() %>%
                dplyr::mutate(
                    category = key,
                    fold_change = fold_change,
                    direction = alternative
                )
        },
        error = function(e) {
            return(NULL)
        }
    )

    return(test_result)
}


wilcox_with_fold_median <- function(data, key, group_var = "group", value_var = "value", groups_order = NULL, paired = T) {
    # Accès dynamique aux colonnes
    group_sym <- sym(group_var)
    value_sym <- sym(value_var)

    # Extraire les 2 groupes
    if (is.null(groups_order)) {
        groups <- unique(data[[group_var]])
    } else {
        groups <- groups_order
    }

    if (length(groups) != 2) {
        return(NULL)
    }

    group1 <- groups[1]
    group2 <- groups[2]

    # Moyennes des deux groupes
    mean1 <- median(data %>% filter(!!group_sym == group1) %>% pull(!!value_sym))
    mean2 <- median(data %>% filter(!!group_sym == group2) %>% pull(!!value_sym))

    # Hypothèse directionnelle
    alternative <- if (mean2 < mean1) "greater" else "less"
    fold_change <- log2(mean2 / mean1)

    # Formule dynamique pour le test
    test_formula <- as.formula(paste(value_var, "~", group_var))

    # Appliquer le test
    test_result <- tryCatch(
        {
            data %>%
                rstatix::wilcox_test(formula = test_formula, alternative = alternative, paired = paired) %>%
                rstatix::add_significance() %>%
                dplyr::mutate(
                    category = key,
                    fold_change = fold_change,
                    direction = alternative,
                    median1 = mean1,
                    median2 = mean2
                )
        },
        error = function(e) {
            return(NULL)
        }
    )

    return(test_result)
}


# --------------------------
### Enrichment functions ----
# --------------------------


.load_MSIGdb <- function(ref, GeneSetClasses) {
    if ((!ref %in% c("hg38", "mm10"))) {
        stop(
            "Reference genome (ref) must be ",
            "'hg38' or 'mm10' if gene sets not specified."
        )
    }
    stopifnot(is.character(GeneSetClasses))
    message(
        paste0(
            "Loading ",
            ref,
            " MSigDB gene sets."
        )
    )
    columns <- c("gs_name", "gs_cat", "gene_symbol")
    if (ref == "hg38") {
        GeneSetsDf <- msigdbr::msigdbr("Homo sapiens")[, columns]
    }
    if (ref == "mm10") {
        GeneSetsDf <- msigdbr::msigdbr("Mus musculus")[, columns]
    }
    colnames(GeneSetsDf) <- c("Gene.Set", "Class", "Genes")
    system.time({
        GeneSetsDf <- GeneSetsDf %>%
            dplyr::group_by(
                .data$Gene.Set, .data$Class
            ) %>%
            dplyr::summarise("Genes" = paste(.data$Genes,
                collapse = ","
            ))
    })
    corres <- data.frame(
        long_name = c(
            "c1_positional", "c2_curated", "c3_motif",
            "c4_computational", "c5_GO", "c6_oncogenic",
            "c7_immunologic", "hallmark"
        ), short_name = c(
            paste0("C", seq_len(7)), "H"
        )
    )
    GeneSetsDf$Class <- corres$long_name[
        match(GeneSetsDf$Class, corres$short_name)
    ]
    GeneSetsDf <- GeneSetsDf[which(GeneSetsDf$Class %in% GeneSetClasses), ]
    GeneSets <- lapply(GeneSetsDf$Gene.Set, function(x) {
        unlist(strsplit(
            GeneSetsDf$Genes[which(GeneSetsDf$Gene.Set == x)],
            split = ","
        ))
    })
    names(GeneSets) <- GeneSetsDf$Gene.Set
    return(GeneSets)
}

enrichment_markers_provided_markers <- function(top_markers, ref = "hg38", GeneSetClasses = c("c2_curated", "hallmark"), qval = 0.05, gene_sets_msigdb = NULL, reference_gene_list = NULL) {
    if (is.null(gene_sets_msigdb)) {
        gene_sets_msigdb <- .load_MSIGdb(ref, GeneSetClasses)
    }

    gene_sets <- gene_sets_msigdb
    sep <- ";"
    if (is.null(reference_gene_list)) {
        possibleIds <- unique(unlist(gene_sets))
    } else {
        possibleIds <- intersect(unique(unlist(gene_sets)), unique(reference_gene_list))
    }


    mylist <- unique(top_markers)
    gene.sets <- lapply(gene_sets, unique)
    nids <- length(possibleIds)
    gene.sets <- lapply(gene.sets, function(x) intersect(x, possibleIds))
    nref <- as.numeric(lapply(gene.sets, length))
    gene.sets <- gene.sets[nref > 0]
    n <- length(mylist)
    fun <- function(x) {
        y <- intersect(x, mylist)
        nx <- length(x)
        ny <- length(y)
        pval <- stats::phyper(ny - 1, nx, nids - nx, n, lower.tail = FALSE)
        c(nx, ny, pval, paste(y, collapse = sep))
    }
    tmp <- as.data.frame(t(as.matrix(vapply(gene.sets, fun, FUN.VALUE = c(
        "Nb_of_genes" = 0, "Nb_of_deregulated_genes" = 0,
        "p-value" = 0, "Deregulated_genes" = ""
    )))))
    rownames(tmp) <- names(gene.sets)
    for (i in seq_len(3)) {
        tmp[, i] <- as.numeric(
            as.character(tmp[, i])
        )
    }
    tmp <- data.frame(
        tmp[, seq_len(3)], p.adjust(tmp[, 3], method = "BH"), tmp[, 4]
    )
    names(tmp) <- c(
        "Nb_of_genes", "Nb_of_deregulated_genes",
        "p-value", "q-value", "Deregulated_genes"
    )
    tmp
    tmp <- tmp[which(tmp$`q-value` < qval), ]
    tmp <- tmp[order(tmp$`q-value`), ]

    return(tmp)
}


get_enrichment_metaprograms <- function(state_genes, GeneSetClasses = c("hallmark", "c2_curated"), qval = 0.05, reference_gene_list = NULL, gene_sets_msigdb = NULL) {
    if (is.null(gene_sets_msigdb)) {
        gene_sets_msigdb <- .load_MSIGdb(ref = "hg38", GeneSetClasses = GeneSetClasses)
    }

    list_enrichment <- list()
    program <- c()
    pathway <- c()
    qvalue <- c()
    for (prog in names(state_genes)) {
        tmp <- enrichment_markers_provided_markers(ref = "hg38", GeneSetClasses = GeneSetClasses, top_markers = state_genes[[prog]], qval = qval, gene_sets_msigdb = gene_sets_msigdb, reference_gene_list = reference_gene_list)
        list_enrichment[[prog]] <- tmp
        if (nrow(tmp) != 0) {
            program <- c(program, rep(prog, nrow(tmp)))
            pathway <- c(pathway, rownames(tmp))
            qvalue <- c(qvalue, tmp$`q-value`)
        }
    }
    df <- data.frame(program, pathway, qvalue)
    if (nrow(df) > 0) {
        df$log10qval <- (-log10(df$qvalue))
    }

    return(list(enrichment_df = df, enrichment_list = list_enrichment))
}





# --------------------------
### Palette colors ----
# --------------------------


palette_samples_H3K4me1 <- c(
    "P1_NT" = "#1B998B",
    "P2_RD" = "navyblue",
    "P3_NT" = "#F8BBD0",
    "P3_RD" = "#F06292",
    "P4_NT" = "#DECBB7",
    "P5_NT" = "#F08080",
    "P5_RD" = "#BF3100",
    "P6_NT" = "#FF7F51",
    "P7_NT" = "#FFFD82",
    "P7_RD" = "#FFB74D",
    "P8_NT" = "#88D498",
    "P8_RD" = "#1A936F",
    "P9_NT" = "#B6A3E0",
    "P9_RD" = "#BA68C8",
    "P10_NT" = "#B5D1EC",
    "P10_P" = "#0a2472",
    "P11_NT" = "#4DD0E1",
    "P12_RD" = "#893168",
    "P13_NT" = "#ACF39D",
    "P13_RD" = "#3DDC97",
    "P14_NT" = "#F7AC8A",
    "P14_RD" = "#EE6C4D",
    "P15_NT" = "#B3E5FC",
    "P15_RD" = "#64B5F6"
)

palette_samples_frozen <- c(
    "P1_NT" = "#1B998B",
    "P1_RD" = "#0B7A75",
    "P2_NT" = "#7986CB",
    "P2_RD" = "navyblue",
    "P3_NT" = "#F8BBD0",
    "P3_RD" = "#F06292",
    "P4_NT" = "#DECBB7",
    "P4_RD" = "#CF9893",
    "P5_NT" = "#F08080",
    "P5_RD" = "#BF3100",
    "P6_NT" = "#FF7F51",
    "P6_P" = "#813405",
    "P7_NT" = "#FFFD82",
    "P7_RD" = "#FFB74D",
    "P8_NT" = "#88D498",
    "P8_RD" = "#1A936F",
    "P8_P" = "#132a13",
    "P9_NT" = "#B6A3E0",
    "P9_RD" = "#BA68C8",
    "P9_P" = "#4d194d",
    "P10_NT" = "#B5D1EC",
    "P10_RD" = "#3B58A7",
    "P10_P" = "#0a2472",
    "P11_NT" = "#4DD0E1",
    "P11_RD" = "#0fa3b1",
    "P12_NT" = "#EB8A90",
    "P12_RD" = "#893168",
    "P13_NT" = "#ACF39D",
    "P13_RD" = "#3DDC97",
    "P14_NT" = "#F7AC8A",
    "P14_RD" = "#EE6C4D",
    "P15_NT" = "#B3E5FC",
    "P15_RD" = "#64B5F6",
    "P16_NT" = "lightgrey",
    "P16_RD" = "#687386"
)

palette_samples_frozen_fig3 <- c(
    "P1_NT" = "#1B998B",
    "P1_RD" = "#0B7A75",
    "P2_NT" = "#7986CB",
    "P2_RD" = "navyblue",
    "P3_NT" = "#F8BBD0",
    "P3_RD" = "#F06292",
    "P4_NT" = "#DECBB7",
    "P4_RD" = "#CF9893",
    "P5_NT" = "#F08080",
    "P5_RD" = "#BF3100",
    "P6_NT" = "#FF7F51",
    "P6_RD" = "#d95221",
    "P6_P" = "#813405",
    "P7_NT" = "#FFFD82",
    "P7_RD" = "#FFB74D",
    "P8_NT" = "#88D498",
    "P8_RD" = "#1A936F",
    "P8_P" = "#132a13",
    "P9_NT" = "#B6A3E0",
    "P9_RD" = "#BA68C8",
    "P9_P" = "#4d194d",
    "P10_NT" = "#B5D1EC",
    "P10_RD" = "#3B58A7",
    "P10_P" = "#0a2472",
    "P11_NT" = "#4DD0E1",
    "P11_RD" = "#0fa3b1",
    "P12_NT" = "#EB8A90",
    "P12_RD" = "#893168",
    "P13_NT" = "#ACF39D",
    "P13_RD" = "#3DDC97",
    "P14_NT" = "#F7AC8A",
    "P14_RD" = "#EE6C4D",
    "P15_NT" = "#B3E5FC",
    "P15_RD" = "#64B5F6",
    "P16_NT" = "lightgrey",
    "P16_RD" = "#687386"
)

palette_patient <- c(
    "P1" = "#0B7A75",
    "P2" = "#7986CB",
    "P3" = "#F06292",
    "P4" = "#DECBB7",
    "P5" = "#BF3100",
    "P6" = "#FF7F51",
    "P7" = "#FFFD82",
    "P8" = "#88D498",
    "P9" = "#BA68C8",
    "P10" = "#3B58A7",
    "P11" = "#4DD0E1",
    "P12" = "#893168",
    "P13" = "#3DDC97",
    "P14" = "#EE6C4D",
    "P15" = "#64B5F6",
    "P16" = "#687386"
)

palette_cell_type <- c(
    "Epithelial_cells" = "#00A08A",
    "Fibroblasts" = "#85D4E3",
    "Macrophages" = "#FFB74D",
    "T_cells" = "#FD6467",
    "Endothelial_cells" = "#AED581"
)

palette_cell_type_broad <- c(
    "Epithelial" = "#00A08A",
    "Mesenchymal" = "#85D4E3",
    "Myeloid" = "#FFB74D",
    "Lymphoid" = "#FD6467",
    "Vascular" = "#AED581",
    "Not_analyzed" = "grey",
    "Nervous" = "grey"
)

palette_cell_type_scRNA <- c(
    "Tumor_cells" = "#00A08A",
    "Fibroblasts" = "#85D4E3",
    "Macrophages" = "#FFB74D",
    "Mast_cells" = "#f39b17",
    "Dendritic_cells" = "#FAD510",
    "T_cells" = "#FD6467",
    "B_cells" = "#E34345",
    "Plasma_cells" = "#A12A2D",
    "NK_cells" = "#FC7E7F",
    "Endothelial_cells" = "#AED581"
)

palette_TME <- c(
    "Not_analyzed" = "grey90",
    "Mye_Macro" = "#FFB74D",
    "Lym_T.CD4" = "#FFD1D1", "Lym_T.CD4_Exhausted" = "#FFA3A3", "Lym_NK_Exhausted" = "#FC7E7F", "Lym_T.CD8" = "#FD6467", "Lym_B" = "#E34345", "Lym_B_Plasmocyte" = "#A12A2D",
    "Mes_Adipocyte" = "#BFEAF0", "Mes_Fibro" = "#85D4E3", "Mes_Fibro_Myofibro" = "#3399AA",
    "Vas_Endo_Veinous" = "#E9F5DB", "Vas_Endo_Capillary" = "#D0EDA8", "Vas_Endo_Arterial" = "#75a156ff", "Vas_Pericyte" = "#C1E69A", "Vas_Fibro" = "#AED581", "Vas_SMC" = "#8BC34A"
)

palette_TME_fig6 <- c(
    "Mye_Macro" = "#FFB74D",
    "Lym_T" = "#FFA3A3", "Lym_NK_Exhausted" = "#FC7E7F", "Lym_B" = "#E34345", "Lym_B_Plasmocyte" = "#A12A2D",
    "Mes_Adipocyte" = "#BFEAF0", "Mes_Fibro" = "#85D4E3", "Mes_Fibro_Myofibro" = "#3399AA",
    "Vas_Endo" = "#D0EDA8", "Vas_Mural" = "#8BC34A"
)

palette_condition <- c(
    "NT" = "#E6A0C4",
    "RD" = "#5785C1",
    "P" = "#C6CDF7"
)

palette_primary_metastasis <- c(
    "Primary" = "#EBA131",
    "Metastasis" = "#625793"
)


palette_states_public_short <- list(
    Cell_Cycle = "#FAD510",
    Cilia = "#85D4E3",
    EMT = "#5785C1",
    Interferons = "#E6A0C4",
    MHCII = "#C6CDF7",
    TNFA = "#FD6467",
    Hypoxia.Glycolysis = "#54D8B1",
    Ribosomal = "#02401B",
    Mitochondrial = "#EAD3BF",
    Metallothioneins = "#B6854D",
    MYC_targets = "#E54E21",
    Heat_Stress = "#6C8645",
    Unfolded_Protein_Response = "#8D8680"
)

palette_states_public <- list(
    Cell_Cycle = "#FAD510",
    Cilia = "#85D4E3",
    EMT = "#5785C1",
    Interferons = "#E6A0C4",
    MHCII = "#C6CDF7",
    TNFA = "#FD6467",
    Hypoxia.Glycolysis = "#54D8B1",
    Ribosomal = "#02401B",
    Mitochondrial = "#EAD3BF",
    Metallothioneins = "#B6854D",
    MYC_targets = "#E54E21",
    Heat_Stress = "#6C8645",
    Unfolded_Protein_Response = "#8D8680",
    Unknown = "#D9D0D3"
)

palette_states_public_all <- list(
    Cell_Cycle_G2M = "#FAD510",
    Cell_Cycle_G1S = "#C93312",
    Cilia = "#85D4E3",
    EMT.ECM = "#0A9F9D",
    EMT.Inflamm.MHCII = "#5785C1",
    Interferons.MHCII = "#E6A0C4",
    MHCII = "#C6CDF7",
    Stress = "#FD6467",
    Hypoxia.Glycolysis = "#54D8B1",
    Ribosomal = "#02401B",
    Mitochondrial = "#EAD3BF",
    Metallothioneins = "#B6854D",
    MYC_targets = "#E54E21",
    Heat_Stress = "#6C8645",
    Unfolded_Protein_Response = "#8D8680",
    Unknown = "#D9D0D3"
)

palette_datasets <- list(SCANDARE = "#E25E27", ZHANG = "#C05950", VASQUEZ = "#F99B86")


# {
#     library(RColorBrewer)
#     n <- length(unique(metadata_cells$patient_id_nd))
#     qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#     col_vector = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
#     set.seed(123)
#     palette_patient_RNA=sample(col_vector, n,replace = F)
#     #saveRDS(palette_patient_RNA,"./palette_patient_3datasets.rds")
#     #palette_samples_all=col_vector[1:n]
# }
palette_patient_scRNA <- readRDS("./palette_patient_3datasets.rds")



# --------------------------
### Significance symbols ----
# --------------------------

symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", ""))


# --------------------------
### List cell type markers ----
# --------------------------

list_cell_type_markers_snRNA <- list(
    Epi = c("KRT7", "EPCAM", "MUC1", "MUC16", "PAX8"),
    Lym_B = c("BANK1", "BLK", "MS4A1"),
    Lym_B_Plasmocyte = c("IGHG1", "IGHM", "JCHAIN"),
    Lym_T.CD4 = c("IL7R", "CD4", "CCR7", "SELL"),
    Lym_T.CD4_Exhausted = c("IL7R", "CD4", "FOXP3", "ICOS", "CTLA4", "TIGIT"),
    Lym_T.CD8 = c("CD8A", "CD8B", "GZMK", "THEMIS"),
    Lym_NK_Exhausted = c("NCAM1", "KLRD1", "PRF1", "HAVCR2"),
    Mye = c("PTPRC", "CD68", "MRC1", "CD163"),
    Mye_Macro = c("MSR1", "AREG", "MRC1", "CD163", "EREG"),
    Mye_Macro_M2 = c("APOE", "MARCO", "MRC1", "CD163", "SELENOP"),
    Mye_Proliferation = c("APOE", "MKI67", "MELK", "CD163", "AURKB"),
    Mes_Fibro = c("DCN", "PDGFRA", "LUM"),
    Mes_Fibro_Myofibro = c("POSTN", "ACTA2", "ACTG2", "CNN1", "MYL9"),
    Mes_Adipocyte = c("ADIPOQ", "CIDEC", "FABP4"),
    Vas_SMC = c("CSPG4", "PLN", "PDGFRB", "MCAM", "MYH11"),
    Vas_Pericyte = c("CSPG4", "ACTA2", "PDGFRB", "RGS5"),
    Vas_Fibro = c("COL1A1", "DCN", "FBLN1", "LUM", "SFRP2", "RGS5"),
    Vas_Endo_Arterial = c("ACE", "EMCN", "PRDM16", "VWF", "GJA5", "HEY1"),
    Vas_Endo_Veinous = c("AQP1", "EMCN", "ACKR1", "VWF", "NR2F2"),
    Vas_Endo_Capillary = c("INSR", "EMCN", "ETS1", "VWF", "ANGPT2")
)

list_cell_type_markers_scRNA <- list(
    Tumor_cells = c("EPCAM", "MUC1", "KRT7", "WFDC2"),
    B_cells = c("CD79A", "MS4A1"),
    Plasma_cells = c("IGHG1", "IGHM", "JCHAIN"),
    T_cells = c("IL7R", "CD2"),
    NK_cells = c("NKG7", "GZMA"),
    Macrophages = c("MSR1", "C1QA", "MRC1", "CD163", "SELENOP"),
    Dendritic_cells = c("CPVL"),
    Mast_cells = c("TPSB2", "CPA3", "KIT", "MS4A2"),
    Fibroblasts = c("COL1A1", "DCN", "PDGFRA", "LUM"),
    Endothelial_cells = c("VWF", "CDH5")
)
