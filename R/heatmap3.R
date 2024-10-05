## Acquired 3/17/2021 from
## https://github.com/obigriffith/biostar-tutorials/raw/da01d867ec229117da0f5a2f1d8e0e3b1ff5b4ec/Heatmaps/heatmap.3.R
## GPLv2 - see https://github.com/obigriffith/biostar-tutorials/blob/da01d867ec229117da0f5a2f1d8e0e3b1ff5b4ec/LICENSE
## Various typos, assumed imports, and code style issues have been fixed

## Internal functions used by heatmap.3

.TENETHeatmap3Scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low) / (high - low)
    x
}

.TENETHeatmap3Invalid <- function(x) {
    if (missing(x) || is.null(x) || length(x) == 0) {
        return(TRUE)
    }
    if (is.list(x)) {
        return(all(vapply(x, .TENETHeatmap3Invalid, logical(1))))
    } else if (is.vector(x)) {
        return(all(is.na(x)))
    } else {
        return(FALSE)
    }
}

## Main heatmap.3 function

.TENETInternalHeatmap3 <- function(
    x,
    Rowv = TRUE,
    Colv = if (symm) "Rowv" else TRUE,
    distfun = stats::dist,
    hclustfun = stats::hclust,
    dendrogram = c("both", "row", "column", "none"),
    symm = FALSE,
    scale = c("none", "row", "column"),
    na.rm = TRUE,
    revC = identical(Colv, "Rowv"),
    addExpr,
    breaks,
    symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
    col = "heat.colors",
    colsep,
    rowsep,
    sepcolor = "white",
    sepwidth = c(0.05, 0.05),
    cellnote,
    notecex = 1,
    notecol = "cyan",
    naColor = graphics::par("bg"),
    trace = c("none", "column", "row", "both"),
    tracecol = "cyan",
    hline = stats::median(breaks),
    vline = stats::median(breaks),
    linecol = tracecol,
    margins = c(5, 5),
    ColSideColors,
    RowSideColors,
    sideHeightFraction = 0.3,
    cexRow = 0.2 + 1 / log10(nr),
    cexCol = 0.2 + 1 / log10(nc),
    labRow = NULL,
    labCol = NULL,
    key = TRUE,
    keysize = 1.5,
    densityInfo = c("none", "histogram", "density"),
    denscol = tracecol,
    symkey = max(x < 0, na.rm = TRUE) || symbreaks,
    densadj = 0.25,
    main = NULL,
    xlab = NULL,
    ylab = NULL,
    lmat = NULL,
    lhei = NULL,
    lwid = NULL,
    ColSideColorsSize = 1,
    RowSideColorsSize = 1,
    KeyValueName = "Value",
    ...) {
    x <- as.matrix(x)
    retval <- list()
    scale <- if (symm && missing(scale)) {
        "none"
    } else {
        match.arg(scale)
    }
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    densityInfo <- match.arg(densityInfo)
    if (length(col) == 1 && is.character(col)) {
        col <- get(col, mode = "function")
    }
    if (!missing(breaks) && (scale != "none")) {
        warning(
            "Using scale=\"row\" or scale=\"column\" when breaks are ",
            "specified can produce unpredictable results. ",
            "Please consider using only one or the other."
        )
    }
    if (is.null(Rowv) || .isSingleNA(Rowv)) {
        Rowv <- FALSE
    }
    if (is.null(Colv) || .isSingleNA(Colv)) {
        Colv <- FALSE
    } else if (length(Colv) == 1 && Colv == "Rowv" && !isTRUE(Rowv)) {
        Colv <- FALSE
    }
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) {
        stop("`x' must be a numeric matrix")
    }
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) {
        stop("`x' must have at least 2 rows and 2 columns")
    }
    if (!is.numeric(margins) || length(margins) != 2) {
        stop("`margins' must be a numeric vector of length 2")
    }
    if (missing(cellnote)) {
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    }
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv)) {
                dendrogram <- "column"
            } else {
                dendrogram <- "none"
            }
            warning(
                "Discrepancy: Rowv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting row dendrogram."
            )
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv)) {
                dendrogram <- "row"
            } else {
                dendrogram <- "none"
            }
            warning(
                "Discrepancy: Colv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting column dendrogram."
            )
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- stats::order.dendrogram(ddr)
    } else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- stats::as.dendrogram(hcr)
        ddr <- stats::reorder(ddr, Rowv)
        rowInd <- stats::order.dendrogram(ddr)
        if (nr != length(rowInd)) {
            stop("row dendrogram ordering gave index of wrong length")
        }
    } else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- stats::as.dendrogram(hcr)
        ddr <- stats::reorder(ddr, Rowv)
        rowInd <- stats::order.dendrogram(ddr)
        if (nr != length(rowInd)) {
            stop("row dendrogram ordering gave index of wrong length")
        }
    } else {
        rowInd <- rev(seq_len(nr))
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- stats::order.dendrogram(ddc)
    } else if (identical(Colv, "Rowv")) {
        if (nr != nc) {
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        }
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- stats::order.dendrogram(ddc)
        } else {
            colInd <- rowInd
        }
    } else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm) {
            x
        } else {
            t(x)
        }))
        ddc <- stats::as.dendrogram(hcc)
        ddc <- stats::reorder(ddc, Colv)
        colInd <- stats::order.dendrogram(ddc)
        if (nc != length(colInd)) {
            stop("column dendrogram ordering gave index of wrong length")
        }
    } else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm) {
            x
        } else {
            t(x)
        }))
        ddc <- stats::as.dendrogram(hcc)
        ddc <- stats::reorder(ddc, Colv)
        colInd <- stats::order.dendrogram(ddc)
        if (nc != length(colInd)) {
            stop("column dendrogram ordering gave index of wrong length")
        }
    } else {
        colInd <- seq_len(nc)
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow)) {
        labRow <- if (is.null(rownames(x))) {
            seq_len(nr)[rowInd]
        } else {
            rownames(x)
        }
    } else {
        labRow <- labRow[rowInd]
    }
    if (is.null(labCol)) {
        labCol <- if (is.null(colnames(x))) {
            seq_len(nc)[colInd]
        } else {
            colnames(x)
        }
    } else {
        labCol <- labCol[colInd]
    }
    if (scale == "row") {
        retval$rowMeans <- rmns <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rmns)
        retval$rowSDs <- sx <- apply(x, 1, stats::sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    } else if (scale == "column") {
        retval$colMeans <- rmns <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rmns)
        retval$colSDs <- sx <- apply(x, 2, stats::sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
        if (missing(col) || is.function(col)) {
            breaks <- 16
        } else {
            breaks <- length(col) + 1
        }
    }
    if (length(breaks) == 1) {
        if (!symbreaks) {
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                length = breaks
            )
        } else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    ncol <- length(breaks) - 1
    if (is.function(col)) {
        col <- col(ncol)
    }
    minBreaks <- min(breaks)
    maxBreaks <- max(breaks)
    x[x < minBreaks] <- minBreaks
    x[x > maxBreaks] <- maxBreaks
    if (missing(lhei) || is.null(lhei)) {
        lhei <- c(keysize, 4)
    }
    if (missing(lwid) || is.null(lwid)) {
        lwid <- c(keysize, 4)
    }
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)

        if (!missing(ColSideColors)) {
            if (!is.character(ColSideColors) || nrow(ColSideColors) != nc) {
                stop("'ColSideColors' must be a matrix of nrow(x) rows")
            }
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
            lhei <- c(
                lhei[1], sideHeightFraction * ColSideColorsSize / 2, lhei[2]
            )
        }

        if (!missing(RowSideColors)) {
            if (!is.character(RowSideColors) || ncol(RowSideColors) != nr) {
                stop("'RowSideColors' must be a matrix of ncol(x) columns")
            }
            lmat <- cbind(
                lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[, 2] + 1
            )
            lwid <- c(
                lwid[1], sideHeightFraction * RowSideColorsSize / 2, lwid[2]
            )
        }
        lmat[is.na(lmat)] <- 0
    }

    if (length(lhei) != nrow(lmat)) {
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    }
    if (length(lwid) != ncol(lmat)) {
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    }
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op))

    graphics::layout(lmat, widths = lwid, heights = lhei, respect = FALSE)

    if (!missing(RowSideColors)) {
        if (!is.matrix(RowSideColors)) {
            graphics::par(mar = c(margins[1], 0, 0, 0.5))
            graphics::image(
                rbind(seq_len(nr)),
                col = RowSideColors[rowInd], axes = FALSE
            )
        } else {
            graphics::par(mar = c(margins[1], 0, 0, 0.5))
            rsc <- t(RowSideColors[, rowInd, drop = FALSE])
            rscColors <- matrix()
            rscNames <- names(table(rsc))
            rscIndex <- 1
            for (rscName in rscNames) {
                rscColors[rscIndex] <- rscName
                rsc[rsc == rscName] <- rscIndex
                rscIndex <- rscIndex + 1
            }
            rsc <- matrix(as.numeric(rsc), nrow = nrow(rsc))
            graphics::image(t(rsc), col = as.vector(rscColors), axes = FALSE)
            if (length(rownames(RowSideColors)) > 0) {
                graphics::axis(1,
                    (seq_len(ncol(rsc)) - 1) / max(1, (ncol(rsc) - 1)),
                    rownames(RowSideColors),
                    las = 2,
                    tick = FALSE
                )
            }
        }
    }

    if (!missing(ColSideColors)) {
        if (!is.matrix(ColSideColors)) {
            graphics::par(mar = c(0.5, 0, 0, margins[2]))
            graphics::image(
                cbind(seq_len(nc)),
                col = ColSideColors[colInd], axes = FALSE
            )
        } else {
            graphics::par(mar = c(0.5, 0, 0, margins[2]))
            csc <- ColSideColors[colInd, , drop = FALSE]
            cscColors <- matrix()
            cscNames <- names(table(csc))
            cscIndex <- 1
            for (cscName in cscNames) {
                cscColors[cscIndex] <- cscName
                csc[csc == cscName] <- cscIndex
                cscIndex <- cscIndex + 1
            }
            csc <- matrix(as.numeric(csc), nrow = nrow(csc))
            graphics::image(csc, col = as.vector(cscColors), axes = FALSE)
            if (length(colnames(ColSideColors)) > 0) {
                graphics::axis(
                    2,
                    (seq_len(ncol(csc)) - 1) / max(1, (ncol(csc) - 1)),
                    colnames(ColSideColors),
                    las = 2,
                    tick = FALSE
                )
            }
        }
    }

    graphics::par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- rev(seq_len(nr))
        if (exists("ddr")) {
            ddr <- rev(ddr)
        }
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    } else {
        iy <- seq_len(nr)
    }
    graphics::image(
        seq_len(nc),
        seq_len(nr),
        x,
        xlim = 0.5 + c(
            0,
            nc
        ),
        ylim = 0.5 + c(
            0,
            nr
        ),
        axes = FALSE,
        xlab = "",
        ylab = "",
        col = col,
        breaks = breaks,
        ...
    )
    retval$carpet <- x
    if (exists("ddr")) {
        retval$rowDendrogram <- ddr
    }
    if (exists("ddc")) {
        retval$colDendrogram <- ddc
    }
    retval$breaks <- breaks
    retval$col <- col
    if (!.TENETHeatmap3Invalid(naColor) & any(is.na(x))) {
        mmat <- ifelse(is.na(x), 1, NA)
        graphics::image(seq_len(nc), seq_len(nr), mmat,
            axes = FALSE, xlab = "", ylab = "",
            col = naColor, add = TRUE
        )
    }
    graphics::axis(1, seq_len(nc),
        labels = labCol, las = 2, line = -0.5, tick = 0,
        cex.axis = cexCol
    )
    if (!is.null(xlab)) {
        graphics::mtext(xlab, side = 1, line = margins[1] - 1.25)
    }
    graphics::axis(4, iy,
        labels = labRow, las = 2, line = -0.5, tick = 0,
        cex.axis = cexRow
    )
    if (!is.null(ylab)) {
        graphics::mtext(ylab, side = 4, line = margins[2] - 1.25)
    }
    if (!missing(addExpr)) {
        eval(substitute(addExpr))
    }
    if (!missing(colsep)) {
        for (csep in colsep) {
            graphics::rect(
                xleft = csep + 0.5,
                ybottom = rep(0, length(csep)),
                xright = csep + 0.5 + sepwidth[1],
                ytop = rep(ncol(x) + 1, csep),
                lty = 1,
                lwd = 1,
                col = sepcolor,
                border = sepcolor
            )
        }
    }
    if (!missing(rowsep)) {
        for (rsep in rowsep) {
            graphics::rect(
                xleft = 0,
                ybottom = (ncol(x) + 1 - rsep) - 0.5,
                xright = nrow(x) + 1,
                ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2],
                lty = 1,
                lwd = 1,
                col = sepcolor,
                border = sepcolor
            )
        }
    }
    minScale <- min(breaks)
    maxScale <- max(breaks)
    xScaled <- .TENETHeatmap3Scale01(t(x), minScale, maxScale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vlineVals <- .TENETHeatmap3Scale01(vline, minScale, maxScale)
        for (i in colInd) {
            if (!is.null(vline)) {
                graphics::abline(
                    v = i - 0.5 + vlineVals, col = linecol,
                    lty = 2
                )
            }
            xv <- rep(i, nrow(xScaled)) + xScaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- seq_along(xv) - 0.5
            graphics::lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hlineVals <- .TENETHeatmap3Scale01(hline, minScale, maxScale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                graphics::abline(
                    h = i - 0.5 + hlineVals, col = linecol,
                    lty = 2
                )
            }
            yv <- rep(i, ncol(xScaled)) + xScaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- rev(seq_along(yv)) - 0.5
            graphics::lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote)) {
        graphics::text(
            x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
            col = notecol, cex = notecex
        )
    }
    graphics::par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    } else {
        graphics::plot.new()
    }
    graphics::par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    } else {
        graphics::plot.new()
    }
    if (!is.null(main)) {
        graphics::title(main, cex.main = 1.5 * op$cex.main)
    }
    if (key) {
        graphics::par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            maxRaw <- max(abs(c(x, breaks)), na.rm = TRUE)
            minRaw <- -maxRaw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        } else {
            minRaw <- min(x, na.rm = TRUE)
            maxRaw <- max(x, na.rm = TRUE)
        }

        z <- seq(minRaw, maxRaw, length = length(col))
        graphics::image(
            z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
            xaxt = "n", yaxt = "n"
        )
        graphics::par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- .TENETHeatmap3Scale01(as.numeric(lv), minRaw, maxRaw)
        graphics::axis(1, at = xv, labels = lv)
        if (scale == "row") {
            graphics::mtext(side = 1, "Row Z-Score", line = 2)
        } else if (scale == "column") {
            graphics::mtext(side = 1, "Column Z-Score", line = 2)
        } else {
            graphics::mtext(side = 1, KeyValueName, line = 2)
        }
        if (densityInfo == "density") {
            dens <- stats::density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- .TENETHeatmap3Scale01(dens$x, minRaw, maxRaw)
            graphics::lines(dens$x, dens$y / max(dens$y) * 0.95,
                col = denscol,
                lwd = 1
            )
            graphics::axis(
                2,
                at = pretty(dens$y) / max(dens$y) * 0.95,
                pretty(dens$y)
            )
            graphics::title("Color Key\nand Density Plot")
            graphics::par(cex = 0.5)
            graphics::mtext(side = 2, "Density", line = 2)
        } else if (densityInfo == "histogram") {
            h <- graphics::hist(x, plot = FALSE, breaks = breaks)
            hx <- .TENETHeatmap3Scale01(breaks, minRaw, maxRaw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            graphics::lines(hx, hy / max(hy) * 0.95,
                lwd = 1, type = "s",
                col = denscol
            )
            graphics::axis(2, at = pretty(hy) / max(hy) * 0.95, pretty(hy))
            graphics::title("Color Key\nand Histogram")
            graphics::par(cex = 0.5)
            graphics::mtext(side = 2, "Count", line = 2)
        } else {
            graphics::title("Color Key")
        }
    } else {
        graphics::plot.new()
    }
    retval$colorTable <- data.frame(
        low = retval$breaks[-length(retval$breaks)],
        high = retval$breaks[-1], color = retval$col
    )
    invisible(retval)
}
