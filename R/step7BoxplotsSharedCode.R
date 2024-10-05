## Internal functions used by the step 7 boxplot functions

## Internal function that, when given a list of genes or RE DNA methylation
## sites, will plot a boxplot showing the expression/methylation level of that
## gene/RE DNA methylation site in the case vs control samples
.quadrantBoxplotFunction <- function(
    geneOrMethSiteID,
    expOrMet,
    expOrMetData,
    geneIDNameDF = NA,
    groupInfo) {
    ## Convert the gene ID into the gene name, if genes are being analyzed,
    ## and set the plot title and result name accordingly
    if (expOrMet == "expression") {
        displayName <- geneIDNameDF[geneOrMethSiteID, "geneName"]
        plotTitle <- paste(displayName, "-", geneOrMethSiteID)
    } else {
        displayName <- geneOrMethSiteID
        plotTitle <- displayName
    }

    ## Get expression/methylation values
    unlistedExpOrMetData <- c(
        unlist(expOrMetData[geneOrMethSiteID, groupInfo$group])
    )

    ## Assemble a new data frame with expression/methylation and sample type
    boxplotDF <- data.frame(
        "expOrMetValues" = unlistedExpOrMetData,
        "sampleType" = groupInfo$cluster
    )

    ## Manually coloring samples - blue for control, red for case data points
    groupColors <- c("Control" = "dodgerblue3", "Case" = "red3")

    ## Student's t-test
    if (
        sum(!is.na(boxplotDF[
            boxplotDF$sampleType == "Case", "expOrMetValues"
        ])) > 1 &
            sum(!is.na(boxplotDF[
                boxplotDF$sampleType == "Control", "expOrMetValues"
            ])) > 1
    ) {
        ## Do the t-test properly since there are enough samples present
        tResults <- stats::t.test(
            boxplotDF[
                boxplotDF$sampleType == "Case", "expOrMetValues"
            ],
            boxplotDF[
                boxplotDF$sampleType == "Control", "expOrMetValues"
            ]
        )

        ## Truncating the t-test p-value for display in plot (rounding to 3
        ## digits)
        simpleTPValue <- formatC(tResults$p.value, format = "e", digits = 3)
    } else {
        ## t-test can't be done, as one group has only 0 or 1 valid samples
        simpleTPValue <- NA
    }

    ## Create basic boxplot. Have to use get because R CMD check does not
    ## understand lazy evaluation
    boxplot <- ggplot2::ggplot(
        boxplotDF,
        ggplot2::aes(x = get("sampleType"), y = get("expOrMetValues"))
    ) +
        ggplot2::geom_boxplot(ggplot2::aes(fill = get("sampleType"))) +
        ggplot2::ggtitle(
            paste0(plotTitle, "\nt-test p = ", simpleTPValue)
        ) +
        ggplot2::ylab(paste(displayName, expOrMet)) +
        ggplot2::xlab("Sample Grouping") +
        ggplot2::guides(fill = "none") +
        ggplot2::theme_bw() +
        ggplot2::scale_fill_manual(values = groupColors) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
            panel.border = ggplot2::element_rect(
                color = "black", fill = NA, linewidth = 1
            ),
            plot.background = ggplot2::element_rect(fill = "white"),
            axis.title.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(size = 18, color = "black"),
            axis.text.y = ggplot2::element_text(size = 16, color = "black"),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank()
        )

    ## To display this plot correctly, the plot area should be at least 10
    ## inches wide by 7 inches tall
    boxplot <- TENETSavedSizePlot(boxplot, width = 10, height = 7)

    ## To name the returned plot, return it as a 1-element named list
    ## The name of the plot in the metadata list will just be the ID
    namedResult <- list("boxplot" = boxplot)
    names(namedResult) <- geneOrMethSiteID
    return(namedResult)
}
