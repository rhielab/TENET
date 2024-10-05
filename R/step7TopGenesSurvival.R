## Internal functions used by step7TopGenesSurvival

## Internal function to generate KM survival statistics
## for gene expression of genes or DNA methylation of sites
## and create plots if specified
.survivalFunction <- function(
    inputID, ## The gene or RE DNA methylation site's ID
    expressionOrMethylation, ## "Expression" or "Methylation"
    geneIDdf = NULL, ## Specify when expressionOrMethylation is "Expression"
    clinicalObject,
    TENETMultiAssayExperiment,
    highThreshold,
    lowThreshold,
    createPlot ## TRUE or FALSE - affects plots for KM only
    ) {
    ## If a gene was specified, get the gene name corresponding to the gene ID
    if (expressionOrMethylation == "Expression") {
        inputName <- geneIDdf[
            inputID,
            "geneName"
        ]
    } else if (expressionOrMethylation == "Methylation") {
        ## The RE DNA methylation site's name is its ID
        inputName <- inputID
    }

    ## Get the expression/methylation of the gene/RE DNA methylation site and
    ## add it to the clinical object
    if (expressionOrMethylation == "Expression") {
        clinicalObject$inputValue <- TENETMultiAssayExperiment@
        ExperimentList$expression@assays@data$expression[
            inputID,
            clinicalObject$expressionSampleNames
        ]
    } else if (expressionOrMethylation == "Methylation") {
        clinicalObject$inputValue <- TENETMultiAssayExperiment@
        ExperimentList$methylation@assays@data$methylation[
            inputID,
            clinicalObject$methylationSampleNames
        ]
    }

    ## Calculate some basic data
    controlSampleN <- nrow(
        clinicalObject[clinicalObject$sampleType == "Control", ]
    )

    caseSampleN <- nrow(
        clinicalObject[clinicalObject$sampleType == "Case", ]
    )

    ## Count the number samples with NA expression/methylation
    NACountControlInputValue <- sum(
        is.na(
            clinicalObject[
                clinicalObject$sampleType == "Control",
                "inputValue"
            ]
        )
    )

    NACountCaseInputValue <- sum(
        is.na(
            clinicalObject[
                clinicalObject$sampleType == "Case",
                "inputValue"
            ]
        )
    )

    ## Calculate mean expression/methylation for all samples
    ## that have expression/methylation
    controlMeanInputValue <- mean(
        clinicalObject[
            clinicalObject$sampleType == "Control",
            "inputValue"
        ],
        na.rm = TRUE
    )

    caseMeanInputValue <- mean(
        clinicalObject[
            clinicalObject$sampleType == "Case",
            "inputValue"
        ],
        na.rm = TRUE
    )

    ## Calculate the number of case samples remaining with missing
    ## clinical data
    NACountCaseClinical <- sum(
        !stats::complete.cases(
            clinicalObject[
                clinicalObject$sampleType == "Case",
                c("vitalStatus", "time")
            ]
        )
    )

    ## Calculate the number of case/control samples considered
    ## Case samples considered have complete expression/methylation + survival
    ## clinical data.
    ## Control samples only need expression/methylation (they aren't included
    ## in the survival analyses, only for expression/methylation statistics)
    controlPresentSampleN <- sum(
        !is.na(
            clinicalObject[
                clinicalObject$sampleType == "Control",
                "inputValue"
            ]
        )
    )

    casePresentSampleN <- sum(
        stats::complete.cases(
            clinicalObject[
                clinicalObject$sampleType == "Case",
                c("vitalStatus", "time", "inputValue")
            ]
        )
    )

    ## Create a subsetted dataset with only the case samples
    ## that have complete information for the vitalStatus, time, and
    ## expression/methylation variables
    completeCasesClinicalObject <- clinicalObject[
        clinicalObject$sampleType == "Case",
    ]

    completeCasesClinicalObject <- completeCasesClinicalObject[
        stats::complete.cases(
            completeCasesClinicalObject[
                , c("vitalStatus", "time", "inputValue")
            ]
        ),
    ]

    ## Calculate quantiles
    highCutoffQuantile <- unname(
        stats::quantile(
            completeCasesClinicalObject[, "inputValue"],
            highThreshold,
            na.rm = TRUE
        )[1]
    )

    lowCutoffQuantile <- unname(
        stats::quantile(
            completeCasesClinicalObject[, "inputValue"],
            lowThreshold,
            na.rm = TRUE
        )[1]
    )

    ## Determine if each sample is in the high, low, or intermediate quartiles
    completeCasesClinicalObject$grouping <- ifelse(
        completeCasesClinicalObject$inputValue > highCutoffQuantile,
        "High",
        ifelse(
            completeCasesClinicalObject$inputValue <= lowCutoffQuantile,
            "Low",
            "Intermediate"
        )
    )

    ## Get the counts of case samples in the three categories of
    ## expression/methylation level
    caseSampleHighN <- sum(
        completeCasesClinicalObject[, "grouping"] == "High"
    )

    caseSampleIntermediateN <- sum(
        completeCasesClinicalObject[, "grouping"] == "intermediate"
    )

    caseSampleLowN <- sum(
        completeCasesClinicalObject[, "grouping"] == "Low"
    )

    ## Get the mean expression/methylation in each group
    ## and convert the value for the Intermediate group to NA if it is NaN
    ## i.e. that group has no members (which can happen)
    caseSampleHighMean <- mean(
        completeCasesClinicalObject[
            completeCasesClinicalObject$grouping == "High",
            "inputValue"
        ],
        na.rm = TRUE
    )

    caseSampleIntermediateMean <- mean(
        completeCasesClinicalObject[
            completeCasesClinicalObject$grouping == "Intermediate",
            "inputValue"
        ],
        na.rm = TRUE
    )

    caseSampleLowMean <- mean(
        completeCasesClinicalObject[
            completeCasesClinicalObject$grouping == "Low",
            "inputValue"
        ],
        na.rm = TRUE
    )

    ## Calculate the proportion of samples which have reached the event
    ## rather than not/censored in the groupings
    proportionEventHigh <- as.numeric(
        nrow(
            completeCasesClinicalObject[
                completeCasesClinicalObject$grouping == "High" &
                    completeCasesClinicalObject$vitalStatus == 2,
            ]
        ) /
            nrow(
                completeCasesClinicalObject[
                    completeCasesClinicalObject$grouping == "High",
                ]
            )
    )

    proportionEventIntermediate <- as.numeric(
        nrow(
            completeCasesClinicalObject[
                completeCasesClinicalObject$grouping == "Intermediate" &
                    completeCasesClinicalObject$vitalStatus == 2,
            ]
        ) /
            nrow(
                completeCasesClinicalObject[
                    completeCasesClinicalObject$grouping == "Intermediate",
                ]
            )
    )

    proportionEventLow <- as.numeric(
        nrow(
            completeCasesClinicalObject[
                completeCasesClinicalObject$grouping == "Low" &
                    completeCasesClinicalObject$vitalStatus == 2,
            ]
        ) /
            nrow(
                completeCasesClinicalObject[
                    completeCasesClinicalObject$grouping == "Low",
                ]
            )
    )

    ## Convert NaN values for intermediate group to NA
    ## (Which can happen as group can be empty)
    if (is.nan(caseSampleIntermediateMean)) {
        caseSampleIntermediateMean <- NA
        proportionEventIntermediate <- NA
    }

    ## Determine which high vs low group had greater proportion of reaching
    ## the event of interest
    if (expressionOrMethylation == "Expression") {
        highestEventProportionGroup <- ifelse(
            proportionEventHigh > proportionEventLow,
            "highExpressionLowSurvival",
            ifelse(
                proportionEventHigh < proportionEventLow,
                "lowExpressionLowSurvival",
                "unclear"
            )
        )
    } else if (expressionOrMethylation == "Methylation") {
        highestEventProportionGroup <- ifelse(
            proportionEventHigh > proportionEventLow,
            "highMethylationLowSurvival",
            ifelse(
                proportionEventHigh < proportionEventLow,
                "lowMethylationLowSurvival",
                "unclear"
            )
        )
    }

    ## Start with KM statistics first

    ## Create survival objects for KM analyses
    ## KM will include only High/Low samples
    KMSurvivalObject <- survival::Surv(
        completeCasesClinicalObject[
            completeCasesClinicalObject$grouping %in% c("High", "Low"),
            "time"
        ],
        completeCasesClinicalObject[
            completeCasesClinicalObject$grouping %in% c("High", "Low"),
            "vitalStatus"
        ]
    )

    rownames(KMSurvivalObject) <- rownames(
        completeCasesClinicalObject[
            completeCasesClinicalObject$grouping %in% c("High", "Low"),
        ]
    )

    ## Get the grouping of the samples without any
    ## Intermediate samples (to go with our KM survival object)
    inputValueGroup <- completeCasesClinicalObject[
        !completeCasesClinicalObject$grouping == "Intermediate",
        "grouping"
    ]

    ## Perform the survival analysis.
    ## This uses the expression grouping as the x variable
    ## and the KM survival object as the y variable to create
    ## a table with information about the test
    ## including chi-squared p-value
    KMSurvivalTable <- survival::survdiff(KMSurvivalObject ~ inputValueGroup)

    ## Get the chi-squared test statistic from the analysis above
    KMChiSquared <- KMSurvivalTable$chisq

    ## Calculate a p-value based on the test statistic to get
    ## a precise p-value for KM analysis
    KMSurvivalPvalue <- as.numeric(
        1 - stats::pchisq(abs(KMChiSquared), df = 1)
    )

    ## Create and return KM plots if createPlot is TRUE
    ## Otherwise, return vector of statistics for both KM and Cox
    ## to later combine into a data frame
    if (createPlot) {
        ## Create legend names for high/low expression groups
        legendNameHigh <- paste(inputName, "high")
        legendNameLow <- paste(inputName, "low")

        ## Round the p-value displayed the graph to 3 digits
        KMSurvivalPvalueFormatted <- formatC(
            KMSurvivalPvalue,
            format = "e",
            digits = 3
        )

        ## Create the plot title
        ## with gene name and p-value included
        ## If expression is specified, include the gene name and ESNG.
        ## If it's methylation, include only the methylation site ID
        if (expressionOrMethylation == "Expression") {
            KMSurvivalTitle <- paste0(
                inputName,
                " - ",
                inputID,
                "\nKaplan-Meier Survival analysis\np = ",
                KMSurvivalPvalueFormatted
            )
        } else if (expressionOrMethylation == "Methylation") {
            KMSurvivalTitle <- paste0(
                inputName,
                "\nKaplan-Meier Survival analysis\np = ",
                KMSurvivalPvalueFormatted
            )
        }

        ## Change margins to increase size at top of plot for title
        ## Mar: bottom, left, top and right
        graphics::par("mar" = c(5, 4, 4.75, 2))

        ## Actually create the survival plot
        ## Using a similar structure to what was used to generate the p-value
        .newInvisibleRecordablePlot()
        plot(
            survival::survfit(KMSurvivalObject ~ inputValueGroup),

            ## Color the lines (high in red first!)
            col = c("red", "black"),

            ## Add thickness to the lines
            lwd = 3,

            ## Use the title that was created earlier as the title of the plot
            main = KMSurvivalTitle,

            ## Set titles of the x and y axis
            ## Note: TCGA measures survival in days as noted
            xlab = "Time",
            ylab = "Survival Proportion",

            ## Change axis size
            cex.axis = 1,
            cex.main = 1.5,
            cex.lab = 1.25
        )

        ## Add a legend to the plot
        graphics::legend(

            ## Set X position of legend in graph
            x = max(completeCasesClinicalObject$time) * (2 / 3),

            ## Set Y position of legend in graph
            y = 1,

            ## Use the legend titles that were created earlier
            legend = c(legendNameHigh, legendNameLow),

            ## As above, use black for low and red for high
            col = c("red", "black"),

            ## Coloring the text in the legend as well
            text.col = c("red", "black"),
            cex = 1.5,

            ## Change the shape of the labels in the legend
            pch = 15
        )

        ## Save the plot to an object
        KMSurvivalPlotObject <- .recordTENETSavedSizePlot()

        ## Close the plot
        grDevices::dev.off()

        ## Return the plot
        return(KMSurvivalPlotObject)
    } else {
        ## Do the Cox analysis
        ## Note: This will only be done if createPlot is set to FALSE
        ## since if createPlot is set to TRUE only the KM plots will be output
        ## without any statistics, and no plots are generated from Cox analyses

        ## Create survival objects for Cox analyses. Cox uses all samples.
        CoxSurvivalObject <- survival::Surv(
            completeCasesClinicalObject$time,
            completeCasesClinicalObject$vitalStatus
        )

        rownames(CoxSurvivalObject) <- rownames(
            completeCasesClinicalObject
        )

        ## We can use the CoxSurvivalObject along with the
        ## expression/methylation of our gene/RE DNA methylation site of
        ## interest to do the Cox regression survival analysis
        coxRegressionResults <- survival::coxph(
            CoxSurvivalObject ~ inputValue,
            data = completeCasesClinicalObject
        )

        ## Assemble a vector of results with information relevant to
        ## the given gene/RE DNA methylation site

        namesTemplate <- c(
            "controlSampleCount",
            "caseSampleCount",
            "controlSampleCount@TYPE@Missing",
            "caseSampleCount@TYPE@Missing",
            "controlMean@TYPE@",
            "caseMean@TYPE@",
            "caseSampleCountClinicalMissing",
            "controlSampleCountWithData",
            "caseSampleCountWithData",
            "caseSampleCountHigh@TYPE@Group",
            "caseSampleCountIntermediate@TYPE@Group",
            "caseSampleCountLow@TYPE@Group",
            "caseMean@TYPE@High@TYPE@Group",
            "caseMean@TYPE@Intermediate@TYPE@Group",
            "caseMean@TYPE@Low@TYPE@Group",
            "caseProportionEventHigh@TYPE@Group",
            "caseProportionEventIntermediate@TYPE@Group",
            "caseProportionEventLow@TYPE@Group",
            "KMSurvivalDirectionOfEffect",
            "KMSurvivalPValue",
            "CoxRegressionCoefficient",
            "CoxHazardRatio",
            "CoxSurvivalDirectionOfEffect",
            "CoxSurvivalPValue"
        )

        survivalReturnVector <- c(
            controlSampleN,
            caseSampleN,
            NACountControlInputValue,
            NACountCaseInputValue,
            controlMeanInputValue,
            caseMeanInputValue,
            NACountCaseClinical,
            controlPresentSampleN,
            casePresentSampleN,
            caseSampleHighN,
            caseSampleIntermediateN,
            caseSampleLowN,
            caseSampleHighMean,
            caseSampleIntermediateMean,
            caseSampleLowMean,
            proportionEventHigh,
            proportionEventIntermediate,
            proportionEventLow,
            highestEventProportionGroup,
            as.numeric(KMSurvivalPvalue),
            as.numeric(unname(coxRegressionResults$coefficients)),
            as.numeric(unname(exp(coxRegressionResults$coefficients))),
            ifelse(
                as.numeric(
                    unname(exp(coxRegressionResults$coefficients))
                ) > 1,
                paste0("high", expressionOrMethylation, "LowSurvival"),
                ifelse(
                    as.numeric(
                        unname(
                            exp(coxRegressionResults$coefficients)
                        )
                    ) < 1,
                    paste0("low", expressionOrMethylation, "LowSurvival"),
                    "unclear"
                )
            ),
            as.numeric(summary(coxRegressionResults)$coefficients[, 5])
        )

        if (expressionOrMethylation == "Expression") {
            ## For expression, add gene ID and gene name columns
            survivalReturnVector <- c(
                inputID,
                inputName,
                survivalReturnVector
            )

            namesTemplate <- c(
                "geneID",
                "geneName",
                namesTemplate
            )
        }

        ## Add names to the return values
        names(survivalReturnVector) <- gsub(
            "@TYPE@", expressionOrMethylation, namesTemplate
        )

        ## Return the statistics
        return(survivalReturnVector)
    }
}

## Internal function to return survival statistics or graphs for a given
## quadrant
.returnSurvivalStatisticsOrGraphs <- function(
    hyperHypo,
    geneIDdf,
    clinicalObject,
    TENETMultiAssayExperiment,
    topGeneNumber,
    highThreshold,
    lowThreshold,
    geneOrTF, ## Return info for top genes ("Gene") or TFs ("TF")
    ## Return results for genes ("Genes") or RE DNA methylation sites linked to
    ## genes ("DNAMethylationSites")
    genesOrMethSites,
    statsOrPlots, ## Return stats ("Stats") or plots ("Plots")
    coreCount) {
    ## Generate the quadrant result name to grab data for
    quadrantResultsName <- paste0(hyperHypo, "methGplusResults")

    ## Ensure the quadrant's results are present in step 5
    .ensureStepPresent(
        TENETMultiAssayExperiment,
        stepName = "step5OptimizeLinks",
        substepName = quadrantResultsName
    )

    ## Get the IDs of the top genes/TFs. If there are fewer genes/TFs than the
    ## topGeneNumber specified by the user, get all the genes/TFs available.
    topQuadrantGeneOrTFIDs <- .getQuadrantTopGenesOrTFs(
        TENETMultiAssayExperiment = TENETMultiAssayExperiment,
        geneOrTF = geneOrTF,
        hyperHypo = hyperHypo,
        topGeneNumber = topGeneNumber
    )$geneID
    if (.isSingleNA(topQuadrantGeneOrTFIDs)) {
        return(NA)
    }

    ## Convert the gene IDs to gene names
    topQuadrantGeneName <- geneIDdf[topQuadrantGeneOrTFIDs, "geneName"]

    ## Get all unique RE DNA methylation sites linked to at least one of the top
    ## genes selected
    quadrantMethSitesLinkedToSignificantGenes <- unique(
        TENETMultiAssayExperiment@metadata$step5OptimizeLinks[[
            quadrantResultsName
        ]][
            TENETMultiAssayExperiment@metadata$step5OptimizeLinks[[
                quadrantResultsName
            ]]$geneID %in% topQuadrantGeneOrTFIDs,
            "DNAMethylationSiteID"
        ]
    )

    quadrantMethSitesLinkedToSignificantGenes <- sort(
        quadrantMethSitesLinkedToSignificantGenes
    )

    ## Perform and return analyses for genes
    if (genesOrMethSites == "Genes") {
        if (statsOrPlots == "Stats") {
            ## Return survival statistics for genes
            returnValue <- as.data.frame(
                do.call(
                    rbind,
                    parallel::mclapply(
                        X = topQuadrantGeneOrTFIDs,
                        FUN = .survivalFunction,
                        expressionOrMethylation = "Expression",
                        geneIDdf = geneIDdf,
                        clinicalObject = clinicalObject,
                        TENETMultiAssayExperiment = TENETMultiAssayExperiment,
                        highThreshold = highThreshold,
                        lowThreshold = lowThreshold,
                        createPlot = FALSE,
                        mc.cores = coreCount
                    )
                )
            )

            rownames(returnValue) <- topQuadrantGeneOrTFIDs

            return(returnValue)
        } else {
            ## Return survival plots for genes
            returnValue <- parallel::mclapply(
                X = topQuadrantGeneOrTFIDs,
                FUN = .survivalFunction,
                expressionOrMethylation = "Expression",
                geneIDdf = geneIDdf,
                clinicalObject = clinicalObject,
                TENETMultiAssayExperiment = TENETMultiAssayExperiment,
                highThreshold = highThreshold,
                lowThreshold = lowThreshold,
                createPlot = TRUE,
                mc.cores = coreCount
            )

            names(returnValue) <- topQuadrantGeneOrTFIDs

            return(returnValue)
        }
    } else {
        if (statsOrPlots == "Stats") {
            ## For RE DNA methylation sites, initialize a data frame with the
            ## methylation site IDs
            returnValue <- data.frame(
                "DNAMethylationSiteID" =
                    quadrantMethSitesLinkedToSignificantGenes,
                stringsAsFactors = FALSE
            )

            ## Add columns to that data frame indicating which of the RE DNA
            ## methylation sites is linked to each of the top genes
            for (i in seq_along(topQuadrantGeneOrTFIDs)) {
                ## Identify if the quadrantMethSitesLinkedToSignificantGenes
                ## are among RE DNA methylation sites linked to the specific
                ## gene of interest
                TFVector <- quadrantMethSitesLinkedToSignificantGenes %in%
                    TENETMultiAssayExperiment@metadata$step5OptimizeLinks[[
                        quadrantResultsName
                    ]][
                        TENETMultiAssayExperiment@
                        metadata$step5OptimizeLinks[[
                            quadrantResultsName
                        ]]$geneID %in% topQuadrantGeneOrTFIDs[i],
                        "DNAMethylationSiteID"
                    ]

                returnValue[i + 1] <- TFVector
            }

            ## Reset the colnames of the DF as the methylation site IDs,
            ## then the combined gene names and IDs
            colnames(returnValue) <- c(
                "DNAMethylationSiteID",
                paste(
                    topQuadrantGeneName,
                    topQuadrantGeneOrTFIDs,
                    "linked",
                    sep = "_"
                )
            )

            ## Return survival statistics for RE DNA methylation sites linked to
            ## top genes
            returnValue2 <- as.data.frame(
                do.call(
                    rbind,
                    parallel::mclapply(
                        X = quadrantMethSitesLinkedToSignificantGenes,
                        FUN = .survivalFunction,
                        expressionOrMethylation = "Methylation",
                        clinicalObject = clinicalObject,
                        TENETMultiAssayExperiment = TENETMultiAssayExperiment,
                        highThreshold = highThreshold,
                        lowThreshold = lowThreshold,
                        createPlot = FALSE,
                        mc.cores = coreCount
                    )
                )
            )

            ## Combine the two return data frames before doing the final return
            returnValue <- cbind(returnValue, returnValue2)

            rownames(returnValue) <-
                quadrantMethSitesLinkedToSignificantGenes

            return(returnValue)
        } else {
            ## Return survival plots for RE DNA methylation sites linked to top
            ## genes
            returnValue <- parallel::mclapply(
                X = quadrantMethSitesLinkedToSignificantGenes,
                FUN = .survivalFunction,
                expressionOrMethylation = "Methylation",
                clinicalObject = clinicalObject,
                TENETMultiAssayExperiment = TENETMultiAssayExperiment,
                highThreshold = highThreshold,
                lowThreshold = lowThreshold,
                createPlot = TRUE,
                mc.cores = coreCount
            )

            names(returnValue) <- quadrantMethSitesLinkedToSignificantGenes

            return(returnValue)
        }
    }
}

#' Perform Kaplan-Meier and Cox regression analyses to assess the association
#' of top gene expression and linked RE DNA methylation site methylation with
#' patient survival
#'
#' This function takes the top genes/transcription factors (TFs) by number of
#' linked RE DNA methylation sites identified by the
#' `step6DNAMethylationSitesPerGeneTabulation` function up to the number
#' specified by the user and generates survival plots and tables with statistics
#' from survival analyses assessing the survival association of the expression
#' level of each gene as well as the methylation level of each RE DNA
#' methylation site linked to them, using percentile cutoffs as specified by the
#' user for Kaplan-Meier analyses.
#'
#' @param TENETMultiAssayExperiment Specify a MultiAssayExperiment object
#' containing expression and methylation SummarizedExperiment objects, such as
#' one created by the TCGADownloader function. This MultiAssayExperiment object
#' should also contain the results from the
#' `step2GetDifferentiallyMethylatedSites`, `step5OptimizeLinks`, and
#' `step6DNAMethylationSitesPerGeneTabulation` functions in
#' its metadata. Additionally, the colData of the MultiAssay object must contain
#' a 'vital_status' and 'time' column, containing data on the patients' survival
#' status and time to event/censorship, respectively.
#' @param geneAnnotationDataset Specify a gene annotation dataset which is
#' used to identify names for genes by their Ensembl IDs. The argument must be
#' either a GRanges object (such as one imported via `rtracklayer::import`) or a
#' path to a GFF3 or GTF file. Both GENCODE and Ensembl annotations are
#' supported. Other annotation datasets may work, but have not been tested.
#' See the "Input data" section of the vignette for information on the required
#' dataset format.
#' Specify NA to use the names for genes listed in the "geneName" column of the
#' elementMetadata of the rowRanges of the "expression" SummarizedExperiment
#' object within the TENETMultiAssayExperiment object. Defaults to NA.
#' @param hypermethGplusAnalysis Set to TRUE to perform survival analyses on
#' top genes/TFs by most hypermethylated RE DNA methylation sites with G+ links,
#' as well as their linked RE DNA methylation sites.
#' @param hypomethGplusAnalysis Set to TRUE to perform survival analyses on
#' top genes/TFs by most hypomethylated RE DNA methylation sites with G+ links,
#' as well as their linked RE DNA methylation sites.
#' @param topGeneNumber Specify the number of top genes/TFs, based on the most
#' linked RE DNA methylation sites of a given analysis type, for which to
#' perform survival analyses. Defaults to 10.
#' @param vitalStatusData Specify the vital status data for samples in the
#' TENETMultiAssayExperiment. Vital status should be given in the form of
#' either "alive" or "dead" (case-insensitive), or 1 or 2, indicating that the
#' sample was collected from a patient who was alive/censored or dead/reached
#' the outcome of interest, respectively. These
#' data can be given in a variety of forms, including a vector, data frame/
#' matrix, or a path to a file that contains the vital status data. If a vector
#' is given, the names of the vector elements must correspond to the names of
#' the samples in the rownames of the colData of the TENETMultiAssayExperiment
#' object. If no names are provided for the vector, then the number of elements
#' in the vector must equal the number of samples in the colData, and are
#' assumed to align with the samples as they are ordered in the colData. If a
#' data frame or matrix is given, its rownames must include the sample names as
#' they appear in the colData of the TENETMultiAssayExperiment object, and its
#' first column must include the vital status data. If a single string
#' is provided, then it is assumed to be a path to a tab-delimited file
#' containing vital status data in the second column, and the names of the
#' samples, again corresponding with the sample names in the colData, in the
#' first column, which will be loaded as the row names. The first row of the
#' file must contain column names. If this variable is set to NA, then
#' the vital status data will be assumed to already be contained in the colData
#' of the TENETMultiAssayExperiment under a column titled "vital_status".
#' Defaults to NA.
#' @param survivalTimeData Specify the survival time data for samples in the
#' TENETMultiAssayExperiment. Survival time should be given in the form of a
#' numeric variable. These data can be given in a variety of forms, including a
#' vector, data frame/matrix, or a path to a file that contains the survival
#' time data. See the documentation for the vitalStatusData argument for more
#' information. If this variable is set to NA, then the survival time
#' data will be assumed to already be contained in the colData of the
#' TENETMultiAssayExperiment under a column titled "time". Defaults to NA.
#' @param highProportion Set a number ranging from 0 to 1, indicating the
#' proportion of all samples to include in the high expression/methylation
#' group for Kaplan-Meier survival analyses. Defaults to 0.5.
#' @param lowProportion Set a number ranging from 0 to 1, indicating the
#' proportion of all samples to include in the low expression/methylation
#' group for Kaplan-Meier survival analyses. The total value of the
#' highProportion and lowProportion arguments should not exceed 1. **Note:**
#' If both this value and the highProportion value are set to 0.5, samples at
#' exactly the 50th percentile will be assigned to the "Low" group.
#' Defaults to 0.5.
#' @param generatePlots Set to TRUE to create and save plots displaying the
#' Kaplan-Meier survival results for the genes/TFs of interest, as well as the
#' RE DNA methylation sites linked to them. Defaults to TRUE.
#' @param coreCount Argument passed as the mc.cores argument to mclapply. See
#' `?parallel::mclapply` for more details. Defaults to 1.
#' @return Returns the MultiAssayExperiment object given as the
#' TENETMultiAssayExperiment argument with an additional list of information
#' named 'step7TopGenesSurvival' in its metadata with the output of this
#' function. This list is subdivided into hypermethGplus or hypomethGplus
#' results as selected by the user, which are further subdivided into lists with
#' data for the top overall genes, and for top TF genes only. Each contains a
#' list with data frames containing survival statistics for the top genes/TFs as
#' well as their linked RE DNA methylation sites, from both Kaplan-Meier and Cox
#' regression analyses. Additionally, they will also include a list of
#' Kaplan-Meier plots if generatePlots is TRUE.
#' @export
#'
#' @examplesIf interactive()
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to perform Kaplan-Meier and Cox regression
#' ## survival analyses for the top 10 genes/TFs, by number of linked hyper- or
#' ## hypomethylated RE DNA methylation sites, as well as for all unique RE DNA
#' ## methylation sites linked to those
#' ## 10 genes/TFs. The vital status and survival time of patients will be
#' ## taken from the "vital_status" and "time" columns present in the colData of
#' ## the example MultiAssayExperiment. Gene names will be retrieved from the
#' ## rowRanges of the 'expression' SummarizedExperiment object in the example
#' ## MultiAssayExperiment. For Kaplan-Meier analyses, the patient samples with
#' ## complete clinical information in the highest half of
#' ## expression/methylation will be compared to the patient samples with
#' ## complete clinical information in the lowest half. Kaplan-Meier plots will
#' ## be saved for the genes and RE DNA methylation sites, and the analysis will
#' ## be performed using one CPU core.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Use the example dataset to do the survival analysis
#' returnValue <- step7TopGenesSurvival(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment
#' )
#'
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to perform Kaplan-Meier and Cox regression
#' ## survival analyses for only the top 5 genes/TFs, by number of linked
#' ## hypomethylated RE DNA methylation sites, as well as for all unique
#' ## RE DNA methylation sites linked to those 5 genes/TFs only. The vital
#' ## status and survival time of patients will be
#' ## taken from specific columns in a separate data frame with example patient
#' ## data from the TENET.ExperimentHub package. Gene names will be retrieved
#' ## from the rowRanges of the 'expression' SummarizedExperiment object in the
#' ## example MultiAssayExperiment. For Kaplan-Meier analyses, the patient
#' ## samples with complete clinical information in the highest quartile of
#' ## expression/methylation will be compared to the patient samples with
#' ## complete clinical information in the lowest quartile. Kaplan-Meier plots
#' ## will not be saved for the genes and RE DNA methylation sites, and the
#' ## analysis will be performed using 8 CPU cores.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Again, this loads the data frame with example clinical data for patients
#' ## in the TENET MultiAssayExperiment object from the TENET.ExperimentHub
#' ## package
#' exampleTENETClinicalDataFrame <-
#'     TENET.ExperimentHub::exampleTENETClinicalDataFrame()
#'
#' ## Use the example datasets to do the survival analysis
#' returnValue <- step7TopGenesSurvival(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     hypermethGplusAnalysis = FALSE,
#'     topGeneNumber = 5,
#'     vitalStatusData = exampleTENETClinicalDataFrame["vital_status"],
#'     survivalTimeData = exampleTENETClinicalDataFrame["time"],
#'     highProportion = 0.25,
#'     lowProportion = 0.25,
#'     generatePlots = FALSE,
#'     coreCount = 8
#' )
step7TopGenesSurvival <- function(
    TENETMultiAssayExperiment,
    geneAnnotationDataset = NA,
    hypermethGplusAnalysis = TRUE,
    hypomethGplusAnalysis = TRUE,
    topGeneNumber = 10,
    vitalStatusData = NA,
    survivalTimeData = NA,
    highProportion = 0.5,
    lowProportion = 0.5,
    generatePlots = TRUE,
    coreCount = 1) {
    ## Validate the analysis types and get a vector of the ones selected
    analysisTypes <- .validateAnalysisTypes(
        hypermethGplusAnalysis, hypomethGplusAnalysis
    )

    ## Return an error message if the input MultiAssayExperiment is invalid
    .validateMultiAssayExperiment(
        TENETMultiAssayExperiment,
        needGeneName = is.na(geneAnnotationDataset)
    )

    ## Check for nonsensical high/lowProportion values which would cause invalid
    ## results
    if (!is.numeric(highProportion) || !is.numeric(lowProportion)) {
        .stopNoCall(
            "Invalid highProportion and/or lowProportion specified. ",
            "Both must be numeric values between 0 and 1."
        )
    }

    if ((highProportion + lowProportion) > 1 ||
        highProportion <= 0 || lowProportion <= 0
    ) {
        .stopNoCall(
            "Invalid highProportion and/or lowProportion specified. ",
            "Both must be positive, and their sum may not be greater than 1."
        )
    }

    ## Process the status data of the samples if vitalStatusData is not NA.
    ## Otherwise, look into the supplied colData of the MultiAssayExperiment
    vitalStatusResults <- .importClinicalData(
        userInput = vitalStatusData,
        argumentName = "vitalStatusData",
        clinicalDataColumn = "vital_status",
        returnType = "single",
        TENETMultiAssayExperiment = TENETMultiAssayExperiment
    )

    ## Process the survival time data of the samples if vitalStatusData is
    ## not NA. Otherwise, look into the supplied colData of the
    ## MultiAssayExperiment
    survivalTimeResults <- .importClinicalData(
        userInput = survivalTimeData,
        argumentName = "survivalTimeData",
        clinicalDataColumn = "time",
        returnType = "single",
        TENETMultiAssayExperiment = TENETMultiAssayExperiment
    )

    ## Get gene IDs and names from the MAE, or gene annotation dataset if
    ## provided
    geneIDdf <- .getGeneIDsAndNames(
        TENETMultiAssayExperiment, geneAnnotationDataset
    )

    ## Set thresholds based on specified fractions.
    ## 1/3 is specified, the top 1/3 of the samples should be returned
    highThresh <- 1 - highProportion

    ## Get the names of the control and case samples in the
    ## methylation data first
    methylationSampleNames <- .getExpOrMetSamplesOfType(
        TENETMultiAssayExperiment,
        "methylation",
        namesOnly = TRUE
    )

    ## Get the methylation values that match with expression values
    ## using the mapping data. This assumes the methylation and expression
    ## values share a clinical data match within the mapping.
    metToExpSampleConversion <- .createMetToExpSampleConversionVector(
        TENETMultiAssayExperiment
    )

    ## Use the conversion vector to get the names of the control and case
    ## samples in the expression data which pair with the
    ## methylationSampleNames
    expressionSampleNames <- metToExpSampleConversion[
        methylationSampleNames
    ]

    ## Start creating a clinical data frame by matching the methylation
    ## and expression sample names with their primary names from the sampleMap
    ## of the TENETMultiAssayExperiment
    clinicalDF <- data.frame(
        "methylationSampleNames" = methylationSampleNames,
        "expressionSampleNames" = expressionSampleNames,
        "primarySampleNames" = TENETMultiAssayExperiment@sampleMap[
            match(
                methylationSampleNames,
                TENETMultiAssayExperiment@sampleMap$colname
            ),
            "primary"
        ],
        "sampleType" = TENETMultiAssayExperiment@sampleMap[
            match(
                methylationSampleNames,
                TENETMultiAssayExperiment@sampleMap$colname
            ),
            "sampleType"
        ],
        stringsAsFactors = FALSE
    )

    ## Add the vitalStatus and time variables that were imported.
    ## Convert vitalStatus to lowercase to ease matching below.
    clinicalDF$vitalStatus <- tolower(vitalStatusResults[
        clinicalDF$primarySampleNames,
        "vital_status"
    ])

    clinicalDF$time <- survivalTimeResults[
        clinicalDF$primarySampleNames,
        "time"
    ]

    ## Format the vitalStatus and time columns, inducing NA values where
    ## samples were not input properly

    ## Format vitalStatus. There are three acceptable formats for the values -
    ## "alive"/"dead", 1/2, and "1"/"2" - which we will harmonize to 1/2
    ## (numeric, for alive/censored and dead, respectively).
    ## Note: 1 (numeric) is equal to 1 (as character) - 1=="1" returns TRUE.
    vitalStatusConversionDF <- data.frame(
        "values" = c("alive", "1", "dead", "2"),
        "return" = c(1, 1, 2, 2),
        stringsAsFactors = FALSE
    )
    rownames(vitalStatusConversionDF) <- vitalStatusConversionDF$values

    clinicalDF$vitalStatus <- ifelse(
        clinicalDF$vitalStatus %in% vitalStatusConversionDF$values,
        vitalStatusConversionDF[
            as.character(clinicalDF$vitalStatus),
            "return"
        ],
        NA
    )

    ## Time should be a numeric variable; other values should be NA.
    ## suppressWarnings() is used here, as it will return a warning
    ## of "NAs introduced by coercion" when non-numeric values are
    ## converted to NA, but that is the intended use of this function.
    clinicalDF$time <- suppressWarnings(as.numeric(clinicalDF$time))

    ## Create an empty list to hold the results from this step 7 function
    resultsList <- list()

    ## Do the analysis for the analysis types selected by the user
    for (hyperHypo in analysisTypes) {
        ## Return results for all genes then TFs for each analysis type
        for (geneOrTF in c("Gene", "TF")) {
            ## Return results for genes then RE DNA methylation sites for all
            ## genes and TFs
            for (genesOrMethSites in c("Genes", "DNAMethylationSites")) {
                ## Finally return statistics then plots for genes and RE DNA
                ## methylation sites
                for (statsOrPlots in c("Stats", "Plots")) {
                    resultsList[[
                        paste0(hyperHypo, "methGplusResults")
                    ]][[
                        paste0("top", geneOrTF, "s")
                    ]][[
                        paste0(
                            "top",
                            genesOrMethSites,
                            "Survival",
                            statsOrPlots
                        )
                    ]] <- .returnSurvivalStatisticsOrGraphs(
                        hyperHypo = hyperHypo,
                        geneIDdf = geneIDdf,
                        clinicalObject = clinicalDF,
                        TENETMultiAssayExperiment = TENETMultiAssayExperiment,
                        topGeneNumber = topGeneNumber,
                        highThreshold = highThresh,
                        lowThreshold = lowProportion,
                        geneOrTF = geneOrTF,
                        genesOrMethSites = genesOrMethSites,
                        statsOrPlots = statsOrPlots,
                        coreCount = coreCount
                    )
                }
            }
        }
    }

    ## Add the results list to the MultiAssayExperiment object
    TENETMultiAssayExperiment@metadata$step7TopGenesSurvival <- resultsList

    return(TENETMultiAssayExperiment)
}
