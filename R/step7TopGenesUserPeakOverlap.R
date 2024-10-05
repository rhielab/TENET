## Internal functions used by step7TopGenesUserPeakOverlap

## GRanges overlap function
## Assumes methSiteGRanges has methylation site IDs in the names()
.methSitePeakGRangesOverlapTFFunction <- function(
    methSiteGRanges,
    peakFileGRanges) {
    overlapNames <- names(
        methSiteGRanges[
            unique(
                S4Vectors::queryHits(
                    GenomicRanges::findOverlaps(
                        methSiteGRanges,
                        peakFileGRanges
                    )
                )
            ),
        ]
    )

    return(names(methSiteGRanges) %in% overlapNames)
}

## For the top genes/TFs, get the RE DNA methylation sites linked to them,
## overlap them with each of the files, then create files noting which datasets
## the RE DNA methylation sites overlapped with
.userPeakOverlapInternalDatasetOutputFunction <- function(
    TENETMultiAssayExperiment,
    hyperHypo,
    geneOrTF,
    geneIDdf,
    methSiteIDdf,
    peakFileGRangesList,
    topGeneNumber,
    distanceFromREDNAMethylationSites,
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
    topQuadrantGeneName <- geneIDdf[
        topQuadrantGeneOrTFIDs,
        "geneName"
    ]

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

    ## Initialize a data frame with the methylation site IDs, info about
    ## the RE DNA methylation sites, and the search start and end coordinates
    methSitesLinkedToGenesDF <- data.frame(
        "DNAMethylationSiteID" = quadrantMethSitesLinkedToSignificantGenes,
        "chromosome" = methSiteIDdf[
            quadrantMethSitesLinkedToSignificantGenes,
            "chromosome"
        ],
        "start" = methSiteIDdf[
            quadrantMethSitesLinkedToSignificantGenes,
            "start"
        ],
        "end" = methSiteIDdf[
            quadrantMethSitesLinkedToSignificantGenes,
            "end"
        ],
        "searchStart" = methSiteIDdf[
            quadrantMethSitesLinkedToSignificantGenes,
            "start"
        ] - distanceFromREDNAMethylationSites,
        "searchEnd" = methSiteIDdf[
            quadrantMethSitesLinkedToSignificantGenes,
            "end"
        ] + distanceFromREDNAMethylationSites,
        stringsAsFactors = FALSE
    )

    ## Add columns to that data frame indicating
    ## which of the RE DNA methylation sites is linked to each of the top genes
    for (i in seq_along(topQuadrantGeneOrTFIDs)) {
        ## Identify if the quadrantMethSitesLinkedToSignificantGenes are among
        ## RE DNA methylation sites linked to the specific gene of interest
        TFVector <- quadrantMethSitesLinkedToSignificantGenes %in%
            TENETMultiAssayExperiment@metadata$step5OptimizeLinks[[
                quadrantResultsName
            ]][
                TENETMultiAssayExperiment@metadata$step5OptimizeLinks[[
                    quadrantResultsName
                ]]$geneID %in% topQuadrantGeneOrTFIDs[i],
                "DNAMethylationSiteID"
            ]

        methSitesLinkedToGenesDF[i + 6] <- TFVector
    }

    ## Reset the colnames and rownames of the DF
    colnames(methSitesLinkedToGenesDF) <- c(
        c(
            "DNAMethylationSiteID",
            "chromosome",
            "start",
            "end",
            "searchStart",
            "searchEnd"
        ),
        paste(
            topQuadrantGeneName,
            topQuadrantGeneOrTFIDs,
            "linked",
            sep = "_"
        )
    )
    rownames(methSitesLinkedToGenesDF) <- methSitesLinkedToGenesDF$
        DNAMethylationSiteID

    ## Create a GRanges copy of this dataset. This maintains the methylation
    ## site IDs in the names() of the GRanges object.
    methSitesLinkedToGenesDFGRanges <- GenomicRanges::makeGRangesFromDataFrame(
        df = methSitesLinkedToGenesDF[
            ,
            c("chromosome", "searchStart", "searchEnd")
        ],
        keep.extra.columns = FALSE,
        starts.in.df.are.0based = FALSE
    )

    ## For each of the peak files, identify if each RE DNA methylation site is
    ## found in the vicinity of at least one peak in each file, and add that
    ## information to methSitesLinkedToGenesDF
    methSitesLinkedToGenesDF <- cbind(
        methSitesLinkedToGenesDF,
        do.call(
            "cbind",
            parallel::mclapply(
                peakFileGRangesList,
                .methSitePeakGRangesOverlapTFFunction,
                methSiteGRanges = methSitesLinkedToGenesDFGRanges,
                mc.cores = coreCount
            )
        )
    )

    ## Return the data frame
    return(methSitesLinkedToGenesDF)
}

## Main step7TopGenesUserPeakOverlap function

#' Identify if RE DNA methylation sites linked to top genes/transcription
#' factors are located within a specific distance of specified genomic regions
#'
#' This function takes the top genes/transcription factors (TFs) by number of
#' linked RE DNA methylation sites identified by the
#' `step6DNAMethylationSitesPerGeneTabulation` function, up to the number
#' specified by the user, and identifies if the RE DNA methylation sites linked
#'  to those genes/TFs from the hyper- or hypomethylated G+ analysis quadrants
#' are found in the vicinity of genomic regions (peaks) of interest, supplied
#' by the user in the form of .bed, .narrowPeak, .broadPeak, and/or gappedPeak
#' files in a specified directory or given as a data frame or GRanges object.
#'
#' @param TENETMultiAssayExperiment Specify a MultiAssayExperiment object
#' containing expression and methylation SummarizedExperiment objects, such as
#' one created by the TCGADownloader function. This MultiAssayExperiment object
#' should also contain the results from the `step5OptimizeLinks` and
#' `step6DNAMethylationSitesPerGeneTabulation` functions in its metadata.
#' @param peakData Specify a data frame, matrix, or GRanges object with
#' genomic regions (peaks) of interest, organized in a bed-like manner (see
#' <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>), or a path to a
#' directory containing .bed, .narrowPeak, .broadPeak, and/or .gappedPeak files
#' with peaks of interest. The files may optionally be compressed
#' (.gz/.bz2/.xz).
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
#' @param DNAMethylationArray Specify the name of a DNA methylation probe
#' array supported by the sesameData package (see
#' `?sesameData::sesameData_getManifestGRanges`). If an array is specified,
#' RE DNA methylation sites and their locations in that array's manifest are
#' cross-referenced with RE DNA methylation site IDs included in the rownames
#' of the methylation dataset provided in the "methylation"
#' SummarizedExperiment object within the TENETMultiAssayExperiment object, and
#' only those overlapping will be considered for analysis. If set to NA, all RE
#' DNA methylation sites with locations listed in the rowRanges of the
#' "methylation" SummarizedExperiment object are used. Defaults to NA.
#' @param hypermethGplusAnalysis Set to TRUE to create data frames with the
#' peak overlap information for the unique hypermethylated RE DNA methylation
#' sites linked to the top genes/TFs by most hypermethylated RE DNA methylation
#' sites with G+ links. Defaults to TRUE.
#' @param hypomethGplusAnalysis Set to TRUE to create data frames with the peak
#' overlap information for the unique hypomethylated RE DNA methylation sites
#' linked to the top genes/TFs by most hypomethylated RE DNA methylation sites
#' with G+ links. Defaults to TRUE.
#' @param topGeneNumber Specify the number of top genes/TFs, based on the most
#' linked RE DNA methylation sites of a given analysis type, for which to
#' generate data showing overlap with user peak files for the RE DNA
#' methylation sites linked to those genes. Defaults to 10.
#' @param distanceFromREDNAMethylationSites Specify either 0 or a positive
#' integer in base pairs to be the distance from the linked RE DNA methylation
#' sites within which an RE DNA methylation site will be considered to overlap
#' a peak. Defaults to 100.
#' @param coreCount Argument passed as the mc.cores argument to mclapply. See
#' `?parallel::mclapply` for more details. Defaults to 1.
#' @return Returns the MultiAssayExperiment object given as the
#' TENETMultiAssayExperiment argument with an additional list of information
#' named 'step7TopGenesUserPeakOverlap' in its metadata with the output of this
#' function. This list is subdivided into hypermethGplus or hypomethGplus
#' results as selected by the user, which are further subdivided into data
#' frames with data for the unique RE DNA methylation sites linked to the top
#' overall genes, and for top TF genes only. Each of these data frames contain
#' a row for each of the unique RE DNA methylation sites linked to the top
#' genes/TFs for the specified analysis types. These note the locations of the
#' RE DNA methylation sites and the specified search windows, whether the RE
#' DNA methylation site is linked to each of the top genes/TFs, and whether
#' each RE DNA methylation site was found within the specified distance to any
#' peak in each of the user's peak files.
#' @export
#'
#' @examplesIf interactive()
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to do peak overlapping for all unique RE DNA
#' ## methylation sites linked to the top 10 genes, by number of linked hyper-
#' ## or hypomethylated RE DNA methylation sites, using a GRanges object
#' ## containing the genomic coordinates of peaks of interest. Gene names, and
#' ## the locations of RE DNA methylation sites, will be retrieved from the
#' ## rowRanges of the 'expression' and 'methylation' SummarizedExperiment
#' ## objects in the example MultiAssayExperiment. The analysis will be
#' ## performed using one CPU core.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Load the example peak GRanges object from the TENET.ExperimentHub package
#' exampleTENETPeakRegions <- TENET.ExperimentHub::exampleTENETPeakRegions()
#'
#' ## Use the example datasets to do the peak overlapping
#' returnValue <- step7TopGenesUserPeakOverlap(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     peakData = exampleTENETPeakRegions
#' )
#'
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to do peak overlapping for all unique RE DNA
#' ## methylation sites linked to the top 5 genes only, by number of linked
#' ## hypomethylated RE DNA methylation sites only, using bed-like files
#' ## containing the genomic coordinates of peaks of interest in the user's R
#' ## working directory. Gene names will be retrieved from the rowRanges of the
#' ## 'expression' SummarizedExperiment object in the example
#' ## MultiAssayExperiment, and RE DNA methylation sites and their locations
#' ## retrieved from the HM450 array via the sesameData package. A window of
#' ## 500 base pairs will be used to identify if the RE DNA methylation sites
#' ## lie within the vicinity of peaks. The analysis will be performed using 8
#' ## CPU cores.
#'
#' ## Load the example TENET MultiAssayExperiment object from the
#' ## TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Use the example dataset to do the peak overlapping
#' returnValue <- step7TopGenesUserPeakOverlap(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     peakData = ".",
#'     DNAMethylationArray = "HM450",
#'     hypermethGplusAnalysis = FALSE,
#'     topGeneNumber = 5,
#'     distanceFromREDNAMethylationSites = 500,
#'     coreCount = 8
#' )
step7TopGenesUserPeakOverlap <- function(
    TENETMultiAssayExperiment,
    peakData,
    geneAnnotationDataset = NA,
    DNAMethylationArray = NA,
    hypermethGplusAnalysis = TRUE,
    hypomethGplusAnalysis = TRUE,
    topGeneNumber = 10,
    distanceFromREDNAMethylationSites = 100,
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

    ## Check the type of information user has provided for peakData,
    ## whether it is a single peak file, or a directory with peak files.
    if (!is.character(peakData)) {
        ## Check that the file is a matrix or data frame and return an error
        ## if it is not
        if (!inherits(peakData, "GRanges")) {
            if (!is.data.frame(peakData)) {
                if (!is.matrix(peakData)) {
                    ## Return an error given it isn't a matrix or data frame
                    .stopNoCall(
                        "Please give a data frame, matrix, or GRanges object ",
                        "with information for peaks of interest, organized in ",
                        "a bed-like manner, as the peakData argument, or ",
                        "a path to a directory containing such files."
                    )
                }
            }

            ## Since the file is not a GRanges object, convert it to one

            ## Change the column names for the first three columns in
            ## the object
            colnames(peakData)[seq_len(3)] <- c("chr", "start", "end")

            ## Create a GRanges object. Assume starts are 0-based
            peakData <- GenomicRanges::makeGRangesFromDataFrame(
                df = peakData,
                keep.extra.columns = FALSE,
                starts.in.df.are.0based = TRUE
            )
        }

        ## To make this functionality consistent with functionality written to
        ## load files from a given directory, create a one-item peak file list
        ## containing the provided object
        peakFileGRangesList <- list("peakFile" = peakData)
    } else {
        ## Ensure that the supplied directory exists. If it does, load any .bed,
        ## .narrowPeak, .broadPeak, and/or .gappedPeak files inside.
        peakFileList <- .listExtBedFiles(
            extDir = peakData,
            paramName = "peakData",
            paramDescription = paste(
                "peaks for factors of interest to overlap with",
                "linked RE DNA methylation sites"
            )
        )

        ## Create a list to store the loaded peak files as GRanges objects
        peakFileGRangesList <- list()

        for (i in peakFileList) {
            ## Load the first three columns of the file as a GRanges object
            peaksGRanges <- rtracklayer::import.bed(
                i,
                colnames = c("chrom", "start", "end")
            )

            ## Add that GRanges object to the list
            peakFileGRangesList <- c(peakFileGRangesList, peaksGRanges)
        }

        ## Set the names to the peaks GRanges list to the names of the files
        names(peakFileGRangesList) <- basename(unlist(c(peakFileList)))
    }

    ## Get gene IDs and names from the MAE, or gene annotation dataset if
    ## provided
    geneIDdf <- .getGeneIDsAndNames(
        TENETMultiAssayExperiment, geneAnnotationDataset
    )

    ## Get methylation site IDs and names from the MAE, or methylation array if
    ## provided
    methSiteIDdf <- .getMethSiteIDsAndLocations(
        TENETMultiAssayExperiment, DNAMethylationArray
    )

    ## Create an empty list to hold the results from this step 7 function
    resultsList <- list()

    ## Generate results for the selected analysis types
    for (hyperHypo in analysisTypes) {
        ## Return results for all genes then TFs for each analysis type
        for (geneOrTF in c("Gene", "TF")) {
            ## Return DFs for each analysis type and top genes/TFs
            resultsList[[
                paste0(hyperHypo, "methGplusResults")
            ]][[
                paste0("top", geneOrTF, "s")
            ]] <- .userPeakOverlapInternalDatasetOutputFunction(
                TENETMultiAssayExperiment = TENETMultiAssayExperiment,
                hyperHypo = hyperHypo,
                geneOrTF = geneOrTF,
                geneIDdf = geneIDdf,
                methSiteIDdf = methSiteIDdf,
                peakFileGRangesList = peakFileGRangesList,
                topGeneNumber = topGeneNumber,
                distanceFromREDNAMethylationSites =
                    distanceFromREDNAMethylationSites,
                coreCount = coreCount
            )
        }
    }

    ## Add the results list to the MultiAssayExperiment
    TENETMultiAssayExperiment@metadata$
        step7TopGenesUserPeakOverlap <- resultsList

    return(TENETMultiAssayExperiment)
}
