## Internal functions used by step7TopGenesTADTables

## Internal function to identify RE DNA methylation sites or genes that overlap
## each TAD in a given file
.TADOverlapper <- function(chromosome, start, stop, referencePeaks, buffer) {
    paste(
        names(
            referencePeaks[(
                (GenomicRanges::seqnames(referencePeaks) == chromosome) &
                    ((GenomicRanges::start(referencePeaks) - buffer) <= stop) &
                    ((GenomicRanges::end(referencePeaks) + buffer) >= start)
            ), ]
        ),
        collapse = ","
    )
}

## Return genes in TAD function
.geneRowNumbersInTAD <- function(
    methSiteTADValue,
    geneIDdfTADColumnName) {
    columnNameSplit <- strsplit(geneIDdfTADColumnName, ",")

    rowNumbers <- paste(
        unique(
            which(
                vapply(
                    columnNameSplit,
                    `%in%`,
                    FUN.VALUE = logical(1),
                    x = methSiteTADValue
                )
            ),
        ),
        collapse = ","
    )

    return(rowNumbers)
}

## Internal function to count the genes in a TAD with a given RE DNA methylation
## site
.geneCountInTAD <- function(geneNumbers) {
    geneNumbersSeparated <- as.numeric(unlist(strsplit(geneNumbers, ",")))
    return(length(geneNumbersSeparated))
}

## Internal function to convert the gene numbers to the gene name
.geneNameFromNumberLister <- function(
    geneNumbers,
    downregulatedGeneDF,
    returnType) {
    geneNumbersSeparated <- as.numeric(unlist(strsplit(geneNumbers, ",")))
    geneNames <- downregulatedGeneDF[
        geneNumbersSeparated,
        returnType
    ]

    uniqueGeneNamesCombined <- paste(unique(geneNames), collapse = ",")

    return(uniqueGeneNamesCombined)
}

## Internal function to identify TADs in a given quadrant
.identifyTADsInQuadrant <- function(
    hyperHypo,
    TENETMultiAssayExperiment,
    topGeneNumber,
    geneOrTF,
    geneIDdf,
    methSiteIDdf,
    TADFileList) {
    quadrantResultsName <- paste0(hyperHypo, "methGplusResults")

    ## Ensure the quadrant's results are present in step 5
    .ensureStepPresent(
        TENETMultiAssayExperiment,
        stepName = "step5OptimizeLinks",
        substepName = quadrantResultsName
    )

    ## Get the IDs of the top genes/TFs in this quadrant
    topQuadrantGeneOrTFIDs <- .getQuadrantTopGenesOrTFs(
        TENETMultiAssayExperiment, geneOrTF, hyperHypo, topGeneNumber
    )$geneID
    if (.isSingleNA(topQuadrantGeneOrTFIDs)) {
        return(NA)
    }

    ## Create a data frame with some basic info about the RE DNA methylation
    ## sites of interest linked to the top genes selected
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

    ## If there are sites which aren't found in the specified methSiteIDdf,
    ## return a warning with those RE DNA methylation sites
    quadrantMethSitesLinkedToGenesNotFound <-
        quadrantMethSitesLinkedToSignificantGenes[
            !(quadrantMethSitesLinkedToSignificantGenes %in%
                methSiteIDdf$DNAMethylationSiteID
            )
        ]

    if (length(quadrantMethSitesLinkedToGenesNotFound) > 0) {
        .warningNoCall(
            "Genomic locations for the following RE DNA methylation sites ",
            "linked to the specified top genes were not identified: ",
            paste(
                quadrantMethSitesLinkedToGenesNotFound,
                collapse = ", "
            ),
            ". These DNA methylation sites have been excluded from this ",
            "analysis."
        )
    }

    ## Get the RE DNA methylation sites which are found in the specified
    ## methSiteIDdf
    quadrantMethSitesLinkedToSignificantGenes <-
        quadrantMethSitesLinkedToSignificantGenes[
            (quadrantMethSitesLinkedToSignificantGenes %in%
                methSiteIDdf$DNAMethylationSiteID
            )
        ]

    ## Create a data frame with the RE DNA methylation sites as the first and
    ## only column for now
    quadrantMethSiteDatasetLinkedToSignificantGenes <- data.frame(
        "DNAMethylationSiteID" = sort(
            quadrantMethSitesLinkedToSignificantGenes
        ),
        stringsAsFactors = FALSE
    )

    ## Add the RE DNA methylation site location info
    quadrantMethSiteDatasetLinkedToSignificantGenes$chromosome <-
        methSiteIDdf[
            quadrantMethSiteDatasetLinkedToSignificantGenes$
                DNAMethylationSiteID,
            "chromosome"
        ]

    quadrantMethSiteDatasetLinkedToSignificantGenes$start <- methSiteIDdf[
        quadrantMethSiteDatasetLinkedToSignificantGenes$DNAMethylationSiteID,
        "start"
    ]

    quadrantMethSiteDatasetLinkedToSignificantGenes$end <- methSiteIDdf[
        quadrantMethSiteDatasetLinkedToSignificantGenes$DNAMethylationSiteID,
        "end"
    ]

    ## For each of the genes/TFs of interest, create a new column in the
    ## DNA methylation site data frame indicating if each RE DNA methylation
    ## site is linked to that gene or not
    for (geneID in topQuadrantGeneOrTFIDs) {
        ## Get the gene name corresponding to the gene ID
        correspondingGeneName <- geneIDdf[
            geneID,
            "geneName"
        ]

        columnName <- paste(
            correspondingGeneName,
            geneID,
            "linked",
            sep = "_"
        )

        ## Get the RE DNA methylation sites linked to that gene specifically
        sitesLinkedToGene <- unique(
            TENETMultiAssayExperiment@metadata$step5OptimizeLinks[[
                quadrantResultsName
            ]][
                TENETMultiAssayExperiment@metadata$step5OptimizeLinks[[
                    quadrantResultsName
                ]]$geneID == geneID,
                "DNAMethylationSiteID"
            ]
        )

        ## Add a TRUE/FALSE vector to tell if each of the unique RE DNA
        ## methylation sites linked to at least one of the top n TFs is linked
        ## to this specific one
        quadrantMethSiteDatasetLinkedToSignificantGenes[columnName] <- (
            quadrantMethSiteDatasetLinkedToSignificantGenes$
                DNAMethylationSiteID %in%
                sitesLinkedToGene
        )
    }

    ## For each of the RE DNA methylation sites, see what genes are in the same
    ## TAD as it. To  do this, we will need to identify the TAD number(s) each
    ## RE DNA methylation site is in  and each gene is in
    for (TADFileIndex in seq_along(TADFileList)) {
        ## Create name for the TAD file overlap columns
        TADFileOverlapColumnName <- paste0(
            names(TADFileList[TADFileIndex]),
            "_overlapNumber"
        )

        ## Create name for the TAD file gene numbers found in the TAD
        geneRowOverlapColumnName <- paste0(
            names(TADFileList[TADFileIndex]),
            "_geneNumbers"
        )

        ## Create a name for the total number of genes in a TAD of the RE DNA
        ## methylation site
        geneCountColumnName <- paste0(
            names(TADFileList[TADFileIndex]),
            "_geneCountInTAD"
        )

        ## Create names for the gene IDs and gene names found in TADs
        geneIDColumnName <- paste0(
            names(TADFileList[TADFileIndex]),
            "_TADGeneIDs"
        )

        geneNameColumnName <- paste0(
            names(TADFileList[TADFileIndex]),
            "_TADGeneNames"
        )

        ## Get overlaps for each of the TADs with the RE DNA methylation sites
        ## of interest
        quadrantMethSiteDatasetLinkedToSignificantGenes[
            TADFileOverlapColumnName
        ] <- unname(
            mapply(
                .TADOverlapper,
                chromosome = quadrantMethSiteDatasetLinkedToSignificantGenes[
                    , 2
                ],
                start = quadrantMethSiteDatasetLinkedToSignificantGenes[
                    , 3
                ],
                stop = quadrantMethSiteDatasetLinkedToSignificantGenes[
                    , 4
                ],
                MoreArgs = list(
                    referencePeaks = TADFileList[[TADFileIndex]],
                    buffer = 0
                )
            )
        )

        ## Get the numbers of the genes that have overlapping TADs with the
        ## RE DNA methylation sites
        quadrantMethSiteDatasetLinkedToSignificantGenes[
            geneRowOverlapColumnName
        ] <- unname(
            vapply(
                quadrantMethSiteDatasetLinkedToSignificantGenes[
                    , TADFileOverlapColumnName # We do need drop = TRUE here
                ],
                .geneRowNumbersInTAD,
                character(1),
                geneIDdfTADColumnName = geneIDdf[
                    , TADFileOverlapColumnName # We do need drop = TRUE here
                ]
            )
        )

        ## Count the number of total genes within a TAD of the RE DNA
        ## methylation site
        quadrantMethSiteDatasetLinkedToSignificantGenes[
            geneCountColumnName
        ] <- unname(
            vapply(
                quadrantMethSiteDatasetLinkedToSignificantGenes[
                    , geneRowOverlapColumnName # We do need drop = TRUE here
                ],
                .geneCountInTAD,
                numeric(1)
            )
        )

        ## Convert the numbers into gene IDs
        quadrantMethSiteDatasetLinkedToSignificantGenes[
            geneIDColumnName
        ] <- unname(
            vapply(
                quadrantMethSiteDatasetLinkedToSignificantGenes[
                    , geneRowOverlapColumnName # We do need drop = TRUE here
                ],
                .geneNameFromNumberLister,
                character(1),
                downregulatedGeneDF = geneIDdf,
                returnType = "geneID"
            )
        )

        ## Convert the numbers into gene names
        quadrantMethSiteDatasetLinkedToSignificantGenes[
            geneNameColumnName
        ] <- unname(
            vapply(
                quadrantMethSiteDatasetLinkedToSignificantGenes[
                    , geneRowOverlapColumnName # We do need drop = TRUE here
                ],
                .geneNameFromNumberLister,
                character(1),
                downregulatedGeneDF = geneIDdf,
                returnType = "geneName"
            )
        )

        ## If no TAD is found, make a note of that in the gene name/ID
        ## columns. Additionally if no genes are found in the TAD, make a
        ## separate note of that
        for (i in seq_len(
            nrow(quadrantMethSiteDatasetLinkedToSignificantGenes)
        )) {
            if (quadrantMethSiteDatasetLinkedToSignificantGenes[
                i, TADFileOverlapColumnName
            ] == "") {
                quadrantMethSiteDatasetLinkedToSignificantGenes[
                    i, geneIDColumnName
                ] <- "No_TAD_identified"
                quadrantMethSiteDatasetLinkedToSignificantGenes[
                    i, geneNameColumnName
                ] <- "No_TAD_identified"
            }

            if (quadrantMethSiteDatasetLinkedToSignificantGenes[
                i, geneIDColumnName
            ] == "") {
                quadrantMethSiteDatasetLinkedToSignificantGenes[
                    i, geneIDColumnName
                ] <- "No_genes_identified_in_TAD"
                quadrantMethSiteDatasetLinkedToSignificantGenes[
                    i, geneNameColumnName
                ] <- "No_genes_identified_in_TAD"
            }
        }

        ## Remove unneeded columns
        quadrantMethSiteDatasetLinkedToSignificantGenes[
            geneRowOverlapColumnName
        ] <- NULL
    }

    return(quadrantMethSiteDatasetLinkedToSignificantGenes)
}

#' Create tables using user-supplied topologically associating domain (TAD)
#' information which identify the topologically associating domains
#' containing each RE DNA methylation site linked to the top genes/transcription
#' factors, as well as other genes in the same topologically associating domain
#' as potential downstream targets
#'
#' This function takes the top genes/transcription factors (TFs) by number of
#' linked RE DNA methylation sites identified by the
#' `step6DNAMethylationSitesPerGeneTabulation` function up to the number
#' specified by the user and generates tables with information for each of the
#' RE DNA methylation sites linked to them for both of the hyper- or
#' hypomethylated G+ analysis quadrants, as selected by the user. These tables
#' note which of the top genes/TFs each RE DNA methylation site is linked to,
#' as well as the total number and names of genes which happen to lie within
#' the same topologically associating domain (TAD) of each RE DNA methylation
#' site in each of the user-supplied TAD files.
#'
#' @param TENETMultiAssayExperiment Specify a MultiAssayExperiment object
#' containing expression and methylation SummarizedExperiment objects, such as
#' one created by the TCGADownloader function. This MultiAssayExperiment object
#' should also contain the results from the `step5OptimizeLinks` and
#' `step6DNAMethylationSitesPerGeneTabulation` functions in its metadata.
#' @param TADFiles Specify a data frame, matrix, or GRanges object with
#' information on the TAD compartments of interest, organized in a bed-like
#' manner (see <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>), or a path
#' to a directory that contains bed-like files that contain such TAD
#' information. If a path is provided, multiple bed-formatted TAD files can be
#' included in the specified directory. The files may optionally be compressed
#' (.gz/.bz2/.xz). **Note:** Data frames and matrices will be assumed to use
#' 1-indexed coordinates; `rtracklayer::import.bed` converts coordinates to
#' 1-indexed upon import.
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
#' @param DNAMethylationArray Specify the name of a DNA methylation probe array
#' supported by the sesameData package (see
#' `?sesameData::sesameData_getManifestGRanges`). If an array is specified, RE
#' DNA methylation sites and their locations in that array's manifest are
#' cross-referenced with RE DNA methylation site IDs included in the rownames
#' of the methylation dataset provided in the "methylation"
#' SummarizedExperiment object within the TENETMultiAssayExperiment object, and
#' only those overlapping will be considered for analysis. If set to NA, all RE
#' DNA methylation sites with locations listed in the rowRanges of the
#' "methylation" SummarizedExperiment object are used. Defaults to NA.
#' @param hypermethGplusAnalysis Set to TRUE to create TAD tables for the RE
#' DNA methylation sites linked to the top genes/TFs by most hypermethylated RE
#' DNA methylation sites with G+ links.
#' @param hypomethGplusAnalysis Set to TRUE to create TAD tables for the RE DNA
#' methylation sites linked to the top genes/TFs by most hypomethylated RE DNA
#' methylation sites with G+ links.
#' @param topGeneNumber Specify the number of top genes/TFs, based on the most
#' linked RE DNA methylation sites of a given analysis type, for which to
#' generate TAD tables for the RE DNA methylation sites linked to those genes.
#' Defaults to 10.
#' @param coreCount Argument passed as the mc.cores argument to mcmapply. See
#' `?parallel::mcmapply` for more details. Defaults to 1.
#' @return Returns the MultiAssayExperiment object given as the
#' TENETMultiAssayExperiment argument with an additional list of information
#' named 'step7TopGenesTADTables' in its metadata with the output of this
#' function. This list is subdivided into hypermethGplus or hypomethGplus
#' results as selected by the user, which are further subdivided into data
#' frames with data for the unique RE DNA methylation sites linked to the top
#' overall genes, and for top TF genes only. This includes the top genes/TFs
#' each RE DNA methylation site is linked to, and, for each TAD file, if an RE
#' DNA methylation site was found in a TAD in that file, as well as the gene
#' count and identities of other genes found in the same TAD as each RE DNA
#' methylation site.
#' @export
#'
#' @examplesIf interactive()
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to do overlapping for all unique RE DNA
#' ## methylation sites linked to the top 10 genes, by number of linked hyper-
#' ## or hypomethylated RE DNA methylation sites, using a GRanges object
#' ## containing topologically associating domain (TAD) data from the
#' ## TENET.ExperimentHub package. Gene names and locations, and the locations
#' ## of RE DNA methylation sites, will be retrieved from the rowRanges of the
#' ## 'expression' and 'methylation' SummarizedExperiment objects in the
#' ## example MultiAssayExperiment. The analysis will be performed using one
#' ## CPU core.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Load the example TAD GRanges object from the TENET.ExperimentHub package
#' exampleTENETTADRegions <- TENET.ExperimentHub::exampleTENETTADRegions()
#'
#' ## Use the example datasets to do the TAD overlapping
#' returnValue <- step7TopGenesTADTables(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     TADFiles = exampleTENETTADRegions
#' )
#'
#' ## This example uses the example MultiAssayExperiment provided in the
#' ## TENET.ExperimentHub package to do overlapping for all unique RE DNA
#' ## methylation sites linked to the top 5 genes only, by number of linked
#' ## hypomethylated RE DNA methylation sites only, with bed-like files
#' ## containing topologically associating domain (TAD) data located in the
#' ## user's R working directory. This analysis will be done using gene names
#' ## and their locations which are included in the "geneName" column of the
#' ## elementMetadata of the rowRanges of the "expression" SummarizedExperiment
#' ## object within the TENETMultiAssayExperiment object, and RE DNA
#' ## methylation sites and their locations which are included in the HM450
#' ## array retrieved via the sesameData package. The analysis will be
#' ## performed using 8 CPU cores.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#'
#' ## Use the example dataset to do the TAD overlapping
#' returnValue <- step7TopGenesTADTables(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     TADFiles = ".",
#'     DNAMethylationArray = "HM450",
#'     hypermethGplusAnalysis = FALSE,
#'     topGeneNumber = 5,
#'     coreCount = 8
#' )
step7TopGenesTADTables <- function(
    TENETMultiAssayExperiment,
    TADFiles,
    geneAnnotationDataset = NA,
    DNAMethylationArray = NA,
    hypermethGplusAnalysis = TRUE,
    hypomethGplusAnalysis = TRUE,
    topGeneNumber = 10,
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

    ## Check the type of information user has provided for TADFiles, whether it
    ## is a single TAD file, or a directory with TAD files.
    if (!(is.character(TADFiles) && length(TADFiles) == 1)) {
        ## Check that the object is a GRanges object, matrix, or data frame and
        ## return an error if it is not
        if (!is.data.frame(TADFiles) && !is.matrix(TADFiles)) {
            if (!inherits(TADFiles, "GRanges")) {
                .stopNoCall(
                    "Please give a data frame, matrix, or GRanges object ",
                    "with information on the TAD compartments of ",
                    "interest, organized in a bed-like manner, as the ",
                    "TADFiles argument, or a path to a directory ",
                    "containing bed-like files, which may optionally be ",
                    "compressed (.gz/.bz2/.xz)."
                )
            }
        } else {
            ## It's a matrix or data frame; convert it into a GRanges object
            TADFiles <- GenomicRanges::makeGRangesFromDataFrame(
                TADFiles,
                starts.in.df.are.0based = FALSE
            )

            ## Add the row numbers as names
            names(TADFiles) <- seq_len(length(TADFiles))
        }

        ## To make this functionality consistent with functionality written to
        ## load files from a given directory, add the data frame or matrix to
        ## a list
        TADFileList <- list("TADFile" = TADFiles)
    } else {
        ## Check to see that there are TAD files included in the specified
        ## folder. Do not filter on file extension because it is not necessarily
        # .bed.
        TADFilePaths <- list.files(path = TADFiles, full.names = TRUE)

        ## If no TAD files are found, return an error
        if (length(TADFilePaths) == 0) {
            ## No TAD files found. Quit the function and note that to the user
            .stopNoCall(
                "No TAD files were found in the directory specified for ",
                "TADFiles. Please place bed-like files containing TAD ",
                "information in the directory specified for TADFiles. ",
                "The files may optionally be compressed (.gz/.bz2/.xz)."
            )
        }

        ## Create an empty list to load each of the files into
        TADFileList <- list()

        ## Prepare and load the TAD files
        ## This assumes the files are at least bed3 formatted with the first 3
        ## columns being the chromosome, start, and end coordinates, and are
        ## 0-indexed
        for (filePathIndex in seq_along(TADFilePaths)) {
            ## Load the first three columns of the file as a GRanges object
            TADFiles <- rtracklayer::import.bed(
                TADFilePaths[filePathIndex],
                colnames = c("chrom", "start", "end")
            )

            ## Add the row numbers as names
            names(TADFiles) <- seq_len(length(TADFiles))

            TADFileList[[filePathIndex]] <- TADFiles
        }

        ## Assign names of the files (without extensions) to the TAD file list
        names(TADFileList) <- sub(
            "\\.[^.]*$",
            "",
            basename(TADFilePaths)
        )
    }

    ## Get entries from the MultiAssayExperiment object which have
    ## data in the geneID and DNAMethylationSiteID DFs
    geneIDdf <- geneIDdf[
        rownames(
            MultiAssayExperiment::assays(
                TENETMultiAssayExperiment@ExperimentList[["expression"]]
            )[[1]]
        ),
    ]

    geneIDdf <- stats::na.omit(geneIDdf)

    methSiteIDdf <- methSiteIDdf[
        rownames(
            MultiAssayExperiment::assays(
                TENETMultiAssayExperiment@ExperimentList[["methylation"]]
            )[[1]]
        ),
    ]

    methSiteIDdf <- stats::na.omit(methSiteIDdf)

    ## Check that there are still genes and RE DNA methylation sites available,
    ## and return an error if not
    if (nrow(geneIDdf) == 0) {
        .stopNoCall(
            "No genes listed in the specified geneAnnotationDataset ",
            "object are present in the rownames of the expression data of the ",
            "'expression' SummarizedExperiment object in the ",
            "TENETMultiAssayExperiment object. Please check the ",
            "geneAnnotationDataset argument and the names of the ",
            "genes in the expression SummarizedExperiment object ",
            "and run this function again."
        )
    }

    if (nrow(methSiteIDdf) == 0) {
        .stopNoCall(
            "No probes in the specified DNA methylation array ",
            "were also found in the rownames of the dataset ",
            "provided in the 'methylation' SummarizedExperiment ",
            "of the TENETMultiAssayExperiment object. Please ",
            "check the specified DNAMethylationArray ",
            "argument and the 'methylation' SummarizedExperiment ",
            "included in the TENETMultiAssayExperiment object ",
            "and run this function again."
        )
    }

    ## Annotate each gene with the TADs it is found in for each file, which will
    ## be overlapped with the information from the RE DNA methylation sites from
    ## each analysis
    for (TADFileIndex in seq_along(TADFileList)) {
        ## Create name for the TAD file overlap columns
        TADFileOverlapColumnName <- paste0(
            names(TADFileList[TADFileIndex]),
            "_overlapNumber"
        )

        ## Get overlaps for each of the TADs with the genes of interest
        geneIDdf[TADFileOverlapColumnName] <- unname(
            parallel::mcmapply(
                .TADOverlapper,
                chromosome = geneIDdf[, "chromosome"],
                start = geneIDdf[, "TSS"],
                stop = geneIDdf[, "TSS"],
                MoreArgs = list(
                    referencePeaks = TADFileList[[TADFileIndex]],
                    buffer = 0
                ),
                mc.cores = coreCount
            )
        )
    }

    ## Create an empty list to hold the results from this step 7 function
    resultsList <- list()

    ## Do the analysis for the analysis types selected by the user
    for (hyperHypo in analysisTypes) {
        for (geneOrTF in c("Gene", "TF")) {
            resultsList[[paste0(hyperHypo, "methGplusResults")]][[
                paste0("top", geneOrTF, "s")
            ]] <- .identifyTADsInQuadrant(
                hyperHypo,
                TENETMultiAssayExperiment,
                topGeneNumber,
                geneOrTF,
                geneIDdf,
                methSiteIDdf,
                TADFileList
            )
        }
    }

    ## Add the results list to the MultiAssayExperiment object
    TENETMultiAssayExperiment@metadata$step7TopGenesTADTables <-
        resultsList

    return(TENETMultiAssayExperiment)
}
