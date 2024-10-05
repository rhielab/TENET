## Internal functions used by step2GetDifferentiallyMethylatedSites

## Internal functions to calculate the densities of all RE DNA methylation sites
## in each case sample, and grab the mean density values for the top n samples

## Get y values for the density of a given plot
.densityY <- function(inputVector, overallDensity) {
    densityValues <- stats::density(
        inputVector,
        from = min(overallDensity$x),
        to = max(overallDensity$x),
        na.rm = TRUE
    )

    return(densityValues$y)
}

## For each of the top n density values for each block within a
## group of samples, get their mean values
.topNMeanGrabber <- function(dataFrame, topN) {
    ## Sort the values in the vector from largest to smallest
    blockValues <- sort(dataFrame, decreasing = TRUE)

    ## Return the mean of the largest topN
    return(mean(blockValues[seq_len(topN)]))
}

#' Identify differentially methylated RE DNA methylation sites
#'
#' This function identifies DNA methylation sites that mark putative regulatory
#' elements (REs), including enhancer and promoter regions. These are sites
#' that lie within regions with specific histone modifications and open
#' chromatin regions, from a user-supplied GRanges object, such as one created
#' by the `step1MakeExternalDatasets` function, and which are located at a
#' user-specified distance relative to the transcription start sites (TSS)
#' listed in either the rowRanges of the elementMetadata of the "expression"
#' SummarizedExperiment in the TENETMultiAssayExperiment object, or the
#' selected `geneAnnotationDataset` (which will be filtered to only genes and
#' transcripts). After identifying DNA methylation sites representing the
#' specified REs, the function classifies the RE DNA methylation sites as
#' methylated, unmethylated, hypermethylated, or hypomethylated based on their
#' differential methylation between the control and case samples supplied by
#' the user, defined by cutoff values which are either automatically based
#' on the mean methylation densities of the identified RE DNA methylation
#' sites, or manually set by the user. **Note:** Using the algorithm to set
#' cutoffs is recommended for use with DNA methylation array data, and may not
#' work for whole-genome DNA methylation data.
#'
#' @param TENETMultiAssayExperiment Specify a MultiAssayExperiment object
#' containing expression and methylation SummarizedExperiment objects,
#' such as one created by the TCGADownloader function. Coordinates for
#' genes and DNA methylation sites should be included in the rowRanges
#' of their respective SummarizedExperiment objects and should be annotated
#' from the same genome build as the regions given in the
#' regulatoryElementGRanges object.
#' @param regulatoryElementGRanges Specify a GRanges object containing genomic
#' regions representing regulatory elements of interest to the user.
#' Coordinates for the regulatory element regions should be annotated to the
#' same genome build as the gene and DNA methylation site coordinates given
#' in the TENETMultiAssayExperiment object. If this argument is set to NA or
#' not specified, this function will use all DNA methylation sites
#' representing regulatory elements of interest as defined by the
#' assessPromoter and TSSDist arguments. Defaults to NA.
#' @param geneAnnotationDataset Specify a gene annotation dataset which is
#' used to identify transcription start sites in order to find DNA methylation
#' sites within regulatory elements of interest (promoters or enhancers) in
#' conjunction with the settings of the assessPromoter and TSSDist arguments.
#' The dataset will be filtered to only genes and transcripts. The argument must
#' be either a GRanges object (such as one imported via `rtracklayer::import`)
#' or a path to a GFF3 or GTF file. Both GENCODE and Ensembl annotations are
#' supported. Other annotation datasets may work, but have not been tested.
#' See the "Input data" section of the vignette for information on the required
#' dataset format.
#' Specify NA to use the start coordinates of all entries in the elementMetadata
#' of the rowRanges of the "expression" SummarizedExperiment object within the
#' TENETMultiAssayExperiment object, in which case no filtering will be done
#' and all entries will be assumed to represent transcripts. Defaults to NA.
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
#' @param assessPromoter Set to TRUE to identify DNA methylation sites
#' that mark promoter regions or FALSE to identify distal enhancer regions.
#' Defaults to FALSE.
#' @param TSSDist Specify a positive integer distance in base pairs to any
#' transcription start site (see `geneAnnotationDataset`) within which DNA
#' methylation sites are considered promoter DNA methylation sites. DNA
#' methylation sites outside of the TSSDist from any transcription start site
#' will be considered enhancer methylation sites. Defaults to 1500.
#' @param purityData Specify a SummarizedExperiment object which contains DNA
#' methylation datasets collected from potential cell types which might affect
#' the purity of the patient samples contained in the
#' TENETMultiAssayExperiment. The coordinates for DNA methylation sites in
#' this dataset should be included in the rowRanges of the purityData
#' SummarizedExperiment object. Additionally, the DNA methylation site IDs
#' in the purityData SummarizedExperiment object should overlap with DNA
#' methylation sites present in the TENETMultiAssayExperiment and only those
#' that do overlap will be considered for analysis. Defaults to NA.
#' @param methCutoff Specify a number from 0 to 1 to be the beta-value cutoff
#' for methylated RE DNA methylation sites. If unspecified or NA, an algorithm
#' will be used to find the optimal cutoff value.
#' @param hypomethCutoff Specify a number from 0 to 1 to be the beta-value
#' cutoff for hypomethylated RE DNA methylation sites. Should be set lower than
#' the methCutoff. If unspecified or NA, an algorithm will be used to find the
#' optimal cutoff value.
#' @param hypermethCutoff Specify a number from 0 to 1 to be the beta-value
#' cutoff for hypermethylated RE DNA methylation sites. Should be set higher
#' than the unmethCutoff. If unspecified or NA, an algorithm will be used to
#' find the optimal cutoff value.
#' @param unmethCutoff Specify a number from 0 to 1 to be the beta-value cutoff
#' for unmethylated RE DNA methylation sites. If unspecified or NA, an
#' algorithm will be used to find the optimal cutoff value.
#' @param methUnmethProportionOffset Specify a number from 0 to 1 to be
#' the proportion of the distance of the region between the first and last
#' local maxima in the density plot of the mean methylation values of the
#' RE DNA methylation sites in the control samples. This value is then added
#' to the position of these local maxima to set the unmethylation and
#' methylation cutoffs if they are not defined by the user. The value ideally
#' should not exceed 0.5. Defaults to 0.2.
#' @param hypomethHypermethProportionOffset Specify a number from 0 to 1
#' to be the proportion of the distance of the region between the first and
#' last local maxima in the density plot of the mean methylation values
#' of the RE DNA methylation sites in the case samples. This value is then
#' added to the calculated unmethylation and methylation cutoffs to then set
#' the hypermethylation and hypomethylation cutoffs if they are not defined by
#' the user. The value ideally should not exceed 0.5. Defaults to 0.1.
#' @param minCaseCount Specify a positive integer to be the
#' minimum number of case samples to be considered for the
#' hyper- or hypomethylated groups. Should be less than the total number
#' of case samples.
#' @param cgDNAMethylationSitesOnly Set to TRUE to include only RE DNA
#' methylation sites with IDs that start with "cg". TRUE means that RE DNA
#' methylation sites whose IDs do not start with "cg" will be removed from
#' TENET analyses. Defaults to TRUE.
#' @return Returns the MultiAssayExperiment object given as the
#' TENETMultiAssayExperiment argument with an additional list of data named
#' "step2GetDifferentiallyMethylatedSites" in its metadata with the output data
#' from this function. These data include the set of calculated cutoff values,
#' the identities and counts of the classified RE DNA methylation sites, as
#' well as plots of the mean methylation distributions of the identified
#' regulatory element DNA methylation sites in the case and control samples and
#' the set cutoff values. Of note for plots, if assessPromoter is TRUE, two
#' distribution plots are saved, one using all promoter DNA methylation sites,
#' and one using DNA methylation sites which are identified to overlap REs.
#' @export
#'
#' @examplesIf interactive()
#' ## This example uses datasets provided in the TENET.ExperimentHub package to
#' ## perform an example analysis, analyzing RE DNA methylation sites in
#' ## potential enhancer elements, located over 1500 bp from transcription
#' ## start sites listed for genes and transcripts in the GENCODE v36 human
#' ## genome annotations, otherwise using default settings and a minimum case
#' ## sample count of 5.
#'
#' ## Load the example TENET MultiAssayExperiment object, and the example
#' ## GRanges object created by the TENET step 1 function, from the
#' ## TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#' exampleStep1MakeExternalDatasetsGRangesObject <-
#'     TENET.ExperimentHub::exampleTENETStep1MakeExternalDatasetsGRanges()
#'
#' ## Use the example datasets to identify differentially methylated
#' ## RE DNA methylation sites
#' returnValue <- step2GetDifferentiallyMethylatedSites(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     regulatoryElementGRanges =
#'         exampleTENETStep1MakeExternalDatasetsGRanges,
#'     minCaseCount = 5
#' )
#'
#' ## This example uses the same datasets, this time analyzing DNA methylation
#' ## sites in promoter elements, considering all RE DNA methylation sites
#' ## found within 2000 bp of all transcription start sites provided in the
#' ## MultiAssayExperiment only. Additionally, the methylation cutoffs are
#' ## manually set to 0.8, 0.7, 0.3, and 0.2 for the `methCutoff`,
#' ## `hypomethCutoff`, `hypermethCutoff`, and `unmethCutoff` respectively. The
#' ## `minCaseCount` is set to 10 samples and all RE DNA methylation sites
#' ## regardless of ID will be considered.
#'
#' ## Load the example TENET MultiAssayExperiment object
#' ## from the TENET.ExperimentHub package as well as the example GRanges object
#' ## created by the TENET step 1 function from the TENET.ExperimentHub package
#' exampleTENETMultiAssayExperiment <-
#'     TENET.ExperimentHub::exampleTENETMultiAssayExperiment()
#' exampleStep1MakeExternalDatasetsGRangesObject <-
#'     TENET.ExperimentHub::exampleTENETStep1MakeExternalDatasetsGRanges()
#'
#' ## Use the example datasets to identify differentially methylated
#' ## RE DNA methylation sites
#' returnValue <- step2GetDifferentiallyMethylatedSites(
#'     TENETMultiAssayExperiment = exampleTENETMultiAssayExperiment,
#'     regulatoryElementGRanges =
#'         exampleTENETStep1MakeExternalDatasetsGRanges,
#'     geneAnnotationDataset = NA,
#'     assessPromoter = TRUE,
#'     TSSDist = 2000,
#'     methCutoff = 0.8,
#'     hypomethCutoff = 0.7,
#'     hypermethCutoff = 0.3,
#'     unmethCutoff = 0.2,
#'     minCaseCount = 10,
#'     cgDNAMethylationSitesOnly = FALSE
#' )
step2GetDifferentiallyMethylatedSites <- function(
    TENETMultiAssayExperiment,
    regulatoryElementGRanges = NA,
    geneAnnotationDataset = NA,
    DNAMethylationArray = NA,
    assessPromoter = FALSE,
    TSSDist = 1500,
    purityData = NA,
    methCutoff = NA,
    hypomethCutoff = NA,
    hypermethCutoff = NA,
    unmethCutoff = NA,
    methUnmethProportionOffset = 0.2,
    hypomethHypermethProportionOffset = 0.1,
    minCaseCount,
    cgDNAMethylationSitesOnly = TRUE) {
    ## Return an error message if the input MultiAssayExperiment is invalid
    .validateMultiAssayExperiment(TENETMultiAssayExperiment)

    ## Bail early if minCaseCount is missing because it takes a long
    ## time to error out
    if (missing(minCaseCount)) {
        .stopNoCall("The minCaseCount parameter must be specified.")
    }

    ## If purity data have been provided, verify the purity dataset provided
    if (!.isSingleNA(purityData)) {
        if (!inherits(purityData, "SummarizedExperiment")) {
            .stopNoCall(
                "The dataset given as the purityData argument does not ",
                "seem to be a SummarizedExperiment object. Please ensure ",
                "the dataset given as the purityData argument is a ",
                "SummarizedExperiment object."
            )
        }
    }

    if (is.na(geneAnnotationDataset)) {
        ## If no geneAnnotationDataset is given, use the GRanges
        ## object in the expression SummarizedExperiment

        ## Extract the TSS locations based on the strand of the listed
        ## genes/transcripts
        TSSLocations <- as.numeric(
            ifelse(
                as.character(
                    GenomicRanges::strand(
                        SummarizedExperiment::rowRanges(
                            TENETMultiAssayExperiment@ExperimentList$expression
                        )
                    )
                ) == "-",
                GenomicRanges::end(
                    SummarizedExperiment::rowRanges(
                        TENETMultiAssayExperiment@ExperimentList$expression
                    )
                ),
                ## Assume that genes with "*" (no strand info) have a positive
                ## strand
                GenomicRanges::start(
                    SummarizedExperiment::rowRanges(
                        TENETMultiAssayExperiment@ExperimentList$expression
                    )
                )
            )
        )

        ## Create a data frame with the TSS locations plus TSSDist buffer
        ## to turn into a GenomicRanges object
        tssDF <- data.frame(
            "chr" = as.character(
                GenomicRanges::seqnames(
                    SummarizedExperiment::rowRanges(
                        TENETMultiAssayExperiment@ExperimentList$expression
                    )
                )
            ),
            "start" = TSSLocations - TSSDist,
            "end" = TSSLocations + TSSDist,
            "name" = names(
                SummarizedExperiment::rowRanges(
                    TENETMultiAssayExperiment@ExperimentList$expression
                )
            ),
            stringsAsFactors = FALSE
        )
        rownames(tssDF) <- make.unique(tssDF$name)

        ## Create a GRanges object from the TSS data frame
        TSSGRanges <- GenomicRanges::makeGRangesFromDataFrame(
            df = tssDF,
            keep.extra.columns = FALSE,
            starts.in.df.are.0based = FALSE
        )
    } else {
        ## Load the gene expression annotations of interest
        tssDF <- .loadGeneAnnotationDataset(
            geneAnnotationDataset,
            featureTypes = c("gene", "transcript")
        )

        ## Add the TSSDist buffer to the TSS coordinates and set those
        ## as the new start and end coordinates
        tssDF$start <- tssDF$TSS - TSSDist
        tssDF$end <- tssDF$TSS + TSSDist

        ## Create a GRanges object from the TSS data frame
        TSSGRanges <- GenomicRanges::makeGRangesFromDataFrame(
            df = tssDF,
            keep.extra.columns = FALSE,
            starts.in.df.are.0based = FALSE
        )
    }

    ## Load the DNA methylation site annotations of interest
    if (is.na(DNAMethylationArray)) {
        ## If no DNAMethylationArray is given, use the GRanges
        ## object in the methylation SummarizedExperiment
        DNAMethylationGRanges <- SummarizedExperiment::rowRanges(
            TENETMultiAssayExperiment@ExperimentList$methylation
        )
    } else {
        ## Create a data frame with relevant info from the
        ## selected methylation array's manifest
        DNAMethylationDF <- .loadMethylationManifest(
            DNAMethylationArray = DNAMethylationArray
        )

        ## Identify the DNA methylation sites in the methylation array
        ## which are present in the rownames of the methylation data
        methSitesInMAE <- rownames(
            MultiAssayExperiment::assays(
                TENETMultiAssayExperiment@ExperimentList$methylation
            )[[1]]
        )

        ## If no DNA methylation sites are found to intersect, return an
        ## error. Otherwise, create a GRanges object of the intersected DNA
        ## methylation sites found in both, with locations provided in the
        ## specified DNA methylation array.
        if (length(intersect(
            DNAMethylationDF$DNAMethylationSiteID, methSitesInMAE
        ))
        == 0
        ) {
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
        } else {
            intersectedMethSites <- intersect(
                methSitesInMAE,
                DNAMethylationDF$DNAMethylationSiteID
            )

            DNAMethylationGRanges <- GenomicRanges::makeGRangesFromDataFrame(
                df = DNAMethylationDF[intersectedMethSites, ],
                keep.extra.columns = TRUE,
                starts.in.df.are.0based = FALSE
            )

            names(DNAMethylationGRanges) <- intersectedMethSites
        }
    }

    ## Make sure to identify sites which lack location information so we can
    ## make sure not to include them later
    methSitesLackingLocation <- unique(
        c(
            names(
                DNAMethylationGRanges[
                    GenomicRanges::seqnames(DNAMethylationGRanges) == "*",
                ]
            ),
            names(
                DNAMethylationGRanges[
                    GenomicRanges::start(DNAMethylationGRanges) == 0,
                ]
            ),
            names(
                DNAMethylationGRanges[
                    GenomicRanges::end(DNAMethylationGRanges) == 0,
                ]
            )
        )
    )

    ## We need to find the DNA methylation sites representing the regulatory
    ## elements of interest to the user, either promoter sites (within TSSDist
    ## of transcription start sites) or enhancer sites (outside of TSSDist of
    ## transcription start sites)

    ## Identify the DNA methylation sites of possible interest to the
    ## analysis based on their distance from the annotated transcription start
    ## sites as set by TSSDist
    if (!assessPromoter) {
        ## Find the DNA methylation sites outside of the TSS regions
        methSitesOfPossibleInterest <- names(
            DNAMethylationGRanges[
                setdiff(
                    seq_along(DNAMethylationGRanges),
                    unique(
                        S4Vectors::queryHits(
                            GenomicRanges::findOverlaps(
                                DNAMethylationGRanges,
                                TSSGRanges
                            )
                        )
                    )
                ),
            ]
        )
    } else {
        ## Find the DNA methylation sites within the TSS regions
        methSitesOfPossibleInterest <- names(
            DNAMethylationGRanges[
                unique(
                    S4Vectors::queryHits(
                        GenomicRanges::findOverlaps(
                            DNAMethylationGRanges,
                            TSSGRanges
                        )
                    )
                ),
            ]
        )
    }

    ## Remove the DNA methylation sites lacking location information from the
    ## IDed sites
    methSitesOfPossibleInterest <- methSitesOfPossibleInterest[
        !(methSitesOfPossibleInterest %in% methSitesLackingLocation)
    ]

    ## Ensure that some sites have been identified and return an error if not
    if (!length(methSitesOfPossibleInterest) > 0) {
        .stopNoCall(
            "No DNA methylation sites in the 'methylation' ",
            "SummarizedExperiment of the TENETMultiAssayExperiment object ",
            "within the specified distance from transcription start sites ",
            "have been identified. Please check the settings of the ",
            "function as well as the regions defined in the 'methylation' ",
            "and 'expression' objects and run this function again."
        )
    }

    ## Get all the DNA methylation sites found in the regulatory elements
    ## supplied by the user if the user has specified a GRanges object as the
    ## regulatoryElementGRanges. Otherwise, if it is NA, we will keep all
    ## enhancer or promoter DNA methylation sites identified.
    if (!.isSingleNA(regulatoryElementGRanges)) {
        ## Warnings are suppressed here because chromosome annotations
        ## found in one GRanges but not the other will return an "Each
        ## of the 2 combined objects has sequence levels not in the other"
        ## warning. It's difficult to guarantee there won't be regions
        ## present in one that aren't in the other, and the presence of
        ## these regions doesn't affect the analysis since they
        ## won't be overlapped anyway. Additionally, the documentation for
        ## the function specifically recommends the use of suppressWarnings:
        ## "(use suppressWarnings() to suppress this warning)"
        regulatoryElementMethSites <- suppressWarnings(
            names(
                DNAMethylationGRanges[
                    unique(
                        S4Vectors::queryHits(
                            GenomicRanges::findOverlaps(
                                DNAMethylationGRanges,
                                regulatoryElementGRanges
                            )
                        )
                    ),
                ]
            )
        )

        ## Identify those RE DNA methylation sites which are also a suitable
        ## distance from TSS as identified earlier
        methSitesOfInterest <- intersect(
            regulatoryElementMethSites,
            methSitesOfPossibleInterest
        )
    } else {
        ## Keep all methSitesOfPossibleInterest
        methSitesOfInterest <- methSitesOfPossibleInterest
    }

    ## Create a data frame for the RE DNA methylation sites
    methSitesDF <- data.frame(
        "DNAMethylationSiteID" = methSitesOfInterest,
        stringsAsFactors = FALSE
    )

    ## Get the names of the control and case samples in the
    ## expression and methylation data
    expressionSampleNamesControl <- .getExpOrMetSamplesOfType(
        TENETMultiAssayExperiment,
        "expression",
        "Control",
        namesOnly = TRUE
    )

    expressionSampleNamesCase <- .getExpOrMetSamplesOfType(
        TENETMultiAssayExperiment,
        "expression",
        "Case",
        namesOnly = TRUE
    )

    methylationSampleNamesControl <- .getExpOrMetSamplesOfType(
        TENETMultiAssayExperiment,
        "methylation",
        "Control",
        namesOnly = TRUE
    )

    methylationSampleNamesCase <- .getExpOrMetSamplesOfType(
        TENETMultiAssayExperiment,
        "methylation",
        "Case",
        namesOnly = TRUE
    )

    ## Calculate the mean methylation of each RE DNA methylation site in the
    ## control and case samples
    methSitesDF$meanControl <- apply(
        MultiAssayExperiment::assays(
            TENETMultiAssayExperiment@ExperimentList$methylation
        )[[1]][
            methSitesOfInterest,
            methylationSampleNamesControl,
            drop = FALSE
        ],
        1,
        mean,
        na.rm = TRUE
    )

    methSitesDF$meanCase <- apply(
        MultiAssayExperiment::assays(
            TENETMultiAssayExperiment@ExperimentList$methylation
        )[[1]][
            methSitesOfInterest,
            methylationSampleNamesCase,
            drop = FALSE
        ],
        1,
        mean,
        na.rm = TRUE
    )

    ## Remove the RE DNA methylation sites with NaN values
    presentMethSites <- methSitesDF[stats::complete.cases(methSitesDF), ]

    ## If purity data have been provided, ensure DNA methylation sites in those
    ## purity datasets are present in the presentMethSites dataset
    if (!.isSingleNA(purityData)) {
        ## Check that the RE DNA methylation sites present in the purityData
        ## SummarizedExperiment object are present in the
        ## methSitesOfPossibleInterest
        purityMethSites <- names(SummarizedExperiment::rowRanges(purityData))

        methSitesOfInterestPurityIntersect <- intersect(
            purityMethSites,
            presentMethSites$DNAMethylationSiteID
        )

        ## If there are no RE DNA methylation sites which intersect between the
        ## methSitesOfInterest and purity dataset, quit the function
        ## and return an error to the user
        if (length(methSitesOfInterestPurityIntersect) == 0) {
            .stopNoCall(
                "No RE DNA methylation sites identified within the user's ",
                "regulatory elements of interest were found in the dataset ",
                "given as the purityData argument. Please check the purity ",
                "data and ensure the RE DNA methylation sites present in the ",
                "rowRanges of the purityData SummarizedExperiment object are ",
                "found in those of the 'methylation' SummarizedExperiment ",
                "object within the TENETMultiAssayExperiment object and are ",
                "also found in the identified RE DNA methylation sites of ",
                "interest."
            )
        }
    }

    ## Use an algorithm to set the cutoffs not selected by the user
    if (any(is.na(c(
        methCutoff, hypomethCutoff, unmethCutoff, hypermethCutoff
    )))) {
        ## The algorithm does not support calculating individual cutoffs
        if (!all(is.na(c(
            methCutoff, hypomethCutoff, unmethCutoff, hypermethCutoff
        )))) {
            .stopNoCall(
                "Some cutoffs are unspecified. The user must specify ",
                "either all or none of the cutoffs."
            )
        }

        ## If all of these cutoffs are set to NA, use the algorithm to
        ## set the cutoffs

        ## Depending on the type of data the user is analyzing, use a
        ## different algorithm for analyzing promoter and distal enhancer
        ## RE DNA methylation sites

        if (assessPromoter) {
            ## To assess promoters, unlike with enhancers where we don't know
            ## where the enhancers are, we need to only work with RE DNA
            ## methylation sites that overlap the HM and NDR datasets. If this
            ## is done for promoters this instead enriches only for highly
            ## active promoters which affects the algorithm since the density
            ## plots involving these promoter DNA methylation sites tend to
            ## show a unimodal rather than bimodal distribution, which affects
            ## how the algorithm works.

            ## Instead, since potential promoter DNA methylation sites are
            ## known, regardless of cell type, we will use the methylation
            ## densities of all RE DNA methylation sites within the specified
            ## TSSDist of potential promoters regardless of if they are found
            ## in histone regions or not.

            ## Identify the promoter DNA methylation sites with non-NA values
            ## in both case and control datasets
            methDataControlPromoter <- MultiAssayExperiment::assays(
                TENETMultiAssayExperiment@ExperimentList$methylation
            )[[1]][
                methSitesOfPossibleInterest,
                methylationSampleNamesControl
            ]

            methDataCasePromoter <- MultiAssayExperiment::assays(
                TENETMultiAssayExperiment@ExperimentList$methylation
            )[[1]][
                methSitesOfPossibleInterest,
                methylationSampleNamesCase
            ]

            ## Get the promoter DNA methylation sites with non-NA values in
            ## both datasets
            methDataControlPromoter <- methDataControlPromoter[
                stats::complete.cases(methDataControlPromoter),
            ]

            methDataCasePromoter <- methDataCasePromoter[
                stats::complete.cases(methDataCasePromoter),
            ]

            ## Identify the RE DNA methylation sites with non-NA values in both
            ## datasets
            nonNAPromoterMethSites <- intersect(
                rownames(methDataControlPromoter),
                rownames(methDataCasePromoter)
            )

            methDataControlPromoter <- methDataControlPromoter[
                nonNAPromoterMethSites,
            ]

            methDataCasePromoter <- methDataCasePromoter[
                nonNAPromoterMethSites,
            ]
        }

        ## If not assessing promoters:
        ## Calculate the density of the mean methylation values for the
        ## enhancer DNA methylation sites in the control and case samples. and
        ## calculate the turnpoints in the density line.

        ## If assessing promoters:
        ## Calculate the density of the mean methylation values for the
        ## promoter DNA methylation sites under HM and NDR datasets or in the
        ## ENCODE datasets in the control and case samples, and calculate the
        ## turnpoints in the density line.

        controlMeanDensity <- stats::density(
            apply(
                MultiAssayExperiment::assays(
                    TENETMultiAssayExperiment@ExperimentList$methylation
                )[[1]][
                    presentMethSites$DNAMethylationSiteID,
                    methylationSampleNamesControl,
                    drop = FALSE
                ],
                1,
                mean,
                na.rm = TRUE
            )
        )
        controlMeanDensityYSeries <- stats::ts(controlMeanDensity$y)
        controlMeanTP <- pastecs::turnpoints(controlMeanDensityYSeries)

        caseMeanDensity <- stats::density(
            apply(
                MultiAssayExperiment::assays(
                    TENETMultiAssayExperiment@ExperimentList$methylation
                )[[1]][
                    presentMethSites$DNAMethylationSiteID,
                    methylationSampleNamesCase,
                    drop = FALSE
                ],
                1,
                mean,
                na.rm = TRUE
            )
        )
        caseMeanDensityYSeries <- stats::ts(caseMeanDensity$y)
        caseMeanTP <- pastecs::turnpoints(caseMeanDensityYSeries)

        ## If assessing promoters, also do the same but for all promoter DNA
        ## methylation sites
        if (assessPromoter) {
            controlMeanDensityAllPromoter <- stats::density(
                apply(methDataControlPromoter, 1, mean, na.rm = TRUE)
            )
            controlMeanDensityYSeriesAllPromoter <- stats::ts(
                controlMeanDensityAllPromoter$y
            )
            controlMeanTPAllPromoter <- pastecs::turnpoints(
                controlMeanDensityYSeriesAllPromoter
            )

            caseMeanDensityAllPromoter <- stats::density(
                apply(methDataCasePromoter, 1, mean, na.rm = TRUE)
            )
            caseMeanDensityYSeriesAllPromoter <- stats::ts(
                caseMeanDensityAllPromoter$y
            )
            caseMeanTPAllPromoter <- pastecs::turnpoints(
                caseMeanDensityYSeriesAllPromoter
            )
        }

        ## Calculate density curves for the enhancer/subset promoter DNA
        ## methylation sites' methylation values in each of the case samples
        caseSampleDensities <- apply(
            MultiAssayExperiment::assays(
                TENETMultiAssayExperiment@ExperimentList$methylation
            )[[1]][
                presentMethSites$DNAMethylationSiteID,
                methylationSampleNamesCase,
                drop = FALSE
            ],
            2,
            .densityY,
            overallDensity = caseMeanDensity
        )

        ## If assessing promoters, also do the same but for all promoter DNA
        ## methylation sites
        if (assessPromoter) {
            caseSampleDensitiesAllPromoter <- apply(
                methDataCasePromoter,
                2,
                .densityY,
                overallDensity = caseMeanDensityAllPromoter
            )
        }

        ## Identify the mean density values for the top n case
        ## samples, where n is the minCaseCount value
        caseTopBoundValues <- apply(
            caseSampleDensities,
            1,
            .topNMeanGrabber,
            topN = (minCaseCount)
        )

        ## If assessing promoters, also identify the mean density values for all
        ## promoter DNA methylation sites
        if (assessPromoter) {
            caseTopBoundValuesAllPromoter <- apply(
                caseSampleDensitiesAllPromoter,
                1,
                .topNMeanGrabber,
                topN = (minCaseCount)
            )
        }

        ## Create a plot with the densities of the mean enhancer/subset promoter
        ## DNA methylation site methylation values in the control samples
        plotTitle <- paste0(
            "TENET Methylation distribution plot\n",
            ifelse(assessPromoter, "promoter", "enhancer"),
            " DNA methylation sites overlapping REs"
        )

        .newInvisibleRecordablePlot()
        plot(
            controlMeanDensity,
            ylim = range(caseTopBoundValues),
            main = plotTitle,
            xlab = "Methylation Beta Value",
            ylab = "Relative Density"
        )

        ## Add the mean density line for the control samples in blue
        graphics::lines(
            controlMeanDensity,
            col = "blue",
            lwd = 2
        )

        ## Add the mean density line for the case samples in red
        graphics::lines(
            caseMeanDensity,
            col = "red",
            lwd = 2
        )

        ## Save the plot to an object
        distributionPlot <- .recordTENETSavedSizePlot()

        ## Close the plot
        grDevices::dev.off()

        if (assessPromoter) {
            ## Create the same plot but for all promoter DNA methylation sites
            ## identified regardless of overlap with promoter elements
            .newInvisibleRecordablePlot()
            plot(
                controlMeanDensityAllPromoter,
                ylim = range(caseTopBoundValuesAllPromoter),
                main = paste(
                    "TENET Methylation distribution plot",
                    "all promoter DNA methylation sites",
                    sep = "\n"
                ),
                xlab = "Methylation Beta Value",
                ylab = "Relative Density"
            )

            ## Add the mean density line for the control samples in blue
            graphics::lines(
                controlMeanDensityAllPromoter,
                col = "blue",
                lwd = 2
            )

            ## Add the mean density line for the case samples in red
            graphics::lines(
                caseMeanDensityAllPromoter,
                col = "red",
                lwd = 2
            )

            ## Save the plot to an object
            distributionPlotAllPromoter <- .recordTENETSavedSizePlot()

            ## Close the plot
            grDevices::dev.off()
        }

        ## If assessing promoters, consider all promoter DNA methylation sites
        ## in the cutoff calculation, not just ones overlapping with the
        ## datasets
        if (assessPromoter) {
            controlMeanTPForCutoff <- controlMeanTPAllPromoter
            caseMeanTPForCutoff <- caseMeanTPAllPromoter
            controlMeanDensityForCutoff <- controlMeanDensityAllPromoter
            caseMeanDensityForCutoff <- caseMeanDensityAllPromoter
        } else {
            controlMeanTPForCutoff <- controlMeanTP
            caseMeanTPForCutoff <- caseMeanTP
            controlMeanDensityForCutoff <- controlMeanDensity
            caseMeanDensityForCutoff <- caseMeanDensity
        }

        ## Ensure that there are at least 3 turnpoints in the control and case
        ## mean density, indicating this is a normal methylation distribution
        ## with at least two peaks and one trough.
        if (
            length(controlMeanTPForCutoff$tppos) < 3 |
                length(caseMeanTPForCutoff$tppos) < 3
        ) {
            ## Output an error regarding the unusual peak distribution
            .displayError(
                "Atypical DNA methylation distribution (non-bimodal)",
                "prevents the algorithm from properly setting cutoffs.",
                "Please plot the returned plot object and use it to",
                "manually set methylation cutoffs instead and",
                "rerun this function."
            )

            ## Return the plot so the user can look at it
            return(cutoffsPlot)
        }
        controlMeanPeak1Block <- controlMeanTPForCutoff$tppos[1]
        controlMeanPeak1X <- controlMeanDensityForCutoff$x[
            controlMeanPeak1Block
        ]

        controlMeanPeak2Block <- controlMeanTPForCutoff$tppos[
            length(controlMeanTPForCutoff$tppos)
        ]
        controlMeanPeak2X <- controlMeanDensityForCutoff$x[
            controlMeanPeak2Block
        ]

        caseMeanPeak1Block <- caseMeanTPForCutoff$tppos[1]
        caseMeanPeak1X <- caseMeanDensityForCutoff$x[
            caseMeanPeak1Block
        ]

        caseMeanPeak2Block <- caseMeanTPForCutoff$tppos[
            length(caseMeanTPForCutoff$tppos)
        ]
        caseMeanPeak2X <- caseMeanDensityForCutoff$x[
            caseMeanPeak2Block
        ]

        ## Calculate the number of these density plots between the
        ## first and last peaks in the control and case samples
        controlPeakDifferenceBlocks <- (
            controlMeanPeak2X - controlMeanPeak1X
        )

        casePeakDifferenceBlocks <- (
            caseMeanPeak2X - caseMeanPeak1X
        )

        ## Identify the block interval value by multiplying the
        ## methUnmethProportionOffset by the
        ## controlPeakDifferenceBlocks and the
        ## hypomethHypermethProportionOffset by the
        ## casePeakDifferenceBlocks
        methUnmethIntervalBlockValue <- round(
            (controlPeakDifferenceBlocks * methUnmethProportionOffset),
            3
        )

        hypomethHypermethIntervalBlockValue <- round(
            (casePeakDifferenceBlocks *
                hypomethHypermethProportionOffset),
            3
        )

        ## Calculate the block values where the cutoffs will be set.
        ## Set the methylation and unmethylation block cutoffs using
        ## the methUnmethIntervalBlockValue.
        methCutoffX <- (
            controlMeanPeak2X - methUnmethIntervalBlockValue
        )

        unmethCutoffX <- (
            controlMeanPeak1X + methUnmethIntervalBlockValue
        )

        ## Then calculate the blocks for the hyper- and hypomethylation
        ## cutoffs by taking the case peak locations, and
        ## subtracting both interval block values
        hypomethCutoffX <- (
            (caseMeanPeak2X - methUnmethIntervalBlockValue) -
                hypomethHypermethIntervalBlockValue
        )

        hypermethCutoffX <- (
            (caseMeanPeak1X + methUnmethIntervalBlockValue) +
                hypomethHypermethIntervalBlockValue
        )

        ## Set cutoffs by rounding the x values calculated
        methCutoff <- round(methCutoffX, 3)
        unmethCutoff <- round(unmethCutoffX, 3)
        hypomethCutoff <- round(hypomethCutoffX, 3)
        hypermethCutoff <- round(hypermethCutoffX, 3)

        ## Plot the mean densities plot with the cutoffs included

        ## Create a plot with the densities of the mean enhancer/promoter
        ## DNA methylation site methylation values in the control samples
        .newInvisibleRecordablePlot()
        plot(
            controlMeanDensityForCutoff,
            ylim = range(caseTopBoundValues),
            main = "TENET Cutoffs Selection Plot",
            xlab = "Methylation Beta Value",
            ylab = "Relative Density"
        )

        ## Add the mean density line for the control samples in blue
        graphics::lines(
            controlMeanDensityForCutoff,
            col = "blue",
            lwd = 2
        )

        ## Add the mean density line for the case samples in red
        graphics::lines(
            caseMeanDensityForCutoff,
            col = "red",
            lwd = 2
        )

        ## Plot the methylation/unmethylation cutoffs on the graph as lines
        graphics::abline(
            v = c(unmethCutoffX, methCutoffX),
            col = c("#CE68FC", "#FFD14F")
        )

        ## Plot the hypermethylation/hypomethylation cutoffs on the graph as
        ## lines
        graphics::abline(
            v = c(hypermethCutoffX, hypomethCutoffX),
            col = c("#7204A4", "#D09A00")
        )

        ## Save the plot to an object
        cutoffsPlot <- .recordTENETSavedSizePlot()

        ## Close the plot
        grDevices::dev.off()
    } else {
        ## If the cutoffs were manually specified, save NAs to the
        ## MultiAssayExperiment instead of the plot
        cutoffsPlot <- NA
        distributionPlot <- NA
        distributionPlotAllPromoter <- NA
    }

    if (!.isSingleNA(purityData)) {
        ## If purity datasets are specified, identify the RE DNA methylation
        ## sites for which purity data are present, then identify
        ## methylated/unmethylated/hypermethylated/hypomethylated RE DNA
        ## methylation sites the same way as without purity, but with
        ## additional exclusion criteria that hypermethylated RE DNA
        ## methylation sites should have their lowest purity dataset mean
        ## methylation value greater than the methylation cutoff, as for
        ## control samples, and their highest purity dataset mean methylation
        ## value less than the unmethylation cutoff, again as for control
        ## samples.

        ## For subsetting purposes later, set the rownames of the
        ## presentMethSites data frame to be the RE DNA methylation sites
        rownames(presentMethSites) <- presentMethSites$DNAMethylationSiteID

        ## For each of the assays with data in the purityData object, calculate
        ## the mean methylation for each of the RE DNA methylation sites in
        ## those datasets across the samples provided in each assay
        puritySiteMethylationMeansList <- list()

        for (i in seq_along(purityData@assays@data)) {
            ## Get the matrix/DF with data from the given assay and
            ## calculate a mean along the rows of it
            siteMethylationMeans <- apply(
                purityData@assays@data[[i]],
                1, ## rows
                FUN = mean
            )

            names(siteMethylationMeans) <- rownames(purityData@assays@data[[i]])

            puritySiteMethylationMeansList[[i]] <- siteMethylationMeans
        }

        ## Assemble a data frame listing all the RE DNA methylation sites with
        ## any data in the purity dataset with the mean methylation of those RE
        ## DNA methylation sites in each of the datasets contained in the
        ## purity dataset
        methSitesOfInterestPurityIntersectMeanDF <- data.frame(
            DNAMethylationSiteIDs = methSitesOfInterestPurityIntersect,
            stringsAsFactors = FALSE
        )

        ## For each of the purity datasets, go through the purity RE DNA
        ## methylation sites and record the mean value for each of the RE DNA
        ## methylation sites in that individual dataset, or otherwise record NA
        ## if that site isn't present in that dataset
        for (i in seq_along(purityData@assays@data)) {
            methSitesOfInterestPurityIntersectMeanDF[
                names(purityData@assays@data)[i]
            ] <- puritySiteMethylationMeansList[[i]][
                methSitesOfInterestPurityIntersectMeanDF$DNAMethylationSiteIDs
            ]
        }

        ## Set the RE DNA methylation site column as rownames and remove the
        ## column with DNA methylation site IDs
        rownames(
            methSitesOfInterestPurityIntersectMeanDF
        ) <- methSitesOfInterestPurityIntersectMeanDF$DNAMethylationSiteIDs

        methSitesOfInterestPurityIntersectMeanDF$DNAMethylationSiteIDs <- NULL

        ## Remove sites with all NA values for every purity dataset
        methSitesOfInterestPurityIntersectMeanDF <-
            methSitesOfInterestPurityIntersectMeanDF[
                !(rowSums(is.na(methSitesOfInterestPurityIntersectMeanDF)) ==
                    ncol(methSitesOfInterestPurityIntersectMeanDF)
                ),
            ]

        ## For each RE DNA methylation site, calculate the lowest and highest
        ## mean value from the purity datasets
        methSitesOfInterestPurityIntersectMeanDF$highestMean <- apply(
            methSitesOfInterestPurityIntersectMeanDF,
            1,
            max,
            na.rm = TRUE
        )

        ## Identify the lowest mean value
        methSitesOfInterestPurityIntersectMeanDF$lowestMean <- apply(
            methSitesOfInterestPurityIntersectMeanDF,
            1,
            min,
            na.rm = TRUE
        )

        ## Assemble a presentMethSites out of the RE DNA methylation sites with
        ## purity info
        presentMethSites <- data.frame(
            "DNAMethylationSiteID" = rownames(
                methSitesOfInterestPurityIntersectMeanDF
            ),
            "meanControl" = presentMethSites[
                rownames(methSitesOfInterestPurityIntersectMeanDF),
                "meanControl"
            ],
            "meanCase" = presentMethSites[
                rownames(methSitesOfInterestPurityIntersectMeanDF),
                "meanCase"
            ],
            "highestMean" = methSitesOfInterestPurityIntersectMeanDF$
                highestMean,
            "lowestMean" = methSitesOfInterestPurityIntersectMeanDF$
                lowestMean,
            stringsAsFactors = FALSE
        )
    }

    ## Get the methylated quadrant RE DNA methylation sites

    ## Get a copy of the case data for RE DNA methylation sites with average
    ## methylation above the methCutoff in control and case samples
    methDataSubset <- MultiAssayExperiment::assays(
        TENETMultiAssayExperiment@ExperimentList$methylation
    )[[1]][
        match(
            presentMethSites[
                (
                    presentMethSites$meanControl > methCutoff &
                        presentMethSites$meanCase > methCutoff
                ),
                "DNAMethylationSiteID"
            ],
            rownames(
                MultiAssayExperiment::assays(
                    TENETMultiAssayExperiment@ExperimentList$methylation
                )[[1]]
            )
        ),
        methylationSampleNamesCase
    ]

    ## For each case sample, identify if it individually is below the
    ## methCutoff
    sampleCutoffMet <- (methDataSubset < methCutoff)

    ## For each of these RE DNA methylation sites, calculate the number of
    ## individual case samples below the methCutoff
    numOfSamplesMeetingCutoff <- apply(sampleCutoffMet, 1, sum, na.rm = TRUE)

    ## Create a new matrix with the case methylation data for RE DNA
    ## methylation sites that have a number of individual samples below the
    ## methCutoff less than the minCaseCount value. These represent the
    ## methylated RE DNA methylation sites.
    methDataCase <- methDataSubset[names(
        numOfSamplesMeetingCutoff[numOfSamplesMeetingCutoff < minCaseCount]
    ), ]

    ## Create a matrix with the matching control data
    methDataControl <- MultiAssayExperiment::assays(
        TENETMultiAssayExperiment@ExperimentList$methylation
    )[[1]][
        rownames(methDataCase),
        methylationSampleNamesControl
    ]

    ## Get the hypomethylated quadrant RE DNA methylation sites

    ## Get a copy of the case data for RE DNA methylation sites with
    ## average methylation above the methCutoff in the control samples
    ## and all purity datasets (if specified)
    methDataSubset <- MultiAssayExperiment::assays(
        TENETMultiAssayExperiment@ExperimentList$methylation
    )[[1]][
        match(
            presentMethSites[
                .ifelseNoIterate(
                    .isSingleNA(purityData),
                    presentMethSites$meanControl > methCutoff,
                    (
                        presentMethSites$meanControl > methCutoff &
                            presentMethSites$lowestMean > methCutoff
                    )
                ),
                "DNAMethylationSiteID"
            ],
            rownames(
                MultiAssayExperiment::assays(
                    TENETMultiAssayExperiment@ExperimentList$methylation
                )[[1]]
            )
        ),
        methylationSampleNamesCase
    ]

    ## For each case sample, identify if it
    ## individually is below the hypomethCutoff
    sampleCutoffMet <- (methDataSubset < hypomethCutoff)

    ## For each of these RE DNA methylation sites, calculate the number of
    ## individual case samples below the hypomethCutoff
    numOfSamplesMeetingCutoff <- apply(sampleCutoffMet, 1, sum, na.rm = TRUE)

    ## Create a new matrix with the case methylation data for RE DNA
    ## methylation sites that have a number of individual samples below the
    ## hypomethCutoff greater than or equal to the minCaseCount value. These
    ## represent the hypomethylated RE DNA methylation sites.
    hypomethDataCase <- methDataSubset[names(
        numOfSamplesMeetingCutoff[numOfSamplesMeetingCutoff >= minCaseCount]
    ), ]

    hypomethDataControl <- MultiAssayExperiment::assays(
        TENETMultiAssayExperiment@ExperimentList$methylation
    )[[1]][
        rownames(hypomethDataCase),
        methylationSampleNamesControl
    ]

    ## Get the unmethylated quadrant RE DNA methylation sites

    ## Get a copy of the case data for RE DNA methylation sites with average
    ## methylation below the unmethCutoff in control and case samples
    methDataSubset <- MultiAssayExperiment::assays(
        TENETMultiAssayExperiment@ExperimentList$methylation
    )[[1]][
        match(
            presentMethSites[
                (
                    presentMethSites$meanControl < unmethCutoff &
                        presentMethSites$meanCase < unmethCutoff
                ),
                "DNAMethylationSiteID"
            ],
            rownames(
                MultiAssayExperiment::assays(
                    TENETMultiAssayExperiment@ExperimentList$methylation
                )[[1]]
            )
        ),
        methylationSampleNamesCase
    ]

    ## For each case sample, identify if it individually is above the
    ## unmethCutoff
    sampleCutoffMet <- (methDataSubset > unmethCutoff)

    ## For each of these RE DNA methylation sites, calculate the number of
    ## individual case samples above the unmethCutoff
    numOfSamplesMeetingCutoff <- apply(sampleCutoffMet, 1, sum, na.rm = TRUE)

    ## Create a new matrix with the case methylation data for RE DNA
    ## methylation sites that have a number of individual samples below the
    ## unmethCutoff less than the minCaseCount value. These represent the
    ## unmethylated RE DNA methylation sites.
    unmethDataCase <- methDataSubset[names(
        numOfSamplesMeetingCutoff[numOfSamplesMeetingCutoff < minCaseCount]
    ), ]

    ## Create a matrix with the matching control data
    unmethDataControl <- MultiAssayExperiment::assays(
        TENETMultiAssayExperiment@ExperimentList$methylation
    )[[1]][
        rownames(unmethDataCase),
        methylationSampleNamesControl
    ]

    ## Get the hypermethylated quadrant RE DNA methylation sites

    ## Get a copy of the case data for RE DNA methylation sites with
    ## average methylation below the unmethCutoff in control samples
    ## and all purity datasets (if specified)
    methDataSubset <- MultiAssayExperiment::assays(
        TENETMultiAssayExperiment@ExperimentList$methylation
    )[[1]][
        match(
            presentMethSites[
                .ifelseNoIterate(
                    .isSingleNA(purityData),
                    presentMethSites$meanControl < unmethCutoff,
                    (
                        presentMethSites$meanControl < unmethCutoff &
                            presentMethSites$highestMean < unmethCutoff
                    )
                ),
                "DNAMethylationSiteID"
            ],
            rownames(
                MultiAssayExperiment::assays(
                    TENETMultiAssayExperiment@ExperimentList$methylation
                )[[1]]
            )
        ),
        methylationSampleNamesCase
    ]

    ## For each case sample, identify if it individually is above the
    ## hypermethCutoff
    sampleCutoffMet <- (methDataSubset > hypermethCutoff)

    ## For each of these RE DNA methylation sites, calculate the number of
    ## individual case samples above the hypermethCutoff
    numOfSamplesMeetingCutoff <- apply(sampleCutoffMet, 1, sum, na.rm = TRUE)

    ## Create a new matrix with the case methylation data for RE DNA
    ## methylation sites that have a number of individual samples above the
    ## hypermethCutoff greater than or equal to the minCaseCount value. These
    ## represent the hypermethylated RE DNA methylation sites.
    hypermethDataCase <- methDataSubset[names(
        numOfSamplesMeetingCutoff[numOfSamplesMeetingCutoff >= minCaseCount]
    ), ]

    ## Create a matrix with the matching control data
    hypermethDataControl <- MultiAssayExperiment::assays(
        TENETMultiAssayExperiment@ExperimentList$methylation
    )[[1]][
        rownames(hypermethDataCase),
        methylationSampleNamesControl
    ]

    ## Create the final dataset as well as step 2 metadata

    ## Identify the genes and RE DNA methylation sites that have all NA values
    ## in either the control or case expression and methylation samples
    expNAGenesControl <- rownames(
        MultiAssayExperiment::assays(
            TENETMultiAssayExperiment@ExperimentList$expression
        )[[1]][
            apply(
                MultiAssayExperiment::assays(
                    TENETMultiAssayExperiment@ExperimentList$expression
                )[[1]][
                    ,
                    expressionSampleNamesControl,
                    drop = FALSE
                ],
                1,
                function(x) {
                    all(is.na(x))
                }
            ),
        ]
    )

    expNAGenesCase <- rownames(
        MultiAssayExperiment::assays(
            TENETMultiAssayExperiment@ExperimentList$expression
        )[[1]][
            apply(
                MultiAssayExperiment::assays(
                    TENETMultiAssayExperiment@ExperimentList$expression
                )[[1]][
                    ,
                    expressionSampleNamesCase,
                    drop = FALSE
                ],
                1,
                function(x) {
                    all(is.na(x))
                }
            ),
        ]
    )

    sitesWithAllNAsControl <- rownames(
        MultiAssayExperiment::assays(
            TENETMultiAssayExperiment@ExperimentList$methylation
        )[[1]][
            apply(
                MultiAssayExperiment::assays(
                    TENETMultiAssayExperiment@ExperimentList$methylation
                )[[1]][
                    ,
                    methylationSampleNamesControl,
                    drop = FALSE
                ],
                1,
                function(x) {
                    all(is.na(x))
                }
            ),
        ]
    )

    sitesWithAllNAsCase <- rownames(
        MultiAssayExperiment::assays(
            TENETMultiAssayExperiment@ExperimentList$methylation
        )[[1]][
            apply(
                MultiAssayExperiment::assays(
                    TENETMultiAssayExperiment@ExperimentList$methylation
                )[[1]][
                    ,
                    methylationSampleNamesCase,
                    drop = FALSE
                ],
                1,
                function(x) {
                    all(is.na(x))
                }
            ),
        ]
    )

    ## Combine the NA gene and methylation site IDs
    expNAGenes <- unique(c(expNAGenesControl, expNAGenesCase))
    sitesWithAllNAs <- unique(c(sitesWithAllNAsControl, sitesWithAllNAsCase))

    ## If the user has chosen to remove them, identify all the RE DNA
    ## methylation sites whose IDs do not start with "cg"
    if (cgDNAMethylationSitesOnly) {
        nonCGMethSiteIDs <- rownames(
            MultiAssayExperiment::assays(
                TENETMultiAssayExperiment@ExperimentList$methylation
            )[[1]]
        )[
            which(
                !startsWith(
                    rownames(
                        MultiAssayExperiment::assays(
                            TENETMultiAssayExperiment@ExperimentList$
                                methylation
                        )[[1]]
                    ),
                    "cg"
                )
            )
        ]
    }

    ## Get the IDs of all the cg RE DNA methylation sites identified
    REMethylationSites <- presentMethSites$DNAMethylationSiteID
    unmethylatedSites <- rownames(unmethDataControl)
    methylatedSites <- rownames(methDataControl)
    hypomethylatedSites <- rownames(hypomethDataControl)
    hypermethylatedSites <- rownames(hypermethDataControl)

    ## Remove the NA methylation sites from the identified site types
    REMethylationSites <- REMethylationSites[
        !REMethylationSites %in% sitesWithAllNAs
    ]
    unmethylatedSites <- unmethylatedSites[
        !unmethylatedSites %in% sitesWithAllNAs
    ]
    methylatedSites <- methylatedSites[!methylatedSites %in% sitesWithAllNAs]
    hypomethylatedSites <- hypomethylatedSites[
        !hypomethylatedSites %in% sitesWithAllNAs
    ]
    hypermethylatedSites <- hypermethylatedSites[
        !hypermethylatedSites %in% sitesWithAllNAs
    ]

    ## If the user has chosen to do so, remove the RE DNA methylation sites
    ## whose IDs do not start with "cg"
    if (cgDNAMethylationSitesOnly) {
        REMethylationSites <- REMethylationSites[
            !REMethylationSites %in% nonCGMethSiteIDs
        ]
        unmethylatedSites <- unmethylatedSites[
            !unmethylatedSites %in% nonCGMethSiteIDs
        ]
        methylatedSites <- methylatedSites[
            !methylatedSites %in% nonCGMethSiteIDs
        ]
        hypomethylatedSites <- hypomethylatedSites[
            !hypomethylatedSites %in% nonCGMethSiteIDs
        ]
        hypermethylatedSites <- hypermethylatedSites[
            !hypermethylatedSites %in% nonCGMethSiteIDs
        ]
        sitesWithAllNAs <- sitesWithAllNAs[
            !sitesWithAllNAs %in% nonCGMethSiteIDs
        ]
    }

    ## Ensure that there are sites in each quadrant of the analysis which are
    ## most relevant to TENET. If there are not, return an error directing the
    ## user to modify the parameters and/or check the dataset and rerun the
    ## analysis
    if (length(hypermethylatedSites) < 1) {
        .stopNoCall(
            "No hypermethylated RE DNA methylation sites were identified. ",
            "Please modify the parameters and/or check the dataset and rerun ",
            "the analysis."
        )
    }

    if (length(hypomethylatedSites) < 1) {
        .stopNoCall(
            "No hypomethylated RE DNA methylation sites were identified. ",
            "Please modify the parameters and/or check the dataset and rerun ",
            "the analysis."
        )
    }

    ## Create a vector with the number of different types of RE DNA methylation
    ## sites found as well as the total number of genes and RE DNA methylation
    ## sites analyzed and create a metadata table for it
    metadataCounts <- list(
        "unmethylatedSites" = length(unmethylatedSites),
        "hypermethylatedSites" = length(hypermethylatedSites),
        "methylatedSites" = length(methylatedSites),
        "hypomethylatedSites" = length(hypomethylatedSites),
        "REMethylationSites" = length(REMethylationSites),
        "NAGeneCount" = length(expNAGenes),
        "NADNAMethylationSites" = length(sitesWithAllNAs),
        "totalGenes" = nrow(
            MultiAssayExperiment::assays(
                TENETMultiAssayExperiment@ExperimentList$expression
            )[[1]]
        ),
        "totalDNAMethylationSites" = nrow(
            MultiAssayExperiment::assays(
                TENETMultiAssayExperiment@ExperimentList$methylation
            )[[1]]
        )
    )

    ## Also create a list with the identified RE DNA methylation sites of each
    ## type
    methSiteTypeList <- list(
        "unmethylatedSites" = unmethylatedSites,
        "hypermethylatedSites" = hypermethylatedSites,
        "methylatedSites" = methylatedSites,
        "hypomethylatedSites" = hypomethylatedSites,
        "REDNAMethylationSites" = REMethylationSites
    )

    ## Combine these elements into a list and save it in the metadata of the
    ## TENET MultiAssayExperiment object. We will also include the
    ## regulatoryElementGRanges object used in this step.
    if (assessPromoter) {
        TENETMultiAssayExperiment@metadata$
            step2GetDifferentiallyMethylatedSites <- list(
            "unmethCutoff" = unmethCutoff,
            "methCutoff" = methCutoff,
            "hypermethCutoff" = hypermethCutoff,
            "hypomethCutoff" = hypomethCutoff,
            "minCaseCount" = minCaseCount,
            "counts" = metadataCounts,
            "siteIdentitiesList" = methSiteTypeList,
            "regulatoryElementGRanges" = regulatoryElementGRanges,
            "methylationDistributionPlotREPromoter" = distributionPlot,
            "methylationDistributionPlotAllPromoter" =
                distributionPlotAllPromoter,
            "methylationCutoffsPlot" = cutoffsPlot
        )
    } else {
        TENETMultiAssayExperiment@metadata$
            step2GetDifferentiallyMethylatedSites <- list(
            "unmethCutoff" = unmethCutoff,
            "methCutoff" = methCutoff,
            "hypermethCutoff" = hypermethCutoff,
            "hypomethCutoff" = hypomethCutoff,
            "minCaseCount" = minCaseCount,
            "counts" = metadataCounts,
            "siteIdentitiesList" = methSiteTypeList,
            "regulatoryElementGRanges" = regulatoryElementGRanges,
            "methylationDistributionPlot" = distributionPlot,
            "methylationCutoffsPlot" = cutoffsPlot
        )
    }

    return(TENETMultiAssayExperiment)
}
