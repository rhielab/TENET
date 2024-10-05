## Internal functions used by step 1

## Internal function to load bed-like files from a directory and
## return it as a GRanges object
.loadExtBedFiles <- function(extDir, paramName, paramDescription) {
    ## Ensure that the supplied directory exists. If it does, load any .bed,
    ## .narrowPeak, .broadPeak, and/or .gappedPeak files inside.
    extFileList <- .listExtBedFiles(extDir, paramName, paramDescription)

    ## For each of the files, load it then combine the peaks into an extGRanges
    ## object
    extGRanges <- NULL

    for (i in extFileList) {
        ## Load the first three columns of the user's bed-like file
        bedlikeFileGRanges <- rtracklayer::import.bed(
            i,
            colnames = c("chrom", "start", "end")
        )

        ## Add the loaded object to the extGRanges object
        extGRanges <- c(extGRanges, bedlikeFileGRanges)
    }

    ## Return the combined extGRanges object
    return(bedlikeFileGRanges)
}

## Main step 1 function

#' Create a GRanges object representing putative regulatory element regions,
#' based on the data sources selected for inclusion, to be used in later
#' TENET steps
#'
#' This function creates a GRanges object containing regions representing
#' putative regulatory elements, either enhancers or promoters, of
#' interest to the user, based on the presence of specific histone marks
#' and open chromatin/nucleosome-depleted regions. This function can take input
#' from user-specified bed-like files (see
#' <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>) containing regions with
#' histone modification (via the `extHM` argument) and/or open
#' chromatin/nucleosome-depleted regions (via the `extNDR` argument), as well
#' as preprocessed enhancer, promoter, and open chromatin datasets from many
#' cell/tissue types included in the TENET.AnnotationHub repository. The
#' resulting GRanges object will be returned. GRanges objects created by this
#' function can be used by the `step2GetDifferentiallyMethylatedSites` function
#' or other downstream functions. **Note:** Using datasets from
#' TENET.AnnotationHub requires an internet connection, as those datasets are
#' hosted in the Bioconductor AnnotationHub Data Lake.
#'
#' @param extHM To use custom histone modification datasets, specify a path
#' to a directory containing .bed, .narrowPeak, .broadPeak, and/or
#' .gappedPeak files with these datasets. The files may optionally be compressed
#' (.gz/.bz2/.xz). Otherwise, specify NA or do not specify this argument.
#' @param extNDR To use custom open chromatin or NDR datasets, specify a
#' path to a directory containing .bed, .narrowPeak, .broadPeak, and/or
#' .gappedPeak files with these datasets. The files may optionally be compressed
#' (.gz/.bz2/.xz). Otherwise, specify NA or do not specify this argument.
#' @param consensusEnhancer Set to TRUE to use the consensus enhancer data
#' included in TENET.AnnotationHub. Defaults to TRUE.
#' @param consensusPromoter Set to TRUE to use the consensus promoter data
#' included in TENET.AnnotationHub. Defaults to FALSE.
#' @param consensusNDR Set to TRUE to use the consensus open chromatin
#' data included in TENET.AnnotationHub. Defaults to TRUE.
#' @param publicEnhancer Set to TRUE to use the preprocessed publicly available
#' enhancer (H3K27ac) datasets included in TENET.AnnotationHub. If set to TRUE,
#' `cancerType` must be specified. Defaults to FALSE.
#' @param publicPromoter Set to TRUE to use the preprocessed publicly available
#' promoter (H3K4me3) datasets included in TENET.AnnotationHub. If set to TRUE,
#' `cancerType` must be specified. Defaults to FALSE.
#' @param publicNDR Set to TRUE to use the preprocessed publicly available
#' open chromatin (ATAC-seq, DNase-seq) datasets included in
#' TENET.AnnotationHub. If set to TRUE, `cancerType` must be specified.
#' Defaults to FALSE.
#' @param cancerType If `publicEnhancer`, `publicPromoter`, and/or `publicNDR`
#' is TRUE, specify a vector of cancer types from 'BLCA', 'BRCA', 'COAD',
#' 'ESCA', 'HNSC', 'KIRP', 'LIHC', 'LUAD', 'LUSC', and 'THCA' to include the
#' public data relevant to those cancer types. Defaults to NA.
#' @param ENCODEPLS Set to TRUE to use the ENCODE promoter-like elements
#' dataset included in TENET.AnnotationHub. Defaults to FALSE.
#' @param ENCODEpELS Set to TRUE to use the ENCODE proximal enhancer-like
#' elements dataset included in TENET.AnnotationHub. Defaults to FALSE.
#' @param ENCODEdELS Set to TRUE to use the ENCODE distal enhancer-like
#' elements dataset included in TENET.AnnotationHub. Defaults to FALSE.
#' @return Returns the created regulatory element GRanges object.
#' @export
#'
#' @examplesIf interactive()
#' ## This example creates a dataset of putative enhancer regulatory elements
#' ## from consensus datasets and breast invasive carcinoma-relevant sources
#' ## collected in the TENET.AnnotationHub package.
#' returnGRanges <- step1MakeExternalDatasets(
#'     extHM = NA,
#'     extNDR = NA,
#'     publicEnhancer = TRUE,
#'     publicNDR = TRUE,
#'     cancerType = "BRCA",
#'     ENCODEdELS = TRUE
#' )
#'
#' ## This example creates a dataset of putative promoter regulatory elements
#' ## using user provided bed-like files contained in the working
#' ## directory, consensus NDR and promoter regions, and regions with
#' ## promoter-like signatures from the ENCODE SCREEN project. This excludes any
#' ## cancer type-specific public datasets.
#' returnGRanges <- step1MakeExternalDatasets(
#'     extHM = ".",
#'     extNDR = ".",
#'     consensusEnhancer = FALSE,
#'     consensusPromoter = TRUE,
#'     ENCODEPLS = TRUE
#' )
step1MakeExternalDatasets <- function(
    extHM = NA,
    extNDR = NA,
    consensusEnhancer = TRUE,
    consensusPromoter = FALSE,
    consensusNDR = TRUE,
    publicEnhancer = FALSE,
    publicPromoter = FALSE,
    publicNDR = FALSE,
    cancerType = NA,
    ENCODEPLS = FALSE,
    ENCODEpELS = FALSE,
    ENCODEdELS = FALSE) {
    ## If publicEnhancer, publicPromoter, or publicNDR is TRUE, ensure that we
    ## have data for all selected cancer types
    if (any(publicEnhancer, publicPromoter, publicNDR)) {
        if (.isSingleNA(cancerType)) {
            .stopNoCall(
                "The cancerType argument must be specified when using public ",
                "datasets."
            )
        }

        allCancers <- c(
            "BLCA", "BRCA", "COAD", "ESCA", "HNSC", "KIRP", "LIHC", "LUAD",
            "LUSC", "THCA"
        )

        if (!all(cancerType %in% allCancers)) {
            .stopNoCall(
                "At least one of the cancer types specified by the ",
                "cancerType argument is not one of the cancer types for ",
                "which we have created datasets. Please ensure the cancer ",
                "types specified by the cancerType argument are among the ",
                "following: ",
                paste(allCancers, collapse = ", ")
            )
        }
    }

    ## Create a list to store GR objects from all bed files
    TENETGRHMList <- list()
    TENETGRNDRList <- list()
    ENCODEList <- list()

    ## If external HM files are used, validate the directory and load all
    ## bed-like files from it, adding the resulting GRanges object to
    ## TENETGRHMList
    if (!is.na(extHM)) {
        TENETGRHMList <- c(
            TENETGRHMList,
            extHMFiles = .loadExtBedFiles(
                extDir = extHM,
                paramName = "extHM",
                paramDescription = "regions with histone modification"
            )
        )
    }

    ## If external NDR files are used, validate the directory and load all
    ## bed-like files from it, adding the resulting GRanges object to
    ## TENETGRNDRList
    if (!is.na(extNDR)) {
        TENETGRNDRList <- c(
            TENETGRNDRList,
            extNDRFiles = .loadExtBedFiles(
                extDir = extNDR,
                paramName = "extNDR",
                paramDescription = "NDR/open chromatin regions"
            )
        )
    }

    ## Check to see if the user has specified to use any of the datasets that
    ## need to be downloaded from our AnnotationHub package, and if so, load
    ## those datasets
    if (
        any(
            consensusEnhancer, consensusPromoter, consensusNDR,
            publicEnhancer, publicPromoter, publicNDR,
            ENCODEPLS, ENCODEpELS, ENCODEdELS
        )
    ) {
        ## Create an AnnotationHub instance to pull data from
        ah <- AnnotationHub::AnnotationHub()

        ## For each of the datasets that the user might specify, load it from
        ## the AnnotationHub

        ## Consensus enhancer
        if (consensusEnhancer) {
            consensusEnhancerRegions <- .loadFromAnnotationHub(
                ah, "AH116724"
            )

            TENETGRHMList <- c(
                TENETGRHMList,
                consensusEnhancerRegions = consensusEnhancerRegions
            )
        }

        ## Consensus open chromatin
        if (consensusNDR) {
            consensusOpenChromatinRegions <- .loadFromAnnotationHub(
                ah, "AH116725"
            )

            TENETGRNDRList <- c(
                TENETGRNDRList,
                consensusOpenChromatinRegions =
                    consensusOpenChromatinRegions
            )
        }

        ## Consensus promoter
        if (consensusPromoter) {
            consensusPromoterRegions <- .loadFromAnnotationHub(
                ah, "AH116726"
            )

            TENETGRHMList <- c(
                TENETGRHMList,
                consensusPromoterRegions = consensusPromoterRegions
            )
        }

        ## Public (cancer-type specific) enhancers
        if (publicEnhancer) {
            publicEnhancerRegions <- .loadFromAnnotationHub(ah, "AH116721")

            TENETGRHMList <- c(
                TENETGRHMList,
                publicEnhancerRegions = subset(
                    publicEnhancerRegions,
                    get("TYPE") %in% cancerType
                )
            )
        }

        ## Public (cancer-type specific) open chromatin
        if (publicNDR) {
            publicOpenChromatinRegions <- .loadFromAnnotationHub(
                ah, "AH116722"
            )

            TENETGRNDRList <- c(
                TENETGRNDRList,
                publicOpenChromatinRegions = subset(
                    publicOpenChromatinRegions,
                    get("TYPE") %in% cancerType
                )
            )
        }

        ## Public (cancer-type specific) promoter
        if (publicPromoter) {
            publicPromoterRegions <- .loadFromAnnotationHub(ah, "AH116723")

            TENETGRHMList <- c(
                TENETGRHMList,
                publicPromoterRegions = subset(
                    publicPromoterRegions,
                    get("TYPE") %in% cancerType
                )
            )
        }

        ## ENCODE dELS
        if (ENCODEdELS) {
            ENCODEdELSRegions <- .loadFromAnnotationHub(ah, "AH116727")

            ENCODEList <- c(
                ENCODEList,
                ENCODEdELSRegions = ENCODEdELSRegions
            )
        }

        ## ENCODE pELS
        if (ENCODEpELS) {
            ENCODEpELSRegions <- .loadFromAnnotationHub(ah, "AH116728")

            ENCODEList <- c(
                ENCODEList,
                ENCODEpELSRegions = ENCODEpELSRegions
            )
        }

        ## ENCODE PLS
        if (ENCODEPLS) {
            ENCODEPLSRegions <- .loadFromAnnotationHub(ah, "AH116729")

            ENCODEList <- c(
                ENCODEList,
                ENCODEPLSRegions = ENCODEPLSRegions
            )
        }
    }

    ## Put all GR objects together
    ## Warnings are suppressed here because chromosome annotations found in
    ## one GRanges but not the other will return an "Each of the 2 combined
    ## objects has sequence levels not in the other" warning. It's hard to
    ## guarantee there won't be regions present in one that aren't in
    ## the other, and the presence of these regions doesn't affect the
    ## analysis since they won't be overlapped anyway.
    ## Additionally, the documentation for the
    ## function specifically recommends the use of suppressWarnings:
    ## "(use suppressWarnings() to suppress this warning)"
    TENETGRHM <- suppressWarnings(unique(unlist(methods::as(
        TENETGRHMList, "GRangesList"
    ))))
    TENETGRNDR <- suppressWarnings(unique(unlist(methods::as(
        TENETGRNDRList, "GRangesList"
    ))))

    ## Merge peaks using reduce() function
    TENETGRHM <- GenomicRanges::reduce(TENETGRHM)
    TENETGRNDR <- GenomicRanges::reduce(TENETGRNDR)

    ## Intersect NDR with HM
    ## Warnings are suppressed here because chromosome annotations found in
    ## one GRanges but not the other will return an "Each of the 2 combined
    ## objects has sequence levels not in the other" warning. It's hard to
    ## guarantee there won't be regions present in one that aren't in
    ## the other, and the presence of these regions doesn't affect the
    ## analysis since they won't be overlapped anyway.
    ## Additionally, the documentation for the
    ## function specifically recommends the use of suppressWarnings:
    ## "(use suppressWarnings() to suppress this warning)"
    TENETGRNDRHMIntersect <- suppressWarnings(
        GenomicRanges::intersect(TENETGRNDR, TENETGRHM)
    )

    ## Include ENCODEList in the total list and reduced with intersect list
    ## if it is not empty
    if (length(ENCODEList) > 0) {
        ENCODEGR <- unique(unlist(methods::as(ENCODEList, "GRangesList")))
        TENETGRNDRHMIntersect <- GenomicRanges::reduce(
            c(TENETGRNDRHMIntersect, ENCODEGR)
        )
    }

    return(TENETGRNDRHMIntersect)
}
