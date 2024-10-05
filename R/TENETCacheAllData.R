## Based on sesameDataCache.R from the sesameData package by Zhou et al
## (https://github.com/zwdzwd/sesameData, commit
## efac93fedb4d6e74998897bdb3a3391f7f7199ed). sesameData is licensed under the
## Artistic License 2.0, which is GPL-compatible due to its relicensing clause.

#' Cache all online datasets required by TENET examples and optional features
#'
#' This function locally caches all online TENET and SeSAMe datasets required
#' by TENET examples and optional features (TENET.ExperimentHub objects used in
#' examples, TENET.AnnotationHub datasets used in step 1, and SeSAMe datasets
#' loaded via the `DNAMethylationArray` argument). The main purpose of this
#' function is to enable the use of TENET in an HPC cluster environment where
#' compute nodes do not have internet access. In this case, you must run
#' `TENETCacheAllData()` once while connected to the internet before using
#' TENET examples or these optional features.
#'
#' @return Returns NULL.
#' @examples
#' TENETCacheAllData()
#' @export
TENETCacheAllData <- function() {
    message("Caching sesameData datasets if needed")
    tryCatch(
        {
            sesameData::sesameDataCacheAll()
        },
        error = function(cond) {
            .stopNoCall(
                "Failed to cache sesameData datasets. ",
                "Please ensure you have internet access and the latest ",
                "versions of R, sesameData, Bioconductor, and ExperimentHub."
            )
        }
    )

    ## Defined in utils.R
    ahIDs <- .TENETAnnotationHubIDs
    ehIDs <- .TENETExperimentHubIDs

    suppressMessages(try(
        {
            ahIDs <- ahIDs[!(ahIDs %in% names(
                AnnotationHub::AnnotationHub(localHub = TRUE)
            ))]
        },
        silent = TRUE
    ))

    suppressMessages(try(
        {
            ehIDs <- ehIDs[!(ehIDs %in% names(
                ExperimentHub::ExperimentHub(localHub = TRUE)
            ))]
        },
        silent = TRUE
    ))

    ahCount <- length(ahIDs)
    ehCount <- length(ehIDs)

    if (ahCount != 0) {
        ## Temporarily raise the number of maximum downloads so the user is not
        ## asked for confirmation
        ahOldMaxDLs <- AnnotationHub::getAnnotationHubOption("MAX_DOWNLOADS")
        AnnotationHub::setAnnotationHubOption("MAX_DOWNLOADS", ahCount)

        tryCatch(
            {
                ## Cache metadata
                message(
                    "Caching TENET.AnnotationHub metadata (N=", ahCount, ")"
                )
                suppressMessages(ah <- AnnotationHub::query(
                    AnnotationHub::AnnotationHub(), "TENET.AnnotationHub"
                )[ahIDs])

                ## Cache the actual data
                lapply(seq_along(ahIDs), function(i) {
                    id <- ahIDs[i]
                    message(sprintf(
                        "Caching TENET.AnnotationHub datasets (%d/%d): %s (%s)",
                        i, ahCount, names(id), id
                    ))
                    suppressMessages(AnnotationHub::cache(ah[id]))
                })
            },
            error = function(cond) {
                .stopNoCall(
                    "Failed to cache TENET.AnnotationHub datasets. ",
                    "Please ensure you have internet access and the latest ",
                    "versions of R, TENET, Bioconductor, and AnnotationHub."
                )
            },
            finally = AnnotationHub::setAnnotationHubOption(
                "MAX_DOWNLOADS", ahOldMaxDLs
            )
        )
    }

    if (ehCount != 0) {
        ## Temporarily raise the number of maximum downloads so the user is not
        ## asked for confirmation
        ehOldMaxDLs <- ExperimentHub::getExperimentHubOption("MAX_DOWNLOADS")
        ExperimentHub::setExperimentHubOption("MAX_DOWNLOADS", ehCount)

        tryCatch(
            {
                ## Cache metadata
                message(
                    "Caching TENET.ExperimentHub metadata (N=", ehCount, ")"
                )
                ## Testing this. The devel version of ExperimentHub doesn't seem
                ## to have the query function
                suppressMessages(eh <- ExperimentHub::ExperimentHub()[ehIDs])

                ## Cache the actual data
                lapply(seq_along(ehIDs), function(i) {
                    id <- ehIDs[i]
                    message(sprintf(
                        "Caching TENET.ExperimentHub datasets (%d/%d): %s (%s)",
                        i, ehCount, names(id), id
                    ))
                    suppressMessages(ExperimentHub::cache(eh[id]))
                })
            },
            error = function(cond) {
                .stopNoCall(
                    "Failed to cache TENET.ExperimentHub datasets. ",
                    "Please ensure you have internet access and the latest ",
                    "versions of R, TENET, Bioconductor, and ExperimentHub."
                )
            },
            finally = ExperimentHub::setExperimentHubOption(
                "MAX_DOWNLOADS", ehOldMaxDLs
            )
        )
    }
    invisible(NULL)
}
