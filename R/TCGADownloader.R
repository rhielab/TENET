## Note: TCGAbiolinks requires sesameData to download methylation data, but
## does not depend on it for installation

## Internal functions used by TCGADownloader

## See https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/

## Internal function to extract the sample type from one or more TCGA sample IDs
.getSampleType <- function(sampleIDs) {
    return(as.numeric(substring(sampleIDs, 14, 15)))
}

## Internal function to get the TCGA sample IDs in a list that correspond to
## tumor samples
.getTumorSamples <- function(sampleIDs) {
    return(sampleIDs[.getSampleType(sampleIDs) < 10])
}

## Internal function to get the TCGA sample IDs in a list that correspond to
## normal or control samples
.getNormalAndControlSamples <- function(sampleIDs) {
    return(sampleIDs[.getSampleType(sampleIDs) >= 10])
}

## Internal function to extract the patient ID from one or more TCGA sample IDs
.getPatientIDs <- function(sampleIDs) {
    return(substring(sampleIDs, 1, 12))
}

## Internal function to strip the assay and center info from one or more TCGA
## sample IDs, to allow matching them between expression and methylation samples
.stripAssayAndCenterInfo <- function(sampleIDs) {
    return(substring(sampleIDs, 1, 19))
}

#' Download TCGA gene expression, DNA methylation, and clinical datasets
#' and compile them into a MultiAssayExperiment object
#'
#' This function downloads and compiles TCGA gene expression and DNA
#' methylation datasets, as well as clinical data primarily intended for use
#' with the TENET package. This simplifies the TCGAbiolinks download functions,
#' identifies samples with matching gene expression and DNA methylation data,
#' and can also remove duplicate tumor samples taken from the same patient
#' donor. Data are compiled into a MultiAssayExperiment object, which is
#' returned and optionally saved in an .rda file at the path specified by the
#' `outputFile` argument.
#'
#' @param rawDataDownloadDirectory Specify the path to the directory where
#' TCGAbiolinks should download data. Note that this dataset can be very
#' sizable.
#' @param TCGAStudyAbbreviation Input a four-letter code for a TCGA dataset
#' for which to download data. See
#' <https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations>
#' for more information and a complete list of options.
#' @param RNASeqWorkflow Select the type of RNA-seq data to download. For
#' TENET purposes, choose either "STAR - FPKM", "STAR - FPKM-UQ",
#' "STAR - FPKM-UQ - old formula", or "STAR - TPM". "STAR - Counts" can also
#' be used but is not recommended for TENET analyses. See
#' <https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/>
#' for more information on these.
#' @param RNASeqLog2Normalization Set to TRUE to do log2 normalization of
#' RNA-seq expression values. Defaults to TRUE.
#' @param removeDupTumor Set to TRUE to remove duplicate tumor samples
#' taken from the same subject, leaving only one sample per subject in
#' alphanumeric order. **Note:** To properly create a dataset for use with
#' TENET, both the removeDupTumor and matchingExpAndMetSamples arguments must be
#' set to TRUE. Defaults to TRUE.
#' @param matchingExpAndMetSamples Select the type of expression and
#' methylation sample data matching to perform. If set to TRUE, only samples
#' with at least one methylation and expression sample annotated to their
#' patient, will be kept. If set to FALSE, all samples will be kept, including
#' those without matching expression and methylation data. **Note:** To properly
#' create a dataset for use with TENET, both the removeDupTumor and
#' matchingExpAndMetSamples arguments must be set to TRUE. Defaults to TRUE.
#' @param clinicalSurvivalData Select how clinical data should be prepared
#' from the TCGA data, with respect to patient vital status and survival time.
#' Valid options include "bcrBiotabPatient" to use survival data contained
#' only in the 'patient' data in the BCR Biotab files downloaded using
#' TCGAbiolinks, or "combined", which uses clinical information from the
#' 'patient' and 'follow_up' datasets in the BCR Biotab files, as well as data
#' from the BCR XML files. Data from the same patient in each of the datasets
#' are combined, and the data with the most recent (highest patient survival
#' time) entry for each patient are kept. Additionally, for both options, the
#' 'days_to_last_followup' and 'days_to_death' variables are collapsed into a
#' single time variable, which is combined with the other patient clinical data
#' in the 'patient' BCR Biotab data. See
#' <https://bioconductor.org/packages/devel/bioc/vignettes/TCGAbiolinks/inst/doc/clinical.html>
#' for more information on how TCGAbiolinks prepares different clinical
#' datasets. Defaults to "combined".
#' @param outputFile Specify the path to an .rda file in which to save the
#' MultiAssayExperiment object with downloaded datasets. If set to NA or
#' undefined, this results in the function only returning the
#' MultiAssayExperiment object and not saving it. Defaults to NA.
#' @return Returns and/or saves to an .rda file a MultiAssayExperiment
#' object with expression and methylation data included SummarizedExperiment
#' objects within the MultiAssayExperiment object, as well as clinical data
#' included in the colData of the MultiAssayExperiment object.
#' @export
#'
#' @examplesIf interactive()
#' ## Download a TCGA LUAD dataset with log2-normalized
#' ## FPKM-UQ expression values from tumor and adjacent normal tissue samples
#' ## with matching expression and methylation data and keeping only one tumor
#' ## sample from each patient. Additionally, survival data will be combined
#' ## from three clinical datasets downloaded by TCGAbiolinks. Raw data files
#' ## will be saved to the working directory, and the processed dataset will
#' ## be returned as a variable.
#' TCGADataset <- TCGADownloader(
#'     rawDataDownloadDirectory = ".",
#'     TCGAStudyAbbreviation = "LUAD",
#'     RNASeqWorkflow = "STAR - FPKM-UQ"
#' )
#'
#' ## Another example, which downloads a TCGA BRCA dataset with FPKM expression
#' ## values with no normalization and no duplicate samples removed. Survival
#' ## data are derived from just the patient BCR Biotab file downloaded by
#' ## TCGAbiolinks. Both raw data files and an .rda file containing the data
#' ## as a MultiAssayExperiment object will be saved to the working directory.
#' ## Note: This functionality is useful for downloading samples from
#' ## TCGA but will *not* work for a TENET assay due to the lack of sample
#' ## matching and duplicate tumor sample removal.
#' TCGADownloader(
#'     rawDataDownloadDirectory = ".",
#'     TCGAStudyAbbreviation = "BRCA",
#'     RNASeqWorkflow = "STAR - FPKM",
#'     RNASeqLog2Normalization = FALSE,
#'     removeDupTumor = FALSE,
#'     matchingExpAndMetSamples = FALSE,
#'     clinicalSurvivalData = "bcrBiotabPatient",
#'     outputFile = "BRCAMultiAssayExperimentObject.rda"
#' )
TCGADownloader <- function(
    rawDataDownloadDirectory,
    TCGAStudyAbbreviation,
    RNASeqWorkflow,
    RNASeqLog2Normalization = TRUE,
    removeDupTumor = TRUE,
    matchingExpAndMetSamples = TRUE,
    clinicalSurvivalData = "combined",
    outputFile = NA) {
    ## Ensure that rawDataDownloadDirectory is defined
    if (missing(rawDataDownloadDirectory)) {
        .stopNoCall(
            "rawDataDownloadDirectory is not defined. Please set it ",
            "to the directory where raw files from TCGA should be ",
            "downloaded by TCGAbiolinks."
        )
    }

    ## Ensure that a study abbreviation has been provided
    if (missing(TCGAStudyAbbreviation)) {
        .stopNoCall(
            "TCGAStudyAbbreviation is not defined. Please set it to the ",
            "TCGA study abbreviation for which to download data (as listed ",
            "in https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/",
            "tcga-study-abbreviations)."
        )
    }

    ## Ensure that a proper RNA-seq workflow type has been specified
    validWorkflows <- c(
        "STAR - FPKM",
        "STAR - FPKM-UQ",
        "STAR - FPKM-UQ - old formula",
        "STAR - TPM",
        "STAR - Counts"
    )

    if (missing(RNASeqWorkflow) || !(RNASeqWorkflow %in% validWorkflows)) {
        .stopNoCall(
            "RNASeqWorkflow is not properly defined. Please set it to ",
            'one of the following: "', paste(validWorkflows, collapse = '", "'),
            '"'
        )
    }

    ## By default, leave the expression values alone
    expressionMultiplier <- 1
    if (RNASeqWorkflow == "STAR - FPKM") {
        expressionColumn <- "fpkm_unstrand"
    } else if (RNASeqWorkflow == "STAR - FPKM-UQ") {
        expressionColumn <- "fpkm_uq_unstrand"
    } else if (RNASeqWorkflow == "STAR - FPKM-UQ - old formula") {
        expressionColumn <- "fpkm_uq_unstrand"
        ## Multiply the expressionData values by 19,029 to
        ## correct for FPKM UQ values similar to ones TCGA used to use
        ## prior to Data Release 37.0 - March 29, 2023
        ## See
        ## https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
        ## "Number of protein coding genes on autosomes: 19,029"
        expressionMultiplier <- 19029
    } else if (RNASeqWorkflow == "STAR - TPM") {
        expressionColumn <- "tpm_unstrand"
    } else {
        ## STAR - Counts
        expressionColumn <- "unstrand"
    }

    ## Ensure that a proper clinical survival data type has been specified
    if (!(clinicalSurvivalData %in% c("bcrBiotabPatient", "combined"))) {
        .stopNoCall(
            "clinicalSurvivalData is not properly defined. Please set it to ",
            'either "bcrBiotabPatient" or "combined"'
        )
    }

    ## Create the all-caps TCGA abbreviation used for downloads
    TCGAStudyAbbreviationDownload <- paste0(
        "TCGA-", toupper(TCGAStudyAbbreviation)
    )

    ## Set up the expression query
    expressionQuery <- TCGAbiolinks::GDCquery(
        project = TCGAStudyAbbreviationDownload,
        data.category = "Transcriptome Profiling",
        data.type = "Gene Expression Quantification",
        experimental.strategy = "RNA-Seq",
        workflow.type = "STAR - Counts"
    )

    ## Download the expression data
    TCGAbiolinks::GDCdownload(
        expressionQuery,
        directory = rawDataDownloadDirectory
    )

    ## Assemble the expression SummarizedExperiment object
    expressionExperiment <- TCGAbiolinks::GDCprepare(
        expressionQuery,
        directory = rawDataDownloadDirectory
    )

    message("Expression SummarizedExperiment assembled.")

    ## Set up the methylation query
    methylationQuery <- TCGAbiolinks::GDCquery(
        project = TCGAStudyAbbreviationDownload,
        data.category = "DNA Methylation",
        data.type = "Methylation Beta Value",
        platform = "Illumina Human Methylation 450"
    )

    ## Download the methylation data
    TCGAbiolinks::GDCdownload(
        methylationQuery,
        directory = rawDataDownloadDirectory
    )

    ## Assemble the methylation SummarizedExperiment object
    methylationExperiment <- TCGAbiolinks::GDCprepare(
        methylationQuery,
        directory = rawDataDownloadDirectory
    )

    message("Methylation SummarizedExperiment assembled.")

    ## Set up the clinical query for the BCR Biotab data
    clinicalQueryBcrBiotab <- TCGAbiolinks::GDCquery(
        project = TCGAStudyAbbreviationDownload,
        data.category = "Clinical",
        data.type = "Clinical Supplement",
        data.format = "BCR Biotab"
    )

    ## Download the BCR Biotab data
    TCGAbiolinks::GDCdownload(
        clinicalQueryBcrBiotab,
        directory = rawDataDownloadDirectory
    )

    ## Load the BCR Biotab clinical data into a data frame
    clinicalDataBcrBiotab <- TCGAbiolinks::GDCprepare(
        clinicalQueryBcrBiotab,
        directory = rawDataDownloadDirectory
    )

    ## The BCR Biotab data have a few different sub-datasets, so we will
    ## get only the patient datasets
    clinicalDataBcrBiotabPatient <- clinicalDataBcrBiotab[[
        grep("patient", names(clinicalDataBcrBiotab))
    ]]

    ## Organize the datasets properly by extracting the true column
    ## names from the first row, then removing the first two rows of
    ## now irrelevant data
    colnames(clinicalDataBcrBiotabPatient) <- unname(
        unlist(clinicalDataBcrBiotabPatient[1, ])
    )

    clinicalDataBcrBiotabPatient <- clinicalDataBcrBiotabPatient[
        seq(from = 3, to = nrow(clinicalDataBcrBiotabPatient), by = 1),
    ]

    ## Add a time variable
    clinicalDataBcrBiotabPatient$time <- ifelse(
        clinicalDataBcrBiotabPatient$vital_status == "Alive",
        clinicalDataBcrBiotabPatient$days_to_last_followup,
        clinicalDataBcrBiotabPatient$days_to_death
    )

    ## Remove the tribble class from the data frame
    clinicalDataBcrBiotabPatient <- as.data.frame(clinicalDataBcrBiotabPatient)

    ## Add rownames to the dataset
    rownames(clinicalDataBcrBiotabPatient) <- clinicalDataBcrBiotabPatient$
        bcr_patient_barcode

    ## Convert factors to characters
    for (i in seq_len(ncol(clinicalDataBcrBiotabPatient))) {
        if (is.factor(clinicalDataBcrBiotabPatient[, i])) {
            clinicalDataBcrBiotabPatient[, i] <- as.character(
                clinicalDataBcrBiotabPatient[, i]
            )
        }
    }

    if (clinicalSurvivalData == "combined") {
        ## TCGAbiolinks downloads various types of clinical data for patients.
        ## Some of these datasets have more updated patient survival
        ## information, but don't contain data for all patients. Thus, the
        ## strategy here is to grab all the relevant datasets, then extract the
        ## most recent data for each patient where they are available.

        ## Get the BCR XML data, which are useful if the 'combined' survival
        ## option is selected
        clinicalQueryBcrXml <- TCGAbiolinks::GDCquery(
            project = TCGAStudyAbbreviationDownload,
            data.category = "Clinical",
            data.format = "bcr xml"
        )

        ## Download the BCR XML data
        TCGAbiolinks::GDCdownload(
            clinicalQueryBcrXml,
            directory = rawDataDownloadDirectory
        )

        ## Load the BCR XML data as a data frame
        clinicalDataBcrXml <- TCGAbiolinks::GDCprepare_clinic(
            clinicalQueryBcrXml,
            clinical.info = "patient",
            directory = rawDataDownloadDirectory
        )

        ## The BCR Biotab data have a few different sub-datasets, so we will
        ## get only the follow_up datasets
        clinicalDataFollowUpNames <- grep(
            "follow_up", names(clinicalDataBcrBiotab),
            value = TRUE
        )

        ## We don't want any _nte_ datasets, as those seem to lack clinical data
        ## in the follow_up datasets
        clinicalDataFollowUpNames <- grep(
            "_nte_", clinicalDataFollowUpNames,
            value = TRUE, invert = TRUE
        )

        ## Sort the names so the most recent of the non-nte datasets is last;
        ## it will be the BCR Biotab follow-up dataset
        clinicalDataFollowUpNames <- sort(clinicalDataFollowUpNames)

        clinicalDataBcrBiotabFollowUp <- clinicalDataBcrBiotab[[
            clinicalDataFollowUpNames[length(clinicalDataFollowUpNames)]
        ]]

        ## Organize the BCR Biotab follow-up data properly by extracting
        ## the true column names from the first row, then removing the first
        ## two rows of now irrelevant data
        colnames(clinicalDataBcrBiotabFollowUp) <- unname(
            unlist(clinicalDataBcrBiotabFollowUp[1, ])
        )

        clinicalDataBcrBiotabFollowUp <- clinicalDataBcrBiotabFollowUp[
            seq(from = 3, to = nrow(clinicalDataBcrBiotabFollowUp), by = 1),
        ]

        ## Add a time variable to the datasets
        clinicalDataBcrXml$time <- ifelse(
            clinicalDataBcrXml$vital_status == "Alive",
            clinicalDataBcrXml$days_to_last_followup,
            clinicalDataBcrXml$days_to_death
        )

        clinicalDataBcrBiotabFollowUp$time <- ifelse(
            clinicalDataBcrBiotabFollowUp$vital_status == "Alive",
            clinicalDataBcrBiotabFollowUp$days_to_last_followup,
            clinicalDataBcrBiotabFollowUp$days_to_death
        )

        ## Remove duplicates from the BCR Biotab follow-up data.
        ## suppressWarnings() is here due to use of as.numeric, which
        ## intentionally induces NAs to non-numeric values of the time
        ## variable.
        suppressWarnings(
            clinicalDataBcrBiotabFollowUp <- clinicalDataBcrBiotabFollowUp[
                order(
                    as.numeric(clinicalDataBcrBiotabFollowUp$time),
                    decreasing = TRUE
                ),
            ]
        )

        clinicalDataBcrBiotabFollowUpDedup <- clinicalDataBcrBiotabFollowUp[
            !duplicated(clinicalDataBcrBiotabFollowUp$bcr_patient_barcode),
        ]

        ## Remove the tribble class from the BiotabFollowUpDedup data frame (the
        ## BcrXml dataset is not a tribble)
        clinicalDataBcrBiotabFollowUpDedup <- as.data.frame(
            clinicalDataBcrBiotabFollowUpDedup
        )

        ## Add rownames to the datasets
        rownames(clinicalDataBcrXml) <- clinicalDataBcrXml$bcr_patient_barcode
        rownames(clinicalDataBcrBiotabFollowUpDedup) <-
            clinicalDataBcrBiotabFollowUpDedup$bcr_patient_barcode

        ## Convert factors to characters
        for (i in seq_len(ncol(clinicalDataBcrXml))) {
            if (is.factor(clinicalDataBcrXml[, i])) {
                clinicalDataBcrXml[, i] <- as.character(clinicalDataBcrXml[, i])
            }
        }

        for (i in seq_len(ncol(clinicalDataBcrBiotabFollowUpDedup))) {
            if (is.factor(clinicalDataBcrBiotabFollowUpDedup[, i])) {
                clinicalDataBcrBiotabFollowUpDedup[, i] <- as.character(
                    clinicalDataBcrBiotabFollowUpDedup[, i]
                )
            }
        }

        ## Create a combined list of patient samples from all the relevant
        ## clinical datasets
        clinicalSamplesCombined <- unique(c(
            rownames(clinicalDataBcrBiotabPatient),
            rownames(clinicalDataBcrBiotabFollowUpDedup),
            rownames(clinicalDataBcrXml)
        ))

        ## For each clinical sample, go through the different
        ## clinical datasets and pull info in a prioritized manner
        ## i.e. the dataset which has the most recent (latest time)
        ## for that patient
        vitalStatusVector <- NULL
        timeVector <- NULL

        for (i in clinicalSamplesCombined) {
            ## Get the time and vital status from the BCR Biotab patient data
            if (i %in% rownames(clinicalDataBcrBiotabPatient)) {
                bcrBiotabPatientStatus <-
                    clinicalDataBcrBiotabPatient[i, "vital_status"]
                bcrBiotabPatientTime <- clinicalDataBcrBiotabPatient[i, "time"]
            } else {
                bcrBiotabPatientStatus <- NA
                bcrBiotabPatientTime <- NA
            }

            ## Get the time and vital status from the BCR Biotab follow-up data
            if (i %in% rownames(clinicalDataBcrBiotabFollowUpDedup)) {
                bcrBiotabFollowUpStatus <-
                    clinicalDataBcrBiotabFollowUpDedup[i, "vital_status"]
                bcrBiotabFollowUpTime <-
                    clinicalDataBcrBiotabFollowUpDedup[i, "time"]
            } else {
                bcrBiotabFollowUpStatus <- NA
                bcrBiotabFollowUpTime <- NA
            }

            ## Get the time and vital status from the BCR XML data
            if (i %in% rownames(clinicalDataBcrXml)) {
                bcrXmlStatus <- clinicalDataBcrXml[i, "vital_status"]
                bcrXmlTime <- clinicalDataBcrXml[i, "time"]
            } else {
                bcrXmlStatus <- NA
                bcrXmlTime <- NA
            }

            ## Combine the values with names
            vitalStatusCombined <- c(
                bcrBiotabFollowUpStatus, bcrBiotabPatientStatus, bcrXmlStatus
            )

            ## Combine the time values. Again, suppressWarnings is here as we
            ## want to convert non-numeric values to NA.
            timeCombined <- suppressWarnings(
                as.numeric(
                    c(bcrBiotabFollowUpTime, bcrBiotabPatientTime, bcrXmlTime)
                )
            )

            ## Get the data with the largest time value
            if (all(is.na(timeCombined))) {
                maxTime <- NA
                maxTimeStatus <- NA
            } else {
                maxTime <- max(timeCombined, na.rm = TRUE)
                maxTimeStatus <- vitalStatusCombined[
                    which(timeCombined == max(timeCombined, na.rm = TRUE))
                ]
            }

            ## Reduce the max- vectors to a single entry
            ## unless there is conflict on whether the
            ## patient is alive or dead with the same time,
            ## in which case, return NA due to ambiguity
            if (length(maxTimeStatus) > 1) {
                if (length(unique(maxTimeStatus)) > 1) {
                    maxTime <- NA
                    maxTimeStatus <- NA
                } else {
                    maxTimeStatus <- unique(maxTimeStatus)
                }
            }

            ## Add the resulting values to the vectors
            vitalStatusVector <- c(vitalStatusVector, maxTimeStatus)
            timeVector <- c(timeVector, maxTime)
        }

        ## Create a data frame of these new values.
        ## We will keep snake_case variable names here so they match with the
        ## variable names we would see directly out of TCGAbiolinks functions.
        survivalDF <- data.frame(
            "vital_status" = vitalStatusVector,
            "time" = timeVector,
            stringsAsFactors = FALSE
        )
        rownames(survivalDF) <- clinicalSamplesCombined

        ## Add the values to the clinicalDataBcrBiotabPatient DF
        clinicalDataBcrBiotabPatient$vital_status <- survivalDF[
            rownames(clinicalDataBcrBiotabPatient),
            "vital_status"
        ]

        clinicalDataBcrBiotabPatient$time <- survivalDF[
            rownames(clinicalDataBcrBiotabPatient),
            "time"
        ]
    }

    ## Use the final clinicalDataBcrBiotabPatient dataset as the clinicalData
    clinicalData <- clinicalDataBcrBiotabPatient

    message("Clinical dataset assembled.")

    ## Multiply the expressionData values by the expressionMultiplier value
    SummarizedExperiment::assays(
        expressionExperiment
    )[[expressionColumn]] <- SummarizedExperiment::assays(
        expressionExperiment
    )[[expressionColumn]] * expressionMultiplier

    ## Perform log2 normalization if selected
    if (RNASeqLog2Normalization) {
        SummarizedExperiment::assays(expressionExperiment)[[
            expressionColumn
        ]] <- log2(
            SummarizedExperiment::assays(
                expressionExperiment
            )[[expressionColumn]] + 1 ## + 1 to avoid log2(0) which is -Inf
        )

        message("Expression dataset log2 transformed.")
    }

    ## Strip the period and trailing number(s) annotation from the gene
    ## Ensembl IDs in the rowRanges
    names(SummarizedExperiment::rowRanges(expressionExperiment)) <- sub(
        "_[0-9]+$", "", sub(
            "\\.[0-9]+", "",
            names(SummarizedExperiment::rowRanges(expressionExperiment))
        )
    )

    ## Strip the period and trailing number(s) annotation from the gene
    ## Ensembl IDs in the main body of expression data
    rownames(
        SummarizedExperiment::assays(expressionExperiment)[[expressionColumn]]
    ) <- sub("_[0-9]+$", "", sub(
        "\\.[0-9]+", "",
        rownames(
            SummarizedExperiment::assays(
                expressionExperiment
            )[[expressionColumn]]
        )
    ))

    ## To expedite later steps, isolate the expression and methylation sample
    ## names
    expNames <- colnames(
        SummarizedExperiment::assays(expressionExperiment)[[expressionColumn]]
    )

    metNames <- colnames(
        SummarizedExperiment::assays(methylationExperiment)[[1]]
    )

    ## If removeDupTumor is set to TRUE, remove the tumor samples that are
    ## duplicates from the same patient (leaving one per patient)
    if (removeDupTumor) {
        ## Identify the matched expression samples which are from tumor and
        ## normal samples
        expNamesTumor <- .getTumorSamples(expNames)
        expNamesNormal <- .getNormalAndControlSamples(expNames)

        ## Do the same for the matched methylation samples which are from tumor
        ## and normal samples
        metNamesTumor <- .getTumorSamples(metNames)
        metNamesNormal <- .getNormalAndControlSamples(metNames)

        ## Create data frames from the tumor IDs to help in identifying the
        ## tumor samples that come from duplicate patients
        dupTumorExpressionDF <- data.frame(
            "fullID" = expNamesTumor,
            stringsAsFactors = FALSE
        )

        rownames(dupTumorExpressionDF) <- expNamesTumor

        dupTumorMethylationDF <- data.frame(
            "fullID" = metNamesTumor,
            stringsAsFactors = FALSE
        )

        rownames(dupTumorMethylationDF) <- metNamesTumor

        ## Extract the patient IDs from the sample IDs
        dupTumorExpressionDF$patientID <- .getPatientIDs(
            dupTumorExpressionDF$fullID
        )

        dupTumorMethylationDF$patientID <- .getPatientIDs(
            dupTumorMethylationDF$fullID
        )

        ## Sort the data frames alphanumerically to facilitate in keeping the
        ## first sample from each patient based on its alphanumeric ID
        dupTumorExpressionDF <- dupTumorExpressionDF[sort(expNamesTumor), ]
        dupTumorMethylationDF <- dupTumorMethylationDF[sort(metNamesTumor), ]

        ## Identify which of the patient samples are duplicates
        dupTumorExpressionDF$duplicatedPatientID <- duplicated(
            dupTumorExpressionDF$patientID
        )

        dupTumorMethylationDF$duplicatedPatientID <- duplicated(
            dupTumorMethylationDF$patientID
        )

        ## Remove all the duplicated tumor samples, isolate the IDs of the
        ## remaining tumor samples, then add the normal sample names back to
        ## those as ones to keep
        expNames <- c(
            rownames(dupTumorExpressionDF[
                dupTumorExpressionDF$duplicatedPatientID == FALSE,
            ]),
            expNamesNormal
        )

        metNames <- c(
            rownames(
                dupTumorMethylationDF[
                    dupTumorMethylationDF$duplicatedPatientID == FALSE,
                ]
            ),
            metNamesNormal
        )
    }

    ## If matchingExpAndMetSamples is set to TRUE, find both
    ## tumor and normal samples with just matched methylation and expression
    ## data
    if (matchingExpAndMetSamples) {
        ## Find the sample names that are present in both datasets
        matchedExpMetNamesSub <- intersect(
            .stripAssayAndCenterInfo(expNames),
            .stripAssayAndCenterInfo(metNames)
        )

        matchedExpNames <- expNames[
            .stripAssayAndCenterInfo(expNames) %in% matchedExpMetNamesSub
        ]

        matchedMetNames <- metNames[
            .stripAssayAndCenterInfo(metNames) %in% matchedExpMetNamesSub
        ]
    } else {
        ## Find the names of all the samples with either expression or
        ## methylation data
        matchedExpNames <- expNames
        matchedMetNames <- metNames
    }

    ## Sort the expression and methylation names alphanumerically (for
    ## later use)
    matchedExpNames <- sort(matchedExpNames)
    matchedMetNames <- sort(matchedMetNames)

    ## Process the clinical data so it can be added as the colData of the
    ## MultiAssayExperiment object

    ## Combine the exp and met names and sort those as well
    matchedExpMetNames <- sort(c(matchedExpNames, matchedMetNames))

    ## Get the unique samples then add in the patient IDs as the row names
    clinical <- unique(clinicalData)
    rownames(clinical) <- clinical$bcr_patient_barcode

    ## Get the patient IDs from the exp and met sample names, which will
    ## match the patient barcode in the clinical data
    expMetClinicalBarcodeMatch <- unique(.getPatientIDs(matchedExpMetNames))

    ## Get a data frame of the clinical data for the patients in the dataset.
    ## Rownames are reset so missing expMetClinicalBarcodeMatch values are
    ## still reflected in the rownames of the multiAssayClinical object with NA
    ## values, instead of the missing values also having NA rownames.
    multiAssayClinical <- clinical[expMetClinicalBarcodeMatch, ]
    rownames(multiAssayClinical) <- expMetClinicalBarcodeMatch

    ## Change the time variable back to numeric in the clinical data
    multiAssayClinical$time <- as.numeric(multiAssayClinical$time)

    ## Set up the mapping list we will need to create the
    ## MultiAssayExperiment object and link expression/methylation samples
    ## to the clinical colData, along with the sampleType information

    ## Create the data frames for the expression/methylation maps
    expressionMap <- data.frame(
        "primary" = .getPatientIDs(matchedExpNames),
        "colname" = matchedExpNames,
        "sampleType" = ifelse(
            .getSampleType(matchedExpNames) < 10, "Case", "Control"
        ),
        stringsAsFactors = FALSE
    )

    methylationMap <- data.frame(
        "primary" = .getPatientIDs(matchedMetNames),
        "colname" = matchedMetNames,
        "sampleType" = ifelse(
            .getSampleType(matchedMetNames) < 10, "Case", "Control"
        ),
        stringsAsFactors = FALSE
    )

    ## Create a list of those data frames with the name of the eventual
    ## components of the MultiAssayExperiment object
    mappingList <- list(expressionMap, methylationMap)
    names(mappingList) <- c("expression", "methylation")

    ## Create the MultiAssayExperiment object
    TENETMultiAssayExperiment <- MultiAssayExperiment::MultiAssayExperiment(
        experiments = list(
            "expression" = SummarizedExperiment::SummarizedExperiment(
                assays = list(
                    "expression" = SummarizedExperiment::assays(
                        expressionExperiment
                    )[[expressionColumn]][
                        ,
                        matchedExpNames
                    ]
                ),
                rowRanges = SummarizedExperiment::rowRanges(
                    expressionExperiment
                )
            ),
            "methylation" = SummarizedExperiment::SummarizedExperiment(
                assays = list(
                    "methylation" = SummarizedExperiment::assays(
                        methylationExperiment
                    )[[1]][
                        ,
                        matchedMetNames
                    ]
                ),
                rowRanges = SummarizedExperiment::rowRanges(
                    methylationExperiment
                )
            )
        ),
        colData = multiAssayClinical,
        sampleMap = MultiAssayExperiment::listToMap(mappingList)
    )

    ## For whatever reason, when setting the sampleMap of the
    ## TENETMultiAssayExperiment above, it has lost the "sampleType"
    ## column that was added, so we manually replace the sampleMap in the object
    TENETMultiAssayExperiment@sampleMap <-
        MultiAssayExperiment::listToMap(mappingList)

    ## Create a metadata data frame with info about how the samples were
    ## processed
    metadataDF <- data.frame(
        "dataset" = TCGAStudyAbbreviationDownload,
        "expressionDataRelease" = unname(
            unlist(expressionExperiment@metadata)
        ),
        "methylationDataRelease" = unname(
            unlist(methylationExperiment@metadata)
        ),
        "RNASeqWorkflow" = RNASeqWorkflow,
        "log2ExpressionNormalization" = RNASeqLog2Normalization,
        "duplicateTumorSamplesFromSamePatientRemoved" = removeDupTumor,
        "expressionAndMethylationSampleMatching" = matchingExpAndMetSamples,
        "clinicalSurvivalData" = clinicalSurvivalData,
        "DNAMethylationDataType" = "DNAMethylationArray",
        stringsAsFactors = FALSE
    )

    ## Add that metadata data frame to the MultiAssayExperiment as a list so
    ## more data can be added in later TENET steps
    TENETMultiAssayExperiment@metadata <- list(
        "TCGADownloaderFunction" = metadataDF
    )

    ## If outputFile is specified, save the object to the .rda file
    if (!is.na(outputFile)) {
        save(TENETMultiAssayExperiment, file = outputFile)
    }

    ## Also return the object
    return(TENETMultiAssayExperiment)
}
