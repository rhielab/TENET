## Internal function to validate an external bed file directory and return a
## list of the files in it
.listExtBedFiles <- function(extDir, paramName, paramDescription) {
    ## Ensure that the supplied directory exists. If it does, load any .bed,
    ## .narrowPeak, .broadPeak, and/or .gappedPeak files inside.
    ## Otherwise, return an error.
    if (!dir.exists(extDir)) {
        .stopNoCall(
            "The directory path given for ", paramName, " was not found. ",
            "Please ensure that the argument specifies a path to a ",
            "directory containing .bed, .narrowPeak, .broadPeak, and/or ",
            ".gappedPeak files containing ", paramDescription, ". ",
            "The files may optionally be compressed (.gz/.bz2/.xz)."
        )
    }

    ## List all the external files found
    extFileList <- list.files(
        path = extDir,
        pattern = "\\.(bed|(narrow|broad|gapped)Peak)(\\.gz|bz2|xz)?$",
        ignore.case = TRUE,
        full.names = TRUE
    )

    ## If there are no files return an error
    if (length(extFileList) == 0) {
        .stopNoCall(
            "The directory path given for ", paramName, " does not ",
            "contain any .bed, .narrowPeak, .broadPeak, and/or .gappedPeak ",
            "files. Please ensure that the argument specifies a path ",
            "to a directory containing .bed, .narrowPeak, .broadPeak, ",
            "and/or .gappedPeak files containing ", paramDescription, ". ",
            "The files may optionally be compressed (.gz/.bz2/.xz)."
        )
    }

    return(extFileList)
}

## Internal function to validate the analysis type arguments and return a
## vector of the selected analysis types
.validateAnalysisTypes <- function(hypermethAnalysis, hypomethAnalysis) {
    analysisTypes <- NULL
    if (!any(hypermethAnalysis, hypomethAnalysis)) {
        .stopNoCall(
            "All analysis types have been set to FALSE. Set at least one ",
            "analysis type to TRUE."
        )
    }
    if (hypermethAnalysis) {
        analysisTypes <- c(analysisTypes, "hyper")
    }
    if (hypomethAnalysis) {
        analysisTypes <- c(analysisTypes, "hypo")
    }
    return(analysisTypes)
}

## Internal function to validate the input MultiAssayExperiment and return an
## error message if it is invalid
.validateMultiAssayExperiment <- function(
    MAE,
    needGeneName = FALSE) {
    ## Ensure that the user has provided a properly formatted
    ## MultiAssayExperiment object
    if (!inherits(MAE, "MultiAssayExperiment")) {
        .stopNoCall(
            "The object given as the TENETMultiAssayExperiment does not ",
            "seem to be a MultiAssayExperiment object. Please ensure the ",
            "TENETMultiAssayExperiment is a MultiAssayExperiment object."
        )
    }

    ## Ensure that there are expression and methylation objects present; if not,
    ## stop the function and return an error message
    if (!"expression" %in% names(MAE@ExperimentList)) {
        .stopNoCall(
            "There is no SummarizedExperiment named 'expression' in the ",
            "TENETMultiAssayExperiment object. Please ensure a ",
            "SummarizedExperiment object named 'expression' is included in ",
            "the TENETMultiAssayExperiment object with expression data from ",
            "matched samples with methylation data in the 'methylation' ",
            "SummarizedExperiment object."
        )
    }
    if (!"methylation" %in% names(MAE@ExperimentList)) {
        .stopNoCall(
            "There is no SummarizedExperiment named 'methylation' in the ",
            "TENETMultiAssayExperiment object. Please ensure a ",
            "SummarizedExperiment object named 'methylation' is included in ",
            "the TENETMultiAssayExperiment object with methylation data from ",
            "matched samples with expression data in the 'expression' ",
            "SummarizedExperiment object."
        )
    }

    ## Ensure that the rowRanges are provided in the expression and methylation
    ## objects as a GRanges object; if not, stop the function and return an
    ## error message
    if (!inherits(
        SummarizedExperiment::rowRanges(MAE@ExperimentList$expression),
        "GRanges"
    )) {
        .stopNoCall(
            "A GRanges object with coordinates of genes has not been given as ",
            "the rowRanges for the 'expression' SummarizedExperiment in the ",
            "TENETMultiAssayExperiment object."
        )
    }
    if (!inherits(
        SummarizedExperiment::rowRanges(MAE@ExperimentList$methylation),
        "GRanges"
    )) {
        .stopNoCall(
            "A GRanges object with coordinates of RE DNA methylation sites ",
            "has not been given as the rowRanges for the 'methylation' ",
            "SummarizedExperiment in the TENETMultiAssayExperiment object."
        )
    }

    ## Ensure that "gene_name" is present in the expression metadata and return
    ## an error if not
    if (needGeneName) {
        ExpressionGRangesMetadata <- GenomicRanges::elementMetadata(
            SummarizedExperiment::rowRanges(MAE@ExperimentList$expression)
        )
        if (!("gene_name" %in% colnames(ExpressionGRangesMetadata))) {
            .stopNoCall(
                "A column named 'gene_name' with the names of genes in the ",
                "TENETMultiAssayExperiment object was not found. Because ",
                "geneAnnotationDataset is not selected, please ensure this ",
                "column is present so gene names can be extracted alongside ",
                "the gene IDs."
            )
        }
    }

    ## Ensure that the same number of expression and methylation samples are
    ## present; if not, stop the function and return an error message implying
    ## improper sample matching
    if (ncol(MAE@ExperimentList$expression@assays@data[[1]]) !=
        ncol(MAE@ExperimentList$methylation@assays@data[[1]])
    ) {
        .stopNoCall(
            "The number of samples in the expression and methylation ",
            "SummarizedExperiments differs. Please ensure that these datasets ",
            "contain the same number of samples and that the samples are ",
            "matched on a 1:1 basis."
        )
    }

    ## Ensure that the mapping object of the MultiAssayExperiment object
    ## contains a "sampleType" column, and that the column contains only the
    ## values "Case" and "Control"
    if (!"sampleType" %in% colnames(MAE@sampleMap)) {
        .stopNoCall(
            "The 'sampleType' column was not found in the mapping object of ",
            "the specified TENETMultiAssayExperiment object. Please ensure ",
            "this column is present in the mapping object and it notes ",
            "whether each sample is a 'Case' or 'Control' sample, using those ",
            "values only."
        )
    } else if (
        all(unique(sort(MAE@sampleMap$sampleType)) != c("Case", "Control"))
    ) {
        .stopNoCall(
            "The 'sampleType' column in the mapping object of the specified ",
            "TENETMultiAssayExperiment object contains values other than ",
            "'Case' and 'Control', or is missing one of those two values. ",
            "Please ensure that this column notes whether each sample is a ",
            "'Case' or 'Control' sample, using those values only."
        )
    } else if (
        any(duplicated(
            paste(
                MAE@sampleMap$primary,
                MAE@sampleMap$assay,
                MAE@sampleMap$sampleType,
                sep = "\b" ## Backspace; should never occur in input data
            )
        ))
    ) {
        .stopNoCall(
            "There are multiple samples from the same source, and ",
            "of the same data and sample type, which will affect TENET ",
            "analyses. Please ensure only a single unique sample from each ",
            "source, data type, and sample type have been provided."
        )
    } else if (
        !all(as.vector(table(
            paste(
                MAE@sampleMap$primary,
                MAE@sampleMap$sampleType,
                sep = "\b" ## Backspace; should never occur in input data
            )
        )) == 2)
    ) {
        .stopNoCall(
            "Data samples in the provided MultiAssayExperiment object do not ",
            "all appear to be paired. Please ensure that matched methylation ",
            "and expression objects have been provided for each source and ",
            "sample type."
        )
    }
}

## Internal function to ensure that data from a certain TENET step are present
## in the MultiAssayExperiment and return an error message if not
.ensureStepPresent <- function(
    MAE,
    stepName,
    substepName = NA,
    substepDescription = NA,
    substepParamDescription = NA) {
    if (!stepName %in% names(MAE@metadata)) {
        .stopNoCall(
            "Data output by the ", stepName, " function were not ",
            "found in the supplied TENETMultiAssayExperiment object. Please ",
            "ensure the data are present by running the ", stepName,
            " function and supplying the output MultiAssayExperiment object."
        )
    }

    if (!is.na(substepName)) {
        if (!substepName %in% names(MAE@metadata[[stepName]])) {
            .stopNoCall(
                "Data output by the ", stepName, " function for ",
                substepDescription, " are not found in the supplied ",
                "TENETMultiAssayExperiment object. Please ensure the data ",
                "are present by running the ", stepName, " function with ",
                substepParamDescription,
                "and supplying the output MultiAssayExperiment object."
            )
        }
    }
}

## Internal function to validate and load the requested gene annotation
## dataset, returning it as a variable. Can take a GRanges object directly, or
## a path to a GFF3 or GTF file.
.loadGeneAnnotationDataset <- function(
    geneAnnotationDataset,
    featureTypes = "gene") {
    ## If the argument is a character string, it must be a path to a GFF3 or
    ## GTF file
    if (length(geneAnnotationDataset) == 1 &&
        is.character(geneAnnotationDataset) &&
        grepl(
            "\\.(gtf|gff3)(\\.(gz|bz2|xz))?$",
            geneAnnotationDataset,
            ignore.case = TRUE
        )) {
        isGFF3 <- grepl(
            "\\.gff3(\\.(gz|bz2|xz))?$",
            geneAnnotationDataset,
            ignore.case = TRUE
        )

        ## Import it using rtracklayer (which handles GFF and GTF, including
        ## gzip, bzip2, and xz compression)
        geneAnnotationDataset <- rtracklayer::import(geneAnnotationDataset)
    } else if (!inherits(geneAnnotationDataset, "GRanges")) {
        ## Otherwise, it must be a GRanges object
        .stopNoCall(
            "The input for geneAnnotationDataset is incorrect. ",
            "Please specify a GRanges object (such as one imported via ",
            "rtracklayer::import) or a path to a GFF3 or GTF file."
        )
    } else {
        ## Only GFF3 files should have an ID column
        isGFF3 <- ("ID" %in%
            colnames(GenomicRanges::elementMetadata(geneAnnotationDataset)))
    }

    ## Save the current elementMetadata into a variable to avoid making the code
    ## too long
    elementMetadata <- GenomicRanges::elementMetadata(geneAnnotationDataset)

    ## Ensembl GFF3 files call the "gene_name" column "Name"
    nameColumn <- "gene_name"
    if (!(nameColumn %in% colnames(elementMetadata))) {
        nameColumn <- "Name"
        ## Temporarily disallow Ensembl GFF3 files until we fix the issue
        if (nameColumn %in% colnames(elementMetadata)) {
            .stopNoCall(
                "This looks like an Ensembl GFF3 file. Ensembl GFF3 files are ",
                "not currently supported due to an unresolved issue. We ",
                "apologize for the inconvenience."
            )
        }
    }

    ## Ensure all requested transcript types are present
    typesNotPresent <- featureTypes[
        !(featureTypes %in% elementMetadata$type)
    ]

    if (isGFF3) {
        if (length(typesNotPresent) > 0) {
            ## For GFF3 files, if the type column does not contain the requested
            ## types, synthesize a new type column from the ID column. Ensembl
            ## GFF3 files have GENCODE-style types there. GENCODE GFF3 files
            ## have a different set of types there, which we can't use. The ID
            ## column may contain multiple pieces of information separated by
            ## colons, in which case the type is the first field.
            elementMetadata$type <- sub(":.*", "", elementMetadata$ID)

            ## Copy the type column back into the GRanges object for later use
            GenomicRanges::elementMetadata(geneAnnotationDataset)$type <-
                elementMetadata$type

            ## Recheck for the requested types using the new type column
            typesNotPresent <- featureTypes[
                !(featureTypes %in% elementMetadata$type)
            ]
        }
    }

    if (length(typesNotPresent) > 0) {
        ## Ensembl releases < 75 have no "transcript", so they can't be
        ## supported easily. See
        ## https://groups.google.com/g/rsem-users/c/aKmnk6_LnMU
        .stopNoCall(
            "Failed to retrieve data of the following type(s) from the ",
            "gene annotation dataset: ",
            paste(typesNotPresent, collapse = ", "), ". ",
            "TENET has only been tested with GENCODE and Ensembl gene ",
            "annotation datasets. If you are using another dataset, please ",
            "ensure that it uses the values \"gene\" and \"transcript\" for ",
            "feature types, which must be stored in a column named \"type\". ",
            "In GFF3 files, feature types may alternatively be stored in the ",
            "first colon-separated field of the \"ID\" column, the second ",
            "field of which must be the ID itself. Types stored there will ",
            "only be used if the \"type\" column does not contain the ",
            "required types. Gene names must be stored in a column named ",
            "\"gene_name\" or \"Name\". GTF files must contain a \"gene_id\" ",
            "column, and GFF3 files must contain an \"ID\" column. An ",
            "annotation dataset specified as a GRanges object will be assumed ",
            "to be derived from a GFF3 file if it contains an \"ID\" column, ",
            "and from a GTF file otherwise. Ensembl GTF files older than ",
            "release 75 are not supported."
        )
    }

    ## For GFF3 files, we need to synthesize a new gene_id column from the ID
    ## column rather than using the ID column directly, because in GENCODE
    ## GFF3 files, the ID column has the same ID for the two copies of a
    ## gene in the X and Y chromosome's pseudoautosomal regions, while the ID
    ## column differentiates them using a _PAR_Y suffix. The ID column may
    ## contain multiple pieces of information separated by colons, in which case
    ## the ID itself is the second field.
    if (isGFF3) {
        GenomicRanges::elementMetadata(geneAnnotationDataset)$gene_id <- sub(
            ":.*", "", sub("[^:]*:", "", elementMetadata$ID)
        )
    }

    ## Note that there are duplicate gene IDs, since transcripts for the
    ## same gene are given the same ID. We want to prioritize the entries
    ## for the gene itself over its transcripts and other entries. Create a
    ## temporary column which will be used for this prioritization. It must
    ## be set directly; modifying the elementMetadata variable would just
    ## modify the previous copy of the elementMetadata.
    GenomicRanges::elementMetadata(geneAnnotationDataset)$isGene <-
        (elementMetadata$type == "gene")

    ## Filter the dataset to the selected transcript types before sorting it,
    ## because the Boolean vector of the criteria will no longer be valid after
    ## sorting
    geneAnnotationDataset <- geneAnnotationDataset[
        elementMetadata$type %in% featureTypes
    ]

    ## Sort the dataset by the variable added above in descending order, which
    ## will put the genes first. It must be accessed directly; accessing the
    ## elementMetadata variable would just access the previous copy of the
    ## elementMetadata.
    geneAnnotationDataset <- geneAnnotationDataset[order(
        GenomicRanges::elementMetadata(geneAnnotationDataset)$isGene,
        decreasing = TRUE
    )]

    ## Create an info data frame with gene IDs, gene names, and coordinates
    ## from the gene annotation dataset. Add a TSS location for each transcript
    ## depending on its strand and start/end positions. Remove GENCODE mapping
    ## versions and ENSG ID versions, but keep anything else extra like _PAR_Y.
    ## Convert the factors into character strings (stringsAsFactors does not
    ## help here, since they are already factors).
    geneIDdf <- data.frame(
        "geneID" = sub("_[0-9]+$", "", sub(
            "\\.[0-9]+", "",
            GenomicRanges::elementMetadata(geneAnnotationDataset)$gene_id
        )),
        "geneName" = GenomicRanges::elementMetadata(geneAnnotationDataset)[[
            nameColumn
        ]],
        "type" = as.character(
            GenomicRanges::elementMetadata(geneAnnotationDataset)$type
        ),
        "chromosome" = as.character(GenomicRanges::seqnames(
            geneAnnotationDataset
        )),
        "start" = GenomicRanges::start(geneAnnotationDataset),
        "end" = GenomicRanges::end(geneAnnotationDataset),
        "TSS" = as.numeric(
            ifelse(
                as.character(
                    GenomicRanges::strand(geneAnnotationDataset)
                ) == "-",
                GenomicRanges::end(geneAnnotationDataset),
                ## Assume that genes with "*" (no strand info) have a positive
                ## strand
                GenomicRanges::start(geneAnnotationDataset)
            )
        ),
        stringsAsFactors = FALSE
    )

    ## Add unique gene IDs as the row names. This will ensure the gene
    ## IDs for the gene listings themselves stay the same, but the listings
    ## for all other transcripts will change.
    rownames(geneIDdf) <- make.unique(geneIDdf$geneID)

    ## Return the final data frame
    return(geneIDdf)
}

## Internal function to get gene IDs and names from the MAE rowRanges, or the
## gene annotation dataset if provided
.getGeneIDsAndNames <- function(MAE, geneAnnotationDataset) {
    if (!is.na(geneAnnotationDataset)) {
        return(.loadGeneAnnotationDataset(geneAnnotationDataset))
    } else {
        ## Extract the metadata from the expression GRanges object in the MAE,
        ## and create a data frame which pairs gene IDs to gene names
        return(.loadGeneAnnotationDataset(SummarizedExperiment::rowRanges(
            MAE@ExperimentList$expression
        )))
    }
}

## Internal function to get the expression or methylation samples of the given
## type from the MAE, or optionally only their names. If sampleType is NA (the
## default), get all expression or methylation samples. This makes use of the
## mapping object, which must include the sampleType column.
.getExpOrMetSamplesOfType <- function(
    MAE,
    expOrMet,
    sampleType = NA,
    namesOnly = FALSE) {
    ## Get the sample mapping, and sort it by sample type then the
    ## primary name. This is done to ensure that when sample names are grabbed
    ## for various expression or methylation sample subsets, they come out in
    ## the same order and can be assumed to be paired already.
    MAEMap <- as.data.frame(MAE@sampleMap)
    MAEMap <- MAEMap[with(MAEMap, order(sampleType, primary)), ]

    ## Get the names of the samples in the expression and methylation data
    if (is.na(sampleType)) {
        sampleNamesSubset <- MAEMap[MAEMap$assay == expOrMet, "colname"]
    } else {
        sampleNamesSubset <- MAEMap[
            (MAEMap$assay == expOrMet & MAEMap$sampleType == sampleType),
            "colname"
        ]
    }

    if (namesOnly) {
        return(sampleNamesSubset)
    }

    ## Otherwise, return the expression or methylation dataset for the given
    ## sample type
    dataSubset <- MultiAssayExperiment::assays(
        MAE@ExperimentList[[expOrMet]]
    )[[1]][
        ,
        sampleNamesSubset
    ]

    return(dataSubset)
}

## Internal function to ensure that all specified RE DNA methylation sites are
## present in the methylation data of interest. If any RE DNA methylation sites
## are not found, it will warn the user and exclude them.
.excludeMissingMethylationSites <- function(methSiteList, methylationData) {
    missingMethSites <- methSiteList[
        which(!(methSiteList %in% rownames(methylationData)))
    ]
    if (length(missingMethSites) > 0) {
        .warningNoCall(
            "The following RE DNA methylation sites were not found in the ",
            "methylation data and will be ignored: ",
            paste(missingMethSites, collapse = ", ")
        )
        methSiteList <- methSiteList[
            which(methSiteList %in% rownames(methylationData))
        ]
    }

    return(methSiteList)
}

## Internal function to validate and load the requested DNA methylation array's
## probe manifest, returning it as a variable
.loadMethylationManifest <- function(DNAMethylationArray) {
    ## Load the manifest via sesameData (it will validate the argument)
    DNAMethylationArrayManifest <- suppressMessages(
        sesameData::sesameData_getManifestGRanges(
            DNAMethylationArray,
            genome = "hg38"
        )
    )

    if (is.null(DNAMethylationArrayManifest)) {
        .stopNoCall(
            "The probe manifest for the methylation array \"",
            DNAMethylationArray,
            "\" could not be loaded via the sesameData package. Please ",
            "ensure that the DNAMethylationArray argument refers to a ",
            "methylation array supported by sesameData (see ",
            "?sesameData::sesameData_getManifestGRanges). If you are in ",
            "a computing cluster environment where compute nodes do not have ",
            "internet access, you must run TENETCacheAllData() once while ",
            "connected to the internet before using the DNAMethylationArray ",
            "argument."
        )
    }

    ## Remove probes with missing chromosomes
    DNAMethylationArrayManifest <- DNAMethylationArrayManifest[
        GenomicRanges::seqnames(DNAMethylationArrayManifest) != "*"
    ]

    ## Remove probes with no coordinates
    DNAMethylationArrayManifest <- DNAMethylationArrayManifest[!(
        GenomicRanges::end(DNAMethylationArrayManifest) == 0 &
            GenomicRanges::start(DNAMethylationArrayManifest) == 0
    )]

    ## Create an info data frame with methylation site IDs, chromosomes, and
    ## coordinates from the RE DNA methylation site manifest
    methSiteIDdf <- data.frame(
        "DNAMethylationSiteID" = names(DNAMethylationArrayManifest),
        "chromosome" = as.character(
            GenomicRanges::seqnames(DNAMethylationArrayManifest)
        ),
        "start" = GenomicRanges::start(DNAMethylationArrayManifest),
        "end" = GenomicRanges::end(DNAMethylationArrayManifest),
        stringsAsFactors = FALSE
    )

    ## Sort the dataset by increasing methylation site ID
    methSiteIDdf <- methSiteIDdf[order(methSiteIDdf$DNAMethylationSiteID), ]

    ## Remove any extra "_TC21"/"_BC21"/etc. from the end of the IDs (found in
    ## EPIC v2). These are mostly replicates of the same probe at the same
    ## location (see the Infinium MethylationEPIC v2.0 Manifest File Release
    ## Notes). Some have different locations, but there are so few (45) that it
    ## isn't really an issue.
    methSiteIDdf$DNAMethylationSiteID <- sub(
        "_[TB][CO][12][0-9]+$", "", methSiteIDdf$DNAMethylationSiteID
    )

    ## Remove any duplicates
    methSiteIDdf <- methSiteIDdf[
        !duplicated(methSiteIDdf$DNAMethylationSiteID),
    ]

    ## Set the rownames to be the methylation site IDs for subsetting
    rownames(methSiteIDdf) <- make.unique(methSiteIDdf$DNAMethylationSiteID)

    return(methSiteIDdf)
}

## Internal function to get methylation site IDs and locations from the MAE, or
## the methylation array's probe manifest if provided
.getMethSiteIDsAndLocations <- function(MAE, DNAMethylationArray) {
    if (!is.na(DNAMethylationArray)) {
        return(.loadMethylationManifest(DNAMethylationArray))
    } else {
        ## Extract the metadata from the methylation GRanges object
        ## in the MAE, and create a data frame with methylation site IDs and
        ## locations with the methylation site IDs as the row names
        returnDF <- data.frame(
            "DNAMethylationSiteID" = names(SummarizedExperiment::rowRanges(
                MAE@ExperimentList$methylation
            )),
            "chromosome" = as.character(
                GenomicRanges::seqnames(SummarizedExperiment::rowRanges(
                    MAE@ExperimentList$methylation
                ))
            ),
            "start" = GenomicRanges::start(SummarizedExperiment::rowRanges(
                MAE@ExperimentList$methylation
            )),
            "end" = GenomicRanges::end(SummarizedExperiment::rowRanges(
                MAE@ExperimentList$methylation
            )),
            stringsAsFactors = FALSE
        )
        rownames(returnDF) <- returnDF$DNAMethylationSiteID
        return(returnDF)
    }
}

## Internal clustering functions for step 7 heatmaps
.step7HeatmapDistf <- function(d) {
    stats::dist(d, method = "euclidean")
}

.step7HeatmapClustf <- function(e) {
    stats::hclust(e, method = "ward.D2")
}

## Internal function to get the IDs and link counts of all/top genes/TFs
## in the given quadrant. If there are fewer genes/TFs than the topGeneNumber
## specified, warn the user and return all the genes/TFs available. If
## topGeneNumber is NA, get all genes/TFs.
.getQuadrantTopGenesOrTFs <- function(
    TENETMultiAssayExperiment,
    geneOrTF,
    hyperHypo,
    topGeneNumber) {
    quadrantResultsName <- paste0(hyperHypo, "methGplusResults")

    ## Ensure the quadrant's results are present in step 6
    .ensureStepPresent(
        TENETMultiAssayExperiment,
        stepName = "step6DNAMethylationSitesPerGeneTabulation",
        substepName = quadrantResultsName
    )

    ## Load the quadrant's link counts from step 6
    quadrantAllGeneOrTFFreq <- TENETMultiAssayExperiment@metadata$
        step6DNAMethylationSitesPerGeneTabulation[[quadrantResultsName]]

    ## Subset the genes to TFs only, if requested
    if (geneOrTF == "TF") {
        ## Create a temporary environment to store data from the TENET package
        dataEnvironment <- new.env(parent = emptyenv())

        ## Get the IDs of the human transcription factors
        utils::data(
            "humanTranscriptionFactorList",
            package = "TENET",
            envir = dataEnvironment
        )
        humanTFs <- dataEnvironment$humanTranscriptionFactorList
        genesAreTFs <- quadrantAllGeneOrTFFreq$geneID %in% humanTFs

        quadrantAllGeneOrTFFreq <- quadrantAllGeneOrTFFreq[
            genesAreTFs,
        ]

        ## Check whether all genes are TFs. If so, return NA to avoid creating
        ## two identical outputs when the user used TFOnly.
        if (all(genesAreTFs)) {
            message(
                "All genes with ", hyperHypo, "methylated G+ links are TFs, ",
                "so the separate output for TFs will be skipped."
            )
            ## Return a data frame so that the subset code works unmodified
            return(data.frame("geneID" = NA))
        }
    }

    ## If topGeneNumber is NA, get all genes/TFs
    if (is.na(topGeneNumber)) {
        return(quadrantAllGeneOrTFFreq)
    }

    ## If there are fewer top genes/TFs than the topGeneNumber
    ## specified by the user, warn the user of this
    if (nrow(quadrantAllGeneOrTFFreq) < topGeneNumber) {
        .warningNoCall(
            "Fewer ", geneOrTF, "s with any significant ", hyperHypo,
            "methylated G+ RE DNA methylation sites linked to them than the ",
            "specified 'topGeneNumber' were identified, so only those genes ",
            "will be included in this analysis."
        )
    }

    ## Return the names of the top genes/TFs. If there are fewer genes than the
    ## topGeneNumber specified by the user, return all the genes available.
    return(quadrantAllGeneOrTFFreq[
        seq_len(min(nrow(quadrantAllGeneOrTFFreq), topGeneNumber)),
    ])
}

## Internal function to check if an object is a single NA (since R 4.3 changed
## the behavior of if(!is.na(...))). Based on
## https://github.com/GeoBosh/gbutils/blob/16a311ce2e40813911461c6084393c1915568c5b/R/isNA.R
## (license: GPL v2 or later)
.isSingleNA <- function(x) {
    return(is.atomic(x) && length(x) == 1 && is.na(x))
}

## Internal function to emulate ifelse without the undesirable behavior of
## iterating on vector inputs. "if" is internally a 3-argument function and
## therefore can be used for this purpose.
.ifelseNoIterate <- `if`

## Internal function to create an invisible, recordable plot of a specified size
.newInvisibleRecordablePlot <- function(width = 10, height = 7) {
    ## If file is not set to NULL, Rplots.pdf is created, which we don't want
    grDevices::pdf(file = NULL, width = width, height = height)
    ## Enable recording
    grDevices::dev.control(displaylist = "enable")
    ## Set background color to white to fix dark mode issue
    graphics::par(bg = "white")
}

## Internal function to import and verify various forms of clinical data
.importClinicalData <- function(
    userInput,
    argumentName,
    clinicalDataColumn = NA, ## Column(s) in colData (if userInput is NA)
    returnType, ## "single" for 1 column, or "multiple" for multiple columns
    TENETMultiAssayExperiment) {
    if (.isSingleNA(userInput)) {
        ## Since userInput is NA, look in the colData of the
        ## TENETMultiAssayExperiment object
        if (returnType == "single") {
            ## Check that the specified column is present in the colData of the
            ## TENETMultiAssayExperiment object
            if (!(clinicalDataColumn %in% colnames(
                MultiAssayExperiment::colData(TENETMultiAssayExperiment)
            ))) {
                .stopNoCall(
                    "Properly formatted data were not found in the ",
                    "colData of the TENETMultiAssayExperiment ",
                    "object. Please ensure a '", clinicalDataColumn,
                    "' column is included in the colData of the ",
                    "TENETMultiAssayExperiment object."
                )
            }

            ## This is to harmonize with functionality when returnType
            ## is "multiple" below
            columnsPresent <- clinicalDataColumn
        } else if (returnType == "multiple") {
            ## Identify which of the specified clinical data columns are
            ## present in the colData of the TENETMultiAssayExperiment object
            columnsNotPresent <- clinicalDataColumn[
                !(clinicalDataColumn %in% colnames(
                    MultiAssayExperiment::colData(TENETMultiAssayExperiment)
                ))
            ]

            columnsPresent <- setdiff(
                clinicalDataColumn, columnsNotPresent
            )

            ## If some columns aren't present, return a warning to the user
            ## alerting them to this. If none of the columns are present, return
            ## an error instead.
            if (length(columnsNotPresent) == length(clinicalDataColumn)) {
                .stopNoCall(
                    "The ", argumentName,
                    " argument is NA, but the corresponding variables in the ",
                    "colData of the TENETMultiAssayExperiment object ",
                    "were not found (",
                    paste(columnsNotPresent, collapse = ", "),
                    "). Please ensure at least some of these columns are ",
                    "present in the colData of the TENETMultiAssayExperiment ",
                    "object and run this function again. Any missing ",
                    "columns will result in missing results from this function."
                )
            } else if (length(columnsNotPresent) > 0) {
                .warningNoCall(
                    "The ", argumentName, " argument is NA, but some of the ",
                    "corresponding clinical variables in the colData of the ",
                    "TENETMultiAssayExperiment object were not found (",
                    paste(columnsNotPresent, collapse = ", "),
                    "). Results will still be generated using variables which ",
                    "are present, but the missing columns will result in ",
                    "missing results from this function."
                )
            }
        }

        ## Get the columns that were present from the colData of the
        ## TENETMultiAssayExperiment object
        import <- as.data.frame(
            TENETMultiAssayExperiment@colData[columnsPresent]
        )
    } else if (!is.null(ncol(userInput))) {
        ## This means the userInput argument is a matrix or data frame

        ## Check that at least a subset of the userInput object's rownames
        ## are present in the names of the colData of the
        ## TENETMultiAssayExperiment object
        if (!any(rownames(userInput) %in%
            rownames(TENETMultiAssayExperiment@colData))) {
            .stopNoCall(
                "None of the rownames of the object given as the ",
                argumentName,
                " argument match the IDs of samples within the ",
                "TENETMultiAssayExperiment object. ",
                "Sample IDs must be included as the rownames of the ",
                "object, and they must match with sample IDs given as the ",
                "rownames of the colData of the ",
                "TENETMultiAssayExperiment object."
            )
        }

        ## If there are multiple columns required for output, we will need to
        ## check that the needed columns are present in the user data
        if (returnType == "multiple") {
            ## Identify which of the specified clinical data columns are
            ## present in the specified userInput object
            columnsNotPresent <- clinicalDataColumn[
                !(clinicalDataColumn %in% colnames(userInput))
            ]

            columnsPresent <- setdiff(
                clinicalDataColumn, columnsNotPresent
            )

            ## If there are some columns that aren't present, return a warning
            ## to the user alerting them as such.
            if (length(columnsNotPresent) > 0) {
                .warningNoCall(
                    "The following columns containing necessary information ",
                    "are missing from the object given as the ", argumentName,
                    " argument: ",
                    paste(columnsNotPresent, collapse = ", "),
                    ". Results will still be generated using variables which ",
                    "are present, but the missing columns will result in ",
                    "missing results from this function."
                )
            }

            ## Organize the userInput object so the samples provided in it
            ## are in the same order as the data for the samples in the
            ## TENETMultiAssayExperiment colData.
            ## This will induce NAs for the data provided in the
            ## userInput object if a given sample in the colData isn't
            ## present in the userInput object.
            import <- as.data.frame(userInput[
                rownames(TENETMultiAssayExperiment@colData),
                columnsPresent
            ])
        } else {
            ## Organize the userInput object so the samples provided in it are
            ## in the same order as the data for the samples in the
            ## TENETMultiAssayExperiment colData.
            ## This will induce NAs for the data provided in the userInput
            ## object if a given sample in the colData isn't present in the
            ## userInput object.

            ## Note: This code looks a little awkward, but it is done in order
            ## to ensure a data frame with a single column isn't converted into
            ## a vector yet.
            import <- as.data.frame(userInput[
                rownames(TENETMultiAssayExperiment@colData),
            ])

            colnames(import) <- colnames(userInput)

            ## Additionally, if dealing with a single column, we need to adjust
            ## the rownames
            if (ncol(import) == 1) {
                rownames(import) <- rownames(TENETMultiAssayExperiment@colData)
            }
        }
    } else {
        ## If the input is a character string, assume it is a path to a file to
        ## load
        if (length(userInput) == 1 && is.character(userInput)) {
            ## Read in the file. Assume patient IDs have been provided in the
            ## first column, and there are headers in the first row
            userInput <- utils::read.delim(
                userInput,
                header = TRUE,
                row.names = 1,
                stringsAsFactors = FALSE
            )

            ## Check that at least a subset of the userInput file's rownames
            ## are present in the names of the colData of the
            ## TENETMultiAssayExperiment object
            if (!any(rownames(userInput) %in%
                rownames(TENETMultiAssayExperiment@colData))) {
                .stopNoCall(
                    "None of the rownames of the file loaded from the ",
                    "path given as the ", argumentName,
                    " argument match the IDs of samples within the ",
                    "TENETMultiAssayExperiment object. ",
                    "Sample IDs must be included as the first column of the ",
                    "file, and they must match with sample IDs given as the ",
                    "rownames of the colData of the ",
                    "TENETMultiAssayExperiment object."
                )
            }

            ## If there are multiple columns required for output, we will need
            ## to check the needed columns are present in the user data
            if (returnType == "multiple") {
                ## Identify which of the specified clinical data columns are
                ## present in the specified userInput object
                columnsNotPresent <- clinicalDataColumn[
                    !(clinicalDataColumn %in% colnames(userInput))
                ]

                columnsPresent <- setdiff(
                    clinicalDataColumn, columnsNotPresent
                )

                ## If there are some columns that aren't present, return a
                ## warning to the user alerting them as such.
                if (length(columnsNotPresent) > 0) {
                    .warningNoCall(
                        "The following columns containing necessary ",
                        "information are missing ",
                        "from the object loaded from the path given as the ",
                        argumentName, " argument: ",
                        paste(columnsNotPresent, collapse = ", "),
                        ". Results will still be generated using variables ",
                        "which are present, but the missing columns will ",
                        "result in missing results from this function."
                    )
                }

                ## Organize the userInput object so the samples provided in it
                ## are in the same order as the data for the samples in the
                ## TENETMultiAssayExperiment colData.
                ## This will induce NAs for the data provided in the
                ## userInput object if a given sample in the colData isn't
                ## present in the userInput object.
                import <- as.data.frame(userInput[
                    rownames(TENETMultiAssayExperiment@colData),
                    columnsPresent
                ])
            } else {
                ## Organize the loaded data so the samples provided in it are in
                ## the same order as the data for the samples in the
                ## TENETMultiAssayExperiment colData.
                ## This will induce NAs for the data provided in the
                ## userInput object if a given sample in the colData isn't
                ## present in the userInput object.

                ## Note: This code looks a little awkward, but it ensures that
                ## a data frame with a single column isn't converted into a
                ## vector.
                import <- as.data.frame(userInput[
                    rownames(TENETMultiAssayExperiment@colData),
                ])

                colnames(import) <- colnames(userInput)
            }
        } else {
            ## Otherwise, assume the input is a vector

            ## Check to see if there are names given in the vector. If there
            ## are, assume they are the sample IDs which are also present in the
            ## rownames of the TENETMultiAssayExperiment colData
            if (!is.null(names(userInput))) {
                ## Check that at least a subset of the userInput object's names
                ## are present in the names of the colData of the
                ## TENETMultiAssayExperiment object
                if (!any(names(userInput) %in%
                    rownames(TENETMultiAssayExperiment@colData))) {
                    .stopNoCall(
                        "None of the names of the values given as the ",
                        argumentName, " argument ",
                        "match the IDs of samples within the ",
                        "TENETMultiAssayExperiment object. The names of the ",
                        "given values must match with sample IDs given as ",
                        "the rownames of the colData of the ",
                        "TENETMultiAssayExperiment object."
                    )
                }

                ## Organize the userInput object so the samples provided in it
                ## are in the same order as the data for the samples in the
                ## TENETMultiAssayExperiment colData.
                ## This will induce NAs if a given sample in the colData isn't
                ## present in the names of the userInput
                userInput <- userInput[
                    rownames(TENETMultiAssayExperiment@colData)
                ]
            } else {
                ## If there are no names, check to ensure the same
                ## number of values were given as there are samples in the
                ## TENETMultiAssayExperiment colData. If this is the case, we
                ## will assume data are given in the same order as the samples
                ## in the TENETMultiAssayExperiment colData. Otherwise return
                ## an error
                if (length(userInput) !=
                    nrow(TENETMultiAssayExperiment@colData)) {
                    .stopNoCall(
                        "There are no sample IDs given as names for the ",
                        "values given as the ", argumentName,
                        " argument, and the number of values given as the ",
                        argumentName,
                        " does not match the number of samples present in the ",
                        "colData of the TENETMultiAssayExperiment object, ",
                        "making accurate data matching difficult. Please ",
                        "check the dataset supplied as the ",
                        argumentName, " argument. The names of the ",
                        "given values must match with sample IDs given as ",
                        "the rownames of the colData of the ",
                        "TENETMultiAssayExperiment object, or there must be ",
                        "the same number of values as there are samples in ",
                        "the TENETMultiAssayExperiment colData and they must ",
                        "be arranged in the same order."
                    )
                } else {
                    ## Issue a warning to the user, but assume the samples are
                    ## given in the correct order and proceed
                    .warningNoCall(
                        "There are no sample IDs given as names for the ",
                        "values given as the ", argumentName,
                        " argument. They will be assumed to derive from ",
                        "samples in the same order the samples are listed in ",
                        "the colData of the TENETMultiAssayExperiment ",
                        "object. If this is not the case, please check the ",
                        "data provided and rerun this function."
                    )
                }
            }

            ## Assemble a single column data frame with the data from the
            ## userInput, with the sample IDs in the colData of the
            ## TENETMultiAssayExperiment
            import <- data.frame(
                "placeholder" = userInput,
                stringsAsFactors = FALSE
            )
            rownames(import) <- rownames(TENETMultiAssayExperiment@colData)
            colnames(import) <- argumentName
        }
    }

    if (returnType == "single") {
        ## If there are multiple columns in the imported data, return a
        ## warning to the user that only the first column will be output
        if (ncol(import) > 1) {
            .warningNoCall(
                "Multiple columns of data have been provided as the ",
                argumentName,
                " argument, but only one column of data is expected. Only ",
                "the first column of data in the object will be used."
            )
        }
    } else if (returnType == "multiple") {
        ## If only a single column is present in the imported data, but multiple
        ## columns are expected, return an error and stop
        if (ncol(import) == 1) {
            .stopNoCall(
                "Only a single column of data has been provided as the ",
                argumentName,
                " argument, but multiple columns of data are expected. ",
                "Please check the data provided and rerun this function."
            )
        }
    }

    return(import)
}

## Internal list of all of our AnnotationHub datasets and their names
.TENETAnnotationHubIDs <- c(
    "Public enhancer regions" = "AH116721",
    "Public open chromatin regions" = "AH116722",
    "Public promoter regions" = "AH116723",
    "Consensus enhancer regions" = "AH116724",
    "Consensus open chromatin regions" = "AH116725",
    "Consensus promoter regions" = "AH116726",
    "ENCODE dELS regions" = "AH116727",
    "ENCODE pELS regions" = "AH116728",
    "ENCODE PLS regions" = "AH116729"
)

## Internal list of all our ExperimentHub datasets and their names
.TENETExperimentHubIDs <- c(
    "Example TENET MultiAssayExperiment" = "EH9587",
    "Example TENET clinical dataframe" = "EH9588",
    "Example TENET step1MakeExternalDatasets GRanges" = "EH9589",
    "Example TENET step2GetDifferentiallyMethylatedSites purity SummarizedExperiment" = "EH9590",
    "Example TENET peak regions GRanges" = "EH9591",
    "Example TENET TAD regions GRanges" = "EH9592"
)

## Internal function to attempt to load a dataset from AnnotationHub and stop
## with a user-friendly message if it fails
.loadFromAnnotationHub <- function(ah, id) {
    name <- names(.TENETAnnotationHubIDs[.TENETAnnotationHubIDs == id])
    message("Loading ", name, " (", id, ") from AnnotationHub")
    tryCatch(
        {
            return(ah[[id, verbose = FALSE]])
        },
        error = function(cond) {
            .stopNoCall(
                "Failed to load TENET.AnnotationHub data. Please ",
                "ensure you have internet access and the latest versions of ",
                "R, TENET, Bioconductor, and AnnotationHub. If you are in a ",
                "computing cluster environment where compute nodes do not ",
                "have internet access, you must run TENETCacheAllData() ",
                "once while connected to the internet before using ",
                "TENET.AnnotationHub datasets."
            )
        }
    )
}

## Internal function to display a warning without showing the function call
.warningNoCall <- function(...) {
    warning(..., call. = FALSE)
}

## Internal function to stop execution without showing the function call
.stopNoCall <- function(...) {
    stop(..., call. = FALSE)
}

## Internal function to display an error like stop() would, but without stopping
## execution or showing the function call
.displayError <- function(...) {
    try(stop(..., call. = FALSE))
}

## Internal function to create the metToExpSampleConversion vector
## which maps methylation sample names to their matched expression sample
## names, assuming the methylation and expression values share a clinical
## data match within the mapping object of the MAE.
.createMetToExpSampleConversionVector <- function(
    MAE) {
    ## Get the methylation values that match with expression values
    ## using the mapping data.
    ## This assumes the methylation and expression values share a clinical
    ## data match within the mapping.

    ## Split the mapping data based on the expression/methylation samples
    ## as well as the sampleType - note: the order will likely be the same in
    ## each but in case it isn't, we can use the similarity in primary names
    ## to identify matching expression/methylation sample pairs. The split
    ## between Case and Control is done here as there is ambiguity in mapping
    ## samples with both Case and Control included, as there may be a Case and
    ## Control sample which both map to the same primary name.
    expressionMappingControl <- MAE@sampleMap[
        (MAE@sampleMap$assay == "expression" &
            MAE@sampleMap$sampleType == "Control"),
    ]

    expressionMappingCase <- MAE@sampleMap[
        (MAE@sampleMap$assay == "expression" &
            MAE@sampleMap$sampleType == "Case"),
    ]

    methylationMappingControl <- MAE@sampleMap[
        (MAE@sampleMap$assay == "methylation" &
            MAE@sampleMap$sampleType == "Control"),
    ]

    methylationMappingCase <- MAE@sampleMap[
        (MAE@sampleMap$assay == "methylation" &
            MAE@sampleMap$sampleType == "Case"),
    ]

    ## Set up initial conversion vectors with the expression names
    ## this will later allow us to convert methylation sample names to their
    ## corresponding expression ones
    metToExpSampleConversionControl <- expressionMappingControl$colname
    metToExpSampleConversionCase <- expressionMappingCase$colname

    ## Get the names of the clinical samples that go with the expression samples
    ## These should also match with the ones for the corresponding methylation
    ## sample
    expressionMappedClinicalMatchControl <- expressionMappingControl$primary
    expressionMappedClinicalMatchCase <- expressionMappingCase$primary

    ## Find the names of the methylation samples in the same order as the
    ## clinical values (and likewise, the expression values)
    clinicalNumericalPositionForMethylationControl <- match(
        expressionMappedClinicalMatchControl,
        methylationMappingControl$primary
    )

    clinicalNumericalPositionForMethylationCase <- match(
        expressionMappedClinicalMatchCase,
        methylationMappingCase$primary
    )

    ## Set the names of the initial conversion vectors to be
    ## the methylation samples corresponding to the expression ones
    names(metToExpSampleConversionControl) <- methylationMappingControl[
        clinicalNumericalPositionForMethylationControl,
        "colname"
    ]

    names(metToExpSampleConversionCase) <- methylationMappingCase[
        clinicalNumericalPositionForMethylationCase,
        "colname"
    ]

    ## Combine the two conversion vectors
    metToExpSampleConversion <- c(
        metToExpSampleConversionControl, metToExpSampleConversionCase
    )

    ## Return the conversion vector
    return(metToExpSampleConversion)
}
