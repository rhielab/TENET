## Code to prepare the `humanTranscriptionFactorList` object

## Successfully accessed 2024-09-25

## Annotations come from:
## Lambert SA, Jolma A, Campitelli LF, et al. The
## Human Transcription Factors. Cell. 2018 Feb 8;172(4):650-665. doi:
## 10.1016/j.cell.2018.01.029. Erratum in: Cell. 2018 Oct 4;175(2):598-599.
## PMID: 29425488.
## http://humantfs2.ccbr.utoronto.ca/download.php
## File listed under: Current lists of human TFs and their motifs (v1.01) -
## 1. Human TFs - TF Ensembl IDs (.txt)

## Define the URL for the TF list file
tfFileURL <- "http://humantfs2.ccbr.utoronto.ca/download/v_1.01/TFs_Ensembl_v_1.01.txt"

## Read the lines of the file into a vector
humanTranscriptionFactorList <- readLines(tfFileURL)

## Check if the data subdirectory exists, and if it doesn't, create it
if (!dir.exists("data")) {
    dir.create("data")
}

save(
    humanTranscriptionFactorList,
    file = "data/humanTranscriptionFactorList.rda",
    compress = "xz",
    compression_level = -9 ## Level -n = level n but with xz -e
)
