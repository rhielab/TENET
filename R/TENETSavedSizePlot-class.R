## This class enables us to store the plot width and height when recording
## plots. Doing so prevents errors or incorrect plot appearance when the user's
## plot window is the wrong size when a recorded plot is later recalled.

## Based on https://www.datamentor.io/r-programming/s4-class

## Define the TENETSavedSizePlot class.
## The methods::new import exists to silence an R CMD check warning.
#' @importFrom methods new
TENETSavedSizePlot <- methods::setClass(
    "TENETSavedSizePlot",
    slots = list(plot = "ANY", width = "numeric", height = "numeric")
)

## Define a function that will be called when show or print are passed a
## TENETSavedSizePlot object
methods::setMethod(
    "show", "TENETSavedSizePlot",
    function(object) {
        ## Create a new plot with the specified width and height.
        ## Only do this if the user is in an interactive session, to avoid
        ## interfering with vignette builds.
        ## noRStudioGD: "Do not use the RStudio graphics device even if
        ## specified as the default device: it does not accept arguments such as
        ## width and height."
        if (interactive()) {
            grDevices::dev.new(
                width = object@width,
                height = object@height,
                unit = "in",
                noRStudioGD = TRUE
            )
        }

        ## Get the inner plot object
        tempPlot <- object@plot

        ## Set its R version to the current R version to avoid a warning when
        ## they don't match
        attr(tempPlot, "Rversion") <- getRversion()

        ## Show the inner plot object
        methods::show(tempPlot)
    }
)

## Define a function that will be called when a TENETSavedSizePlot object is
## created. This allows us to save the current plot width and height
## automatically.
## Based on https://stackoverflow.com/questions/61047760/what-is-distinction-between-initialize-method-vs-prototype-in-setclass-to-set-de
setMethod(
    "initialize", "TENETSavedSizePlot",
    function(.Object, plot, width = NA, height = NA) {
        if (.isSingleNA(width) || .isSingleNA(height)) {
            ds <- grDevices::dev.size(units = "in")
            if (.isSingleNA(width)) .Object@width <- ds[[1]]
            if (.isSingleNA(height)) .Object@height <- ds[[2]]
        }
        .Object@plot <- plot
        methods::validObject(.Object)
        return(.Object)
    }
)

## Define a utility function to record a plot and wrap it in a
## TENETSavedSizePlot
.recordTENETSavedSizePlot <- function() {
    return(TENETSavedSizePlot(grDevices::recordPlot()))
}
