## =============================================================================
## Classes

setClassUnion("MaybeEnvironment", c('environment', 'NULL'))
setClass("Cached",
         representation=representation("VIRTUAL",
           cache="MaybeEnvironment",
           range.strict="logical"),
         prototype=prototype(cache=NULL, range.strict=FALSE))

## -----------------------------------------------------------------------------
## Rig up an AlignedRead track of sorts
##
## TODO: have setCoverage honor gaps in alignments
setClass("BamTrack", contains=c("AlignedReadTrack", "Cached"),
         representation=representation(
           bam.file="BamFile",
           params="ScanBamParam"),
         prototype=prototype(
           bam.file=BamFile(""),
           params=ScanBamParam()))

setClass("BigWigTrack", contains=c("DataTrack", "Cached"),
         representation=representation(bw.file="BigWigFile"),
         prototype=prototype(bw.file=BigWigFile("")))

##' Works as a GeneRegionTrack but stores the "unwound/flattend" TranscriptDb
##' GRanges object for quick reuse via subsetting inside the
##' \code{@cache[['annotation]]} environment variable.
##'
##' The \code{annotation} object is hardocded to mimic the structure of the
##' current `range` objects being built via a all to .buildRange,TranscriptDb.
##'
##' The GenomicRanges object in `range` is used as the "active" region of the
##' cache If @range.strict == TRUE, these range fenceposts are used when
##' enforcing plotting coordinates when plotTracks is called.
setClass("TxDbTrack", contains=c("GeneRegionTrack", "Cached"))

## =============================================================================
## Methods
setGeneric("cached", function(x, ...) standardGeneric("cached"))
setMethod("cached", c(x="GdObject"), function(x, ...) {
  ## length(grep('cache', class(x), ignore.case=TRUE)) > 0
  inherits(x, "Cached")
})

setGeneric("rangeStrict", function(x, ...) standardGeneric("rangeStrict"))
setMethod("rangeStrict", c(x="GdObject"), function(x, ...) {
  if (!inherits(x, 'Cached')) TRUE else x@range.strict
})

setGeneric("rangeStrict<-", function(x, value) standardGeneric("rangeStrict<-"))
setReplaceMethod("rangeStrict", "GdObject", function(x, value) {
  if (inherits(x, "Cached")) {
    x@range.strict <- value
  }
  x
})
