BigWigTrack <- function(x, start=NULL, end=NULL, chromosome='chrNA', strand="*",
                        range.strict=FALSE, genome='??', name="BigWigTrack",
                        window=-1, windowSize=5, type="hist", ...) {
  if (missing(genome) || genome == '??') warning("No genome is provided")
  range <- GRanges()
  chromosome <- .chrName(chromosome)[1]

  ## Hopefully you're trying to pass in a BAM file
  if (is.character(x)) {
    x <- BigWigFile(x)
  }
  stopifnot(inherits(x, "BigWigFile"))
  stopifnot(file.exists(path(x)))
  range <- GRanges()
  new("BigWigTrack", bw.file=x, chromosome=chromosome, range=range, name=name,
      genome=genome, window=window, windowSize=windowSize, type=type, ...)
}

bigWigFile <- function(x) {
  stopifnot(is(x, "BigWigTrack"))
  x@bw.file
}

setMethod("initialize", "BigWigTrack",
function(.Object, ..., bw.file=BigWigFile(""), cache=new.env(), range.strict=FALSE) {
  if (missing(bw.file) || !is(bw.file, "BigWigFile") || !file.exists(path(bw.file))) {
    stop("bam required during initialize,BigWigTrack")
  }
  .Object@cache <- cache
  .Object@bw.file <- bw.file
  .Object@range.strict <- range.strict
  callNextMethod(.Object=.Object, ...)
})


## Gviz-methods
setMethod("subset", c(x="BigWigTrack"), function(x, from=NULL, to=NULL, sort=FALSE, ...) {
  ## Assuming you came here from plotTracks where the chromosome would have
  ## been set in the consolidateTrack method
  ## This is called three times
  ##   * directly from plotTracks;
  ##   * then during some setup;
  ##   * then from a scores() call -- this time from/to are NULL and we get hosed
  if (is.null(from) && is.null(to)) {
    return(x)
  }
  args <- list(...)
  si <- seqinfo(bigWigFile(x))
  if (inherits(from, "GenomicRanges")) {
    chr <- as.character(seqnames(from)[1])
    to <- max(end(from))
    from <- min(start(from))
  } else {
    chr <- if (is.character(args$chromosome)) args$chromosome else chromosome(x)
  }
  invalid <- !chr %in% names(seqlengths(si))
  if (invalid) {
    x@range <- GRanges()
    x@data <- matrix(numeric(), nrow=1)
  } else {
    if (chromosome(x) != chr) chromosome(x) <- chr
    if (is.null(from)) from <- 1L
    if (is.null(to)) to <- seqlengths(si)[chr]
    rng <- GRanges(chr, IRanges(from, to))
    bw.sel <- BigWigSelection(rng)
    res <- import(bigWigFile(x), selection=bw.sel, asRangedData=FALSE)
    x@data <- matrix(values(res)$score, nrow=1)
    x@range <- local({values(res)$score <- NULL; res})
  }

  is.strict <- rangeStrict(x)
  if (!rangeStrict(x)) rangeStrict(x) <- TRUE
  x <- callNextMethod(x, from=from, to=to, sort=sort, ...)
  if (!is.strict) rangeStrict(x) <- FALSE

  x
})
