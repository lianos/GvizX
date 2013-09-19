## filter.fn is a function that takes a GRanges result (the reads) from the
## scanBam query and returns the ones you want to keep -- use in conjunction
## with `tags` or other columns returned from scanBam. In addition to
## seqnames,start,end,strand, `scan.bam.what` are columns that are added
## to the query
BamTrack <- function(x, start=NULL, end=NULL, width=NULL, chromosome='chrNA',
                     strand=NULL, genome='??', stacking="squish", name="BamTrack",
                     coverageOnly=FALSE, range.strict=FALSE,
                     window=-1, windowSize=5, type="hist",
                     scan.bam.what=c('mapq', 'cigar'),
                     scan.bam.flag=scanBamFlag(isUnmappedQuery=FALSE),
                     scan.bam.tags=NULL, filter.fn=NULL, ...) {
    if (missing(genome)) warning("No genome is provided")
    chromosome <- .chrName(chromosome)[1]

    ## Hopefully you're trying to pass in a BAM file
    if (is.character(x)) {
      x <- BamFile(x)
    }
    stopifnot(file.exists(path(x)))
    stopifnot(file.exists(paste(path(x), 'bai', sep='.')))

    ## Set range to be originally empty and fill it during a call to
    ## `subset` when plotTracks is triggered
    ## range <- .buildRange(range=range, start=start, end=end, width=width,
    ##                      args=list(strand=strand),
    ##                      defaults=list(strand=strand), chromosome=chromosome)
    ## str <- unique(as.character(GenomicRanges::strand(range)))
    ## if("*" %in% str)
    ##     stop("Only '+' and '-' strand information is allowed for AlignedReadTrack objects.")
    range <- GRanges()

    params <- ScanBamParam(what=scan.bam.what, flag=scan.bam.flag)
    if (is.character(scan.bam.tags)) {
      bamTag(params) <- scan.bam.tags
    }
    cache <- new.env()
    cache$filter.fn <- filter.fn
    ## And finally the object instantiation
    new("BamTrack", bam.file=x, cache=cache, params=params,
        chromosome=chromosome, range=range, name=name,
        genome=genome, stacking=stacking, coverageOnly=coverageOnly,
        range.strict=range.strict, window=window, windowSize=windowSize,
        type=type, ...)
}

bamFile <- function(x) {
  stopifnot(is(x, "BamTrack"))
  x@bam.file
}


setMethod("initialize", "BamTrack",
function(.Object, ..., bam.file=BamFile(""), params=ScanBamParam(),
         cache=new.env(), range.strict=FALSE) {
  .Object@bam.file <- bam.file
  .Object@params <- params
  .Object@cache <- cache
  .Object@range.strict <- range.strict
  callNextMethod(.Object=.Object, ...)
})

setValidity("BamFile", function(object) {
  file.exists(path(object))
  file.exists(paste(path(object), "bai", sep="."))
  TRUE
})


## Gviz-methods
setReplaceMethod("chromosome", "BamTrack", function(GdObject, value) {
    chr <- .chrName(value[1L])
    invalid <- !chr %in% names(seqlengths(bamFile(GdObject)))
    if (invalid) {
      warning("Invalid chromosome for BamFile")
    } else {
      GdObject@chromosome <- chr
    }
    GdObject
})

## Gviz-methods
setMethod("subset", c(x="BamTrack"), function(x, from=NULL, to=NULL, sort=FALSE, stacks=FALSE, ...) {
  ## Assuming you came here from plotTracks where the chromosome would have
  ## been set in the consolidateTrack method
  ## cat("subset,BamTrack :: from", from, "to", to, "\n")
  args <- list(...)
  si <- seqinfo(bamFile(x))

  if (inherits(from, "GenomicRanges")) {
    chr <- as.character(seqnames(from)[1])
    to <- max(end(from))
    from <- min(start(from))
  } else {
    chr <- if (is.character(args$chromosome)) args$chromosome else chromosome(x)
  }
  invalid <- !chr %in% seqlevels(si)
  if (invalid) {
    x@range <- GRanges()
  } else {
    if (chromosome(x) != chr) chromosome(x) <- chr
    if (is.null(from)) from <- 1L
    if (is.null(to)) to <- seqlengths(si)[chr]

    params <- x@params
    bamWhich(params) <- GRanges(chr, IRanges(from, to))
    res <- readGappedAlignments(path(bamFile(x)), param=params)
    meta <- values(res)
    meta$ngap <- ngap(res)
    meta$cigar <- cigar(res)

    gr <- as(res, "GRanges")
    rm.meta <- c('qname', 'rname', 'pos', 'qwidth', 'strand')
    for (wut in intersect(rm.meta, names(meta))) {
      meta[[wut]] <- NULL
    }
    values(gr) <- meta

    if (is.function(args$filter.fn)) {
      filter.fn <- args$filter.fn
    } else {
      filter.fn <- x@cache$filter.fn
    }

    if (length(gr) && is.function(filter.fn)) {
      gr <- filter.fn(gr)
    }
    x@range <- gr
  }

  if (length(x@range)) {
      x <- Gviz:::setCoverage(x)
  }
  is.strict <- rangeStrict(x)
  if (!rangeStrict(x)) rangeStrict(x) <- TRUE
  x <- callNextMethod(x, from=from, to=to, sort=sort, stacks=stacks)
  if (!is.strict) rangeStrict(x) <- FALSE

  return(x)
})

setAs("BamTrack", "DataTrack", function(from) {
  x <- ranges(from)
  if (length(x) == 0) {
    return(DataTrack())
  }
  cov <- coverage(from, strand="*")
  sel <- suppressWarnings(runValue(cov) != 0)
  ans <- DataTrack(start=start(cov)[sel], end=end(cov)[sel],
                   data=runValue(cov)[sel], name=names(from),
                   genome=genome(from), chromosome=chromosome(from))
  displayPars(ans) <- displayPars(from)
  ans
})


