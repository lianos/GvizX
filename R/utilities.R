## There are several scenarios where we want to convert an Rle into pieces
## amenable to pass into a DataTrack. This does that:
##
##   bits <- Rle2DataBits(some.Rle)
##   dt <- DataTrack(start=bits$start, end=bits$end, data=bits$data, ...)
Rle2DataBits <- function(x) {
  sel <- suppressWarnings(runValue(x) != 0)
  list(start=start(x)[sel], end=end(x)[sel], data=runValue(x)[sel])
}


getAnnoPackageName <- function (from, package = NULL) {
  is.anno.package <- length(grep("^org\\..*\\.db$", from) == 1L)
  if (is.anno.package) {
    if (!require(from, character.only = TRUE)) {
      stop("Unknown package: ", from)
    }
    from
  } else {
    annotationPackage(from, package = package)
  }
}

