apa.file <- function(..., .installed=FALSE) {
  ## TODO: Flip installed=TRUE when we're all said and done
  args <- list(...)
  if (.installed) {
    args <- c(list('extdata'), args, list(package="ApAcalypse"))
    fn < do.call(system.file, args)
  } else {
    args <- c(list("~/Documents/Thesis/ApAcalypse/inst/extdata"), args)
    fn <- do.call(file.path, args)
  }
  fn
}

library(Gviz)
library(GvizX)
source("~/Documents/Thesis/ApAcalypse/R/gviz.R")

gm <- readRDS(apa.file("gviz-anno-table.refGene.rds"))
gene.track <- TxDbTrack(gm, background.title="#FFFFFF")

bams <- c(bcells="/Users/stavros/cBio/papers/apa-1/data/bams/unique/clean-final/bcells-1.bam",
          brain="/Users/stavros/cBio/papers/apa-1/data/bams/unique/clean-final/brain.bam",
          breast="/Users/stavros/cBio/papers/apa-1/data/bams/unique/clean-final/breast.bam",
          ovary="/Users/stavros/cBio/papers/apa-1/data/bams/unique/clean-final/ovary.bam",
          skmuscle="/Users/stavros/cBio/papers/apa-1/data/bams/unique/clean-final/skmuscle.bam",
          testis="/Users/stavros/cBio/papers/apa-1/data/bams/unique/clean-final/testis.bam")
btrax <- sapply(names(bams), function(bf) {
    BamTrack(bams[[bf]], name=bf, genome='hg19', window=-1, windowSize=5,
             col='#0080ff', col.histogram='#0080ff', fill.histogram='#0080ff',
             background.title="#FFFFFF", col.title="black", col.axis="black")
}, simplify=FALSE)

bcells <- btrax$bcells
# IL6
bcells <- subset(bcells, 22771021, 22771848, chromosome='chr7')

plotPeaks("chr17", 8062012, 8064598, btrax[1:3], gene.track, same.y=FALSE)
plotPeaks("chr3",  169710731, 169716238, btrax, gene.track, same.y=FALSE)
plotPeaks("chr11", 62505657,  62507801, btrax, gene.track, same.y=FALSE)#TTC9C
