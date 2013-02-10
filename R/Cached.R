# This is necessary so that we return NULL in the inner Gviz calls to
# .defaultRange for our cached tracks. When `rangeStrict(x)` is FALSE, we
# don't want any cached ranges to set the bounds on a plot

setMethod("start", "Cached", function(x) {
  if (rangeStrict(x)) callNextMethod(x) else NULL
})

setMethod("end", "Cached", function(x) {
  if (rangeStrict(x)) callNextMethod(x) else NULL
})
