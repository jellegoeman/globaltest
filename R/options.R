################################
# setting, storing and viewing options
################################
gt.opt.env <- new.env()
gt.opt.env$trace <- interactive()
gt.opt.env$trim <- FALSE
gt.opt.env$transpose <- FALSE
gt.opt.env$max.print <- Inf
gt.opt.env$warn.deprecated <- TRUE


gt.options <- function(trace, trim, transpose, max.print, warn.deprecated) {
  change <- FALSE
  if (!missing(trace)) {
    gt.opt.env$trace <- trace
    change <- TRUE
  }
  if (!missing(trim)) {
    gt.opt.env$trim <- trim
    change <- TRUE
  }
  if (!missing(transpose)) {
    gt.opt.env$transpose <- transpose
    change <- TRUE
  }
  if (!missing(max.print)) {
    gt.opt.env$max.print <- max.print
    change <- TRUE
  }
  if (!missing(warn.deprecated)) {
    gt.opt.env$warn.deprecated <- warn.deprecated
    change <- TRUE
  }
  if (change)
    invisible(as.list(gt.opt.env))
  else
    as.list(gt.opt.env)
}
