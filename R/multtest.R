#==========================================================
# Adds multiple-testing-corrected p-values
#==========================================================
multtest <- function(gt, proc = c("FDR", "FWER")) {
  if (!is(gt, "gt.result")) stop("gt must be a gt.result object")
  proc <- match.arg(proc)
  if (missing(proc) || (proc == "FDR")) {
    mtproc <- "BH"
    colname <- "FDR.adjusted"
  } else {
    mtproc <- "Holm"
    colname <- "FWER.adjusted"
  }
  mt <- mt.rawp2adjp(p.value(gt))
  mt <- mt$adjp[sort.list(mt$index), mtproc]
  
  if (colname %in% colnames(gt@res)) {
    gt@res[colname] <- mt
  } else {
    gt@res <- cbind(gt@res, mt)
    colnames(gt@res)[ncol(gt@res)] = colname
  }
  gt
}
