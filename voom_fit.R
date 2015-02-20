voom_fit=function (counts, design = NULL, lib.size = NULL, normalize.method = "none", 
          plot = FALSE, span = 0.5, ...) 
{
  out <- list()
  if (is(counts, "DGEList")) {
    out$genes <- counts$genes
    out$targets <- counts$samples
    if (is.null(design) && diff(range(as.numeric(counts$sample$group))) > 
          0) 
      design <- model.matrix(~group, data = counts$samples)
    if (is.null(lib.size)) 
      lib.size <- with(counts$samples, lib.size * norm.factors)
    counts <- counts$counts
  }
  else {
    isExpressionSet <- suppressPackageStartupMessages(is(counts, 
                                                         "ExpressionSet"))
    if (isExpressionSet) {
      if (length(fData(counts))) 
        out$genes <- fData(counts)
      if (length(pData(counts))) 
        out$targets <- pData(counts)
      counts <- exprs(counts)
    }
    else {
      counts <- as.matrix(counts)
    }
  }
  if (is.null(design)) {
    design <- matrix(1, ncol(counts), 1)
    rownames(design) <- colnames(counts)
    colnames(design) <- "GrandMean"
  }
  if (is.null(lib.size)) 
    lib.size <- colSums(counts)
  y <- t(log2(t(counts + 0.5)/(lib.size + 1) * 1e+06))
  y <- normalizeBetweenArrays(y, method = normalize.method)
  fit <- lmFit(y, design, ...)
  if (is.null(fit$Amean)) 
    fit$Amean <- rowMeans(y, na.rm = TRUE)
  sx <- fit$Amean + mean(log2(lib.size + 1)) - log2(1e+06)
  sy <- sqrt(fit$sigma)
  allzero <- rowSums(counts) == 0
  if (any(allzero)) {
    sx <- sx[!allzero]
    sy <- sy[!allzero]
  }
  l <- lowess(sx, sy, f = span)
  plot(sx, sy, xlab = "log2( count size + 0.5 )", ylab = "Sqrt( standard deviation )", 
       pch = 16, cex = 0.25)
  title("voom: Mean-variance trend")
  lines(l, col = "red")

  le = loess(sy~sx,span=span)
#  scatter.smooth(sx,sy,pch = 16, cex = 0.25, lpars = list(col = "cyan", lwd=2))
  scatter.smooth(sx,sy,pch = 16, cex = 0.25, lpars = list(col = "green", lwd = 2),degree=2)
  return(le)
}