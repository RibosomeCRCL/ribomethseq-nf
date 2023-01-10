# -------------------------------------------------
# plot COA columns
# plot points with optional rescaling

coa.plot.col <- function(dudi, xax=1, yax=2, rescale=NULL, margin=10, first=NULL,
                     fore.co=c(),
                     pch.co=17,
                     cex.co=1, col.co=1,
                     lab.co=c(),
                     pos.co=1,
                     tex.co=0.5,
                     txt.co=NULL,
                     fac.co=NULL,
                     fac.co.col=NULL,
                     fac.co.lab=NULL,
                     fac.co.ell=1,
                     col.ch=1,
                     xlim=NULL, ylim=NULL,
                     ...) {

  .expand <- function(r, f) {
    m <- (r[2] - r[1]) * f
    r + c(-m, m)
  }
  
  .plot <- function(xy, pch, cex, col, lab, pos, tex, txt, fac, fac.col, fac.ell, fac.lab) {
    points(xy, pch=pch, cex=cex, col=col, ...)
    if (! is.null(lab)) {
      indx <- which(row.names(xy) %in% lab)
      rxy  <- data.frame(xy)[indx,]
      txt <- if (is.null(txt)) row.names(rxy) else txt
      text(rxy, labels=txt, pos=pos, col=col, cex=tex)
    }
    if (! is.null(fac)) {
      if (is.null(fac.col)) fac.col <- rep('black', length(levels(fac)))
      if (is.null(fac.lab)) fac.lab <- levels(fac)
      s.class(xy, fac=fac, col=fac.col, cellipse=fac.ell, label=fac.lab, add.plot=T)
    }
  }
  
  .fore.plot <- function(xy, fore, pch, cex, col) {
    indx <- which(row.names(xy) %in% fore)
    .restrict <- function(arr) {
      if (length(arr) > 1) arr[indx] else arr
    }
    points(data.frame(xy)[indx,], pch=.restrict(pch), cex=.restrict(cex), col=.restrict(col))
  }
  
  dudi <- if (is.null(rescale)) dudi else coa.rescale(dudi, rescale)
  
  eig <- if (! is.null(dudi$eig.cor)) dudi$eig.cor else dudi$eig
  
  xlab <- paste("Axis ", xax, ': ', round(eig[xax] * 100 / sum(eig)), '%', sep='')
  ylab <- paste("Axis ", yax, ': ', round(eig[yax] * 100 / sum(eig)), '%', sep='')
  
  xy.co <- dudi$co[,c(xax, yax)]
  
  if(is.null(xlim)) xlim <- .expand(range(xy.co[,1]), margin/100)
  if(is.null(ylim)) ylim <- .expand(range(xy.co[,2]), margin/100)
  
  #
  # empty plot for labels
  #
  plot(x=0, y=0, cex=0, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
  
  if (is.null(first))
    first <- 'co'
  
  # actual plots     

    .plot(xy.co, pch=pch.co, cex=cex.co, col=col.co, lab=lab.co, pos=pos.co, tex=tex.co, txt=txt.co,
          fac=fac.co, fac.col=fac.co.col, fac.ell=fac.co.ell, fac.lab=fac.co.lab)

  
  # foreground points
  
    .fore.plot(xy.co, fore=fore.co, pch=pch.co, cex=cex.co, col=col.co)
  
  
  # crosshair at 0,0
  if (! is.null(col.ch)) abline(0,0,0,0, col=col.ch)
  
  #TRUE
}
# -------------------------------------------------
plot_coa <- function(count_matrix,col_for_color = NULL, axis = c(1,2)) {
  coa_calculated <- .compute_coa(count_matrix)
  return(coa.plot.col(coa_calculated,xax = 1,
  yax = 2,
  lab.co = colnames(coa_calculated$tab),
  main = paste0("Correspondence analysis\n of the raw counts on all genomic positions (", nrow(coa_calculated$tab),")")))
}

.compute_coa <- function(counts = NULL){
  counts["position"] <- NULL
  res_coa <- ade4::dudi.coa(counts[complete.cases(counts), ],
                            scannf = FALSE,
                            nf = 5)
  return(res_coa)
}
