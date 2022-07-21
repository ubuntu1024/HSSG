Fstat=function (formula, data = NULL, method = "bray", 
                    strata = NULL, contr.unordered = "contr.sum", contr.ordered = "contr.poly", 
                    parallel = getOption("mc.cores"), ...) {
  EPS <- sqrt(.Machine$double.eps)
  TOL <- 1e-07
  Terms <- terms(formula, data = data)
  lhs <- formula[[2]]
  lhs <- eval(lhs, data, parent.frame())
  formula[[2]] <- NULL
  rhs.frame <- model.frame(formula, data, drop.unused.levels = TRUE)
  op.c <- options()$contrasts
  options(contrasts = c(contr.unordered, contr.ordered))
  rhs <- model.matrix(formula, rhs.frame)
  options(contrasts = op.c)
  grps <- attr(rhs, "assign")
  qrhs <- qr(rhs)
  rhs <- rhs[, qrhs$pivot, drop = FALSE]
  rhs <- rhs[, 1:qrhs$rank, drop = FALSE]
  grps <- grps[qrhs$pivot][1:qrhs$rank]
  u.grps <- unique(grps)
  nterms <- length(u.grps) - 1
  if (nterms < 1) 
    stop("right-hand-side of formula has no usable terms")
  H.s <- lapply(2:length(u.grps), function(j) {
    Xj <- rhs[, grps %in% u.grps[1:j]]
    qrX <- qr(Xj, tol = TOL)
    Q <- qr.Q(qrX)
    tcrossprod(Q[, 1:qrX$rank])
  })
  if (inherits(lhs, "dist")) {
    if (any(lhs < -TOL)) 
      stop("dissimilarities must be non-negative")
    dmat <- as.matrix(lhs^2)
  }
  else if ((is.matrix(lhs) || is.data.frame(lhs)) && isSymmetric(unname(as.matrix(lhs)))) {
    dmat <- as.matrix(lhs^2)
    lhs <- as.dist(lhs)
  }
  else {
    dist.lhs <- as.matrix(vegdist(lhs, method = method, ...))
    dmat <- dist.lhs^2
  }
  n <- nrow(dmat)
  G <- -sweep(dmat, 1, rowMeans(dmat))/2
  SS.Exp.comb <- sapply(H.s, function(hat) sum(G * t(hat)))
  SS.Exp.each <- c(SS.Exp.comb - c(0, SS.Exp.comb[-nterms]))
  H.snterm <- H.s[[nterms]]
  tIH.snterm <- t(diag(n) - H.snterm)
  if (length(H.s) > 1) 
    for (i in length(H.s):2) H.s[[i]] <- H.s[[i]] - H.s[[i -1]]
  SS.Res <- sum(G * tIH.snterm)
  df.Exp <- sapply(u.grps[-1], function(i) sum(grps == i))
  df.Res <- n - qrhs$rank
  F.Mod <- (SS.Exp.each/df.Exp)/(SS.Res/df.Res)
  
  SumsOfSqs = c(SS.Exp.each, SS.Res)
  tab <- data.frame(Df = c(df.Exp, df.Res), SumsOfSqs = SumsOfSqs, 
                    F.Model = c(F.Mod, NA))
  rownames(tab) <- c(attr(attr(rhs.frame, "terms"), "term.labels")[u.grps], "Residuals")
  class(tab) <- c("anova", class(tab))
  tab
}
