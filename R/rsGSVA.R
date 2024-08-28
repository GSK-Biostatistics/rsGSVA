#' Run GSVA algorithm while providing separate distribution for underlying kcdf
#'
#' @param X matrix of expression data. Rows are genes and columns are subjects
#' @param gs list of gene sets
#' @param train_expr matrix of reference expression data. Rows are subjects and columns are genes
#' @param min.sz min.sz from \code{GSVA::filterGeneSets}
#' @param max.sz max.sz from \code{GSVA::filterGeneSets}
#' @param verbose TRUE/FALSE whether log messages should be shown
#'
#' @author William Jordan
#'
#' @description
#' This function estimates GSVA enrichment scores (Hänzelmann et al, 2013) with given
#' test and target distribution samples. Expects desired test matrix, gene set, target matrix and min
#' and max sizing for filtering gene sets.
#'
#' @references
#'
#' Hänzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis for microarray and RNA-seq data. BMC bioinformatics. 2013 Dec;14:1-5.
#'
#' @export
rsGSVA <- function(X, gs, train_expr = NULL, min.sz=1, max.sz=500, verbose=TRUE){
  # expr <- GSVA:::.filterFeatures(X, "gsva")

  train_expr <- GSVA:::.filterFeatures(train_expr, "gsva")

  expr <- X[rownames(train_expr),]

  # common_rows <- intersect(rownames(expr), rownames(train_expr))
  #
  # expr <- expr[common_rows,]
  # train_expr <- train_expr[common_rows,]

  gset.idx.list <- gs
  mapped.gset.idx.list <- gset.idx.list
  mapped.gset.idx.list <- GSVA:::.mapGeneSetsToFeatures(mapped.gset.idx.list, rownames(expr))

  mapped.gset.idx.list <- GSVA::filterGeneSets(mapped.gset.idx.list,
                                               min.sz=max(1, min.sz),
                                               max.sz=max.sz)
  method = "gsva"
  kcdf = "Gussian"
  abs.ranking=FALSE
  min.sz=1
  max.sz=Inf
  parallel.sz=1L
  mx.diff=TRUE
  tau=switch(method, gsva=1, ssgsea=0.25, NA)
  ssgsea.norm=TRUE
  #verbose=TRUE
  BPPARAM=BiocParallel::SerialParam(progressbar=verbose)
  rnaseq <- FALSE
  kernel <- TRUE

  rval <- .rsGSVA(expr, mapped.gset.idx.list, train_expr, method, kcdf, rnaseq, abs.ranking,
                             parallel.sz, mx.diff, tau, kernel, ssgsea.norm,
                             verbose, BPPARAM)

  rval
}

# Non-export
.rsGSVA <- function(expr, gset.idx.list, train_expr,
                               method=c("gsva", "ssgsea", "zscore", "plage"),
                               kcdf=c("Gaussian", "Poisson", "none"),
                               rnaseq=FALSE,
                               abs.ranking=FALSE,
                               parallel.sz=1L,
                               mx.diff=TRUE,
                               tau=1,
                               kernel=TRUE,
                               ssgsea.norm=TRUE,
                               verbose=TRUE,
                               BPPARAM=BiocParallel::SerialParam(progressbar=verbose)){
  if (length(gset.idx.list) == 0)
    stop("The gene set list is empty! Filter may be too stringent.")

  if (any(lengths(gset.idx.list) == 1))
    warning("Some gene sets have size one. Consider setting 'min.sz > 1'.")

  parallel.sz <- as.integer(parallel.sz)
  if (parallel.sz < 1L)
    parallel.sz <- 1L

  ## because we keep the argument 'parallel.sz' for backwards compatibility
  ## we need to harmonize it with the contents of BPPARAM
  if (parallel.sz > 1L && class(BPPARAM) == "SerialParam") {
    BPPARAM=MulticoreParam(progressbar=verbose, workers=parallel.sz, tasks=100)
  } else if (parallel.sz == 1L && class(BPPARAM) != "SerialParam") {
    parallel.sz <- bpnworkers(BPPARAM)
  } else if (parallel.sz > 1L && class(BPPARAM) != "SerialParam") {
    bpworkers(BPPARAM) <- parallel.sz
  }

  if (class(BPPARAM) != "SerialParam" && verbose)
    cat(sprintf("Setting parallel calculations through a %s back-end\nwith workers=%d and tasks=100.\n",
                class(BPPARAM), parallel.sz))

  if(verbose)
    cat("Estimating GSVA scores for", length(gset.idx.list),"gene sets.\n")

  n.samples <- ncol(expr)
  n.genes <- nrow(expr)
  n.gset <- length(gset.idx.list)

  es.obs <- matrix(NaN, n.gset, n.samples, dimnames=list(names(gset.idx.list),colnames(expr)))
  colnames(es.obs) <- colnames(expr)
  rownames(es.obs) <- names(gset.idx.list)

  es.obs <- stable_score_compute.geneset.es(expr, gset.idx.list, 1:n.samples,
                                            train_expr,
                                            rnaseq=rnaseq, abs.ranking=abs.ranking,
                                            parallel.sz=parallel.sz,
                                            mx.diff=mx.diff, tau=tau, kernel=kernel,
                                            verbose=verbose, BPPARAM=BPPARAM)

  colnames(es.obs) <- colnames(expr)
  rownames(es.obs) <- names(gset.idx.list)

  es.obs
}

# Non-export
stable_score_compute.gene.density <- function(expr, sample.idxs, train_expr, rnaseq=FALSE, kernel=TRUE){
  # if no training set use input expression set as normal method
  if(is.null(train_expr)){
    print("no input training dist")
    train_expr <- expr
  }
  n.test.samples <- ncol(expr)
  n.genes <- nrow(expr)
  # n.density.samples <- length(sample.idxs)
  n.density.samples <- ncol(train_expr)

  gene.density <- NA
  if (kernel) {
    A = .C("matrix_density_R",
           as.double(t(train_expr[ ,1:n.density.samples, drop=FALSE])),
           as.double(t(expr)),
           R = double(n.test.samples * n.genes),
           n.density.samples,
           n.test.samples,
           n.genes,
           as.integer(rnaseq))$R

    ranks1 <<- A

    gene.density <- t(matrix(A, n.test.samples, n.genes))
  } else {
    print('else')
    gene.density <- t(apply(expr, 1, function(x, sample.idxs) {
      f <- ecdf(train_expr[sample.idxs])
      f(x)
    }, sample.idxs))
    gene.density <- log(gene.density / (1-gene.density))
  }

  return(gene.density)
}

# Non-export
stable_score_compute.geneset.es <- function(expr, gset.idx.list, sample.idxs, train_expr, rnaseq=FALSE,
                                            abs.ranking, parallel.sz=1L,
                                            mx.diff=TRUE, tau=1, kernel=TRUE,
                                            verbose=TRUE, BPPARAM=SerialParam(progressbar=verbose)) {
  num_genes <- nrow(expr)
  if (verbose) {
    if (kernel) {
      if (rnaseq)
        cat("Estimating ECDFs with Poisson kernels\n")
      else
        cat("Estimating ECDFs with Gaussian kernels\n")
    } else
      cat("Estimating ECDFs directly\n")
  }

  ## open parallelism only if ECDFs have to be estimated for
  ## more than 100 genes on more than 100 samples
  if (parallel.sz > 1 && length(sample.idxs > 100) && nrow(expr) > 100) {
    if (verbose)
      cat(sprintf("Estimating ECDFs in parallel on %d cores\n", as.integer(parallel.sz)))
    iter <- function(Y, n_chunks=BiocParallel::multicoreWorkers()) {
      idx <- splitIndices(nrow(Y), min(nrow(Y), n_chunks))
      i <- 0L
      function() {
        if (i == length(idx))
          return(NULL)
        i <<- i + 1L
        Y[idx[[i]], , drop=FALSE]
      }
    }
    gene.density <- BiocParallel::bpiterate(iter(expr, 100),
                                            stable_score_compute.gene.density,
                                            sample.idxs=sample.idxs,
                                            train_expr=train_expr,
                                            rnaseq=rnaseq, kernel=kernel,
                                            REDUCE=rbind, reduce.in.order=TRUE,
                                            BPPARAM=BPPARAM)
  } else
    gene.density <- stable_score_compute.gene.density(expr, sample.idxs, train_expr, rnaseq, kernel)

  stable_score_compute_rank_score <- function(sort_idx_vec){
    tmp <- rep(0, num_genes)
    tmp[sort_idx_vec] <- abs(seq(from=num_genes,to=1) - num_genes/2)
    return (tmp)
  }

  rank.scores <- rep(0, num_genes)
  sort.sgn.idxs <- apply(gene.density, 2, order, decreasing=TRUE) # n.genes * n.samples

  rank.scores <- apply(sort.sgn.idxs, 2, stable_score_compute_rank_score)

  m <- BiocParallel::bplapply(gset.idx.list, GSVA:::ks_test_m,
                              gene.density=rank.scores,
                              sort.idxs=sort.sgn.idxs,
                              mx.diff=mx.diff, abs.ranking=abs.ranking,
                              tau=tau, verbose=verbose,
                              BPPARAM=BPPARAM)
  m <- do.call("rbind", m)
  colnames(m) <- colnames(expr)

  return (m)
}

#' @useDynLib rsGSVA
NULL
