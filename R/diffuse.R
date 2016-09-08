#' Extend the results from a <code>svd.prioritized</code> object by network diffusion.
#'
#' @export
#'
#' @import data.table
#'
#' @param obj  an svd.data object
#' @param method  method that should be used for diffusion
#' \itemize{
#'   \item{neighbors }{ just looks at the neighbors :)}
#'   \item{mrw }{ do a Markov random walk}
#' }
#' @param path  path to the network
#' @param ...  additional parameters
diffuse <- function(obj, method=c("neighbors", "mrw"), path, ...)
  UseMethod("diffuse")

#' @noRd
#' @export
#' @import data.table igraph
diffuse.svd.prioritized.pmm <- function(obj,
                                        method=c("neighbors", "mrw"),
                                        path, ...)
{
  if (!file.exists(path)) stop(paste("Can't find: ", path, "! Yieks!", sep=""))
  graph <- .read.graph(path)
  res  <- .diffuse.lmm(obj, match.arg(method), graph, ...)
  class(res) <- c("svd.diffused.pmm", "svd.diffused", class(res))
  invisible(res)
}

#' @noRd
#' @import data.table igraph
.diffuse.lmm <- function(obj, method, graph, ...)
{
  switch(method,
         "neighbors"= .diffuse.lmm.knn(obj, graph, ...),
         "mrw"      = .diffuse.lmm.mrw(obj, graph, ...),
         stop("No suitable method found"))
}

#' @import data.table
.diffuse.lmm.mrw <- function(obj, graph, ...)
{
  .diffuse.lmm.mrw.pathogen.wise(obj, graph, ...)
}

#' @import data.table igraph
#' @import foreach parallel doParallel
#' @import Matrix
#' @importFrom iterators iter
#' @importFrom dplyr filter
#' @importFrom assertthat assert_that
.diffuse.lmm.mrw.pathogen.wise <- function(obj, gra, ...)
{
  phs <- obj$gene.pathogen.effect.hits
  vir <- unique(phs$Virus)
  adj.mat <-  igraph::get.adjacency(gra, attr="weight")
  W <- .stoch.col.norm(adj.mat)
  len <- nrow(W)
  adj.mat.genes <- colnames(adj.mat)
  neighbors <- foreach::foreach(v=iterators::iter(vir), .combine=cbind) %do% {
    flt <- dplyr::filter(phs, Virus==v)
    vir.gen <- flt$GeneSymbol
    vir.eff <- abs(flt$Effect)
    intr.gen <- intersect(adj.mat.genes, vir.gen)
    idxs.com   <- which(adj.mat.genes %in% intr.gen)
    idxs.set   <- which(vir.gen %in% intr.gen)
    p0 <- rep(0, len)
    p0[idxs.com] <- vir.eff[idxs.set]/sum(vir.eff[idxs.set])
    p.inf <- solve(diag(len) - .5 * W) %*% (.5 * p0)
    p.inf
  }
  mat <- .rankaggreg(neighbors, adj.mat.genes, k=25)
  colnames(mat) <- adj.mat.genes
  nei.tab <- data.table::data.table(Gene=adj.mat.genes,
                                    Mean=rowMeans(neighbors)) %>%
      .[order(Mean, decreasing=T)] %>% .[1:25]

  new.gen.idx <- which(adj.mat.genes %in% nei.tab$Gene)
  new.adj.mat <- adj.mat[new.gen.idx,new.gen.idx]
  nodes <- colnames(new.adj.mat)
  res.gr <- igraph::graph.adjacency(new.adj.mat,mode="undirected")
  igraph::V(res.gr)$color <- igraph::V(res.gr)$Color
  graph.info <- list(graph=res.gr,
                     colors=c("lightblue", "orange"),
                     legend=c("LMM", "'Diffusion'"),
                     type="1-NN",
                     tresh=2)
  list(hits=res, graph.info=graph.info)
  invisible(nei.tab)
}



#' @import data.table
.diffuse.lmm.knn <- function(obj, graph, ...)
{
  .diffuse.lmm.knn.pathogen.wise(obj, graph, ...)
}

#' @import data.table igraph
#' @import foreach parallel doParallel
#' @importFrom iterators iter
#' @importFrom dplyr filter
.diffuse.lmm.knn.pathogen.wise <- function(obj, gra, ...)
{
  phs <- obj$gene.pathogen.effect.hits
  vir <- unique(phs$Virus)
  all.edges  <- igraph::get.edgelist(gra)
  #cl <- parallel::makeCluster(parallel::detectCores() - 2)
  #doParallel::registerDoParallel(cl)
  neighbors <- foreach::foreach(v=iterators::iter(vir), .combine=rbind) %do%
  {
      virgen <- dplyr::filter(phs, Virus==v)$GeneSymbol
      idxs   <- which(all.edges[,1] %in% virgen | all.edges[,2] %in% virgen)
      fr <- t(apply(all.edges[idxs, ], 1, function(e) sort(e)))
      chosen.edges <- data.table(Virus=v,
                                    Gene1=c(fr[, 1]),
                                    Gene2=c(fr[, 2])) %>%
        unique
      chosen.edges
  }
 # parallel::stopCluster(cl)
  rel <- dplyr::group_by(neighbors, Gene1, Gene2) %>%
    dplyr::summarise(Count=n()) %>%
    dplyr::filter(Count > 1) %>%
    as.data.frame
  li <- unique(c(rel$Gene1, rel$Gene2))
  res <- dplyr::filter(neighbors, Gene1 %in% li | Gene2 %in% li ) %>%
    as.data.frame
  nodes <- data.table(Node=unique(c(res[, 2], res[, 3]))) %>%
    dplyr::mutate(Color=ifelse(Node %in% phs$GeneSymbol, "lightblue", "orange")) %>%
    dplyr::mutate(FromLMM=ifelse(Node %in% phs$GeneSymbol, 1, 0)) %>%
    as.data.frame
  edges <- res[, c(2, 3)]
  res.gr <- igraph::graph.data.frame(edges, directed=F, vertices=nodes)
  igraph::V(res.gr)$color <- igraph::V(res.gr)$Color
  graph.info <- list(graph=res.gr,
                     colors=c("lightblue", "orange"),
                     legend=c("Linear mixed model", "Diffusion"),
                     type="1-NN",
                     tresh=2)
  list(hits=res, graph.info=graph.info)
}
