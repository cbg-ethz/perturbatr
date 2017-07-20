library(data.table)
library(dtplyr)
library(dplyr)
library(knockout)
library(optparse)
library(igraph)
library(mvtnorm)
library(ggplot2)
library(uuid)

.create.noiseless.data <- function(rep.cnt, virus.cnt, genes.cnt, gene.vcov)
{
  screens.cnt    <- 5

  viruses      <- paste0("V", 1:virus.cnt)
  virus.effects <- rnorm(virus.cnt, 0, 1)
  if (virus.cnt == 2) { virus.effects[1] <- virus.effects[2] * -1 }
  virus.table <- data.table(Virus=viruses, VirusEffect=virus.effects)

  genes         <- rownames(gene.vcov)
  gene.effects  <- as.vector(mvtnorm::rmvnorm(1, sigma=gene.vcov))
  gene.table    <- data.table(GeneSymbol=genes, GeneEffect=gene.effects)

  screens         <- paste0("S", 1:screens.cnt)
  screen.effects  <- rnorm(screens.cnt, 0, 1)
  if (screens.cnt == 2) screen.effects[1] <- screen.effects[2] * -1
  screen.table    <- data.table(ScreenType=screens, ScreenEffect=screen.effects)

  VG  <- paste(sep=":", viruses, rep(genes, each=virus.cnt))
  vg.effects <- rnorm(genes.cnt*virus.cnt, 0, 1)
  vg.table <-  data.table(VG=VG, VirusGeneEffect=vg.effects)

  VS  <- paste(sep=":", viruses, rep(screens, each=virus.cnt))
  vs.effects <- rnorm(virus.cnt*screens.cnt, 0, 1)
  vs.table <-  data.table(VS=VS, VirusScreenEffect=vs.effects)

  # combine effect names and add the respective sum of effects
  effect.data <- as.data.table(expand.grid(genes, screens, viruses)) %>%
    dplyr::rename(GeneSymbol=Var1, ScreenType=Var2, Virus=Var3) %>%
    dplyr::group_by(GeneSymbol, ScreenType, Virus) %>%
    dplyr::mutate(VG=paste(Virus, GeneSymbol, sep=":"),
                  VS=paste(Virus, ScreenType, sep=":")) %>%
    ungroup

  effect.data <- dplyr::left_join(effect.data, virus.table, by="Virus")
  effect.data <- dplyr::left_join(effect.data, gene.table, by="GeneSymbol")
  effect.data <- dplyr::left_join(effect.data, screen.table, by="ScreenType")
  effect.data <- dplyr::left_join(effect.data, vg.table, by="VG")
  effect.data <-  dplyr::left_join(effect.data, vs.table, by="VS")

  effect.data <- effect.data %>%
    dplyr::group_by(VS, VG, ScreenType, GeneSymbol, Virus) %>%
    dplyr::mutate(Effect=sum(VirusEffect, GeneEffect, ScreenEffect, VirusGeneEffect, VirusScreenEffect)) %>%
    ungroup

  # repeat every line x times to have replicates
  effect.data <- effect.data[rep(1:.N,each=rep.cnt)]
  effect.data$Weight <- 1
  effect.data$Control <- 0
  effect.data
}

.diffuse <- function(effects, file=NULL, graph=NULL)
{
  res <- NA
  if (!is.null(file))
  {
    res <- knockout::diffuse(effects, method="mrw", path=file,
                             r=.25,delete.nodes.on.degree=3)
  }  else if (!is.null(graph))
  {
    res <- knockout::diffuse(effects, method="mrw", graph=abs(graph),
                             r=.25, delete.nodes.on.degree=3)
  }
  else stop("No correct graph given")
  return(list(diffusion=res$diffusion))
}

.effects <- function(ran)
{
  data.table(GeneSymbol=rownames(ran),Effect=ran[,1])
}

.graph <- function(n)
{
  repeat {
    tryCatch({
      adj <- sapply(1:n, function(e) rnorm(n*10, 0, 1))
      g   <- barabasi.game(n, directed=F)
      edj <- igraph::as_edgelist(g)
      for (i in nrow(edj))
      {
        adj[,edj[i, ]] <- adj[,edj[i, ]] + sign(rnorm(1)) * rbeta(1, 2, 1)
      }
      vcov <- cov(adj)
      invisible(mvtnorm::rmvnorm(1, sigma=vcov))
      break
    }, error=function(r){print(r)}, warning=function(r){print(r)})
  }
  vcov
}

.boot <- function(dat, rep.cnt, vir.cnt, v, graph)
{
  bench.list <- list()
  i <- 1

  repeat
  {
    s <- paste0("var:", v, ",vir:", vir.cnt, ",rep:", rep.cnt, ",bootstrap:", i)
    print(s)
    tryCatch({
      subs <- knockout:::bootstrap.svd.lmm.model.data(dat)
      bench.list[[s]] <- list(Var=v,
                              Rep=rep.cnt,
                              Vir=vir.cnt,
                              Bootstrap=i,
                              model=analyse(subs, graph=graph))
      i <- i + 1
    }, error=function(e)   { print(paste0("Didnt fit ", i, ": ", e)); i <- 10000 })
    if (i >= 101) break
  }

  bench.list
}

analyse <- function(md, file=NULL, graph=NULL)
{

  md <- dplyr::filter(md, Virus != "CVB")

  pmm.fit     <- lme4::lmer(Readout ~ Virus + (1 | GeneSymbol) + (1 | Virus:GeneSymbol),
                            data = md, weights = md$Weight, verbose = F)
  lmm.fit     <- lme4::lmer(Readout ~ Virus + (1 | GeneSymbol) + (1 | Virus:GeneSymbol) + (1 | ScreenType) + (1 | Virus:ScreenType),
                            data = md, weights = md$Weight, verbose = F)

  pmm.effects <- .effects(lme4::ranef(pmm.fit)[["GeneSymbol"]])
  pmm.mrw     <- .diffuse(pmm.effects, file, graph)

  lmm.effects <- .effects(lme4::ranef(lmm.fit)[["GeneSymbol"]])
  lmm.mrw     <- .diffuse(lmm.effects, file, graph)

  list(lmm.mrw=lmm.mrw, lmm.fit=lmm.effects,
       pmm.mrw=pmm.mrw, pmm.fit=pmm.effects)
}

ranking.stability.sythetic <- function(output.path, uid, virs.cnt, rep.cnt, var, gene.cnt=100)
{

  graph           <- .graph(gene.cnt)
  rownames(graph) <- colnames(graph) <- paste0("G", 1:gene.cnt)
  noiseless.data  <- .create.noiseless.data(rep.cnt, virs.cnt, gene.cnt, graph)
  noisy.data <- dplyr::mutate(noiseless.data, Readout=Effect+rnorm(nrow(noiseless.data), 0, var))

  bench.list            <- list()
  bench.list[["graph"]] <- graph
  bench.list[["data"]]  <- list(Var=0, Rep.cnt=rep.cnt, Vir=virs.cnt, Bootstrap=0, model=noiseless.data)

  s <- paste0("var:", var, ",vir:", virs.cnt, ",rep:", rep.cnt, ",bootstrap:0")
  m <- analyse(noisy.data, graph=graph)
  bench.list[[s]] <- list(Var=var, Rep=rep.cnt, Vir=virs.cnt, Bootstrap=0, model=m)

  for (vir.cnt in seq(2, virs.cnt))
  {
    print(vir.cnt)
    viruses         <- paste0("V", 1:vir.cnt)
    sample.data     <- dplyr::filter(noisy.data, Virus %in% viruses)
    if (rep.cnt >= 7)
    {
      bench.list <- c(bench.list, .boot(dat=sample.data,
                                        rep.cnt=rep.cnt,
                                        vir.cnt=vir.cnt,
                                        v=var,
                                        graph=graph))
    }
    else
    {
      s <- paste0("var:", var, ",vir:", vir.cnt, ",rep:", rep.cnt, ",bootstrap:0")
      m <- analyse(sample.data, graph=graph)
      bench.list[[s]] <- list(Var=var, Rep.cnt=rep.cnt, Vir=vir.cnt, Bootstrap=0, model=m)
    }
  }

  data.path   <- paste0(output.path, "/lmm_stability_",
                        "_synthetic_data",
                        "_viruscnt_", virs.cnt,
                        "_repcnt_", rep.cnt,
                        "_var_", var,
                        "_", uid,
                        ".rds")
  saveRDS(bench.list, data.path)
}

ranking.stability.bio <- function(model.data, graph.file, output.path, uid)
{
  bench.list              <- list()
  rank.all.data           <- analyse(model.data, file=graph.file)
  bench.list[["full"]]    <- rank.all.data

  vrs <- c("HCV", "DENV", "CHIKV", "SARS")
  for (idx in seq(2, length(vrs)))
  {
    dat <- dplyr::filter(model.data, Virus %in% vrs[1:idx])
    i   <- 1
    run <- 1
    repeat
    {
      print(paste0("Bootstrap bio: ", i, ",v: ", paste0(vrs[1:idx], collapse="_")))
      tryCatch({
        rnai.screen.sample <- knockout:::bootstrap.svd.lmm.model.data(dat)
        s <- paste0("bootstrap:", i, "virs:", paste0(vrs[1:idx], collapse="_"))
        bench.list[[s]] <- list(Bootstrap=i,
                                Virus= paste0(vrs[1:idx], collapse="_"),
                                analyse(rnai.screen.sample, file=graph.file))

        i <- i + 1
      }, error=function(e){ print(paste0("Didnt fit ", i, ": ", e)); i <- 10000 },
      warning=function(e) { print(paste0("Bootstrap bio warning: ",  99)) })
      if (i >= 101) break
      run <- run + 1
    }
  }

  dat.pth <- paste0(output.path, "/lmm_stability__bio_data_",, ".rds")
  saveRDS(bench.list, dat.pth)
}

run.stability.analysis <- function()
{
  option_list <- list(
    make_option(c("-v", "--virus"), action="store", help="virus count", type="integer"),
    make_option(c("-r", "--replicate"), action="store", help="replicate count", type="integer"),
    make_option(c("-s", "--sig"), action="store", help="standard deviation", type="integer"))
  opt_parser <- OptionParser(option_list=option_list)
  opt        <- parse_args(opt_parser)
  if (is.null(opt$virus) || is.null(opt$sig) || is.null(opt$replicate))
  {
    print_help(opt_parser)
    stop("Please provide correct arguments")
  }

  graph.file <- "./mappings/fi_flat.tsv"
  if (!file.exists(graph.file))
  {
    graph.file <- "/cluster/home/simondi/simondi/svd/data/fi_flat.tsv"
  }
  rna.file <- "./integrated_data_files/rnai_screen_normalized.rds"
  if (!file.exists(rna.file))
  {
    rna.file <- "/cluster/home/simondi/simondi/svd/data/rnai_screen_normalized.rds"
  }
  output.path <- "./results/lmm_stability_selection/"
  if (!file.exists(output.path))
  {
    output.path <- "/cluster/home/simondi/simondi/svd/results/stability_selection"
  }

  rnai.screen <- readRDS(rna.file)
  uid <- uuid::UUIDgenerate()

  model.data <- dplyr::filter(rnai.screen, Virus != "CVB")
  model.data <- model.data.lmm(model.data,
                               weights=list("pooled"=1.5, "single"=1))

  ranking.stability.bio(model.data,
                        graph.file,
                        output.path, uid)

  ranking.stability.sythetic(output.path,
                             uid,
                             virs.cnt=opt$virus,
                             rep.cnt=opt$replicate,
                             var=opt$sig)

  s <- warnings()
  print(s)
}

run.stability.analysis()
