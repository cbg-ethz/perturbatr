require(Matrix)

#' @noRd
#' @import data.table
#' @import dplyr select
#' @import dplyr filter
#' @import tidyr separate
write.target.rel.matrices <-
function()
{
  ma.path <- "/Users/simondi/PHD/data/data/sysvirdrug/target_relation_matrices/raw_matrices/"

  for (f in list.files(ma.path))
  {
    f <- "dharmacon_pooled_kinase"
    lib <- ifelse(!is.na(grep("ambion", f)[1]), "Ambion", "DharmaconSMARTPool")
    screen <- ifelse(!is.na(grep("kinase", f)[1]), "Kinome", "DruggableGenome")
    write.target.rel.matrix(f,
                            paste(ma.path, f, sep="/"),
                            lib, screen)
  }
}

#' @noRd
#' @import Matrix
write.target.rel.matrix <-
function
(
  dir.name,
  f,
  lib,
  screen
)
{
  outpath <- paste("./inst/extdata/", dir.name,
                   "_target_relation_matrix.rds", sep="")
  fls <- paste(f, list.files(f), sep="/")
  mat <- fls[grep("mtx", fls)]
  rownames <- fls[grep("rownames", fls)]
  colnames <- fls[grep("colnames", fls)]
  X <- readMM(mat)
  rn <- read.delim(rownames, stringsAsFactors = FALSE, header = FALSE)
  cn <- read.delim(colnames, stringsAsFactors = FALSE, header = FALSE)
  rownames(X)[rn[,1]] <- rn[,2]
  colnames(X)[cn[,1]] <- cn[,2]
  X <- as(1 - 2^X, "dgTMatrix")
  X[X < 0] <- 0
  saveRDS(X, file = outpath)
}

set.rel.mat.entries <- function()
{
  ma.path <- list.files("/Users/simondi/PROJECTS/sysvirdrug_project/svd/inst/extdata",
                        pattern="relation_matrix", full.names=T)
  load("inst/extdata/rnai.screen.raw.rda")
  for (f in ma.path)
  {
    lib <- ifelse(!is.na(grep("ambion", f)[1]), "Ambion", "DharmaconSMARTPool")
    screen <- ifelse(!is.na(grep("kinase", f)[1]), "Kinome", "DruggableGenome")
    dat <- dplyr::filter(rnai.screen.raw, Library==lib, Screen==screen) %>%
      dplyr::select(Entrez, siRNAIDs) %>% unique %>%
      dplyr::filter(!is.na(Entrez), !is.na(siRNAIDs), siRNAIDs != "empty")
    set.rel.mat.entries.single(f, dat$Entrez, dat$siRNAIDs, lib, screen)
  }
}

#' @noRd
set.rel.mat.entries.single <-
function
(
  f,
  entrez,
  sirnas,
  lib, screen
)
{
  py.path <- "/Users/simondi/PROJECTS/sysvirdrug_project/scripts_cml/set_relmat_entries.py"
  require(MASS)
  require(rPython)
  tmp.f  <- tempfile(pattern = "filem", tmpdir = tempdir(), fileext = "")
  tmp.ro <- tempfile(pattern = "filero", tmpdir = tempdir(), fileext = "")
  tmp.co <- tempfile(pattern = "filecol", tmpdir = tempdir(), fileext = "")
  tmp.fr <- tempfile(pattern = "filefr", tmpdir = tempdir(), fileext = "")
  tmp.X <- tempfile(pattern = "fileX", tmpdir = tempdir(), fileext = "")
  X <- as.matrix(readRDS(f))
  MASS::write.matrix(X, file=tmp.f, sep="\t")
  write(rownames(X), ncolumns=1, file=tmp.ro)
  write(colnames(X), ncolumns=1, file=tmp.co)
  write.table(data.frame(entrez, sirnas), file=tmp.fr, quote=F,
              row.names=F, col.names=T, sep="\t")
  python.load(py.path)
  python.call("run", tmp.f, tmp.ro, tmp.co, tmp.fr, tmp.X)
  X.r <- as.matrix(read.csv(tmp.X, sep="\t", header=F))
  rownames(X.r) <- rownames(X)
  colnames(X.r) <- colnames(X)
  X.r <- X.r[,which(colSums(X.r) != 0)]
  X.r <- as(X.r, "dgTMatrix")
  unlink(tmp.f)
  unlink(tmp.co)
  unlink(tmp.ro)
  unlink(tmp.fr)
  unlink(tmp.X)
  saveRDS(X.r, file = f)
}

drop.rows <- function()
{
  load("./inst/extdata/rnai.screen.raw.rda")
  mat.path <- "./inst/extdata/dharmacon_pooled_kinase_target_relation_matrix.rds"
  X <- readRDS(mat.path)
  dat <- as.data.frame(filter(rnai.screen.raw, Screen=="Kinome",
                              Library=="DharmaconSMARTPool") %>%
                         dplyr::select(Entrez, siRNAIDs) %>%
                         dplyr::filter(!is.na(Entrez), !is.na(siRNAIDs), siRNAIDs != "empty") %>%
                         unique)
  siRNAIDs <- as.character(dat$siRNAIDs)
  entrez <- as.integer(dat$Entrez)
  rel.sirnas <- as.character(rownames(X))
  rel.entrez <- as.integer(colnames(X))

  # remove siRNAs with non-max target relation values
  # otherwise the two matrices have different dimensionalities
  sirnas.to.use <- sapply(1:length(siRNAIDs), function(i)
  {
    to.use <- NA_character_
    sispl <- unname(unlist(strsplit(siRNAIDs[i], ",")))
    ent   <- entrez[i]
    rel.gene.idx <- which(rel.entrez == ent)
    rel.sirna.idxs <- which(rel.sirnas %in% sispl)
    # get the siRNA names that should be removed from the
    if (!is.na(rel.gene.idx[1]) & length(rel.gene.idx) == 1 &
        length(rel.sirna.idxs) >= 1 &
        is.whole(length(rel.sirna.idxs) / length(sispl)))
    {
      val.idx <- which.max(X[rel.sirna.idxs, rel.gene.idx])
      if (length(val.idx) != 1) val.idx <- 1
      to.use <- rel.sirna.idxs[val.idx]
    }
    else if ( length(rel.sirna.idxs) >= 1)
    {
      to.use <- rel.sirna.idxs[1]
      warning(paste("Could not find gene in dharmacon rel.mat
                    deconvolution", i, ent, siRNAIDs[i], ". Returning random index.\n"))
    }
    else
    {
      warning(paste("Could not find gene AND sirna in dharmacon rel.mat
                    deconvolution", i, ent, siRNAIDs[i], ".Returning nothing.\n"))
    }
    to.use
  })
  sirnas.to.use <- as.integer(sirnas.to.use)
  X <- X[sirnas.to.use,]
  saveRDS(X, file = mat.path)
}

#if (lib == "DharmaconSMARTPool") X <- drop.rows(X, sirnas, entrez)
do.mtx <- function()
{
  print("Writing matrices")
  write.target.rel.matrices()
  print("Dropping rows from pooled samples")
  drop.rows()
  print("Setting target")
  set.rel.mat.entries()
}
