#' Off-target correct and rank by hit magnitude
#'
#' @export
#' @import data.table
#'
#' @param obj  a summarized data.table
#' @param path  path (or file) to the target-relation matrices
#' @param do.pooled  boolean flag whether pooled libraries should also be corrected
#' @param ... additional arguments
correct <-
function
(
  obj,
  path,
  single=F,
  do.pooled=F,
  ...
)
{
  UseMethod("correct")
}

#' @noRd
#' @import data.table
correct.svd.data <-
function
(
  obj,
  path,
  single=F,
  do.pooled=F,
  ...
)
{
  if (missing(path)) stop("Please provide the path to your target-relation matrices
                          (or one of the files, wel'll parse it for you)!")
  res <-  .off.target.correct(obj, path, do.pooled,...)
  class(res) <- c("svd.analysed.offtc", "svd.analysed", class(res))
  invisible(res)
}

#' @noRd
#' @import data.table
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @importFrom dplyr select
#' @importFrom dplyr filter
.off.target.correct <-
function
(
 obj,
 path,
 single,
 do.pooled,
 ...
)
{
  #TODO: do all at once!!! dont split up matrices, but combine them to have more samples
  if (single)
  {
    stop("Single analysis not yet supported!")
      cordat <-
        dplyr::select(obj, Virus, Replicate, Plate, GeneSymbol, Entrez,
                      Readout, siRNAIDs,
                      ReadoutType, InfectionType, Library, Screen) %>%
        dplyr::filter(!is.na(Entrez), !is.na(GeneSymbol)) %>%
        dplyr::group_by(Virus, Library, Screen, ReadoutType, InfectionType) %>%
        dplyr::mutate(Grp = .GRP)
      grps <- unique(cordat$Grp)
      # TODO: exclude dharmacon
      new.dat <- do.call(
        "rbind",
        lapply
        (
          grps,
          function(grp)
          {
            print(paste("Doing grp: ", grp))
            curr.dat <- dplyr::filter(cordat, Grp==grp)
            fr <- do.correct(curr.dat)
            fr
          }
        )
      )
  } else {
    do.correct.all(obj, path, single)
  }
  invisible(new.dat)
}

#' @noRd
#' @import data.table
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr left_join
#' @importFrom dplyr summarize
do.correct <-
function
(
  obj,
  ...
)
{
  lib <- obj$Library[1]
  nor <- nrow(obj)
  obj <- dplyr::filter(dat, !is.na(Entrez), !is.na(siRNAIDs)) %>% ungroup
  if (nrow(obj) < nor) message("Rows with is.na(Entrez) have been removed!")
  ofc <- off.target.correct.grp(obj)
  nr <- nrow(ofc)
  fr <- data.table::data.table(Virus=dat$Virus[1],
                               Replicate=dat$Replicate[1],
                               Plate=dat$Plate[1],
                               Entrez=ofc$Entrez,
                               Readout=ofc$Readout,
                               ReadoutType=dat$ReadoutType[1],
                               InfectionType=dat$InfectionType[1],
                               Library=dat$Library[1],
                               Screen=dat$Screen[1])
  past.dat <- dplyr::select(dat, Entrez, GeneSymbol, siRNAIDs) %>%
    dplyr::group_by(Entrez, GeneSymbol ) %>%
    dplyr::summarize(siRNAIDs=concat(siRNAIDs, col="", sep=";"))
  fr <- dplyr::left_join(fr, past.dat, by="Entrez") %>%
    .[order(Virus, Replicate, Plate, GeneSymbol)] %>%
    data.table::setcolorder(c("Virus", "Replicate", "Plate", "GeneSymbol",
                              "Entrez", "Readout",  "ReadoutType",
                              "Library", "siRNAIDs","Screen", "InfectionType"))
  fr
}

#' @noRd
off.target.correct.grp <-
function
(
 dat
)
{
  siRNAIDs <- dat$siRNAIDs
  entrez   <- dat$Entrez
  readout  <- dat$Readout
  lib      <- dat$Library[1]
  screen   <- dat$Screen[1]
  vir      <- dat$Virus[1]
  if (!is.character(lib))
  {
    stop("Lib is not a string!")
  }
  gesp.dat <- switch(lib,
                     "Ambion"=get.data.ambion(siRNAIDs, entrez, readout,
                                              screen, lib, vir),
                     "DharmaconSMARTPool"=get.data.dharmacon(siRNAIDs, entrez,
                                                             readout, screen, lib, vir),
                     stop("Non-defined library used! Please report to maintainer!"))
  gesp <- gespeR(gesp.dat$rel.mat, gesp.dat$pheno.frame, screen, vir)
  gesp$Readout[is.na(gesp$Readout)] <- 0
  invisible(list(Readout=gesp$Readout,
                 Entrez=as.integer(levels(gesp$Entrez))[gesp$Entrez]))
}

#' @noRd
#' @importFrom gespeR Phenotypes
#' @importFrom gespeR gespeR
#' @import Matrix
gespeR <-
function
(
 rel.mat,
 pheno.frame,
 screen,
 vir
)
{
  phenos <- gespeR::Phenotypes(Matrix::Matrix(pheno.frame$readout, ncol=1),
                       ids=as.character(pheno.frame$sirnas),
                       pnames="pheno",
                       type="SSP")
  # fit a model using gesper
  gesper.fit <- gespeR::gespeR(phenotypes=phenos,
                               target.relations=rel.mat,
                               mode = "cv")
  gsps  <- as.vector(gesper.fit@GSP@values)
  genes <- as.vector(gesper.fit@GSP@ids)
  # todo what do which other siRNAS that have not been corrected
  invisible(data.frame(Readout=gsps, Entrez=genes))
}

#' @noRd
get.data.ambion <-
function
(
 siRNAIDs,
 entrez,
 readout,
 screen,
 lib,
 vir
)
{
  # get the relations and parse the cols and rows
  rel.mat <- load.rel.mat(lib, screen)
  # create phenotype frame and remove non-intersecting siRNAs
  pheno.list <- init.frames(siRNAIDs, entrez, readout,
                            screen, vir, rel.mat)
  pheno.list
}

#' @noRd
get.data.dharmacon <-
function
(
 siRNAIDs,
 entrez,
 readout,
 screen,
 lib,
 vir
)
{
  # get the raw target-relation matrix
  rel.mat <- load.rel.mat(lib, screen)
  # add rows with multiple rnas
  fr <- do.call(
    "rbind",
    lapply(
        1:length(siRNAIDs),
        function(i)
        {
          sirnas <- unlist(strsplit(siRNAIDs[i], ","))
          data.frame(sirnas=sirnas, entrez=entrez[i], readout=readout[i])
        }
      )
  )
  pheno.list <- init.frames(fr$sirnas, fr$entrez, fr$readout,
                            screen, vir, rel.mat)
  pheno.list
}

#' @noRd
init.frames <-
  function
(
  siRNAIDs,
  entrez,
  readout,
  screen,
  vir,
  rel.mat
)
{
  pheno.sirnas <- as.character(siRNAIDs)
  rel.sirnas   <- as.character(rel.mat@siRNAs)
  rel.genes    <- as.integer(rel.mat@genes)
  # get the siRNAs that are available in both matrices
  intr.sirnas      <- intersect(pheno.sirnas, rel.sirnas)
  pheno.sirna.idxs <- which(pheno.sirnas %in% intr.sirnas)
  rel.sirna.idxs   <- which(rel.sirnas %in% intr.sirnas)
  # remove indexes that are not found
  pheno.sirnas <- pheno.sirnas[pheno.sirna.idxs]
  readouts <- readout[pheno.sirna.idxs]
  rel.sirnas <- rel.sirnas[rel.sirna.idxs]
  vals <- rel.mat@values[rel.sirna.idxs, ]
  # order the pheno sirnas accoding to the ordering of the rel.matrix sirnas
  ord <- order(match(pheno.sirnas, rel.sirnas))
  pheno.ordered <- pheno.sirnas[ord]
  pheno.count <- table(pheno.ordered)
  # add redundant siRNAs to the target relation matrix
  vals <- vals[rep(1:nrow(vals), pheno.count), ]
  # assign to matrix
  rel.mat@siRNAs <- rownames(vals) <- pheno.ordered
  rel.mat@values <- vals
  # create phenotype frame
  pheno.frame <- data.frame(sirnas=siRNAIDs, readout=readout, entrez=entrez)
  #remove not found ids
  pheno.frame <- pheno.frame[pheno.sirna.idxs, ]
  pheno.frame <- pheno.frame[ord, ]
  #check if dimensions are so are right
  assertthat::assert_that(all(pheno.frame$sirnas == rel.mat@siRNAs))
  assertthat::assert_that(all(pheno.frame$sirnas == rownames(rel.mat@values)))
  list(pheno.frame=pheno.frame,
       rel.mat=rel.mat)
}

#' @noRd
#' @importFrom gespeR TargetRelations
load.rel.mat <-
function
(
  lib,
  screen
)
{
  if(!is.na(grep("druggable", tolower(screen))[1])) screen <- "druggable"
  if(!is.na(grep("dharma", tolower(lib))[1])) lib <- "dharmacon"
  files <- list.files(system.file("extdata", package="svd"),
                      pattern="target_rel", full.names=T)
  files <- files[grep(lib, files, ignore.case=T)]
  res.file <- files[grep(screen, files, ignore.case=T)]
  if (is.na(res.file[1])) stop(paste("No target relation matrix found for ",
                                     lib, screen))
  rel.mat <- gespeR::TargetRelations(res.file)
  invisible(rel.mat)
}
