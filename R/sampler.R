#' @noRd
bootstrap <- function(model.data, level=c("sirna", "pathogen"))
{
  UseMethod("bootstrap")
}

#' @method bootstrap svd.lmm.model.data
#' @import data.table
#' @importFrom dplyr left_join mutate select group_by filter
bootstrap.svd.lmm.model.data <- function(model.data, level=c("sirna", "pathogen"))
{

  dat <-
    model.data %>%
    dplyr::group_by(Virus, InfectionType, GeneSymbol) %>%
    dplyr::mutate(cnt=n(), grp=.GRP) %>%
    ungroup
  grps <- unique(dat$grp)
  res <- do.call(
    "rbind",
    lapply(
      grps,
      function (g)
      {
        grp.dat <- dplyr::filter(dat, grp==g)
        idx <- sample(seq(grp.dat$cnt[1]), replace=T) %>% unique
        grp.dat[idx]
      }
    )
  )
  res <- as.svd.lmm.model.data(res)
  res
}

#' @noRd
loocv <- function(model.data, idx)
{
  UseMethod("loocv")
}

#' @method loocv svd.lmm.model.data
#' @import data.table
#' @importFrom dplyr left_join mutate select group_by filter
loocv <- function(model.data, idx)
{
  dat <-
    model.data %>%
    dplyr::group_by(Virus, InfectionType, GeneSymbol) %>%
    dplyr::mutate(cnt=n(), grp=.GRP) %>%
    ungroup
  grps <- unique(dat$grp)
  res <- do.call(
    "rbind",
    lapply(
      grps,
      function (g)
      {
        grp.dat <- dplyr::filter(dat, grp==g)
        idxs <- seq(grp.dat$cnt[1])[-idx]
        grp.dat[idxs]
      }
    )
  )
  res <- as.svd.lmm.model.data(res)
  res
}
