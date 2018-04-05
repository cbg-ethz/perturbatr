#!/usr/bin/env Rscript

# perturbatr: analysis of high-throughput gene perturbation screens
#
# Copyright (C) 2018 Simon Dirmeier
#
# This file is part of perturbatr
#
# perturbatr is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# perturbatr is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with perturbatr. If not, see <http://www.gnu.org/licenses/>.


library(dplyr)
library(tibble)
library(perturbatr)
library(lme4)
library(utils)
library(tidyr)
library(formula.tools)
options(warn=-1)


do.lmm <- function(model.data, formulae)
{
  lmms <- list()
  for (f in formulae)
  {
      tryCatch(
        lmms[[f]] <- lme4::lmer(as.formula(f), data=model.data,
                                weights=model.data$Weight,
                                verbose = F, REML = F),
        error = function(e){ print(paste("Could not fit", f, e) ) }
      )
  }
  lmms
}


find.best <- function(lmm.fits, method=c("BIC", "logLik"))
{
  method <- match.arg(method)
  tab <- sapply(
    lmm.fits,
    function(e)
    {
      c(Model=as.character(formula(e)),  AIC=AIC(e), BIC=BIC(e), logLik=logLik(e))
    }
  ) %>%
    t() %>%
    as.tibble() %>%
    dplyr::mutate(logLig = as.double(logLik), BIC = as.double(BIC))

  l <- list()
  if (method == "BIC") {
    l <- list(Best=tab$Model[which.min(tab$BIC)], Models=tab)
  } else {
    l <- list(Best=tab$Model[which.max(tab$logLik)], Models=tab)
  }

  l
}


mixed.effects.model.selection.stacking <- function(rnai.screen,
                                                   model.data,
                                                   starting.models)
{
  single.model.strings <- c("Cell", "ScreenType", "ReadoutType", "Design")
  model.formulas  <- c(starting.models,
                       paste0(starting.models,
                              " + (1 | ", single.model.strings, ")"))
  lmm.fits        <- do.lmm(model.data, model.formulas)
  best.lmm.single <- find.best(lmm.fits)

  dual.model.strings <- expand.grid(single.model.strings,
                                    single.model.strings) %>%
    apply(1, sort) %>%
    t() %>%
    as.tibble() %>%
    unique() %>%
    dplyr::filter(V1 != V2)

  model.formulas <- c(best.lmm.single$Best,
                      paste0(starting.models,
                             " + (1 | ",
                             dual.model.strings$V1,
                             ") + (1 | ",
                             dual.model.strings$V2, ")"))

  lmm.fits <- do.lmm(model.data, model.formulas)
  best.lmm.double <- find.best(lmm.fits)

  triple.model.strings <- expand.grid(single.model.strings,
                                      single.model.strings,
                                      single.model.strings)  %>%
    apply(1, sort) %>%
    t() %>%
    as.tibble() %>%
    unique() %>%
    dplyr::filter(V1 != V2 & V1 != V3 & V2 != V3)

  model.formulas <- c(best.lmm.double$Best,
                      paste0(starting.models,
                             " + (1 | ",
                             triple.model.strings$V1,
                             ") + (1 | ",
                             triple.model.strings$V2,
                             ") + (1 | ",
                             triple.model.strings$V3, ")"))
  lmm.fits <- do.lmm(model.data, model.formulas)
  best.lmm.triple <- find.best(lmm.fits)

  model.formulas <- c(
    best.lmm.triple$Best,
    "Readout ~ Condition + (1 | GeneSymbol) + (1 | Cell) +  (1 | Design) +  (1 | ScreenType) + (1 | ReadoutType)")
  lmm.fits <- do.lmm(model.data, model.formulas)
  best.lmm <- find.best(lmm.fits)

  best.lmm
}


mixed.effects.model.selection.aggregating <- function(rnai.screen,
                                                      model.data,
                                                      starting.models)
{
  single.model.strings <- c("Cell", "ScreenType", "ReadoutType", "Design")
  model.formulas <- c(starting.models,
                      paste0(starting.models,
                             " + (1 | ",
                             single.model.strings,
                             ")"))
  lmm.fits        <- do.lmm(model.data, model.formulas)
  best.lmm.single <- find.best(lmm.fits)

  dual.model.strings <- expand.grid(single.model.strings, single.model.strings)  %>%
    apply(1, sort) %>%
    t() %>%
    as.tibble() %>%
    unique() %>%
    dplyr::filter(V1 != V2)

  model.formulas <- c(best.lmm.single$Best,
                      paste0(starting.models,
                             " + (1 | ",
                             dual.model.strings$V1,
                             ") + (1 | ",
                             dual.model.strings$V2,
                             ")"))

  lmm.fits        <- do.lmm(model.data, model.formulas)
  best.lmm.double <- find.best(lmm.fits)

  triple.model.strings <- expand.grid(single.model.strings,
                                      single.model.strings,
                                      single.model.strings)  %>%
    apply(1, sort) %>%
    t() %>%
    as.tibble() %>%
    unique() %>%
    dplyr::filter(V1 != V2 & V1 != V3 & V2 != V3)

  model.formulas <- c(best.lmm.double$Best,
                      paste0(starting.models,
                             " + (1 | ",
                             triple.model.strings$V1,
                             ") + (1 | ",
                             triple.model.strings$V2,
                             ") + (1 | " ,
                             triple.model.strings$V3, ")"))

  lmm.fits        <- do.lmm(model.data, model.formulas)
  best.lmm.triple <- find.best(lmm.fits)

  model.formulas <- c(
    best.lmm.triple$Best,
    paste0(starting.models, " + (1 | Cell) +  (1 | Design) +  (1 | ScreenType) + (1 | ReadoutType)"))

  lmm.fits <- do.lmm(model.data, model.formulas)
  best.lmm <- find.best(lmm.fits)

  best.lmm
}


fixed.effects.model.selection.stacking <- function(rnai.screen,
                                                   model.data,
                                                   starting.models)
{
  single.model.strings <- c("Cell", "ReadoutType", "Design")
  model.formulas       <- c(starting.models)
  for (i in seq(length(single.model.strings)))
  {
    co <- utils::combn(single.model.strings, i) %>% t()
    for (j in 1:nrow(co))
    {
      c <- co[j, ]
      mo  <- paste0("Readout ~ Condition + ", paste(c , collapse=" + "), " + (1 | GeneSymbol)")
      mo2 <- paste0("Readout ~ Condition + ", paste(c , collapse=" + "), " + (1 | GeneSymbol) + (1 | Condition:GeneSymbol)")
      model.formulas <- c(model.formulas, mo, mo2)
    }
  }

  lmm.fits <- do.lmm(model.data, model.formulas)
  best.lmm.single <- find.best(lmm.fits)

  single.model.strings <- c("Cell", "ReadoutType", "Design")
  model.formulas       <- c(best.lmm.single$Best)
  for (i in seq(length(single.model.strings)))
  {
    co <- utils::combn(single.model.strings, i) %>% t
    for (j in 1:nrow(co))
    {
      c   <- co[j, ]
      mo  <- paste0("Readout ~ Condition + ", paste(c , collapse=" + "), " + (1 | GeneSymbol)")
      mo2 <- paste0("Readout ~ Condition + ", paste(c , collapse=" + "), " + (1 | GeneSymbol) + (1 | Condition:GeneSymbol)")
      mo3 <- paste0("Readout ~ Condition + ", paste(c , collapse=" + "), " + (1 | GeneSymbol) + (1 | Condition:GeneSymbol) + (1 | ScreenType)")
      mo4 <- paste0("Readout ~ Condition + ", paste(c , collapse=" + "), " + (1 | GeneSymbol) + (1 | Condition:GeneSymbol) + (1 | ScreenType) + (1 | ScreenType:GeneSymbol)")
      model.formulas <- c(model.formulas, mo, mo2)
    }
  }

  lmm.fits <- do.lmm(model.data, model.formulas)
  best.lmm <- find.best(lmm.fits)

  best.lmm
}

run <- function()
{
  cat(paste("Loading data\n"))
  rna.file    <- "data/rnai_screen_normalized_2.rds"
  rnai.screen <- readRDS(rna.file)
  cat(paste("Setting model data\n"))
  model.data <- methods::as(rnai.screen, "PerturbationData")
  model.data <- perturbatr:::setModelData(model.data, drop=T)

  cat(paste("Doing model selection on non-nested models\n"))
  starting.models       <- c("Readout ~ Condition + (1 | GeneSymbol)", "Readout ~ Condition + (1 | GeneSymbol) + (1 | Condition:GeneSymbol)")
  best.stacked.model    <- mixed.effects.model.selection.stacking(rnai.screen, model.data, starting.models)
  best.aggregated.model <- mixed.effects.model.selection.aggregating(rnai.screen, model.data, best.stacked.model$Best)
  best.stacked.model.2  <- mixed.effects.model.selection.stacking(rnai.screen, model.data, best.aggregated.model$Best)
  cat(paste(best.stacked.model.2, "\n"))

  # the last call gives as result:
  # "Readout ~ Condition + (1 | GeneSymbol) + (1 | Condition:GeneSymbol) + (1 | ScreenType)"
  # from here subgrouping random effects can be tested
  cat(paste("Doing model selection on nested-models\n"))
  formu <- c("Readout ~ Condition + (1 | GeneSymbol) + (1 | Condition:GeneSymbol) + (1 | ScreenType)",
             "Readout ~ Condition + (1 | GeneSymbol) + (1 | Condition:GeneSymbol) + (1 | ScreenType) + (1 | Condition:ScreenType)")
  fin <- do.lmm(model.data, formu)
  # result : take  "Readout ~ Condition + (1 | GeneSymbol) + (1 | Condition:GeneSymbol) + (1 | ScreenType) + (1 | Condition:ScreenType)"
  lmerTest::anova(fin[[1]], fin[[2]])
  starting.models <- c("Readout ~ Condition + (1 | GeneSymbol) + (1 | Condition:GeneSymbol) + (1 | ScreenType) + (1 | Condition:ScreenType)",
                       "Readout ~ Condition + (1 | GeneSymbol) + (1 | Condition:GeneSymbol) + (1 | ScreenType)",
                       "Readout ~ Condition + (1 | GeneSymbol) + (1 | Condition:GeneSymbol)")
  cat(paste("Doing model selection on fixed effets models\n"))
  best.saturated.model <- fixed.effects.model.selection.stacking(rnai.screen, model.data, starting.models)

  cat(paste("Best model:", best.saturated.model$Best, "\n"))
}

run()
