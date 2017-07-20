library(data.table)
library(dtplyr)
library(dplyr)
library(knockout)
library(lme4)
library(utils)
library(tidyr)
library(formula.tools)
rm(list=ls())

do.lmm <- function(model.data, formulae)
{
  lmms <- list()
  for (f in formulae)
  {
      tryCatch(
        lmms[[f]] <- lme4::lmer(as.formula(f), data=model.data,
                                weights=model.data$Weight,
                                verbose = F, REML = F),
        error = function(e){ print(paste("Could not fit", f) ) }
      )
  }
  lmms
}

find.best <- function(lmm.fits, method=c("BIC", "logLik"))
{
  method <- match.arg(method)
  tab <- sapply(lmm.fits, function(e) {  c(Model=as.character(formula(e)),  AIC=AIC(e), BIC=BIC(e), logLik=logLik(e)) }) %>% t %>% as.data.table
  setDT(tab)[, logLik := as.double(logLik)]
  setDT(tab)[, BIC := as.double(BIC)]
  l <- list()
  if (method == "BIC") {
    l <- list(Best=tab$Model[which.min(tab$BIC)], Models=tab)
  } else {
    l <- list(Best=tab$Model[which.max(tab$logLik)], Models=tab)
  }
  l
}

mixed.effects.model.selection.stacking <- function(rnai.screen, model.data, starting.models)
{

  single.model.strings <- c("Cell", "InfectionType", "ReadoutType", "Design")
  model.formulas <- c(starting.models, paste0(starting.models, " + (1 | ", single.model.strings, ")"))
  lmm.fits <- do.lmm(model.data, model.formulas)
  best.lmm.single <- find.best(lmm.fits)

  dual.model.strings <- expand.grid(single.model.strings, single.model.strings)  %>% apply(1, sort) %>% t %>% as.data.table %>% unique %>% filter(V1 != V2)
  model.formulas <- c(best.lmm.single$Best,
                      paste0(starting.models, " + (1 | ", dual.model.strings$V1, ") + (1 | ", dual.model.strings$V2, ")"))
  lmm.fits <- do.lmm(model.data, model.formulas)
  best.lmm.double <- find.best(lmm.fits)


  triple.model.strings <- expand.grid(single.model.strings, single.model.strings, single.model.strings)  %>% apply(1, sort) %>% t %>% as.data.table %>% unique %>% filter(V1 != V2 & V1 != V3 & V2 != V3)
  model.formulas <- c(best.lmm.double$Best,
                      paste0(starting.models, " + (1 | ", triple.model.strings$V1, ") + (1 | ", triple.model.strings$V2, ") + (1 | " ,   triple.model.strings$V3, ")"))
  lmm.fits <- do.lmm(model.data, model.formulas)
  best.lmm.triple <- find.best(lmm.fits)

  model.formulas <- c(best.lmm.triple$Best,
                      "Readout ~ Virus + (1 | GeneSymbol) + (1 | Cell) +  (1 | Design) +  (1 | InfectionType) + (1 | ReadoutType)")
  lmm.fits <- do.lmm(model.data, model.formulas)
  best.lmm <- find.best(lmm.fits)

  best.lmm
}

mixed.effects.model.selection.aggregating <- function(rnai.screen, model.data, starting.models)
{
  single.model.strings <- c("Cell", "InfectionType", "ReadoutType", "Design")
  model.formulas <- c(starting.models,
                      paste0(starting.models,  " + (1 | ", single.model.strings, ")"))
  lmm.fits <- do.lmm(model.data, model.formulas)
  best.lmm.single <- find.best(lmm.fits)

  dual.model.strings <- expand.grid(single.model.strings, single.model.strings)  %>% apply(1, sort) %>% t %>% as.data.table %>% unique %>% filter(V1 != V2)
  model.formulas <- c(best.lmm.single$Best,
                      paste0(starting.models, " + (1 | ", dual.model.strings$V1, ") + (1 | ", dual.model.strings$V2, ")"))
  lmm.fits <- do.lmm(model.data, model.formulas)
  best.lmm.double <- find.best(lmm.fits)

  triple.model.strings <- expand.grid(single.model.strings, single.model.strings, single.model.strings)  %>% apply(1, sort) %>% t %>% as.data.table %>% unique %>% filter(V1 != V2 & V1 != V3 & V2 != V3)
  model.formulas <- c(best.lmm.double$Best,
                      paste0(starting.models, " + (1 | ", triple.model.strings$V1, ") + (1 | ", triple.model.strings$V2, ") + (1 | " ,   triple.model.strings$V3, ")"))
  lmm.fits <- do.lmm(model.data, model.formulas)
  best.lmm.triple <- find.best(lmm.fits)

  model.formulas <- c(best.lmm.triple$Best
                      , paste0(starting.models, " + (1 | Cell) +  (1 | Design) +  (1 | InfectionType) + (1 | ReadoutType)"))
  lmm.fits <- do.lmm(model.data, model.formulas)
  best.lmm <- find.best(lmm.fits)
  best.lmm
}

fixed.effects.model.selection.stacking <- function(rnai.screen, model.data, starting.models)
{

  single.model.strings <- c("Cell", "ReadoutType", "Design")
  model.formulas <- c(starting.models)
  for (i in seq(length(single.model.strings)))
  {
    co <- utils::combn(single.model.strings, i) %>% t
    for (j in 1:nrow(co))
    {
      c <- co[j, ]
      mo <- paste0("Readout ~ Virus + ", paste(c , collapse=" + "), " + (1 | GeneSymbol)")
      mo2 <- paste0("Readout ~ Virus + ", paste(c , collapse=" + "), " + (1 | GeneSymbol) + (1 | Virus:GeneSymbol)")
      model.formulas <- c(model.formulas, mo, mo2)
    }
  }
  lmm.fits <- do.lmm(model.data, model.formulas)
  best.lmm.single <- find.best(lmm.fits)

  single.model.strings <- c("Cell", "ReadoutType", "Design")
  model.formulas <- c(best.lmm.single$Best)
  for (i in seq(length(single.model.strings)))
  {
    co <- utils::combn(single.model.strings, i) %>% t
    for (j in 1:nrow(co))
    {
      c <- co[j, ]
      mo <- paste0("Readout ~ Virus + ", paste(c , collapse=" + "), " + (1 | GeneSymbol)")
      mo2 <- paste0("Readout ~ Virus + ", paste(c , collapse=" + "), " + (1 | GeneSymbol) + (1 | Virus:GeneSymbol)")
      mo3 <- paste0("Readout ~ Virus + ", paste(c , collapse=" + "), " + (1 | GeneSymbol) + (1 | Virus:GeneSymbol) + (1 | InfectionType)")
      mo4 <- paste0("Readout ~ Virus + ", paste(c , collapse=" + "), " + (1 | GeneSymbol) + (1 | Virus:GeneSymbol) + (1 | InfectionType) + (1 | InfectionType )")
      model.formulas <- c(model.formulas, mo, mo2)
    }
  }
  lmm.fits <- do.lmm(model.data, model.formulas)
  best.lmm <- find.best(lmm.fits)

  best.lmm
}

run.model.selection <- function()
{
  rna.file <- "~/PHD/data/data/svd/integrated_data_files/rnai_screen_normalized.rds"
  rnai.screen <- readRDS(rna.file)
  model.data <- knockout::model.data.lmm(rnai.screen, drop=T)

  # find best LMM
  starting.models <- c("Readout ~ Virus + (1 | GeneSymbol)", "Readout ~ Virus + (1 | GeneSymbol) + (1 | Virus:GeneSymbol)")
  best.stacked.model <- mixed.effects.model.selection.stacking(rnai.screen, model.data, starting.models)
  best.aggregated.model <- mixed.effects.model.selection.aggregating(rnai.screen, model.data, best.stacked.model$Best)
  best.stacked.model.2 <- mixed.effects.model.selection.stacking(rnai.screen, model.data, best.aggregated.model$Best)
  # the last call gives as result:
  # "Readout ~ Virus + (1 | GeneSymbol) + (1 | Virus:GeneSymbol) + (1 | InfectionType)"
  # from here subgrouping random effects can be tested
  formu <- c("Readout ~ Virus + (1 | GeneSymbol) + (1 | Virus:GeneSymbol) + (1 | InfectionType)",
             "Readout ~ Virus + (1 | GeneSymbol) + (1 | Virus:GeneSymbol) + (1 | InfectionType) + (1 | Virus:InfectionType)")
  fin <- do.lmm(model.data, formu)
  # result : take  "Readout ~ Virus + (1 | GeneSymbol) + (1 | Virus:GeneSymbol) + (1 | InfectionType) + (1 | Virus:InfectionType)"
  lmerTest::anova(fin[[1]], fin[[2]])
  starting.models <- c("Readout ~ Virus + (1 | GeneSymbol) + (1 | Virus:GeneSymbol) + (1 | InfectionType) + (1 | Virus:InfectionType)",
                       "Readout ~ Virus + (1 | GeneSymbol) + (1 | Virus:GeneSymbol) + (1 | InfectionType)",
                       "Readout ~ Virus + (1 | GeneSymbol) + (1 | Virus:GeneSymbol)")
  # compare best LMM to
  best.saturated.model <- fixed.effects.model.selection.stacking(rnai.screen, model.data, starting.models)
  cat(paste("Best model:", best.saturated.model$Best))
}
