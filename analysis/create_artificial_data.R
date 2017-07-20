
vs     <- paste0("V", 1:5)
wells  <- 384
plates <- 10

group_by(dat, Plate, Replicate) %>% summarize(n=n())

graph <- knockout:::.read.graph(path="~/PROJECTS/sysvirdrug_project/data/svd/mappings/fi_flat.tsv", graph=NULL)
gene.names <- V(graph) %>% unlist %>% sample %>% names

dat <- knockout::filter(rnai.screen.raw, Virus=="HCV", Screen=="Kinome")
dat <- dplyr::filter(dat@.data, Plate %in% 1:5, Replicate %in% 1:6)

dat$Virus <- "V1"
dat$ScreenType[dat$Replicate %% 2 == 0] <- "A/R"
dat$NumCells <- rpois(length(dat$NumCells), 1500)

dat <- group_by(dat, RowIdx, ColIdx) %>%
  dplyr::mutate(GeneSymbol=gene.names[.GRP],
                Entrez=.GRP,
                Readout=rgamma(1, 1000, 1),
                siRNAIDs=stringr::str_sub(uuid::UUIDgenerate(), 1, 5))
dat$Control <- 0
dat$Control[ dat$RowIdx %in% c(3,5) & dat$ColIdx %in% c(7, 10) ] <- 1
dat$Control[ dat$RowIdx %in% c(2,6) & dat$ColIdx %in% c(3, 4) ] <- -1

dat$Control[dat$GeneSymbol %in% unique(dat$GeneSymbol[dat$Control == -1])] <- -1
dat$Control[dat$GeneSymbol %in% unique(dat$GeneSymbol[dat$Control == 1])] <- 1

dat$Readout <- dat$Readout + rnorm(length(dat$Readout), 0, 1000)
dat$Control <- as.integer(dat$Control)

dat <- rbind(dat, dat)
dat$Virus <- rep(c("V1", "V2"), each=nrow(dat)/2)

rnaiscreen<- as(as.data.table(dat), "knockout.data")

devtools::use_data(rnaiscreen, overwrite=T)
