setDT(screen.summ.loess.bsc.zplate)[Screen=="DruggableGenome", Screen := "Genome"]
rep.bar <- group_by(screen.summ.loess.bsc.zplate, Virus, Screen) %>%
  dplyr::summarise(Replicates=length(unique(Replicate)),
                   Genes=length(unique(GeneSymbol)),
                   ReadoutType=) %>%
  gather(Type, Count, Replicates, Genes)
rep.bar


ggplot(rep.bar, aes(x=Virus, y = Count)) +
  geom_bar(aes(fill=Virus) ,stat="identity") +
  facet_grid(Type ~ Screen, scales='free_y') +
  scale_fill_brewer(palette="Spectral") +
  geom_text(aes(label = Count, y = Count), size = 5, vjust=-.25) +
  theme_bw() +
  theme(text=element_text(size=20))


chikv.hits.hstat.norm.bsc.zplate              <- prioritize(chikv.hstat.bsc.zplate,   hit.ratio=.6, readout.thresh=2)
sars.hits.hstat.norm.bsc.zplate               <- prioritize(sars.hstat.norm.bsc.zplate, hit.ratio=.6, readout.thresh=2)
hcv.g.hits.hstat.norm.bsc.zplate              <- prioritize(hcv.g.hstat.norm.bsc.zplate, hit.ratio=.6, readout.thresh=2)
hcv.k.hits.hstat.norm.loess.bsc.zplate        <- prioritize(hcv.k.hstat.norm.loess.bsc.zplate, hit.ratio=.6, readout.thresh=1)
denv.k.hits.hstat.norm.loess.bsc.zplate       <- prioritize(denv.k.hstat.norm.loess.bsc.zplate, hit.ratio=.6, readout.thresh=1)
denv.g.hits.hstat.norm.bsc.zplate             <- prioritize(denv.g.hstat.norm.bsc.zplate,     hit.ratio=.6 ,readout.thresh=2)


chikv.genes <- unique(chikv.hits.hstat.norm.bsc.zplate$Entrez)
sars.genes <- unique(sars.hits.hstat.norm.bsc.zplate$Entrez)
hcv.genes <- unique(c(hcv.k.hits.hstat.norm.loess.bsc.zplate$Entrez, hcv.g.hits.hstat.norm.bsc.zplate$Entrez))
denv.genes <- unique(c(denv.k.hits.hstat.norm.loess.bsc.zplate$Entrez, denv.g.hits.hstat.norm.bsc.zplate$Entrez))

chikv.go <- go.mapping(chikv.genes)
denv.go <- go.mapping(denv.genes)
hcv.go <- go.mapping(hcv.genes)
sars.go <- go.mapping(sars.genes)


library(VennDiagram)
library(grid)
library(gridBase)
library(lattice)
pdf("~/Desktop/pbe.pdf")
grid.layout(nrow=1, ncol=1)
s <- venn.diagram(list(CHIKV = chikv.genes,
                       SARS = sars.genes,
                       HCV=hcv.genes,
                        DENV=denv.genes
                       ),
                  fill = c("red", "green", "blue", "yellow"),
                  alpha = rep(0.5, 4), cex = 2,cat.fontface = 4,lty =2,
                  filename=NULL)
grid.draw(s)
dev.off()
s <- venn.diagram(list(CHIKV = chikv.go,
                       SARS = sars.go,
                       HCV=hcv.go,
                       DENV=denv.go
),
fill = c("red", "green", "blue", "yellow"),
alpha = rep(0.5, 4), cex = 2,cat.fontface = 4,lty =2,
filename=NULL)
grid.draw(s)
dev.off()

cool.hits <- filter(screen.hits.pmm.max.effect, FDR<.005 |
                      GeneSymbol %in% c("CD81", "TSG101", "COPB", "CDK6", "CDK5R2",
                                        "PIP5K1A", "PRKR", "HNRNPK", "APOE", "SPCS3",
                                        "EIF4A3", "EP300", "KNCK6", "CRABP1", "HNRPK") )


ggplot(cool.hits) +
  geom_bar(aes(x=GeneSymbol,y=abs(Effect), fill=Effect), stat="identity") +
  scale_fill_distiller(palette="Spectral") +
  theme_bw() +
  theme(text=element_text(size=22)) + ylab("Effect") +
  coord_polar()


y <- rnorm(100)
dfr <- data.table(X=rnorm(1000, 0,1),
                  Y=c(y*1.5 + 1, y*1.5  + 2 , y*1.5 + 3, y*1.5  + 4, y*1.5  + 5),
                  Gene=factor(rep(1:5, each=100)))
ggplot(dfr, aes(x=X, y=Y, color=Gene)) +
  geom_point() +
  geom_smooth(method="lm") +
  scale_color_brewer(palette="Spectral") +
  theme_bw() +
  theme(text=element_text(size=24))

