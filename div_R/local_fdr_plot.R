cpg <- screen.pmms.loess.bsc.zplate$model$ccg

sss <- locfdr::locfdr(cpg$ccg, plot=1)

s <- localfdr(cpg$ccg)

zz <- s$hist.dat$zzz
lo <- min(zz)
up <- max(zz)
zzz <- pmax(pmin(zz, up), lo)
bre <- 120
breaks <- seq(lo, up, length = bre)
zh <- hist(zzz, breaks = breaks, plot = F)

x <- (breaks[-1] + breaks[-length(breaks)])/2
yall <- y <- zh$counts
K <- length(y)
N <- length(zz)
hits <- s$hist.dat
hist(hits$z, breaks = 100, xlab = " ", main = "")
zs <- hist(hits$, breaks = 100, xlab = " ", main = "")
for (k in 1:length(zs$counts)) {
  lines(c(hits$x[k], hits$x[k]), c(0, hits$yt[k]), lwd = 2, col = 6)
}
lines(hits$x, hits$f, lwd = 3, col = 3)
lines(hits$x, hits$p, lwd = 2, lty = 2, col = 4)


library(ggplot2)

hits <- s$hist.dat
fr <- data.frame(Z=hits$zvalues)
fr3 <- data.frame(X=hits$x, y=hits$yt, f=hits$f, f0=hits$f0)


ggplot(fr) + geom_histogram(aes(Z), bins = 130, alpha = 0.3) +
  geom_bar(data=fr3, aes(x=X, y=y), stat = "identity", color="blue") +
  #geom_line(data=fr3, aes(x=X, y=FA), colour = "red", size = 1, linetype=2) +
  #geom_line(data=fr3, aes(x=X, y=f), color="blue", size=1, linetype=1) +
  geom_line(data=fr3, aes(x=X, y=f0,  colour="p0*f0"), size=1, linetype=1) +
  geom_line(data=fr3, aes(x=X, y=f,  colour="f"), size=1, linetype=2) +
  ylab("Frequency") + xlab("Gene-pathogen effect") +
  theme_bw() +
  theme(text = element_text(size=14)) +
  scale_colour_manual(name="Density", values=c("green", "blue"))


all.res <- screen.pmms.loess.bsc.zplate$res$all.virus.results
single.res <- do.call("rbind", lapply(screen.pmms.loess.bsc.zplate$res$single.virus.results, function(i) i))
single.res$Virus <- unname(unlist(sapply(rownames(single.res), function(e) sub(".[[:digit:]]+", "", e))))
single.res <-
  group_by(single.res, Virus, Sign=sign(Effect)) %>% dplyr::summarize(cnt=n()) %>% ungroup %>%
  mutate(Count=cnt*Sign)

ggplot(single.res, aes(x=Virus, y=Count, fill=Sign)) +
  geom_bar(stat="identity", position="dodge") +
  theme_bw() +
  theme(legend.position="none") + geom_hline(yintercept=0) +
  ylab("Count hits") +
  geom_text(aes(x=Virus, y=ifelse(Count>0, Count+1, Count-1), label=Count), size=4, colour="black")



gene.per.virus <- group_by(all.res )
ggplot(all.res, aes(x=GeneSymbol, y=Effect)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(axis.x.text)






