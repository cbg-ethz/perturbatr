checkOptions <- function( args )
{
  i = grep("--l",args);
  LB <- ifelse( !any(i), 0, as.numeric(strsplit(args[i],"=")[[1]][2]))
  i = grep("--u",args);
  UB = ifelse( !any(i), 1, as.numeric(strsplit(args[i],"=")[[1]][2]))
  i = grep("--o",args);
  outputFile = ifelse( !any(i), NA, strsplit(args[i],"=")[[1]][2])
  i = grep("--i",args);
  inputFile = ifelse( !any(i), NA, strsplit(args[i],"=")[[1]][2])
  i = grep("--r$",args);
  reverse  = ifelse( !any(i), FALSE,TRUE)
  i = grep("-b",args);
  bonferroni <- any(i)
  bb = list(LB=LB,UB=UB,outputFile=outputFile,inputFile=inputFile,reverse=reverse,bonferroni=bonferroni);
  return (bb);
}

## Changed Here.
handleOneGroup <- function(i,dataset, optsb, t.rank = NULL)
{
  # browser()
  ## how many of them are up or low than cut-off.
  if(optsb$reverse)
  {
    i_max = sum(dataset$Score[i]>=optsb$LB)
    i_min = max(1,sum(dataset$Score[i]>=optsb$UB))
  }else
  {
    i_max = sum(dataset$Score[i]<=optsb$UB)
    i_min = max(1,sum(dataset$Score[i]<=optsb$LB))
  }
  ## t.r is the true rank, instead of order.
  if (is.null(t.rank))
  {
    t.r = i
  } else {
    t.r = t.rank[i]
  }
  # r = OPIScore(i,nrow(dataset),i_min,i_max)
  r = OPIScore(t.r,nrow(dataset),i_min,i_max,optsb$bonferroni)
  return ( cbind(
    LogP = r["logp"]
    ,OPI_Hit=as.numeric(seq(length(i))<=r["cutoff"])
    ,"#hitWell"=i_max
    ,"#totalWell"=length(i)
    ,rank = i))
}

OPIScore <- function(I_rank, N, i_min=1, i_max=-1, bonferroni=FALSE)
{
  n_drawn = length(I_rank) # number of black
  if(i_max == -1)
  {
    i_max=n_drawn
  }
  r1 = c(logp=1.0,cutoff=0)
  if( i_max < i_min) return (r1)
  # phyper(x, lower.tail = F), x = x-1, when lower.tail = F
  best_logp=1.0
  cutoff=0
  for (i in 1:n_drawn) {
    logp=phyper(i-1,I_rank[i],N-I_rank[i], n_drawn,lower.tail = F,log.p=F)
    logp = max(logp, 0)
    if (logp<=best_logp) {
      best_logp=logp
      cutoff=i
    }
  }
  if (FALSE) {
    best_logp=best_logp+log(i_max-i_min+1)/log(10)
  }
  return (c(logp=best_logp, cutoff=cutoff))

  #logp =  apply(cbind(seq(i_min,i_max),I_rank[i_min:i_max]),1,function(x) { phyper(x[1]-1,x[2] ,N-x[2], n_drawn,lower.tail = F,log.p=T)})
  #
  #	logp = logp/log(10)
  #	logp[logp<(-100)] = -100
  ##	if(all(is.na(logp))) {
  #		return  (r1)
  #	}else
  #		return   ( c(logp=min(logp),cutoff = i_min-1+which.min(logp)))
}

## Changed Here.
OPI<-function(Groups,Scores,opts,Data=NULL)
{
  t = data.frame(cbind(Gene_ID = Groups, Score= Scores))
  Sorted_Order = order(t$Score,decreasing=opts$reverse);
  Data = Data[Sorted_Order,]
  t = t[Sorted_Order,]

  ## get the ranks, "max" for the tie.
  t.rank <- rank(t$Score, ties.method="max")
  t = do.call("rbind", tapply(seq(nrow(t)), list(t$Gene_ID), handleOneGroup, dataset = t, opts, t.rank))
  # t = do.call("rbind", tapply(seq(nrow(t)), list(t$Gene_ID), handleOneGroup, dataset = t, opts))
  t = cbind(Data, t[order(t[,"rank"]),])

  # add OPI_Rank
  t = t[order(t$LogP,t$Score*ifelse(opts$reverse,-1,1)),]
  t$OPI_Rank = cumsum(t$OPI_Hit)
  t$OPI_Rank[t$OPI_Hit == 0] = 999999
  t
}


#t <- data.frame(Gene_ID=fr$genes, Score=-fr$readout , Well_ID=si)

t <- data.frame(Gene_ID=ret$GeneSymbol, Score=-abs(ret$Readout), Well_ID=ret$siRNAIDs)

t = subset(t, !is.na(Gene_ID) & Gene_ID != "" & !is.na(Score))
r = OPI(t$Gene_ID,t$Score,list(UB=100, LB=-100, reverse=F, bonferroni=F),t)

