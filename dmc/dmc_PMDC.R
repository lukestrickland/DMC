
# load_model ("LBA","lbaN_B.R")
require("gridExtra")
require("lme4")
require("plyr")
require("dplyr")
require("data.table")

theme_set(theme_simple())

singlerep.mean.sd <- function(rsamples, fun) {
  inference <- list()
  for (i in 1:length(rsamples)) {
    thetas <- rsamples[[i]]$theta
    inference [[i]] <- fun (thetas)
  }
  inf2 <- unlist(inference)
  dim3 <- c(dim(inference[[1]]), length(inf2)/prod(dim(inference[[1]])))
  dim(inf2) <- dim3
  test<- apply(inf2, 4, function(x) c(mean(x),sd(x)))
  test<- apply(test,1,mean); names(test) <- c("M", "SD")
  paste(round(test[1],3), " (", round(test[2],3), ")", sep="")
}



get.95.50 <- function (rsamples, p.vector, fun){ 
  test <- c(h.check.function.recovery.dmc(
  rsamples, p.vector, fun), h.check.function.recovery.dmc(
  rsamples, p.vector, fun, CI=c(0.25,0.75)))
  paste(round(test[1],2), "% / ", round(test[2],2), "%", sep="")
  }

h.check.function.recovery.dmc <- function(rsamples, p.vector, fun, CI= c(0.025, 0.975)) {
#Get the posterior mean and credible intervals for each replicate
  inference <- list()
  for (i in 1:length(rsamples)) {
    thetas <- rsamples[[i]]$theta
    inference [[i]] <- fun (thetas)
  }
  inf2 <- unlist(inference)
  dim3 <- c(dim(inference[[1]]), length(inf2)/prod(dim(inference[[1]])))
  dim(inf2) <- dim3
  test<- apply(inf2, 4, function(x) c(mean(x), quantile(x, probs=CI)))
#Hacky way to reshape p.vector same way to get functions in same way
#I should probably make fun work with sim.p.vector rather than thetas for future work,
#simpler, but already too far down the rabbit hole with this way for this paper
  pnames <- names(p.vector)
  dim(p.vector) <- c(1,length(p.vector),1); colnames(p.vector) <- pnames
  actual <- fun(p.vector)
#percentage of time CI contains actual  
  sum(apply(test,2, function(x) {x[2]<actual & x[3] >actual})) / length(test[1,]) *100
}

get.cors <- function(thetas, data) {
  
  CORS <- apply(thetas, c(1,2,3), function(x) cor(x, data))

  RAV <- apply(CORS, c(2), postRav, n=length(data), kappa=1, spacing=.01)
  post_medians <- apply(RAV, c(2), postRav.mean)
  post_LCI <- apply(RAV, c(2), function(x) postRav.ci(x)[1])
  post_HCI <- apply(RAV, c(2), function(x) postRav.ci(x)[2])
  out <- list(post_medians, post_LCI, post_HCI)
  names(out) <- c("medians", "LCI", "HCI")
  out
}

get.subj.effects.m <- function (PPS, fun, names) {
  EFFECTS <- lapply(PPS, get.subj.effects, fun=fun)
  for (i in 1:length(names)) EFFECTS[[i]]$model <- names[i]
  do.call(rbind, EFFECTS)
}

get.subj.effects <- function(PP, fun) {
  for (i in 1:length(PP)){
    effects<-get.subj.pp.MCI (PP[[i]], fun)
    cat(i)
    if (i ==1) out <- effects else out <- rbind(out,effects)
  }
  
  ENS <- rownames(out)
  OUT <- data.frame(out)
  OUT$effect <- ENS
  colnames(OUT) <- c("mean", "lower", "upper", "data", "effect")
  OUT
}


subj.meanthetas <- function (samples){
    samps <- lapply(samples, function(x) x["theta"])
    ## thetas into big array for apply
    samps2<- unlist(samps)
    dim3 <- c(dim(samps[[1]]$theta), length(samps2)/prod(dim(samps[[1]]$theta)))
    dim(samps2) <- dim3
    samps3<- apply(samps2, c(4,2), mean)
    ## back to a theta list after applied
    colnames(samps3)<- colnames(samps[[1]]$theta)
    df <- cbind(names(samples), data.frame(samps3))
    names(df)[1] <- "s"
    df
}

get.subj.pp.MCI <- function(sim, fun) {
 
  data <- attr(sim, "data")
  nreps=max(sim$reps)
  data.effects<- fun(data)
  noutput <- length(data.effects)

  sim.effects <- matrix(NA, nrow= nreps, ncol=noutput+1)
  sim.effects[,noutput+1] <- 1:nreps

  colnames(sim.effects) <- c(names(data.effects), "n.rep")

  for (j in 1:nreps) {

    currentsim.effects <- sim[sim$reps==j,]
    sim.effects[j,1:noutput] <- fun(currentsim.effects)

  }
  out <- apply(sim.effects, 2,function(x) c(mean(x, na.rm=T), quantile(x, probs=c(0.025,0.975), na.rm=T))) 
  cbind(t(out[,(!colnames(out) %in% "n.rep")]), data.effects)
  }

  
  
plot.acc.effects.x2 <- function(plot.df, mnam) {
  plot.df <- plot.df[-c(1:4),]
  plot.df$PM <- NA
  plot.df$PM[grep ("imp", rownames(plot.df))] <- "imp - control"
  plot.df$PM[grep ("unimp", rownames(plot.df))] <- "unimp - control"
  plot.df$PM[grep ("pmeffect", rownames(plot.df))] <- "imp- unimp"
  plot.df$PM[grep ("pmrt", rownames(plot.df))] <- "unimp- imp"
  
  plot.df$S <- NA
  plot.df$S[grep ("n", rownames(plot.df))]<- "Non-word Trial"
  plot.df$S[grep ("w", rownames(plot.df))]<- "Word Trial"
  plot.df$S[grep ("pmrtdiffw", rownames(plot.df))]<- "PM W Trial"
  plot.df$S[grep ("pmrtdiffn", rownames(plot.df))]<- "PM N Trial"
  plot.df$S[grep ("pmeffectw", rownames(plot.df))]<- "PM W Trial"
  plot.df$S[grep ("pmeffectn", rownames(plot.df))]<- "PM N Trial"
  
  
  plot.df$PM <- factor(plot.df$PM)
  plot.df$S <- factor(plot.df$S, levels=c("Non-word Trial", "Word Trial", "PM W Trial", "PM N Trial"))
  
  levels(plot.df$PM)<- c("I - C", "I - U", "U - C", "U - I")
  levels(plot.df$S) <- c("Non-word Trial", "Word Trial", "WPM Trial", "NWPM Trial")

  PLOT.DF <- plot.df
  plot.df <- PLOT.DF[c(1:2,13:16),]
  plot <- ggplot(plot.df, aes(PM, mean)) 
  plot <- plot + facet_grid(. ~ S, scales = "free", space = "free") + geom_point(size=3) + geom_errorbar(aes(ymax = upper, ymin = lower), width= 0.2) + geom_point(aes(PM, data), pch=21, size=4, colour="black")+xlab("Block contrast") +theme(text = element_text(size=24)) +ylab("Accuracy difference") + labs(title = mnam)
  plot+ theme(plot.title = element_text(hjust = 0.5, size=24)) 
}
  
plot.RT.effects.x2 <- function(plot.df, mnam) {
  plot.df <- plot.df[-c(1:4),]
  plot.df$PM <- NA
  plot.df$PM[grep ("imp", rownames(plot.df))] <- "imp - control"
  plot.df$PM[grep ("unimp", rownames(plot.df))] <- "unimp - control"
  plot.df$PM[grep ("pmeffect", rownames(plot.df))] <- "imp- unimp"
  plot.df$PM[grep ("pmrt", rownames(plot.df))] <- "unimp- imp"
  
  plot.df$S <- NA
  plot.df$S[grep ("n", rownames(plot.df))]<- "Non-word Trial"
  plot.df$S[grep ("w", rownames(plot.df))]<- "Word Trial"
  plot.df$S[grep ("pmrtdiffw", rownames(plot.df))]<- "PM W Trial"
  plot.df$S[grep ("pmrtdiffn", rownames(plot.df))]<- "PM N Trial"
  plot.df$S[grep ("pmeffectw", rownames(plot.df))]<- "PM W Trial"
  plot.df$S[grep ("pmeffectn", rownames(plot.df))]<- "PM N Trial"
  
  
  plot.df$PM <- factor(plot.df$PM)
  plot.df$S <- factor(plot.df$S, levels=c("Non-word Trial", "Word Trial", "PM W Trial", "PM N Trial"))
  
  levels(plot.df$PM)<- c("I - C", "I - U", "U - C", "U - I")
  levels(plot.df$S) <- c("Non-word Trial", "Word Trial", "WPM Trial", "NWPM Trial")

  PLOT.DF <- plot.df
  plot.df <- PLOT.DF[3:12,]
  plot.df$R <- NA
  plot.df$R[grep ("P", rownames(plot.df))] <- "PM"
  plot.df$R[grep ("N", rownames(plot.df))] <- "Non-word"
  plot.df$R[grep ("W", rownames(plot.df))] <- "Word"
  plot.df$R[grep ("pmrt", rownames(plot.df))] <- "all"
  plot <- ggplot(plot.df, aes(R, mean)) 
  plot <- plot+ facet_grid(. ~ S + PM, scales = "free", space = "free") + geom_point(size=5) + geom_errorbar(aes(ymax = upper, ymin = lower), width= 0.2) + geom_point(aes(R, data), pch=21, size=7, colour="black") +xlab("Response")+theme(text = element_text(size=24)) +ylab("RT difference")+ labs(title = mnam)
  plot + theme(plot.title = element_text(hjust = 0.5, size=32)) + theme(text = element_text(size = 32))
}

plot.RT.effects.x1 <- function(plot.df, mnam) {
  plot.df <- plot.df[-(1:2),]
  plot.df$PM <- NA
  plot.df$PM[grep ("focal", rownames(plot.df))] <- "F - C"
  plot.df$PM[grep ("nonfocal", rownames(plot.df))] <- "NF - C"
  plot.df$PM[grep ("pmeffect", rownames(plot.df))] <- "F- NF"
  plot.df$PM[grep ("pmrt", rownames(plot.df))] <- "NF - F"
  plot.df$S <- NA
  plot.df$S[grep ("n", rownames(plot.df))]<- "Non-word Trial"
  plot.df$S[grep ("w", rownames(plot.df))]<- "Word Trial"
  plot.df$S[grep ("pm", rownames(plot.df))]<- "PM Trial"
  plot.df$PM <- factor(plot.df$PM)
  plot.df$S <- factor(plot.df$S, levels=c("Non-word Trial", "Word Trial", "PM Trial"))
  PLOT.DF <- plot.df
  plot.df <- PLOT.DF[2:10,]
  plot.df$R <- NA
  plot.df$R[grep ("P", rownames(plot.df))] <- "PM"
  plot.df$R[grep ("N", rownames(plot.df))] <- "Non-word"
  plot.df$R[grep ("W", rownames(plot.df))] <- "Word"
  plot.df$R[grep ("pmrt", rownames(plot.df))] <- "all"
  plot <- ggplot(plot.df, aes(R, mean)) 
  plot<- plot+ facet_grid(. ~ S + PM, scales = "free", space = "free") + geom_point(size=4) + geom_errorbar(aes(ymax = upper, ymin = lower), width= 0.2) + geom_point(aes(R, data), pch=21, size=6, colour="black") +xlab("Response")+theme(text = element_text(size=24)) +ylab("RT difference")+ labs(title = mnam)
  plot + theme(plot.title = element_text(hjust = 0.5, size=26)) + theme(text = element_text(size = 26))
}

plot.acc.effects.x1 <- function (plot.df,mnam) {
  plot.df <- plot.df[-c((1:2),4),]
  plot.df$PM <- NA
  plot.df$PM[grep ("focal", rownames(plot.df))] <- "F - C"
  plot.df$PM[grep ("nonfocal", rownames(plot.df))] <- "NF - C"
  plot.df$PM[grep ("pmeffect", rownames(plot.df))] <- "F - NF"
  plot.df$PM[grep ("pmrt", rownames(plot.df))] <- "NF- F"
  plot.df$PM[grep ("focalPM", rownames(plot.df))] <- "focalPM"
  plot.df$PM[grep ("nonfocalPM", rownames(plot.df))] <- "nonfocalPM"
  

  plot.df$S <- NA
  plot.df$S[grep ("n", rownames(plot.df))]<- "Non-word Trial"
  plot.df$S[grep ("w", rownames(plot.df))]<- "Word Trial"
  plot.df$S[grep ("pm", rownames(plot.df))]<- "PM Trial"
  plot.df$S[grep ("focalPM", rownames(plot.df))] <- "PM Trial"
  plot.df$S[grep ("nonfocalPM", rownames(plot.df))] <- "PM Trial"
  plot.df$PM <- factor(plot.df$PM)
  plot.df$S <- factor(plot.df$S, levels=c("Non-word Trial", "Word Trial", "PM Trial"))
  PLOT.DF <- plot.df
  plot.df <- PLOT.DF[c(1, 10:13),]
  plot <- ggplot(plot.df, aes(PM, mean)) 
  plot <- plot + facet_grid(. ~ S, scales = "free", space = "free") + geom_point(size=3) + geom_errorbar(aes(ymax = upper, ymin = lower), width= 0.2) + geom_point(aes(PM, data), pch=21, size=5, colour="black")+xlab("Block contrast") +theme(text = element_text(size=24)) +ylab("Accuracy difference") + labs(title = mnam)
  plot + theme(plot.title = element_text(hjust = 0.5, size=24))
}


convert.magic <- function(obj,types){
    for (i in 1:length(obj)){
        FUN <- switch(types[i],character = as.character, 
                                   numeric = as.numeric, 
                                   factor = as.factor)
        obj[,i] <- FUN(obj[,i])
    }
    obj
}


tabtoAPA <- function (tab){
  for (i in seq (1, length(tab), 2))  {
   newtab<- paste(tab[,i], " (", tab[,i+1], ")", sep="")
  if (i == 1) { new.tab <- newtab } else {new.tab <- cbind (new.tab, newtab)}
   i = i+1
  }
   rownames(new.tab)<- rownames (tab)
   new.tab
  
}

get.theta.array <- function(hsamples){
  samps <- lapply(hsamples, function(x) x["theta"])
  samps2<- unlist(samps)
  dim3 <- c(dim(samps[[1]]$theta), length(samps2)/prod(dim(samps[[1]]$theta)))
  dim(samps2) <- dim3
  colnames(samps2) <- colnames(hsamples[[1]]$theta)
  samps2
}


GET.fitgglist.dmc <- function (
  PP, factors=NA, noR = FALSE,
  quantiles.to.get = c(0.1, 0.5, 0.9), CI= c(0.025, 0.975),
  acc.fun=function(x){as.numeric(x$S)==as.numeric(x$R)},
  correct.only=FALSE,error.only=FALSE
  
) {
 
  sim <- do.call(rbind, PP)
  # Do the same for the data
  data <- lapply(PP, function(x) attr(x, "data"))
  data <- do.call(rbind, data)
  get.fitgglist.dmc (sim,data, factors=factors, noR=noR, quantiles.to.get=quantiles.to.get,
                     CI=CI, acc.fun=acc.fun, correct.only=correct.only, error.only=
                       error.only)
  
}


get.ns.dmc<- function(samples) {
model <- attributes(samples$data)$model
facs <- names(attr(model,"factors"))
table(samples$data[,facs],dnn=facs)}

ggplot.recov <- function(ggdf, ncol=10) {
  ggplot(ggdf, aes(reps, M))+
    geom_point(size=0.5)+
    geom_hline(aes(yintercept=data), linetype=1, size=1.5)+ylab("")+ 
    geom_ribbon(data=ggdf,aes(ymin=LCI,ymax=HCI), alpha=0.3)+ facet_wrap(~param, ncol=ncol) 
}


get.ggdf.recov <- function(post_summaries, msds, grepchr="B") {
  tmp <- list()
  j=0
  for (i in grep(grepchr, colnames(post_summaries[[1]]))){
    j = j + 1
    tmp [[j]] <- data.frame(cbind(post_summaries[[1]][,i], 
                                   post_summaries[[2]][,i], post_summaries[[3]][,i]))
    name <- rownames(msds)[i]
    colnames(tmp[[j]]) <- c("M", "LCI", "HCI")
    tmp[[j]]$reps <- 1:100
    tmp[[j]]$param <- name
    tmp[[j]]$data <- msds$M[i]
  }
  do.call(rbind, tmp)
}


get.participant.median.CIs <- function(hsamples) {
  samps <- lapply(hsamples, function(x) x["theta"])
  samps2<- unlist(samps)
  dim3 <- c(dim(samps[[1]]$theta), length(samps2)/prod(dim(samps[[1]]$theta)))
  dim(samps2) <- dim3
  post_medians <- apply(samps2, c(4,2), mean)
  post_LCI <- apply(samps2, c(4,2), quantile, prob=0.025)
  post_HCI <- apply(samps2, c(4,2), quantile, prob=0.975)
  colnames(post_medians) <- colnames(hsamples[[1]]$theta)
  colnames(post_LCI) <- colnames(hsamples[[1]]$theta)
  colnames(post_HCI) <- colnames(hsamples[[1]]$theta)
  out <- list(post_medians, post_LCI, post_HCI)
  names(out) <- c("medians", "LCI", "HCI")
  out
}

# A few functions for posterior predicctive p values and z scores.
minp <- function (effect) min(ecdf(effect)(0), 1-ecdf(effect)(0))

zandp <- function(samples, fun){
    effect<- group.inference.dist(samples, fun)
    Z <- mean(effect)/sd(effect)
    p <- minp(effect)
    paste(round(Z,2), " (", round(p,3), ")", sep="")
}

mean.sd <- function(samples, fun){
    effect<- group.inference.dist(samples, fun)
    M <- mean(effect)
    SD <- sd(effect)
    paste(round(M,3), " (", round(SD,3), ")", sep="")
}

  Z.p.acrossexp <- function(samples1,samples2, fun){
    effect1<- group.inference.dist(samples1, fun)
    effect2 <- group.inference.dist(samples2, fun)
    effect<- effect1 - effect2
    Z <- mean(effect)/sd(effect)
    p <- minp(effect)
    paste(round(Z,2), " (", round(p,3), ")", sep="")
  }

##accepts a function and does it to the thetas for each subject then averages after
group.inference.dist <- function (hsamples, fun) {
  inference <- list()
  for (i in 1:length(hsamples)) {
    thetas <- hsamples[[i]]$theta
    inference [[i]] <- fun (thetas)
  }
  inf2 <- unlist(inference)
  dim3 <- c(dim(inference[[1]]), length(inf2)/prod(dim(inference[[1]])))
  dim(inf2) <- dim3
  apply(inf2, c(1,2,3), mean)
}



## Below funciton averages parameters across conditions by grepping out
# from av.posts. So if you set av.posts to match anything containing mean_v,
# it would average all rates together and replace all values with the avg before
# simming. We use it to parse the effects of rates/thresholds on the manifests
# in terms of block and cond.

# samples=samples[[3]]
# probs=c(1:99)/100;random=TRUE
# bw="nrd0";report=10;save.simulation=TRUE;factors=NA
# av.posts<-av.posts.threscond
avps.post.predict.dmc = function(samples,n.post=100,probs=c(1:99)/100,random=TRUE,
                                    bw="nrd0",report=10,save.simulation=TRUE,factors=NA, av.posts=c())
  # make list of posterior preditive density, quantiles and response p(robability)
{


  get.dqp <- function(sim,facs,probs,n.post=NA) {

    quantile.names <- function(x,probs=seq(0, 1, 0.25),na.rm=FALSE,type=7, ...) {
      out <- quantile(x,probs=probs,na.rm=na.rm,type=type,names=FALSE,...)
      names(out) <- probs*100
      if (all(is.na(out))) NULL else out
    }

    qs <- tapply(sim$RT,sim[,c(facs,"R")],quantile.names,probs=probs,na.rm=TRUE)
    #     qs <- apply(qs,1:length(dim(qs)),function(x){
    #       if ( is.null(x[[1]]) || all(is.na(x[[1]]))) NULL else x[[1]]})

    # get probabilities given not na
    simOK <- sim[!is.na(sim$RT),]
    p <- tapply(simOK$RT,simOK[,c(facs,"R")],length)
    p[is.na(p)] <- 0 # In case some cells are empty
    np <- rep(apply(p,1:length(facs),sum),times=length(levels(simOK$R)))
    p <- p/np
    # get p NA
    pNA <- tapply(sim$RT,sim[,c(facs,"R")],length)
    pNA[is.na(pNA)] <- 0 # In case some cells are empty
    npNA <- rep(apply(pNA,1:length(facs),sum),times=length(levels(sim$R)))
    pNA <- tapply(is.na(sim$RT),sim[,c(facs,"R")],sum)
    pNA[is.na(pNA)] <- 0 # In case some cells are empty
    pNA <- pNA/npNA

    # For a simulation get proability replicates
    if (!is.na(n.post)) {
      repfac <- rep(1:n.post,each=sum(ns))
      repfac <- repfac[!is.na(sim$RT)]
      ps <- tapply(simOK$RT,cbind(simOK[,c(facs,"R")],rep=repfac),length)
      ps[is.na(ps)] <- 0 # In case some cells are empty
      ps <- n.post*ps/np
      # and NA replicates
      pNAs <- tapply(is.na(sim$RT),cbind(sim[,c(facs,"R")],rep=rep(1:n.post,each=sum(ns))),sum)
      pNAs[is.na(pNAs)] <- 0 # In case some cells are empty
      pNAs <- n.post*pNAs/npNA
    } else {
      ps=NULL
      pNAs=NULL
    }

    # cell names
    cell.names <- dimnames(qs)[[1]]
    n.cell <- length(facs)+1
    for ( i in 2:n.cell )
      cell.names <- outer(cell.names,dimnames(qs)[[i]],"paste",sep=".")
    cell.names <- as.vector(cell.names)
    # Get density and make defective
    dens <- tapply(sim$RT,sim[,c(facs,"R")],function(x){
      if (all(is.na(x))) NULL else {
        x <- x[x>=quantile(x,.01,na.rm=TRUE) & x<=quantile(x,.99,na.rm=TRUE)]
        if (length(x[!is.na(x)])<2) NULL else density(x[!is.na(x)],bw=bw)
      }
    })
    for (i in 1:length(p)) if ( is.finite(p[i]) && !(p[i]==0) )
    {
      if (!is.null(qs[i][[1]])) {
        names(qs[i][[1]]) <- as.numeric(names(qs[i][[1]]))*p[i]/100
        attr(qs[i][[1]],"cell.name") <- cell.names[i]
      }
      if (!is.null(dens[i][[1]]) ) {
        dens[i][[1]]$y <- dens[i][[1]]$y*p[i]
        attr(dens[i][[1]],"cell.name") <- cell.names[i]
      }
    }
    dnd <- dimnames(dens)
    dnq <- dimnames(qs)
    dens <- apply(dens,1:length(facs),function(x){x})
    qs <- apply(qs,1:length(facs),function(x){x})
    if ( is.null(dim(dens)) ) {
      dens <- array(dens,dim=c(length(dens)))
      dimnames(dens) <- dnd[-length(dnd)]
      qs <- array(qs,dim=c(length(qs)))
      dimnames(qs) <- dnq[-length(dnq)]
    }
    list(pdf=dens,cdf=qs,p=p,ps=ps,pNA=pNA,pNAs=pNAs)
  }

  model <- attributes(samples$data)$model
  facs <- names(attr(model,"factors"))
  if (any(is.na(factors))) factors <- facs
  if (!all(factors %in% facs))
    stop(paste("Factors argument must contain one or more of:",paste(facs,collapse=",")))
  resp <- names(attr(model,"responses"))
  ns <- table(samples$data[,facs],dnn=facs)
  n.par <- dim(samples$theta)[2]
  thetas <- matrix(aperm(samples$theta,c(3,1,2)),ncol=n.par)
  colnames(thetas) <- dimnames(samples$theta)[[2]]
  if (is.na(n.post)) use <- c(1:dim(thetas)[1]) else {
    if (random) use <- sample(c(1:dim(thetas)[1]),n.post,replace=F) else
      use <- round(seq(1,dim(thetas)[1],length.out=n.post))
  }
  n.post <- length(use)
  posts <- thetas[use,]


  cat("Below is how I'm averaging (each row is averaged). If this is wrong, adjust your
      av.posts to grep correctly.")
  for (q in 1:length(av.posts)) print(colnames(posts[, grep(av.posts[q], colnames(posts))]))
  ###tweak to average av.posts
  q=1

  if(length(av.posts)!= 0) {
    ### loop through all the averaged posts
    for (q in 1:length(av.posts)) {

      num.params <- dim(posts[, grep(av.posts[q], colnames(posts))])[2]
      average.params <- rowMeans(posts[, grep(av.posts[q], colnames(posts))])
      posts[, grep(av.posts[q], colnames(posts))] <- matrix(average.params,nrow=length(average.params),ncol=num.params,byrow=F)

    }
  }

  ########
  n.rep <- sum(ns)
  sim <- data.frame(matrix(nrow=n.post*n.rep,ncol=dim(samples$data)[2]))
  names(sim) <- names(samples$data)
  # Tweaks for Stop Signal
  if ( !any(names(samples$data)=="SSD") ) {
    SSD <- rep(Inf,sum(ns))
    leave.out <- -c(1:dim(samples$data)[2])[names(samples$data)=="RT"]
  } else {
    # Assumes last two are SSD and RT! FIX ME.
    if (is.null(facs)) SSD <- samples$data$SSD else
      SSD <- unlist(tapply(samples$data$SSD,samples$data[,facs],identity))
    leave.out <- -c((dim(samples$data)[2]-1):dim(samples$data)[2])
  }
  cat("\n")
  cat(paste("Simulating (\'.\'=",report,"): ",sep=""))
  for (i in names(samples$data)[leave.out])
    sim[[i]] <- factor(rep(NA,n.post*n.rep),levels=levels(samples$data[[i]]))
  for (i in 1:n.post) {
    tmp <- simulate.dmc(posts[i,],model,ns,SSD=SSD)
    sim[(1+(i-1)*n.rep):(i*n.rep),names(tmp)] <- tmp
    if ( (i %% report) == 0) cat(".")
  }
  if ( any(names(sim)=="R2") ) { # MTR model
    sim$R <- paste(as.character(sim$R),as.character(sim$R2),sep="")
    sim$R[sim$R2=="DK"] <- "DK"
    sim$R <- factor(sim$R)
    samples$data$R <- paste(as.character(samples$data$R),as.character(samples$data$R2),sep="")
    samples$data$R[samples$data$R2=="DK"] <- "DK"
    samples$data$R <- factor(samples$data$R)
  }
  if (save.simulation) {
    sim <- cbind(reps=rep(1:n.post,each=dim(samples$data)[1]),sim)
    attr(sim,"data") <- samples$data
    sim
  } else {
    sim.dqp <- get.dqp(sim,factors,probs,n.post)
    dat.dqp <- get.dqp(sim=samples$data,factors,probs)
    names(dat.dqp) <- paste("data",names(dat.dqp),sep=".")
    c(sim.dqp,dat.dqp)
  }
}


#lapplys the above function on everybody
avps.h.post.predict.dmc<- function(samples,n.post=100,probs=c(1:99)/100,
                                          bw="nrd0",
                                     save.simulation=FALSE, av.posts=c())
  # apply lost.predict to each subject
{
  lapply(samples,avps.post.predict.dmc,n.post=n.post,probs=probs,bw=bw,
         save.simulation=save.simulation, av.posts=av.posts)
}

# PPs = E1PP
# fun = block.effects.E1A4
# lower=.025
# upper=.975
get.effects.dmc <- function (PPs, fun = function (x) {mean (x)}, lower=.025, upper=.975) {

  simdata<- do.call(rbind, PPs)
  data <- lapply(PPs, function(x) attr(x, "data"))
  data <- do.call(rbind, data)
  nreps=max(PPs[[1]]$reps)


  data.effects<- fun(data)
  noutput <- length(data.effects)

  sim.effects <- matrix(NA, nrow= nreps, ncol=noutput+1)
  sim.effects[,noutput+1] <- 1:nreps

  colnames(sim.effects) <- c(names(data.effects), "n.rep")
  ######

  ##calculate effects separately for each rep
  for (j in 1:nreps) {

    currentsim.effects <- simdata[simdata$reps==j,]
    sim.effects[j,1:noutput] <- fun(currentsim.effects)

  }

  ##Get a ggplot df with posterior mean, lower, and upper.
  effects.ggdf <-  t(apply(sim.effects, c(2), function(x) c(mean(x),
                                                quantile(x, probs=c(lower,upper)))))
  effects.ggdf <- data.frame(effects.ggdf)
  effects.ggdf <- effects.ggdf[(!rownames(effects.ggdf) %in% "n.rep"),]
  colnames(effects.ggdf) <- c("mean", "lower", "upper")
  contrast <- rownames(effects.ggdf)
  effects.ggdf$data<-as.vector(data.effects)
  attr(effects.ggdf, "post.effects.samples") <- sim.effects
  effects.ggdf
}



#The below function picks certain parameters from a list called pickps_set,
# and replaces them with pickps_other, before performing posterior prediciton.
# We use it to turn control mechanisms off in the model. To turn off proactive
# control, we set the ongoing task thresholds equal in the PM block to the control
#threshols. To turn off reactive, we set the ongiong rates on PM trials (PM block)
# to the ongoing rates on non-PM trials (PM block)

# samples=samples[[3]]
# probs=c(1:99)/100;random=TRUE
# bw="nrd0";report=10;save.simulation=TRUE;factors=NA; n.post=100
# pickps_set <- c("B.A2C",        "B.B2C",        "B.C2C" ,
#                 "B.D2C"  ,             "B.A2N"  ,      "B.B2N"  ,      "B.C2N"  ,
#                 "B.D2N")
#
# pickps_others <- c("B.A3C"   ,     "B.B3C",        "B.C3C" ,
#                   "B.D3C",        "B.A3N"   ,     "B.B3N" ,       "B.C3N" ,
#                   "B.D3N" )
#
#


pickps.post.predict.dmc = function(samples,n.post=100,probs=c(1:99)/100,random=TRUE,
                                   bw="nrd0",report=10,save.simulation=TRUE,factors=NA, pickps_others, pickps_set)
  # make list of posterior preditive density, quantiles and response p(robability)
{
  
  
  get.dqp <- function(sim,facs,probs,n.post=NA) {
    
    quantile.names <- function(x,probs=seq(0, 1, 0.25),na.rm=FALSE,type=7, ...) {
      out <- quantile(x,probs=probs,na.rm=na.rm,type=type,names=FALSE,...)
      names(out) <- probs*100
      if (all(is.na(out))) NULL else out
    }
    
    qs <- tapply(sim$RT,sim[,c(facs,"R")],quantile.names,probs=probs,na.rm=TRUE)
    #     qs <- apply(qs,1:length(dim(qs)),function(x){
    #       if ( is.null(x[[1]]) || all(is.na(x[[1]]))) NULL else x[[1]]})
    
    # get probabilities given not na
    simOK <- sim[!is.na(sim$RT),]
    p <- tapply(simOK$RT,simOK[,c(facs,"R")],length)
    p[is.na(p)] <- 0 # In case some cells are empty
    np <- rep(apply(p,1:length(facs),sum),times=length(levels(simOK$R)))
    p <- p/np
    # get p NA
    pNA <- tapply(sim$RT,sim[,c(facs,"R")],length)
    pNA[is.na(pNA)] <- 0 # In case some cells are empty
    npNA <- rep(apply(pNA,1:length(facs),sum),times=length(levels(sim$R)))
    pNA <- tapply(is.na(sim$RT),sim[,c(facs,"R")],sum)
    pNA[is.na(pNA)] <- 0 # In case some cells are empty
    pNA <- pNA/npNA
    
    # For a simulation get proability replicates
    if (!is.na(n.post)) {
      repfac <- rep(1:n.post,each=sum(ns))
      repfac <- repfac[!is.na(sim$RT)]
      ps <- tapply(simOK$RT,cbind(simOK[,c(facs,"R")],rep=repfac),length)
      ps[is.na(ps)] <- 0 # In case some cells are empty
      ps <- n.post*ps/np
      # and NA replicates
      pNAs <- tapply(is.na(sim$RT),cbind(sim[,c(facs,"R")],rep=rep(1:n.post,each=sum(ns))),sum)
      pNAs[is.na(pNAs)] <- 0 # In case some cells are empty
      pNAs <- n.post*pNAs/npNA
    } else {
      ps=NULL
      pNAs=NULL
    }
    
    # cell names
    cell.names <- dimnames(qs)[[1]]
    n.cell <- length(facs)+1
    for ( i in 2:n.cell )
      cell.names <- outer(cell.names,dimnames(qs)[[i]],"paste",sep=".")
    cell.names <- as.vector(cell.names)
    # Get density and make defective
    dens <- tapply(sim$RT,sim[,c(facs,"R")],function(x){
      if (all(is.na(x))) NULL else {
        x <- x[x>=quantile(x,.01,na.rm=TRUE) & x<=quantile(x,.99,na.rm=TRUE)]
        if (length(x[!is.na(x)])<2) NULL else density(x[!is.na(x)],bw=bw)
      }
    })
    for (i in 1:length(p)) if ( is.finite(p[i]) && !(p[i]==0) )
    {
      if (!is.null(qs[i][[1]])) {
        names(qs[i][[1]]) <- as.numeric(names(qs[i][[1]]))*p[i]/100
        attr(qs[i][[1]],"cell.name") <- cell.names[i]
      }
      if (!is.null(dens[i][[1]]) ) {
        dens[i][[1]]$y <- dens[i][[1]]$y*p[i]
        attr(dens[i][[1]],"cell.name") <- cell.names[i]
      }
    }
    dnd <- dimnames(dens)
    dnq <- dimnames(qs)
    dens <- apply(dens,1:length(facs),function(x){x})
    qs <- apply(qs,1:length(facs),function(x){x})
    if ( is.null(dim(dens)) ) {
      dens <- array(dens,dim=c(length(dens)))
      dimnames(dens) <- dnd[-length(dnd)]
      qs <- array(qs,dim=c(length(qs)))
      dimnames(qs) <- dnq[-length(dnq)]
    }
    list(pdf=dens,cdf=qs,p=p,ps=ps,pNA=pNA,pNAs=pNAs)
  }
  
  model <- attributes(samples$data)$model
  facs <- names(attr(model,"factors"))
  if (any(is.na(factors))) factors <- facs
  if (!all(factors %in% facs))
    stop(paste("Factors argument must contain one or more of:",paste(facs,collapse=",")))
  resp <- names(attr(model,"responses"))
  ns <- table(samples$data[,facs],dnn=facs)
  n.par <- dim(samples$theta)[2]
  thetas <- matrix(aperm(samples$theta,c(3,1,2)),ncol=n.par)
  colnames(thetas) <- dimnames(samples$theta)[[2]]
  if (is.na(n.post)) use <- c(1:dim(thetas)[1]) else {
    if (random) use <- sample(c(1:dim(thetas)[1]),n.post,replace=F) else
      use <- round(seq(1,dim(thetas)[1],length.out=n.post))
  }
  n.post <- length(use)
  posts <- thetas[use,]
  
  #more robust
  
  
  ###Replace some parameter vlaues with others.
  posts[,colnames(posts) %in% pickps_others][,pickps_others] <- 
    posts[,colnames(posts) %in% pickps_set][,pickps_set] 
  
  ########
  n.rep <- sum(ns)
  sim <- data.frame(matrix(nrow=n.post*n.rep,ncol=dim(samples$data)[2]))
  names(sim) <- names(samples$data)
  # Tweaks for Stop Signal
  if ( !any(names(samples$data)=="SSD") ) {
    SSD <- rep(Inf,sum(ns))
    leave.out <- -c(1:dim(samples$data)[2])[names(samples$data)=="RT"]
  } else {
    # Assumes last two are SSD and RT! FIX ME.
    if (is.null(facs)) SSD <- samples$data$SSD else
      SSD <- unlist(tapply(samples$data$SSD,samples$data[,facs],identity))
    leave.out <- -c((dim(samples$data)[2]-1):dim(samples$data)[2])
  }
  cat("\n")
  cat(paste("Simulating (\'.\'=",report,"): ",sep=""))
  for (i in names(samples$data)[leave.out])
    sim[[i]] <- factor(rep(NA,n.post*n.rep),levels=levels(samples$data[[i]]))
  for (i in 1:n.post) {
    tmp <- simulate.dmc(posts[i,],model,ns,SSD=SSD)
    sim[(1+(i-1)*n.rep):(i*n.rep),names(tmp)] <- tmp
    if ( (i %% report) == 0) cat(".")
  }
  if ( any(names(sim)=="R2") ) { # MTR model
    sim$R <- paste(as.character(sim$R),as.character(sim$R2),sep="")
    sim$R[sim$R2=="DK"] <- "DK"
    sim$R <- factor(sim$R)
    samples$data$R <- paste(as.character(samples$data$R),as.character(samples$data$R2),sep="")
    samples$data$R[samples$data$R2=="DK"] <- "DK"
    samples$data$R <- factor(samples$data$R)
  }
  if (save.simulation) {
    sim <- cbind(reps=rep(1:n.post,each=dim(samples$data)[1]),sim)
    attr(sim,"data") <- samples$data
    sim
  } else {
    sim.dqp <- get.dqp(sim,factors,probs,n.post)
    dat.dqp <- get.dqp(sim=samples$data,factors,probs)
    names(dat.dqp) <- paste("data",names(dat.dqp),sep=".")
    c(sim.dqp,dat.dqp)
  }
}

#lapply the above to the whole samples object.
pickps.h.post.predict.dmc<- function(samples,n.post=100,probs=c(1:99)/100,
                                   bw="nrd0",
                                   save.simulation=FALSE, pickps_set, pickps_others)
  # apply lost.predict to each subject
{
  lapply(samples,pickps.post.predict.dmc,n.post=n.post,probs=probs,bw=bw,
         save.simulation=save.simulation, pickps_set=pickps_set, pickps_others=pickps_others)
}



get.hdata.dmc <- function(hsamples){
list.wind<-lapply(seq_along(hsamples), function(samples, n, i) cbind(n[[i]], samples[[i]]$data), 
                  samples= hsamples, n = names(hsamples))
out<-do.call(rbind, list.wind)
names(out)[1] <- "s"
out
}

fixedeffects.meanthetas <- function (samples){
samps <- lapply(samples, function(x) x["theta"])
##
## thetas into big array for apply
samps2<- unlist(samps)
dim3 <- c(dim(samps[[1]]$theta), length(samps2)/prod(dim(samps[[1]]$theta)))
dim(samps2) <- dim3
samps3<- apply(samps2, c(1,2,3), mean)
## back to a theta list after applied
colnames(samps3)<- colnames(samps[[1]]$theta)
samps5<- list(samps3)
attributes(samps5)<- attributes(samps[[1]])
samps5
}

get.pinf.subjects <- function(funs=list(mean), hsamples, eff.names= c()) {
  inference<- list()
   for (i in 1:length(hsamples)) {
     thetas <- hsamples[[i]]$theta
     effects<-list()
     for(j in 1:length(funs)){
       effects [[j]] <- funs[[j]](thetas)
       # names(effects) <- eff.names
     }
     inference[[i]] <- effects
     names(inference[[i]]) <- eff.names
   } 
   final.effects <- list()
   for (k in 1:length(funs)) {
     this.inf <- lapply(inference, function(x) x[[k]])
     inf2 <- unlist(this.inf)
     dim3 <- c(dim(this.inf[[1]]), length(inf2)/prod(dim(this.inf[[1]])))
     dim(inf2) <- dim3
     final.effects[[k]] <- inf2
   }
   final.effects
}


###########################Air Traffic Control Project#############################


#- Accounting for non-responses.
#The below function, post.predict.dmc.MATCHORDER, accepts a samples object and
#an okdats object. The latter must have the FULL original data frame for each
#participant, before cleaning and including non-responses.
# The function simulates the same amount of data as in the data frame with the same
#design. Then it arranges the orders of trials to correspond to the data.
# This will allow subsequent functions both to truncate simulations that would
#be non-responses, and to calculate the non-response rate.

# samples=samples.E1[[1]];n.post=100;probs=c(1:99)/100;random=TRUE
# bw="nrd0";report=10;save.simulation=FALSE;factors=NA
# save.simulation.as.attribute=FALSE;ignore.R2=FALSE
# gglist=FALSE; probs.gglist=c(0.1, 0.5, 0.9);CI.gglist=c(0.025, 0.975)
#assumes samples have an attribute called NRdata. This attribute is the unfiltered data for each
#sample i.e. before we removed non-responses, cleaned etc.
post.predict.dmc.MATCHORDER <- function(samples,n.post=100,probs=c(1:99)/100,random=TRUE,
         bw="nrd0",report=10,save.simulation=FALSE,factors=NA,
         save.simulation.as.attribute=FALSE,ignore.R2=FALSE,
         gglist=FALSE, probs.gglist=c(0.1, 0.5, 0.9),CI.gglist=c(0.025, 0.975))
  # make list of posterior preditive density, quantiles and response p(robability)
  # NB: quantiles only calcualted for 2 or more RTs
{

  model <- attributes(samples$data)$model
  facs <- names(attr(model,"factors"))
  cvs <- samples$data[,attr(model,"cvs")]
  attr(cvs,"row.facs") <- apply(apply(
    samples$data[,facs,drop=FALSE],2,as.character),1,paste,collapse=".")
  if ( ignore.R2 & any(names(samples$data)=="R2") )
    samples$data <- samples$data[,names(samples$data)[names(samples$data)!="R2"]]
  if (!is.null(factors) ) {
    if (any(is.na(factors))) factors <- facs
    if (!all(factors %in% facs))
      stop(paste("Factors argument must contain one or more of:",paste(facs,collapse=",")))
  }
  resp <- names(attr(model,"responses"))
  ##LUKE: plug in data size from okdats.noNR
  data <- attr(samples, "NRdata")
  ns <- table(data[,facs])
  # ns <- table(samples$data[,facs],dnn=facs)
  n.par <- dim(samples$theta)[2]
  thetas <- matrix(aperm(samples$theta,c(3,1,2)),ncol=n.par)
  colnames(thetas) <- dimnames(samples$theta)[[2]]
  if (is.na(n.post)) use <- c(1:dim(thetas)[1]) else {
    if (random) use <- sample(c(1:dim(thetas)[1]),n.post,replace=F) else
      use <- round(seq(1,dim(thetas)[1],length.out=n.post))
  }
  n.post <- length(use)
  posts <- thetas[use,]
  n.rep <- sum(ns)
  sim <- data.frame(matrix(nrow=n.post*n.rep,ncol=dim(data)[2]))
  names(sim) <- names(data)
  # Tweaks for Stop Signal
  if ( !any(names(samples$data)=="SSD") ) {
    SSD <- rep(Inf,sum(ns))
    leave.out <- -c(1:dim(samples$data)[2])[names(samples$data) %in% c("RT",names(cvs),"R2")]
  } else {
    # Assumes last two are SSD and RT! FIX ME. EG WONT WORK IF THERE ARE CVS
    if (is.null(facs)) SSD <- samples$data$SSD else
      SSD <- unlist(tapply(samples$data$SSD,samples$data[,facs],identity))
    leave.out <- -c((dim(samples$data)[2]-1):dim(samples$data)[2])
  }
  cat("\n")
  cat(paste("Simulating (\'.\'=",report,"): ",sep=""))
  for (i in names(samples$data)[leave.out])
    sim[[i]] <- factor(rep(NA,n.post*n.rep),levels=levels(samples$data[[i]]))
  for (i in 1:n.post) {
    tmp <- simulate.dmc(posts[i,],model,n=ns,SSD=SSD,cvs=NULL)

    if ( (i %% report) == 0) cat(".")
    ###Luke plug in order matching

    #getting the factor structure of whatever design you feed in
    callargs.data <- list()
    callargs.sim <- list()
    for (p in 1:length(facs)) {callargs.data[[p]] <- data[,facs[p]]
    callargs.sim[[p]] <- tmp[,facs[p]]
    }

    data.ind <- factor(do.call(paste, callargs.data))
    sim.ind <- factor(do.call(paste, callargs.sim))
    swappedsim <- data
    for(q in levels(data.ind)){swappedsim$R[data.ind==q] <- tmp$R[sim.ind==q]
    swappedsim$RT[data.ind==q] <- tmp$RT[sim.ind==q]
    }
    tmp <- swappedsim

    if (ignore.R2) tmp <- tmp[,names(tmp)[names(tmp)!="R2"]]
    sim[(1+(i-1)*n.rep):(i*n.rep),names(tmp)] <- tmp

    ########
  }
  if ( any(names(sim)=="R2") ) { # MTR model
    sim$R <- paste(as.character(sim$R),as.character(sim$R2),sep="")
    sim$R[sim$R2=="DK"] <- "DK"
    sim$R <- factor(sim$R)
    samples$data$R <- paste(as.character(samples$data$R),as.character(samples$data$R2),sep="")
    samples$data$R[samples$data$R2=="DK"] <- "DK"
    samples$data$R <- factor(samples$data$R)
  }
  reps <- rep(1:n.post,each=dim(data)[1])
  if ( save.simulation ) {
    sim <- cbind(reps,sim)
    attr(sim,"data") <- data
    # attr(sim,"data") <- samples$data
    sim
  } else {
    sim.dqp <- get.dqp(sim,facs=factors,probs,n.post,ns=ns,bw=bw)
    dat.dqp <- get.dqp(sim=samples$data,factors,probs,bw=bw)
    names(dat.dqp) <- paste("data",names(dat.dqp),sep=".")
    out <- c(sim.dqp,dat.dqp)
    dpqs <- vector(mode="list",length=length(n.post))
    for (i in 1:n.post) {
      simi <- sim[reps==i,]
      dpqs[[i]] <- get.dqp(simi,factors,probs,1)
    }
    attr(out,"dpqs") <- dpqs
    if (save.simulation.as.attribute)
      attr(out,"sim") <- cbind(reps,sim)
    if (gglist) attr(out, "gglist") <-
      get.fitgglist.dmc(cbind(reps,sim),samples$data,factors=factors, noR=FALSE,
                        quantiles.to.get= probs.gglist, CI = CI.gglist)
    out
  }
}

#this function merely lapplys the function above.
h.post.predict.dmc.MATCHORDER <- function(hsamples, n.post=100) {
  lapply(hsamples, post.predict.dmc.MATCHORDER, save.simulation=TRUE, n.post=n.post)
}

# Gets the cumulative sum of RTs for data.
# assigns trials to the data
#with a loop.  Note this trial index is used for the sim as well.
add.trial.cumsum.data <- function(df) {
  df$trial <- NA
  df$trial.pos <-as.numeric(df$trial.pos)

  df$trial <- NA
  df$trial.pos <-as.numeric(df$trial.pos)
  g=1
  for (t in 1:length(df$RT))  {
    if (t==1) df$trial[t] <- 1 else if (df$trial.pos[t] ==  df$trial.pos[t-1] +1) df$trial[t] <- g
    else {
      g <- g+1
      df$trial[t] <- g}
  }
  df<-cbind(df,unlist(by(df$RT,df$trial,cumsum)))
  names(df)[length(df)]<- "cumsum"
  df
}



# library(car)
# library(heplots)

wsAnova=function(dat,SStype=3,spss=F) {
  has.car=require(car)  
  if (!has.car) return("No \"car\"package, no ANOVA\n") 
  for (i in 1:(dim(dat)[2]-1)) dat[,i] <- factor(dat[,i])
  dat=dat[do.call(order,dat[,-dim(dat)[2]]),]
  snams=levels(dat[,1]); ns=length(snams)
  dvnam=names(dat)[dim(dat)[2]]
  facn=names(dat)[-c(1,dim(dat)[2])]
  nifac=length(facn)
  idata=data.frame(dat[dat[,1]==snams[1],facn])
  names(idata)=facn
  for (i in facn) 
    if (i==facn[1]) ifacn=as.character(idata[,1]) else 
      ifacn=paste(ifacn,as.character(idata[,i]),sep=".")
  facnr=facn[nifac:1]
  e.mv=matrix(unlist(tapply(dat[,dvnam],dat[,facnr],function(x){x})),
              ncol=length(ifacn),dimnames=list(snams,ifacn))
  print(summary(Anova(lm(e.mv ~ 1),
                idata=idata,type=SStype, 
                idesign=formula(paste("~",paste(facn,collapse="*")))),
          multivariate=FALSE))
#LUKE PARTIAL ETA SQUARED  
  SSH <- (summary(Anova(lm(e.mv ~ 1),
                        idata=idata,type=SStype, 
                        idesign=formula(paste("~",paste(facn,collapse="*")))),
                  multivariate=FALSE))$univariate.tests[,1] 
  
  SST <- SSH + (summary(Anova(lm(e.mv ~ 1),
                              idata=idata,type=SStype, 
                              idesign=formula(paste("~",paste(facn,collapse="*")))),
                        multivariate=FALSE))$univariate.tests[,3] 
  print("Partial Eta Sq:")
  print((SSH/SST)[-1])
###  
  if (spss) {
    e.mv=cbind.data.frame(s=row.names(e.mv),e.mv)
    row.names(e.mv)=NULL
    e.mv
  }
}




mneffects=function(df,elist,digits=3,err=F,vars=F,dvnam="y") {
  df <- df[,c(names(df)[names(df)!=dvnam],dvnam)]
  dvnam=dim(df)[2]
  for (i in 1:(dim(df)[2]-1)) df[,i] <- factor(df[,i])  
  for (i in 1:length(elist)) {
    cat(paste(paste(elist[[i]],collapse=":"),"\n"))
    mns=tapply(df[,dvnam],df[,elist[[i]]],mean,na.rm=T)
    if (err) print(round(plogis(mns),digits)) else
      if (vars) print(round(sqrt(mns),digits))  else
        print(round(mns,digits))    
    cat("\n")
  }  
}


se=function(df,facs,sfac="s",dvnam="y",ws=TRUE,ci="SE") {
  df <- df[,c(names(df)[names(df)!=dvnam],dvnam)]
  dvnam=dim(df)[2]
  for (i in 1:(dim(df)[2]-1)) df[,i] <- factor(df[,i])  
  if (ws) {
    smns <- tapply(df[,dvnam],df[,sfac],mean)
    smn <- df[,sfac]
    levels(smn) <- smns
    df[,dvnam] <- df[,dvnam]-as.numeric(as.character(smn))  
  }
  mn=tapply(df[,dvnam],df[,facs],mean)
  se=tapply(df[,dvnam],df[,facs],sd)
  ns <- length(levels(df[,sfac]))
  if (ws) {
    m <- prod(dim(se))
    ns <- ns*(m-1)/m
  }
  if (is.na(ci)) mn else {
    if (ci=="SE") se/sqrt(ns) else
     qt(1-(100-ci)/200,ns-1)*se/sqrt(ns)
  }
}

se2=function(df,facs,sfac="s",dvnam="y",ws=TRUE,ci="SE") {
  df <- df[,c(names(df)[names(df)!=dvnam],dvnam)]
  dvnam=dim(df)[2]
  for (i in 1:(dim(df)[2]-1)) df[,i] <- factor(df[,i])  
  if (ws) {
    smns <- tapply(df[,dvnam],df[,sfac],mean, na.rm=T)
    smn <- df[,sfac]
    levels(smn) <- smns
    df[,dvnam] <- df[,dvnam]-as.numeric(as.character(smn))  
  }
  mn=tapply(df[,dvnam],df[,facs],mean, na.rm=T)
  se=tapply(df[,dvnam],df[,facs],sd, na.rm=T)
  ns <- length(levels(df[,sfac]))
  if (ws) {
    m <- prod(dim(se))
    ns <- ns*(m-1)/m
  }
  if (is.na(ci)) mn else {
    if (ci=="SE") se/sqrt(ns) else
     qt(1-(100-ci)/200,ns-1)*se/sqrt(ns)
  }
}

add.bars=function(mn,se,xvals=NA,len=.1,antiprobit=FALSE,col="black") {
  
  plotbars <- function(x,m,l,h,len,col="black") {
    for (j in 1:length(x)) arrows(x[j],m[j],x[j],l[j],length=len,angle=90,col=col)
    for (j in 1:length(x)) arrows(x[j],m[j],x[j],h[j],length=len,angle=90,col=col)    
  }
  
  if (any(is.na(xvals))) if (is.matrix(mn)) 
    xvals <- as.numeric(dimnames(mn)[[2]]) else
    xvals <- as.numeric(factor(names(mn)))
  lo <- mn-se
  hi <- mn+se
  if (antiprobit) {
    mn=pnorm(mn)
    lo=pnorm(lo)
    hi=pnorm(hi)
  }
  if (!is.matrix(mn)) 
      plotbars(xvals,mn,lo,hi,col=col,len=len) else
    for (i in 1:dim(mn)[1]) 
      plotbars(x=xvals,m=mn[i,],l=lo[i,],h=hi[i,],len=len,col=col)
}    

arr2df=function(arr) {
  if (is.null(dim(arr))) out=data.frame(y=arr) else {
    dn=dimnames(arr)
    if (length(dn)==1) {
      out=cbind.data.frame(factor(dn[[1]],dn[[1]]),arr)
      names(out)=c(names(dn),"y")
      row.names(out)=NULL
    } else {
      tmp=vector(mode="list",length=length(dn))
      names(tmp)=names(dn)
      k=1
      for (j in names(dn)) {
        n=length(dn[[j]])
        tmp[[j]]=gl(n,k,length(arr),dn[[j]])
        k=k*n
      }
      out=cbind(data.frame(tmp),y=as.vector(arr))
      row.names(out)=NULL
    }
  }
  out
}

# spss=F; dvnam=NA; SStype=3; sfac="s"
mixedAnova=function(dat,bsfacn,wsfacn=NULL,sfac="s",SStype=3,spss=F,dvnam=NA) {
  has.car=require(car)  
  if (!has.car) return("No \"car\"package, no ANOVA\n") 
  if (is.na(dvnam)) dvnam <- names(dat)[dim(dat)[2]]
  dat <- dat[,c(sfac,bsfacn,wsfacn,dvnam)]
  if (length(wsfacn)>0) dat <- dat[do.call(order,dat[,c(sfac,wsfacn)]),]
  for (i in 1:(dim(dat)[2]-1)) dat[,i] <- factor(dat[,i])
  snams=levels(dat[,sfac]); ns=length(snams)
  nifac=length(wsfacn)
  lev1s=unlist(lapply(lapply(dat[,wsfacn],levels),function(x){x[1]}))
  bsfacs=dat[apply(dat[,wsfacn,drop=F],1,function(x){all(x==lev1s)}),bsfacn,drop=F]  
  if ( nifac>0 ) {
    idata=data.frame(dat[dat[,sfac]==snams[1],wsfacn])
    names(idata)=wsfacn
    for (i in wsfacn) 
      if (i==wsfacn[1]) ifacn=as.character(idata[,1]) else 
        ifacn=paste(ifacn,as.character(idata[,i]),sep=".")
    e.mv=matrix(unlist(
      tapply(dat[,dvnam],dat[,wsfacn[length(wsfacn):1]],function(x){x})),
                ncol=length(ifacn),dimnames=list(snams,ifacn))
    summary(Anova(
      lm(formula(paste("e.mv ~",paste(bsfacn,collapse="*"))),bsfacs),
      idata=idata,type=SStype,
      idesign=formula(paste("~",paste(wsfacn,collapse="*")))),multivariate=F)
  } else {
    e.mv <- cbind(y=dat[,dvnam]) 
    print(Anova(lm(formula(paste("e.mv ~",paste(bsfacn,collapse="*"))),
                   bsfacs),type=3))
  }
  if (spss) {
    e.mv=cbind.data.frame(s=row.names(e.mv),bsfacs,e.mv)
    row.names(e.mv)=NULL
    e.mv
  }
}

# scol = subject column (numeric or name, same for following)
# dvcol = dependent variable column
# ws and bs name corresponding factors in dat
# mixedAnova=function(dat,ws,bs=NA,scol=1,dvcol=dim(dat)[2],
#   SStype=3,icontrasts=c("contr.sum", "contr.poly"),
#   savedf=F) {
#   dat=dat[do.call(order,dat[,ws,drop=F]),]
#   snams=levels(dat[,scol]); ns=length(snams)
#   dvnam=names(dat)[dvcol]
#   nifac=length(ws)
#   idata=data.frame(dat[dat[,scol]==snams[1],ws])
#   names(idata)=ws  
#   for (i in ws) 
#     if (i==ws[1]) ifacn=as.character(idata[,1]) else 
#                   ifacn=paste(ifacn,as.character(idata[,i]),sep=".")
#   facnr=ws[nifac:1]
#   e.mv=matrix(unlist(tapply(dat[,dvnam],dat[,facnr],function(x){x})),
#               ncol=length(ifacn),dimnames=list(snams,ifacn))
#   if (length(bs)==1) bsf=bs else bsf=paste(bs,collapse="*")
#   if (any(is.na(bs))) {
#     form=formulua("e.mv~1")
#     tmp <- e.mv
#   } else {
#     form=formula(paste("e.mv ~",bsf))
#     trans=snams
#     names(trans)=snams
#     for (i in bs) {
#       for (j in snams) trans[j] <- as.character(dat[dat[,scol]==j,i][1])
#       if (i==bs[1]) bsfac <- factor(trans[dimnames(e.mv)[[1]]]) else
#         bsfac <- cbind.data.frame(bsfac,factor(trans[dimnames(e.mv)[[1]]]))
#     }
#     tmp <- cbind.data.frame(bsfac,e.mv)
#     names(tmp)[1:length(bs)] <- bs
#   }
#   summary(Anova(lm(form,tmp),
#     idata=idata,type=SStype,icontrasts=icontrasts,
#     idesign=formula(paste("~",paste(ws,collapse="*"))) ),
#     multivariate=FALSE)
#   invisible(tmp)
# }


# get.cors <- function(thetas, data) {
#   
#   CORS <- apply(thetas, c(1,2,3), function(x) cor(x, data)) 
#   post_medians <- apply(CORS, c(2), median)
#   post_LCI <- apply(CORS, c(2), quantile, prob=0.025)
#   post_HCI <- apply(CORS, c(2), quantile, prob=0.975)
#   out <- list(post_medians, post_LCI, post_HCI)
#   names(out) <- c("medians", "LCI", "HCI")
#   out
# }


cors.plot <- function(out) {
  tmp <- data.frame(cbind(out[[1]],
                                   out[[2]], out[[3]]))
  colnames(tmp)<- c("M", "LCI", "HCI")
  tmp$param <- names(out[[1]])
  ggdf<-tmp
  ggplot(ggdf, aes(param, M))+
  geom_point(size=3)+
  ylab("")+ 
  geom_hline(aes(yintercept=0), linetype=2) +
  geom_errorbar(data=ggdf,aes(ymin=LCI,ymax=HCI)) 
  
}



x1pmaccs <- function (currentsim) {

pmaccdiff=NA;nonfocalPM=NA; focalPM=NA
  
focalPM= length(currentsim$RT[currentsim$S=="p" & currentsim$R=="P" & currentsim$E=="F"] )/length(currentsim$RT[currentsim$S=="p" & currentsim$E=="F"])  
nonfocalPM= length(currentsim$RT[currentsim$S=="p" & currentsim$R=="P" & currentsim$E=="H"] )/length(currentsim$RT[currentsim$S=="p" & currentsim$E=="H"])
pmaccdiff <- length(currentsim$RT[currentsim$S=="p" & currentsim$R=="P" & currentsim$E=="F"] )/length(currentsim$RT[currentsim$S=="p" & currentsim$E=="F"] ) -  length(currentsim$RT[currentsim$S=="p" & currentsim$R=="P" & currentsim$E=="H"] )/length(currentsim$RT[currentsim$S=="p" & currentsim$E=="H"] )


out <- c(focalPM, nonfocalPM, pmaccdiff)
names(out) <- c("focalPM", "nonfocalPM","pmeffect")
out


}

x1pmaccdiff <- function (currentsim) {

pmaccdiff=NA
pmaccdiff <- length(currentsim$RT[currentsim$S=="p" & currentsim$R=="P" & currentsim$E=="F"] )/length(currentsim$RT[currentsim$S=="p" & currentsim$E=="F"] ) -  length(currentsim$RT[currentsim$S=="p" & currentsim$R=="P" & currentsim$E=="H"] )/length(currentsim$RT[currentsim$S=="p" & currentsim$E=="H"] )


out <- c(pmaccdiff)
names(out) <- c("pmeffect")
out


}

x1pmeffects <- function (currentsim) {

pmaccdiff=NA;focalcostwW =NA;nonfocalcostwW=NA; focalcostwN=NA; nonfocalcostwN=NA
focalcostnN=NA;nonfocalcostnN=NA;focalcostnW=NA; nonfocalcostnW=NA
focalaccw=NA; nonfocalaccw=NA;focalaccn=NA;nonfocalaccn=NA
pmrtdiff=NA; nonfocalPM=NA; focalPM=NA
  
focalPM= length(currentsim$RT[currentsim$S=="p" & currentsim$R=="P" & currentsim$E=="F"] )/length(currentsim$RT[currentsim$S=="p" & currentsim$E=="F"])  
nonfocalPM= length(currentsim$RT[currentsim$S=="p" & currentsim$R=="P" & currentsim$E=="H"] )/length(currentsim$RT[currentsim$S=="p" & currentsim$E=="H"])
pmaccdiff <- length(currentsim$RT[currentsim$S=="p" & currentsim$R=="P" & currentsim$E=="F"] )/length(currentsim$RT[currentsim$S=="p" & currentsim$E=="F"] ) -  length(currentsim$RT[currentsim$S=="p" & currentsim$R=="P" & currentsim$E=="H"] )/length(currentsim$RT[currentsim$S=="p" & currentsim$E=="H"] )
pmrtdiff <-   mean(currentsim$RT[currentsim$S=="p"  & currentsim$E=="H"] ) - mean(currentsim$RT[currentsim$S=="p" & currentsim$E=="F"] )



focalcostwW <- mean(na.rm=T,currentsim$RT[currentsim$S=="w" & currentsim$R=="W" & currentsim$E=="F"] ) -  mean(na.rm=T,currentsim$RT[currentsim$S=="w" & currentsim$R=="W" & currentsim$E=="2"])
nonfocalcostwW <- mean(na.rm=T,currentsim$RT[currentsim$S=="w" & currentsim$R=="W" & currentsim$E=="H"] ) -  mean(na.rm=T,currentsim$RT[currentsim$S=="w" & currentsim$R=="W" & currentsim$E=="2"])

focalcostwN <- mean(na.rm=T,currentsim$RT[currentsim$S=="w" & currentsim$R=="N" & currentsim$E=="F"] ) -  mean(na.rm=T,currentsim$RT[currentsim$S=="w" & currentsim$R=="N" & currentsim$E=="2"])
nonfocalcostwN <- mean(na.rm=T,currentsim$RT[currentsim$S=="w" & currentsim$R=="N" & currentsim$E=="H"] ) -  mean(na.rm=T,currentsim$RT[currentsim$S=="w" & currentsim$R=="N" & currentsim$E=="2"])

focalcostnN <- mean(na.rm=T,currentsim$RT[currentsim$S=="n" & currentsim$R=="N" & currentsim$E=="F"] ) -  mean(na.rm=T,currentsim$RT[currentsim$S=="n" & currentsim$R=="N" & currentsim$E=="2"])
nonfocalcostnN <- mean(na.rm=T,currentsim$RT[currentsim$S=="n" & currentsim$R=="N" & currentsim$E=="H"] ) -  mean(na.rm=T,currentsim$RT[currentsim$S=="n" & currentsim$R=="N" & currentsim$E=="2"])

focalcostnW <- mean(na.rm=T,currentsim$RT[currentsim$S=="n" & currentsim$R=="W" & currentsim$E=="F"] ) -  mean(na.rm=T,currentsim$RT[currentsim$S=="n" & currentsim$R=="W" & currentsim$E=="2"])
nonfocalcostnW <- mean(na.rm=T,currentsim$RT[currentsim$S=="n" & currentsim$R=="W" & currentsim$E=="H"] ) -  mean(na.rm=T,currentsim$RT[currentsim$S=="n" & currentsim$R=="W" & currentsim$E=="2"])


focalaccw <-length(currentsim$RT[currentsim$S=="w" & currentsim$R=="W" & currentsim$E=="F"] )/length(currentsim$RT[currentsim$S=="w" & currentsim$E=="F"] ) -  length(currentsim$RT[currentsim$S=="w" & currentsim$R=="W" & currentsim$E=="2"] )/length(currentsim$RT[currentsim$S=="w" & currentsim$E=="2"] )
nonfocalaccw <-length(currentsim$RT[currentsim$S=="w" & currentsim$R=="W" & currentsim$E=="H"] )/length(currentsim$RT[currentsim$S=="w" & currentsim$E=="H"] ) -  length(currentsim$RT[currentsim$S=="w" & currentsim$R=="W" & currentsim$E=="2"] )/length(currentsim$RT[currentsim$S=="w" & currentsim$E=="2"] )

focalaccn <-length(currentsim$RT[currentsim$S=="n" & currentsim$R=="N" & currentsim$E=="F"] )/length(currentsim$RT[currentsim$S=="n" & currentsim$E=="F"] ) -  length(currentsim$RT[currentsim$S=="n" & currentsim$R=="N" & currentsim$E=="2"] )/length(currentsim$RT[currentsim$S=="n" & currentsim$E=="2"] )
nonfocalaccn <-length(currentsim$RT[currentsim$S=="n" & currentsim$R=="N" & currentsim$E=="H"] )/length(currentsim$RT[currentsim$S=="n" & currentsim$E=="H"] ) -  length(currentsim$RT[currentsim$S=="n" & currentsim$R=="N" & currentsim$E=="2"] )/length(currentsim$RT[currentsim$S=="n" & currentsim$E=="2"] )



out <- c(focalPM, nonfocalPM, pmaccdiff,pmrtdiff ,focalcostwW, focalcostwN, focalcostnN, focalcostnW, nonfocalcostwW, nonfocalcostwN, nonfocalcostnN, nonfocalcostnW, focalaccw, nonfocalaccw, focalaccn, nonfocalaccn)
names(out) <- c("focalPM", "nonfocalPM","pmeffect", "pmrtdiff", "focalwW", "focalwN","focalnN", "focalnW", "nonfocalwW", "nonfocalwN", 
                "nonfocalnN", "nonfocalnW", "focalaccw", "nonfocalaccw", "focalaccn", "nonfocalaccn")
out


}


x1COSTPM <- function (currentsim) {

 nonfocalPM=NA; focalPM=NA
  
focalPM= length(currentsim$RT[currentsim$S=="p" & currentsim$R=="P" & currentsim$E=="F"] )/length(currentsim$RT[currentsim$S=="p" & currentsim$E=="F"])  
nonfocalPM= length(currentsim$RT[currentsim$S=="p" & currentsim$R=="P" & currentsim$E=="H"] )/length(currentsim$RT[currentsim$S=="p" & currentsim$E=="H"])



focalcostw <- mean(na.rm=T,currentsim$RT[currentsim$S=="w" & currentsim$E=="F"] ) -  mean(na.rm=T,currentsim$RT[currentsim$S=="w" & currentsim$E=="2"])
nonfocalcostw <- mean(na.rm=T,currentsim$RT[currentsim$S=="w" & currentsim$E=="H"] ) -  mean(na.rm=T,currentsim$RT[currentsim$S=="w" &  currentsim$E=="2"])

focalcostn <- mean(na.rm=T,currentsim$RT[currentsim$S=="n" &currentsim$E=="F"] ) -  mean(na.rm=T,currentsim$RT[currentsim$S=="n" &  currentsim$E=="2"])
nonfocalcostn <- mean(na.rm=T,currentsim$RT[currentsim$S=="n" & currentsim$E=="H"] ) -  mean(na.rm=T,currentsim$RT[currentsim$S=="n" & currentsim$E=="2"])

out <- c(focalPM, nonfocalPM, focalcostw, focalcostn, nonfocalcostw, nonfocalcostn)
names(out) <- c("focalPM", "nonfocalPM", "focalcostw", "focalcostn", "nonfocalcostw", "nonfocalcostn")
out


}


x2COSTPM<- function (currentsim) {
  
  IPMW=NA;UPMW=NA ;IPMN=NA;UPMN=NA
  impcostw =NA;unimpcostw=NA
  impcostn=NA;unimpcostnN=NA

    
  IPMW = length(currentsim$RT[currentsim$S=="pw" & currentsim$R=="P" & currentsim$E=="I"] )/length(currentsim$RT[currentsim$S=="pw" & currentsim$E=="I"] )
  UPMW = length(currentsim$RT[currentsim$S=="pw" & currentsim$R=="P" & currentsim$E=="U"] )/length(currentsim$RT[currentsim$S=="pw" & currentsim$E=="U"] )
  IPMN = length(currentsim$RT[currentsim$S=="pn" & currentsim$R=="P" & currentsim$E=="I"] )/length(currentsim$RT[currentsim$S=="pn" & currentsim$E=="I"] )
  UPMN = length(currentsim$RT[currentsim$S=="pn" & currentsim$R=="P" & currentsim$E=="U"] )/length(currentsim$RT[currentsim$S=="pn" & currentsim$E=="U"] )
  
  
  impcostw <- mean(currentsim$RT[currentsim$S=="ww"  & currentsim$E=="I"] ) -  mean(currentsim$RT[currentsim$S=="ww"  & currentsim$E=="2"])
  unimpcostw <- mean(currentsim$RT[currentsim$S=="ww"  & currentsim$E=="U"] ) -  mean(currentsim$RT[currentsim$S=="ww"  & currentsim$E=="2"])
  
  impcostw <- mean(currentsim$RT[currentsim$S=="ww"  & currentsim$E=="I"] ) -  mean(currentsim$RT[currentsim$S=="ww"  & currentsim$E=="2"])
  unimpcostw <- mean(currentsim$RT[currentsim$S=="ww"  & currentsim$E=="U"] ) -  mean(currentsim$RT[currentsim$S=="ww"  & currentsim$E=="2"])
  
  impcostn <- mean(currentsim$RT[currentsim$S=="nn"  & currentsim$E=="I"] ) -  mean(currentsim$RT[currentsim$S=="nn"  & currentsim$E=="2"])
  unimpcostn <- mean(currentsim$RT[currentsim$S=="nn"  & currentsim$E=="U"] ) -  mean(currentsim$RT[currentsim$S=="nn"  & currentsim$E=="2"])
  
  impcostn <- mean(currentsim$RT[currentsim$S=="nn"  & currentsim$E=="I"] ) -  mean(currentsim$RT[currentsim$S=="nn"  & currentsim$E=="2"])
  unimpcostn <- mean(currentsim$RT[currentsim$S=="nn"  & currentsim$E=="U"] ) -  mean(currentsim$RT[currentsim$S=="nn"  & currentsim$E=="2"])
  
  
  out <- c(
    IPMW,UPMW,IPMN,UPMN,
    impcostw,impcostn, unimpcostw, unimpcostn)
  names(out) <- c(
     "IPMW","UPMW","IPMN","UPMN",
    "impw", "impn","unimpw", 
                  "unimpn")
  out


}





x1covars <- function (currentsim) {


focalCVwW <-NA
focalCVpW <- NA
focalCVpP <- NA



nonfocalCVwW <- NA
nonfocalCVpW <-NA
nonfocalCVpP <-NA



focalCVwW <- sd(na.rm=T,currentsim$RT[currentsim$S=="w" & currentsim$R=="W" & currentsim$E=="F"] )/ mean(na.rm=T,currentsim$RT[currentsim$S=="w" & currentsim$R=="W" & currentsim$E=="F"] )
focalCVpW <- sd(na.rm=T,currentsim$RT[currentsim$S=="p" & currentsim$R=="W" & currentsim$E=="F"] )/ mean(na.rm=T,currentsim$RT[currentsim$S=="p" & currentsim$R=="W" & currentsim$E=="F"] )
focalCVpP <- sd(na.rm=T,currentsim$RT[currentsim$S=="p" & currentsim$R=="P" & currentsim$E=="F"] )/ mean(na.rm=T,currentsim$RT[currentsim$S=="p" & currentsim$R=="P" & currentsim$E=="F"] )



nonfocalCVwW <- sd(na.rm=T,currentsim$RT[currentsim$S=="w" & currentsim$R=="W" & currentsim$E=="H"] )/ mean(na.rm=T,currentsim$RT[currentsim$S=="w" & currentsim$R=="W" & currentsim$E=="H"] )
nonfocalCVpW <- sd(na.rm=T,currentsim$RT[currentsim$S=="p" & currentsim$R=="W" & currentsim$E=="H"] )/ mean(na.rm=T,currentsim$RT[currentsim$S=="p" & currentsim$R=="W" & currentsim$E=="H"] )
nonfocalCVpP <- sd(na.rm=T,currentsim$RT[currentsim$S=="p" & currentsim$R=="P" & currentsim$E=="H"] )/ mean(na.rm=T,currentsim$RT[currentsim$S=="p" & currentsim$R=="P" & currentsim$E=="H"] )



out <- c(nonfocalCVwW, nonfocalCVpW, nonfocalCVpP,
         focalCVwW, focalCVpW, focalCVpP)
         
         
names(out) <- c("nonfocalCVwW", "nonfocalCVpW", "nonfocalCVpP",
         "focalCVwW", "focalCVpW", "focalCVpP")
out


}



x2covars <- function (currentsim) {


importantCVwwW <-NA
importantCVpwW <- NA
importantCVpwP <- NA

importantCVnnN <-NA
importantCVpnN <- NA
importantCVpnP <- NA


unimportantCVwwW <-NA
unimportantCVpwW <- NA
unimportantCVpwP <- NA

unimportantCVnnN <-NA
unimportantCVpnN <- NA
unimportantCVpnP <- NA



importantCVwwW  <- sd(na.rm=T,currentsim$RT[currentsim$S=="ww" & currentsim$R=="W" & currentsim$E=="I"] )/ mean(na.rm=T,currentsim$RT[currentsim$S=="ww" & currentsim$R=="W" & currentsim$E=="I"] )
importantCVpwW <- sd(na.rm=T,currentsim$RT[currentsim$S=="pw" & currentsim$R=="W" & currentsim$E=="I"] )/ mean(na.rm=T,currentsim$RT[currentsim$S=="pw" & currentsim$R=="W" & currentsim$E=="I"] )
importantCVpwP  <- sd(na.rm=T,currentsim$RT[currentsim$S=="pw" & currentsim$R=="P" & currentsim$E=="I"] )/ mean(na.rm=T,currentsim$RT[currentsim$S=="pw" & currentsim$R=="P" & currentsim$E=="I"] )

importantCVnnN  <- sd(na.rm=T,currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$E=="I"] )/ mean(na.rm=T,currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$E=="I"] )
importantCVpnN <- sd(na.rm=T,currentsim$RT[currentsim$S=="pn" & currentsim$R=="N" & currentsim$E=="I"] )/ mean(na.rm=T,currentsim$RT[currentsim$S=="pn" & currentsim$R=="N" & currentsim$E=="I"] )
importantCVpnP  <- sd(na.rm=T,currentsim$RT[currentsim$S=="pn" & currentsim$R=="P" & currentsim$E=="I"] )/ mean(na.rm=T,currentsim$RT[currentsim$S=="pn" & currentsim$R=="P" & currentsim$E=="I"] )

unimportantCVwwW  <- sd(na.rm=T,currentsim$RT[currentsim$S=="ww" & currentsim$R=="W" & currentsim$E=="U"] )/ mean(na.rm=T,currentsim$RT[currentsim$S=="ww" & currentsim$R=="W" & currentsim$E=="U"] )
unimportantCVpwW <- sd(na.rm=T,currentsim$RT[currentsim$S=="pw" & currentsim$R=="W" & currentsim$E=="U"] )/ mean(na.rm=T,currentsim$RT[currentsim$S=="pw" & currentsim$R=="W" & currentsim$E=="U"] )
unimportantCVpwP  <- sd(na.rm=T,currentsim$RT[currentsim$S=="pw" & currentsim$R=="P" & currentsim$E=="U"] )/ mean(na.rm=T,currentsim$RT[currentsim$S=="pw" & currentsim$R=="P" & currentsim$E=="U"] )

unimportantCVnnN  <- sd(na.rm=T,currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$E=="U"] )/ mean(na.rm=T,currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$E=="U"] )
unimportantCVpnN <- sd(na.rm=T,currentsim$RT[currentsim$S=="pn" & currentsim$R=="N" & currentsim$E=="U"] )/ mean(na.rm=T,currentsim$RT[currentsim$S=="pn" & currentsim$R=="N" & currentsim$E=="U"] )
unimportantCVpnP  <- sd(na.rm=T,currentsim$RT[currentsim$S=="pn" & currentsim$R=="P" & currentsim$E=="U"] )/ mean(na.rm=T,currentsim$RT[currentsim$S=="pn" & currentsim$R=="P" & currentsim$E=="U"] )





out <- c(
  
  importantCVwwW,
importantCVpwW,
importantCVpwP,

importantCVnnN,
importantCVpnN,
importantCVpnP,


unimportantCVwwW ,
unimportantCVpwW,
unimportantCVpwP,

unimportantCVnnN ,
unimportantCVpnN,
unimportantCVpnP

  
)
         
         
names(out) <- c(
  
"importantCVwwW",
"importantCVpwW",
"importantCVpwP",

"importantCVnnN",
"importantCVpnN",
"importantCVpnP",


"unimportantCVwwW ",
"unimportantCVpwW",
"unimportantCVpwP",

"unimportantCVnnN ",
"unimportantCVpnN",
"unimportantCVpnP"

)
out


}



x2pmeffects <- function (currentsim) {
  
  pmaccdiffw=NA;pmaccdiffn=NA;impcostwW =NA;unimpcostwW=NA; impcostwN=NA; unimpcostwN=NA
  impcostnN=NA;unimpcostnN=NA;impcostnW=NA; unimpcostnW=NA
  impaccw=NA; unimpaccw=NA;impaccn=NA;unimpaccn=NA
  pmrtdiffw=NA; pmrtdiffn=NA; IPMW=NA;UPMW=NA ;IPMN=NA;UPMN=NA
    
  IPMW = length(currentsim$RT[currentsim$S=="pw" & currentsim$R=="P" & currentsim$E=="I"] )/length(currentsim$RT[currentsim$S=="pw" & currentsim$E=="I"] )
  UPMW = length(currentsim$RT[currentsim$S=="pw" & currentsim$R=="P" & currentsim$E=="U"] )/length(currentsim$RT[currentsim$S=="pw" & currentsim$E=="U"] )
  IPMN = length(currentsim$RT[currentsim$S=="pn" & currentsim$R=="P" & currentsim$E=="I"] )/length(currentsim$RT[currentsim$S=="pn" & currentsim$E=="I"] )
  UPMN = length(currentsim$RT[currentsim$S=="pn" & currentsim$R=="P" & currentsim$E=="U"] )/length(currentsim$RT[currentsim$S=="pn" & currentsim$E=="U"] )
  
  pmaccwdiff <- length(currentsim$RT[currentsim$S=="pw" & currentsim$R=="P" & currentsim$E=="I"] )/length(currentsim$RT[currentsim$S=="pw" & currentsim$E=="I"] ) -  length(currentsim$RT[currentsim$S=="pw" & currentsim$R=="P" & currentsim$E=="U"] )/length(currentsim$RT[currentsim$S=="pw" & currentsim$E=="U"] )
  pmaccndiff <- length(currentsim$RT[currentsim$S=="pn" & currentsim$R=="P" & currentsim$E=="I"] )/length(currentsim$RT[currentsim$S=="pn" & currentsim$E=="I"] ) -  length(currentsim$RT[currentsim$S=="pn" & currentsim$R=="P" & currentsim$E=="U"] )/length(currentsim$RT[currentsim$S=="pn" & currentsim$E=="U"] )
  
  
  pmrtwdiff <-   mean(currentsim$RT[currentsim$S=="pw"  & currentsim$E=="U"] ) - mean(currentsim$RT[currentsim$S=="pw" & currentsim$E=="I"] )
  pmrtndiff <-   mean(currentsim$RT[currentsim$S=="pn"  & currentsim$E=="U"] ) - mean(currentsim$RT[currentsim$S=="pn" & currentsim$E=="I"] )
  
  
  impcostwW <- mean(currentsim$RT[currentsim$S=="ww" & currentsim$R=="W" & currentsim$E=="I"] ) -  mean(currentsim$RT[currentsim$S=="ww" & currentsim$R=="W" & currentsim$E=="2"])
  unimpcostwW <- mean(currentsim$RT[currentsim$S=="ww" & currentsim$R=="W" & currentsim$E=="U"] ) -  mean(currentsim$RT[currentsim$S=="ww" & currentsim$R=="W" & currentsim$E=="2"])
  
  impcostwN <- mean(currentsim$RT[currentsim$S=="ww" & currentsim$R=="N" & currentsim$E=="I"] ) -  mean(currentsim$RT[currentsim$S=="ww" & currentsim$R=="N" & currentsim$E=="2"])
  unimpcostwN <- mean(currentsim$RT[currentsim$S=="ww" & currentsim$R=="N" & currentsim$E=="U"] ) -  mean(currentsim$RT[currentsim$S=="ww" & currentsim$R=="N" & currentsim$E=="2"])
  
  impcostnN <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$E=="I"] ) -  mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$E=="2"])
  unimpcostnN <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$E=="U"] ) -  mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$E=="2"])
  
  impcostnW <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="W" & currentsim$E=="I"] ) -  mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="W" & currentsim$E=="2"])
  unimpcostnW <- mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="W" & currentsim$E=="U"] ) -  mean(currentsim$RT[currentsim$S=="nn" & currentsim$R=="W" & currentsim$E=="2"])
  
  
  impaccw <-length(currentsim$RT[currentsim$S=="ww" & currentsim$R=="W" & currentsim$E=="I"] )/length(currentsim$RT[currentsim$S=="ww" & currentsim$E=="I"] ) -  length(currentsim$RT[currentsim$S=="ww" & currentsim$R=="W" & currentsim$E=="2"] )/length(currentsim$RT[currentsim$S=="ww" & currentsim$E=="2"] )
  unimpaccw <-length(currentsim$RT[currentsim$S=="ww" & currentsim$R=="W" & currentsim$E=="U"] )/length(currentsim$RT[currentsim$S=="ww" & currentsim$E=="U"] ) -  length(currentsim$RT[currentsim$S=="ww" & currentsim$R=="W" & currentsim$E=="2"] )/length(currentsim$RT[currentsim$S=="ww" & currentsim$E=="2"] )
  
  impaccn <-length(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$E=="I"] )/length(currentsim$RT[currentsim$S=="nn" & currentsim$E=="I"] ) -  length(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$E=="2"] )/length(currentsim$RT[currentsim$S=="nn" & currentsim$E=="2"] )
  unimpaccn <-length(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$E=="U"] )/length(currentsim$RT[currentsim$S=="nn" & currentsim$E=="U"] ) -  length(currentsim$RT[currentsim$S=="nn" & currentsim$R=="N" & currentsim$E=="2"] )/length(currentsim$RT[currentsim$S=="nn" & currentsim$E=="2"] )
  
  
  
  out <- c(
    IPMW,UPMW,IPMN,UPMN,
    pmaccwdiff,pmaccndiff ,pmrtwdiff,pmrtndiff ,impcostwW, impcostwN, impcostnN, impcostnW, unimpcostwW, unimpcostwN, unimpcostnN, unimpcostnW, impaccw, unimpaccw, impaccn, unimpaccn)
  names(out) <- c(
     "IPMW","UPMW","IPMN","UPMN",
    "pmeffectw", "pmeffectn", "pmrtdiffw","pmrtdiffn", "impwW", "impwN","impnN", "impnW", "unimpwW", "unimpwN", 
                  "unimpnN", "unimpnW", "impaccw", "unimpaccw", "impaccn", "unimpaccn")
  out


}



x2pmaccs <- function (currentsim) {
  
  pmaccdiffw=NA;pmaccdiffn=NA; IPMW=NA;UPMW=NA ;IPMN=NA;UPMN=NA
    
  IPMW = length(currentsim$RT[currentsim$S=="pw" & currentsim$R=="P" & currentsim$E=="I"] )/length(currentsim$RT[currentsim$S=="pw" & currentsim$E=="I"] )
  UPMW = length(currentsim$RT[currentsim$S=="pw" & currentsim$R=="P" & currentsim$E=="U"] )/length(currentsim$RT[currentsim$S=="pw" & currentsim$E=="U"] )
  IPMN = length(currentsim$RT[currentsim$S=="pn" & currentsim$R=="P" & currentsim$E=="I"] )/length(currentsim$RT[currentsim$S=="pn" & currentsim$E=="I"] )
  UPMN = length(currentsim$RT[currentsim$S=="pn" & currentsim$R=="P" & currentsim$E=="U"] )/length(currentsim$RT[currentsim$S=="pn" & currentsim$E=="U"] )
  
  pmaccwdiff <- length(currentsim$RT[currentsim$S=="pw" & currentsim$R=="P" & currentsim$E=="I"] )/length(currentsim$RT[currentsim$S=="pw" & currentsim$E=="I"] ) -  length(currentsim$RT[currentsim$S=="pw" & currentsim$R=="P" & currentsim$E=="U"] )/length(currentsim$RT[currentsim$S=="pw" & currentsim$E=="U"] )
  pmaccndiff <- length(currentsim$RT[currentsim$S=="pn" & currentsim$R=="P" & currentsim$E=="I"] )/length(currentsim$RT[currentsim$S=="pn" & currentsim$E=="I"] ) -  length(currentsim$RT[currentsim$S=="pn" & currentsim$R=="P" & currentsim$E=="U"] )/length(currentsim$RT[currentsim$S=="pn" & currentsim$E=="U"] )
  
  out <- c(
    IPMW,UPMW,IPMN,UPMN,
    pmaccwdiff,pmaccndiff 
    )
  names(out) <- c(
     "IPMW","UPMW","IPMN","UPMN",
    "pmeffectw", "pmeffectn")
  out


}

cors.plot <- function(out) {
  tmp <- data.frame(cbind(out[[1]],
                                   out[[2]], out[[3]]))
  colnames(tmp)<- c("M", "LCI", "HCI")
  tmp$param <- factor(names(out[[1]]), levels= names(out[[1]]))
  ggdf<-tmp
  ggplot(ggdf, aes(param, M))+
  geom_point(size=3)+
    xlab("")+
  geom_hline(aes(yintercept=0), linetype=2) +
  geom_errorbar(data=ggdf,aes(ymin=LCI,ymax=HCI)) 
}



#### Hacky stuff to get different style posterior exploration graphs



hacked.RT.dmc <- function (df, xaxis = 'R', panels.ppage=4, do.quantiles=TRUE,
                           nrow=NULL, ncol=NULL)
 
  # xaxis is the name of the factor to put on the X axis
{
 
  if (!is.null(attr(df, "gglist"))) {
    cat ("Treating as a pp object for a single participant")
    df <- attr(df, "gglist")[[2]]
  }

  if (!is.null(attr(attr(df, "av"), "gglist"))) {
    cat ("Treating as a pp object for a group of participants")
    df <- attr(attr(df, "av"), "gglist")[[2]]
  }

  df <- df[!is.na(df$median),]

  # get factors (other than the xaxis factor) for the grid
  if (do.quantiles) grid <- colnames(df) [!colnames (df) %in%
      c(xaxis, "median", "lower", "upper", "data", "quantile")] else
    grid <- colnames(df) [!colnames (df) %in%
      c(xaxis, "median", "lower", "upper", "data")]

  if (length(grid)==0) {
    n.plots <-1
    plot.ind <- rep(1, length(df$data))
  } else {
    n.plots <- sum (!is.na((tapply(df$data, df[,c(grid)], function (x) 1))))
    plot.ind <- as.numeric(factor(with(df, interaction(df[,grid]))))
  }

  plot.pages <- ceiling(n.plots/panels.ppage)

  # each step of the loop grabs enough data to create number of plots p page
  plots <- list()
  for (j in 1:plot.pages) {

    active.df <- df[plot.ind %in% ((j-1)*panels.ppage+1):(j* panels.ppage),]

    plot <- ggplot(active.df, aes_string(x = xaxis, y= 'median')) + geom_point(size=5) +
        geom_errorbar(aes(ymax = upper, ymin = lower), width= 0.2) +
        geom_point(aes_string(x = xaxis, y= 'data'), pch=21, size=5, colour="black") +
        ylab("RT")
    if (do.quantiles)
      plot <- plot + geom_line(aes(group = quantile, y=data), linetype=2) else
      plot <- plot + geom_line(aes(group=1,y=data), linetype=2)

    if(length(grid)!=0) plot <- plot +
      facet_wrap(grid, labeller = labeller(label_value, .multi_line=F), scales = 'free',
                 nrow=nrow, ncol=ncol)

    plots[[j]] <- plot
  }

  if (length(plots)==1) plots[[1]] + theme_simple() else lapply(plots, function(x) x + theme_simple())
  }


hacked.RP.dmc <- function (df, xaxis = 'R', panels.ppage=4, nrow=NULL, ncol=NULL)
  # xaxis is the name of the factor to put on the X axis
{

  if (!is.null(attr(df, "gglist"))) {
    cat ("Treating as a pp object for a single participant")
    df <- attr(df, "gglist")[[1]]
  }

  if (!is.null(attr(attr(df, "av"), "gglist"))) {
    cat ("Treating as a pp object for a group of participants")
    df <- attr(attr(df, "av"), "gglist")[[1]]
  }

  #remove cases where the combination of factors was not present in design
  df <- df[!is.na(df$data),]

  grid <- colnames(df) [!colnames (df) %in% c(xaxis, "median", "lower", "upper", "data")]

  if (length(grid)==0) {
    n.plots <-1
    plot.ind <- rep(1, length(df$data))
  } else {
    n.plots <- sum (!is.na((tapply(df$data, df[,c(grid)], function (x) 1))))
    plot.ind <- as.numeric(factor(with(df, interaction(df[,grid]))))
  }

  plot.pages <- ceiling(n.plots/panels.ppage)

  grid_labels <- labeller(label_value, .multi_line=F)

  plots <- list()
  for (j in 1:plot.pages) {

    active.df <- df[plot.ind %in% ((j-1)*panels.ppage+1):(j* panels.ppage),]

    plot <- ggplot(active.df, aes_string(x = xaxis, y= 'median')) +
      geom_point(size=5) + geom_errorbar(aes(ymax = upper, ymin = lower), width= 0.2) +
      geom_point(aes_string(x = xaxis, y= 'data'), pch=21, size=5, colour="black") +
      geom_line(aes(group = 1, y=data), linetype=2) + ylab("Response Proportion")

    if (length(grid)!=0)
      plot <- plot + facet_wrap(grid, labeller = grid_labels, scales = 'free_x',
                                nrow=nrow, ncol=ncol)
    plots[[j]] <- plot
  }

  if (length(plots)==1) plots[[1]] + theme_simple() else lapply(plots, function(x) x + theme_simple())
}
plot.cost.x1 <- function (test) {

  names(test)[1] <- "median"
  test$effect  <-  c("F", "NF", 
                                              "words F","non-words F",
                                             "words NF",
                                             "non-words NF")
  
  hacked.RT.dmc(test[3:6,], xaxis="effect", do.quantiles=F) +ylab("Cost to RT") +xlab("")+
  scale_y_continuous(breaks=c(-0.04,
                              -0.02, 0.0, 0.02,0.04,0.06,0.08), limits=c(-0.04,.08), name="Cost to RT") 
}

plot.PM.x1 <- function (test) {
  names(test)[1] <- "median"
  test$effect  <- c("F", "NF", 
                                              "words F","non-words F",
                                             "words NF",
                                             "non-words NF")
  
  hacked.RP.dmc(test[1:2,], xaxis="effect") +ylab("PM Accuracy") +ylim(0.3,1) +xlab("
                                                                                    ")
}


plot.cost.x2 <- function (test) {

  names(test)[1] <- "median"
  test$effect  <-  c("Iw", "Uw", "In", "Un", 
                                              "words I","non-words I",
                                             "words U",
                                             "non-words U")
  
  test$effect <- factor(test$effect, levels=c("Iw", "Uw", "In", "Un", 
                                              "words I","non-words I",
                                             "words U",
                                             "non-words U")
  )

  hacked.RT.dmc(test[5:8,], xaxis="effect", do.quantiles=F) +ylab("Cost to RT") +xlab("")+
  scale_y_continuous(breaks=c(-0.04,
                              -0.02, 0.0, 0.02,0.04,0.06,0.08,0.1,0.12), limits=c(-0.04,.12), name="Cost to RT") 
}

plot.PM.x2 <- function (test) {
  names(test)[1] <- "median"
   test$effect  <-  c("Ipw", "Upw", "Ipn", "Upn", 
                                              "words I","non-words I",
                                             "words U",
                                             "non-words U")
  
  hacked.RP.dmc(test[1:4,], xaxis="effect") +ylab("PM Accuracy") +ylim(0.2,1) +xlab("
                                                                                    ")
}


fast.avgsamples <- function (hsamples) {
  nmcs<- sapply(hsamples, function(x) x$nmc)
  nmc <- min(nmcs)
#Different numbers of nmc for each participant... use the min number and then
  #for participants wtih more randomly sample out that many
  for (i in 1:length(hsamples)) if (nmcs[i] > nmc) hsamples[[i]]$theta <- 
    hsamples[[i]]$theta[,,sample(1:dim(hsamples[[i]]$theta)[3], nmc)]
  avg <- list()
  for (i in 1:length(hsamples)) avg [[i]] <- hsamples[[i]]$theta
  avg2 <- unlist(avg)
  dim3 <- c(dim(avg[[1]]), length(avg2)/prod(dim(avg[[1]])))
  dim(avg2) <- dim3
  out <- apply(avg2, c(1,2,3), mean)
  colnames(out) <- dimnames(avg[[1]])[[2]]
  out
}


get.msds <- function(samples) {
  avg_samples <- fast.avgsamples(samples)
  M <- apply(avg_samples, 2, mean)
  SD <- apply(avg_samples, 2, sd)
  out <- cbind(M,SD)
  colnames(out) <- c("M", "SD")
  data.frame(out)
}





# 
# 
# library("ggplot2")
# library("abind")
# theme_luke <- function (base_size = 12, base_family = "") {
#   theme_gray(base_size = base_size, base_family = base_family) %+replace% 
#     theme(
#       panel.background = element_rect(fill="white"),
#       panel.grid.minor.y = element_blank(),
#       legend.key = element_rect(fill="white", colour= "white"),
#       strip.background = element_rect(fill="white")
#       
#     )   
# }
# 
# 
# ###Basic proportion/quantile predictice function. Sims a df same size as data for each n.post.
# quantile.proportion.predictives.dmc <- function(samples,n.post=100,report=10, probs= c (0.1, 0.5, 0.9)) {
#   model <- attributes(samples$data)$model
#   facs <- names(attr(model,"factors")); nfacs <- length(facs)
#   resp <- names(attr(model,"responses")) ## sometimes do not have names
#   ns <- table(samples$data[,facs])  ## number of level of stimulus factor
#   n.rep <- sum(ns)     ## total number of data points
#   n.par <- dim(samples$theta)[2]    ## dim[1]=chains; dim[2]=list of para; dim[3]=samples
#   
#   ## accumulate chain1, chain2, chain3, ... all on one column for each para
#   thetas <- matrix(aperm(samples$theta,c(1,3,2)),ncol=n.par)
#   colnames(thetas) <- dimnames(samples$theta)[[2]]
#   
#   
#   ## random select n.post samples from thetas for each para
#   posts <- thetas[sample(c(1:nrow(thetas)), n.post, replace=F),]
#   
#   ## Construct a NA container with same dimension as data to store simulation data
#   sim <- data.frame(matrix(nrow=n.post*n.rep, ncol=ncol(samples$data)))
#   names(sim) <- names(samples$data)
#   
#   
#   for (i in names(samples$data)) {
#     sim[[i]] <- factor(rep(NA,n.post*n.rep),levels=levels(samples$data[[i]]))
#   }
#   ## Every "report" step prints one dot
#   message(paste0("\n Simulating (each dot = ",report," steps): "))
#   
#   ## before looping paste a string together to tell the aggregate function about the design   
#   form <- "RT~"
#   form <- paste(form, facs[1], sep="")
#   for (k in 2:length(facs)) form <- paste(form,"+", paste(facs[k]));designform<- form; form <- paste (form, " + R")
#   form <- as.formula(form)
#   
#   ##make data frame design.grid with number of trial for each combination. exclude 0s   
#   full.design <- expand.grid(lapply(samples$data[,1:(length(facs)+1)],unique))   
#   npresent <- aggregate(as.formula(designform), data=samples$data, na.action=NULL, FUN="length")
#   designgrid <- merge(full.design, npresent, by =facs)
#   names(designgrid)[length(designgrid)] <- "npresent"
#   designgrid$nresps <- 0  
#   
#   
#   # data calculations
#   data.qqs <- aggregate(form,  data=samples$data, FUN="quantile", probs= probs)
#   data.respsmade <- aggregate(form, data=samples$data, na.action=NULL, FUN="length")
#   data.pps<-designgrid
#   data.pps$nresps[match(do.call(paste,data.respsmade[,1:(length(facs)+1)]),do.call(paste,data.pps[,1:(length(facs)+1)]))]<-data.respsmade$RT
#   data.pps$proportion <- data.pps$nresp/data.pps$npresent
#   
#   
#   
#   ## Posterior sim loop
#   ## 1:1*nrep, 20001:2*n.rep, 40001:3*n.rep
#   
#   
#   
#   for (i in 1:n.post) {
#     
#     currentsim <- simulate.dmc(posts[i,],model,ns)
#     if ( (i %% report) == 0) message(".", appendLF = FALSE)
#     
#     qqs <-  aggregate(form,  data=currentsim, FUN="quantile", probs= probs)
#     qqs$rep <- i
#     respsmade <- aggregate(form,  data=currentsim, FUN="length")
#     pps<-designgrid
#     pps$nresps[match(do.call(paste,respsmade[,1:(length(facs)+1)]),do.call(paste,pps[,1:(length(facs)+1)]))]<-respsmade$RT
#     pps$proportion <- pps$nresp/pps$npresent
#     
#     
#     ###
#     pps$rep <- i
#     if (i==1) { fullpps<- pps; fullqqs <- qqs } 
#     else {fullpps <- rbind(fullpps, pps); fullqqs <- rbind(fullqqs, qqs)}
#     
#   }
#   out <- list(fullpps, fullqqs, data.pps, data.qqs)
#   names(out) <- c("proportions", "quantiles", "data.proportions", "data.quantiles")
#   out
#   
# }
# 
# 
# ## creates a data frame with the data and sim proportions averaged over subjects
# ## averages the proportions over all subjects into a single df for both data and model
# # takes from pp (.proportions), converts from list (df), averages (meansim, meandata)
# group.pp <- function (pp) {
#   sim.proportions <- lapply (pp, function(x) x[["proportions"]])
#   for (i in 1:length(sim.proportions)) sim.proportions[[i]]$subject <- as.character(i)
#   simdf<- do.call("rbind", sim.proportions)
#   form <- "proportion~"
#   form <- paste(form, colnames(simdf)[1],sep="")
#   for (k in 2:(length(colnames(simdf))-5)) form <- paste(form,"+", paste(colnames(simdf)[k]));dataform <-form; form <- paste(form,"+", "rep")
#   meansim <- aggregate(as.formula(form),  data=simdf, FUN="mean")
#   data.proportions <- lapply (pp, function(x) x[["data.proportions"]])
#   for (i in 1:length(data.proportions)) data.proportions[[i]]$subject <- as.character(i)
#   datadf<- do.call("rbind", data.proportions)
#   meandata <- aggregate(as.formula(dataform),  data=datadf, FUN="mean")
#   out <- list(meansim, meandata)
#   names(out) <- c("sim", "data")
#   out
# }
# 
# #### takes averaged proportions and for the model calculates mean, lower quantile, upper quantile (i.e bayesian error bars)
# ggplot.proportions.dmc <- function (pp, form="S+R", lower=0.025, upper= 0.975) {
#   if (names(pp)[1]== "proportions") {sim <- pp[["proportions"]]; dat <- pp[["data.proportions"]]} else {
#     sim <- group.pp (pp)  [["sim"]]
#     dat <- group.pp (pp)  [["data"]]
#   }
#   form <- paste("proportion~", form)
#   sim.df <- aggregate (as.formula(form), data= sim, FUN= function(x) c(mean(x), quantile(x, probs= lower), quantile(x, probs= upper)))
#   data.df <- aggregate (as.formula(form), data= dat, FUN= "mean");names(data.df)[length(data.df)] <- "data"
#   plot.df <- join(sim.df, data.df)
#   plot.df <- do.call(data.frame, plot.df)
#   names(plot.df)[(length(plot.df)-3):(length(plot.df)-1)]<- c("mean", "lower", "upper")
#   plot.df
# }
# 
# single.qq <- function (x) {
#   out <- do.call(data.frame, x)
#   index <-grep("RT.", colnames(out))
#   names(out)[index] <- paste("q", dimnames(x$RT)[[2]], sep="")
#   names(out) <- strsplit(names(out), "%")
#   out[,c(1:index[1]-1, length(out), index[1]:(length(out)-1))]
# }
# 
# group.qq <- function (pp) {
#   sim.RTs <- lapply (pp, function(x) x[["quantiles"]])
#   for (i in 1:length(sim.RTs)) sim.RTs[[i]]$subject <- as.character(i)
#   simdf<- do.call("rbind", sim.RTs)
#   form <- "RT~"
#   form <- paste(form, colnames(simdf)[1],sep="")
#   for (k in 2:(length(colnames(simdf))-3)) form <- paste(form,"+", paste(colnames(simdf)[k]));dataform <-form; form <- paste(form,"+", "rep")
#   # averages the quantiles but note also turns them into columns of the df
#   meansim <- aggregate(as.formula(form),  data=simdf, FUN="mean")
#   ## fix issue with quantile names starting with numbers
#   index <-grep("rep", colnames(meansim))
#   names(meansim)[index+1:(length(meansim)-index)] <- paste("q", names(meansim)[index+1:(length(meansim)-index)], sep="")
#   names(meansim) <- strsplit(names(meansim), "%")
#   
#   data.RTs <- lapply (pp, function(x) x[["data.quantiles"]])
#   for (i in 1:length(data.RTs)) data.RTs[[i]]$subject <- as.character(i)
#   datadf<- do.call("rbind", data.RTs)
#   meandata <- aggregate(as.formula(dataform),  data=datadf, FUN="mean")
#   ##fix naming for data as well 
#   names(meandata)[index:length(meandata)] <- paste("q", names(meandata)[index:length(meandata)], sep="")
#   names(meandata) <- strsplit(names(meandata), "%")
#   out <- list(meansim, meandata)
#   names(out) <- c("sim", "data")
#   out
# }
# ggplot.quantiles.dmc <- function (pp, aform="S+R", lower=0.025, upper= 0.975) {
#   if (names(pp)[1]== "proportions") {sim <- single.qq(pp[["quantiles"]]); dat<- single.qq(pp[["data.quantiles"]])} else {
#     sim <- group.qq (pp)  [["sim"]]
#     dat <- group.qq (pp)  [["data"]]
#   }
#   quantile.list <- list()
#   ind<-  grep("rep", names(sim))
#   n.quants <- length(sim) - ind
#   for (i in 1:n.quants) {
#     form <- as.formula(paste(paste(names(sim)[ind+i], "~"), aform))
#     mean.df <- aggregate (form, data= sim, FUN= function(x) c(mean(x), quantile(x, probs= lower), quantile(x, probs= upper)))
#     data.df <- aggregate (form, data= dat, FUN= "mean");names(data.df)[length(data.df)] <- "data"
#     plot.df <- do.call(data.frame, join(mean.df, data.df))
#     names(plot.df)[(length(plot.df)-n.quants ): (length(plot.df)-1)] <- c("mean", "lower", "upper")
#     plot.df$quantile <- names(sim)[ind+i]
#     quantile.list[[i]] <- plot.df ; names(quantile.list) [i] <- names(sim)[ind+i]
#   }
#   do.call(rbind, quantile.list)
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
