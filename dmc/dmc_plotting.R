# System functions for the DMC (Dynamic Models of Choice)
#    Functions specific to generating graphical output
#    Usually user does not need to edit


plot.cell.density <- function(data.cell,C=NA,xlim=c(0,Inf),ymax=NA,
                              save.density=FALSE,digits=2,main="",show.mean=FALSE)
  # If !is.na(C) plots density for correct and error responses for a data frame
  # with columns R (a factor) and RT, adding a boolean score column for R=C.
  # Otherwise plots each response. Can deal with NA in the RT column, in which
  # case it provides a summary of p(NA)
{
  if (!is.factor(data.cell$R)) data.cell$R <- factor(data.cell$R)
  if (length(C)==1) C <- rep(C,dim(data.cell)[1])
  p.na <- mean(is.na(data.cell$RT))
  is.in <- !is.na(data.cell$RT)
  is.in[is.in] <- data.cell$RT[is.in]>xlim[1] & data.cell$RT[is.in]<xlim[2]
  if ( !any(is.na(C)) ) {
    if ( is.logical(C) & length(C)==dim(data.cell)[1] )
      dat$C <- C[is.in] else dat$C <- dat$R==C[is.in]
    if (length(dat$RT[dat$C])>2)
      dns.correct <- density(dat$RT[dat$C]) else dns.correct <- NULL
    if (length(dat$RT[!dat$C])>2)
      dns.error <- density(dat$RT[!dat$C]) else dns.error <- NULL
    if (is.null(dns.error) & is.null(dns.correct))
      stop("There are no densities to plot")
    acc <- mean(dat$C)
    if (!is.null(dns.correct))
      dns.correct$y <- dns.correct$y*acc*(1-p.na)
    if (!is.null(dns.error))
      dns.error$y <- dns.error$y*(1-acc)*(1-p.na)
    if (is.na(ymax)) ymax <- max(c(dns.correct$y,dns.error$y))
    if (!is.null(dns.correct)) {
      plot(dns.correct,xlab="RT",ylab="density",ylim=c(0,ymax),main=main)
      if (!is.null(dns.error)) lines(dns.error,col="red")
    } else {
      plot(dns.error,xlab="RT",ylab="density",ylim=c(0,ymax),main=main,col="red")
      if (!is.null(dns.correct)) lines(dns.correct)
    }
    nams <- "Accuracy ="
    ps <- round(acc,digits)
    if (p.na!=0) {
      nams <- c(nams,"p(NA) =")
      ps <- c(ps,round(p.na,digits))
    }
    legend("topright",paste(nams,ps),bty="n")
    legend("topright",xjust=0, inset=c(0,0.1), c("correct","error"),bty="n", lty=c(1,1), col=c("black","red"))
    if ( save.density ) list(correct=dns.correct,error=dns.error)
  } else {
    rs <- levels(dat$R)
    dns <- vector(mode="list",length=length(rs))
    names(dns) <- rs
    ps <- table(dat$R)/dim(dat)[1]
    for (i in rs) {
      rt <- dat$RT[dat$R==i]
      if (length(rt)>2) {
        dns[[i]] <- density(rt)
        dns[[i]]$y <- ps[i]*dns[[i]]$y*(1-p.na)
      }
    }
    ymax <- suppressWarnings(max(unlist(lapply(dns,function(x){max(x$y)}))))
    no.dns <- unlist(lapply(dns,is.null))
    if (all(no.dns))
      stop("There are no densities to plot!")
    dns1 <- dns[!no.dns]
    ltys <- c(1:length(dns1))
    plot(dns1[[1]],xlab="RT",ylab="density",ylim=c(0,ymax),lty=ltys[1],
         main=main)
    if (length(dns1)>1) for (i in 2:length(dns1)) lines(dns1[[i]],lty=ltys[i])
    nams <- paste("p(",names(dns1),") =",sep="")
    ps <- round(ps[!no.dns]*(1-p.na),digits)
    if ( p.na!=0 ) {
      nams <- c(nams,"p(NA) =")
      ps <- c(ps,round(p.na,2))
      lty <- c(ltys,NA)
    }
    legend("topright",paste(nams,ps),lty=ltys,bty="n")
    if ( save.density ) dns
  }
}


profile.dmc <- function(p.name,min.p,max.p,p.vector,data,
                        n.point=100,digits=2,ylim=NA)
  # for parameter p.name in p.vector draws the profile likelihood for data and
  # returns the maximum (on a grid of resolution n.point)
{
  if (!(p.name %in% names(p.vector)))
    stop("p.name not in p.vector")
  p <- p.vector
  ps <- seq(min.p,max.p,length.out=n.point)
  ll <- numeric(n.point)
  for (i in 1:n.point)
  {
    p[p.name] <- ps[i]
    ll[i] <- sum(log(likelihood.dmc(p,data)))
  }
  names(ll) <- round(ps,digits)
  if (any(is.na(ylim)))
    plot(ps,ll,type="l",xlab=p.name,ylab="log-likelihood") else
      plot(ps,ll,type="l",xlab=p.name,ylab="log-likelihood",ylim=ylim)
  ll[ll==max(ll)]
}



# xlim=NA;natural=TRUE;n.point=1e3;trans=NA;main=""; trans=NA
plot.prior <- function(i,p.prior,xlim=NA,natural=TRUE,n.point=1e3,trans=NA,
  main="",ylim=NA, line = FALSE, ...)
  # Plot the i'th member of list created by p.prior. If trans then plot
  # on natural scale using transform specified in attr(p.prior[[i]],"trans")
{
  if ( attr(p.prior[[i]],"dist")!="constant" ) {
    if ( any(is.na(trans)) )
      if (natural) trans <- attr(p.prior[[i]],"untrans") else
        trans <- "identity"
    if ( is.numeric(i) ) i <- names(p.prior)[i]
    if ( !(i %in% names(p.prior)) ) stop("Parameter not in prior")
    p <- p.prior[[i]]
    if ( any(is.na(xlim)) ) {
      xlim <- switch(attr(p,"dist"),
                     tnorm=c(pmax(p$lower,p[[1]]-3*p[[2]]),
                             pmin(p$upper,p[[1]]+3*p[[2]])),
                     beta_lu=c(p$lower,p$upper),
                     gamma_l=c(p$lower,p[[1]]*p[[2]]+3*sqrt(p[[1]])*p[[2]]),
                     lnorm_l=c(p$lower,exp(p[[1]]+2*p[[2]]))
      )
    }

    x <- seq(xlim[1], xlim[2], length.out = n.point)

    ### Code to calcualte Jacobian of transform by Quentin Gronau

    # Tranform function
    eval(parse(text = paste0("f <- function(x)", trans, "(x)")))
    xx <- f(x) # Trasformed x

    # Redefine x on finite support
    finite_index <- is.finite(xx)
    x_ok <- x[finite_index]
    x <- seq(x_ok[1], x_ok[length(x_ok)], length.out = n.point)
    xx <- f(x)

    # Make p
    p$x <- x
    p$log <- FALSE

    # Get derivative
    fprime <- try(deriv(parse(text = paste0(trans, "(x)")), "x", func = TRUE),
      silent = TRUE)

    # Calcualte density including normalizing Jacobian if transformed.
    if ( class(fprime) == "try-error" )
      y <- 1/abs( sapply(x, function(x) numDeriv::grad(f, x)) ) *
        do.call(paste("d",attr(p,"dist"),sep=""),p) else
      y <- 1/abs( attr(fprime(x), "gradient") ) *
            do.call(paste("d",attr(p,"dist"),sep=""),p)

    # Plot
    if (any(is.na(ylim))) ylim <- c(0,max(y))

    if (line) lines(xx,y,...) else
      plot(xx,y,type="l",xlab=i,ylab="density",main=main,ylim=ylim,...)

    invisible(cbind(x=xx,y=y))
  }
}



# hyper=FALSE;start=1;end=NA;thin=1;save.ll=FALSE
# main.pll=NULL;pll.chain=FALSE;pll.together=TRUE;pll.barplot=FALSE
# only.prior=FALSE;only.like=FALSE;subject=1;layout=NA
# smooth=FALSE;density=FALSE
plot.dmc <- function(samples,hyper=FALSE,start=1,end=NA,thin=1,save.ll=FALSE,
                     main.pll=NULL,pll.chain=FALSE,pll.together=TRUE,pll.barplot=FALSE,
                     only.prior=FALSE,only.like=FALSE,subject=1,layout=NA,
                     smooth=FALSE,density=FALSE,...)
  # If not hyper, by default changes samples$theta to a mcmc.list object and plots
  # If pll.chain=TRUE changes samples$pll to an mcmc object and plots posterior
  # log-likelihood of chains. If pll.barplot=TRUE plots their means as a barplot.
  # only.prior and only.like pick out one or other component of pll.
  # If hyper same but at the the hyper  (phi) level.
{

  my.chain.plot <- function(chain.pll,thin=1,main.pll=NULL) {
    x <- 1:dim(chain.pll)[1]
    x <- x[x %% thin == 0]
    ylim <- c(min(chain.pll),max(chain.pll))
    plot(x,chain.pll[x,1],xlab="Iterations",ylab="",ylim=ylim,type="l",main=main.pll)
    for (i in 2:dim(chain.pll)[2]) lines(x,chain.pll[x,i],col=i,lty=i)
  }


  auto.layout <- any(is.na(layout))
  if ( hyper ) {
    hyper <- attr(samples,"hyper")
    if (is.null(hyper))
      stop("There are no hyper-parameters to plot.")
    if ( is.na(end) ) end <- hyper$nmc
    if ( end <= start )
      stop("End must be greater than start")
    if ( pll.chain | pll.barplot ) {
      chain.pll <- hyper$h_summed_log_prior[start:end,] +
                   hyper$h_log_likelihoods[start:end,]
      colnames(chain.pll) <- 1:dim(chain.pll)[2]
      if (pll.barplot) {
        mean.ll <- apply(chain.pll,2,mean)
        names(mean.ll) <- 1:length(mean.ll)
        barplot(mean.ll,ylab="Mean Post Like",main=main.pll)
        if (save.ll) mean.ll
      } else {
        if (!auto.layout) par(mfrow=layout)
        if (!pll.together)
          plot(mcmc(chain.pll,thin=thin),auto.layout=auto.layout,
            density=density,smooth=smooth,...) else
            plot(mcmc.list(lapply(data.frame(chain.pll),function(x){
              mcmc(x,thin=thin)})),auto.layout=auto.layout,
              density=density,smooth=smooth,...)
      }
    } else {
      if (!auto.layout) par(mfrow=layout)
      plot(window(phi.as.mcmc.list(hyper,start=start,end=end,thin=1)),
           auto.layout=auto.layout,density=density,smooth=smooth,...)
    }
  } else {
    if (is.null(samples$theta)) samples <- samples[[subject]]
    if ( is.na(end) ) end <- samples$nmc
    if ( end <= start ) stop("End must be greater than start")

    if ( pll.chain | pll.barplot ) {
      if ( only.prior ) chain.pll <- samples$summed_log_prior[start:end,] else
        if ( only.like ) chain.pll <- samples$log_likelihoods[start:end,] else
          chain.pll <- samples$summed_log_prior[start:end,] +
                       samples$log_likelihoods[start:end,]
      colnames(chain.pll) <- 1:dim(chain.pll)[2]

      if ( pll.barplot ) {
        mean.ll <- apply(chain.pll,2,mean)
        names(mean.ll) <- 1:length(mean.ll)
        barplot(mean.ll,ylab="Mean Post Like",main=main.pll)
        if (save.ll) mean.ll
      } else {
        if ( !auto.layout ) par(mfrow=layout)
        if ( !pll.together )
          plot(mcmc(chain.pll,thin=thin),auto.layout=auto.layout,
            density=density,smooth=smooth,main.pll=main.pll) else {
            if (!density & !smooth) my.chain.plot(chain.pll,thin=thin,main.pll=main.pll,...) else
            plot(mcmc.list(lapply(data.frame(chain.pll),function(x){
              mcmc(x)})),auto.layout=auto.layout,
              density=density,smooth=smooth,...)
          }
      }
    } else {
      if (!auto.layout) par(mfrow=layout)
      plot(window(theta.as.mcmc.list(samples,thin=thin),start=start,end=end),
           auto.layout=auto.layout,density=density,smooth=smooth,...)
    }
  }
}



ppl.barplots.dmc <- function(samples,start=1,end=NA,layout=c(5,2))
  # Grid of barplots of pll for set of subjects
{
  if (!is.null(samples$theta))
    stop("This function cannot be applied to a single subject.")
  par(mfrow=layout)
  for (i in 1:length(samples))
    plot.dmc(samples,subject=i,pll.barplot=TRUE,start=start,end=end,
             main.pll=paste("Subject",i))
}



# style="pdf"; layout=NULL; pos=NULL; dname="Data";mname="Model"; model.legend=TRUE; x.min.max=NA
# percentiles=c(10,30,50,70,90); ylim=NA; show.fit=TRUE;show.fits=TRUE; split.responses=FALSE
# show.response=NA; ltys=NA;lwds=NA
# model.legend=FALSE; no.layout=FALSE
plot.pp.dmc <- function(pp,style="pdf",layout=NULL,no.layout=FALSE,
  pos="topleft",percentiles=c(10,30,50,70,90),
                        dname="Data",mname="Model",model.legend=TRUE,
                        show.fit=TRUE,show.fits=TRUE,x.min.max=NA,ylim=NA,
                        show.response=NA,
                        ltys=NA,lwds=NA)
  # pdf or cdf of data and posterior predictive fits
  # Show response is a vector of integers in 1:n.resp of which resposnes to show,
  # if NA shows all
{

  plot.one <- function(pp,style,data.style,pos,x.min.max=NA,ylim,
                       show.fits,dname,mname,model.legend,lgnd,
                       n.resp,show.response,ltys,lwds)
  {

    no.facs <- class(try(pp[[data.style]][1][[1]][[2]],silent=TRUE))=="try-error"
    if (no.facs) ni <- 1 else ni <- length(pp[[1]])
    for ( i in 1:ni ) if ( !all(unlist(lapply(pp[[1]][[i]],is.null))) ) {
      if ( no.facs )
        lgnd <- names(lapply(pp[[style]],function(x){out <- attr(x,"cell.name")})) else
        lgnd <- lapply(pp[[style]][[i]],function(x){out <- attr(x,"cell.name")})
      ok.n <- c(1:n.resp)[!unlist(lapply(pp[[1]][[i]],is.null))]
      is.in <- rep(FALSE,n.resp)
      is.in[show.response] <- TRUE
      is.in <- is.in & ((1:n.resp)==ok.n)
      inj <- c(1:n.resp)[is.in]
      if ( style=="pdf" ) {
        if (is.null(pos)) pos <- "topright"
        if (any(is.na(ylim))) ylim <- c(0,max(
          c(unlist(lapply(pp[[style]][[i]],function(x){
              if (length(x$y)==0) NA else max(x$y)})),
            unlist(lapply(pp[[data.style]][[i]],function(x){
              if (length(x$y)==0) NA else max(x$y)}))),na.rm=TRUE))
        xlim <- c(Inf,-Inf)
        for (j in 1:length(ok.n)) {
          if (no.facs)
            xij <- pp[[style]][[ok.n[j]]][[1]]$x else
            xij <- pp[[style]][i][[1]][[ok.n[j]]]$x
          xlim[1] <- ifelse(min(xij)<xlim[1],
                            min(xij),xlim[1])
          xlim[2] <- ifelse(max(xij)>xlim[2],
                            max(xij),xlim[2])
        }
        if ( !any(is.na(x.min.max)) ) {
          if (x.min.max[1]>xlim[1]) xlim[1] <- x.min.max[1]
          if (x.min.max[2]<xlim[2]) xlim[2] <- x.min.max[2]
        }
        for (j in show.response) {
          if (no.facs) {
            xy <- pp[[style]][[ok.n[j]]][[1]]
            xy.dat <- pp[[data.style]][[ok.n[j]]][[1]]
          } else {
            xy <- pp[[style]][i][[1]][[ok.n[j]]]
            xy.dat <- pp[[data.style]][i][[1]][[ok.n[j]]]
          }
          if ( j == show.response[1] )
            plot(xy,main="",ylim=ylim,xlim=xlim) else
              lines(xy,lty=ltys[j],lwd=lwds[j])
          lines(xy.dat,lty=ltys[j],lwd=lwds[j],col="red")
        }
      } else {
        if (is.null(pos)) pos <- "topleft"
        if (any(is.na(ylim))) ylim <- c(0,1)
        # Get limits
        xlim <- c(Inf,-Inf)
        for ( j in inj ) {
          if (no.facs) ppj <- pp[[data.style]][[ j ]][[1]] else
            ppj <- pp[[data.style]][i][[1]][[ j ]]
          xlim[1] <- ifelse(suppressWarnings(min(ppj))<xlim[1],min(ppj),xlim[1])
          xlim[2] <- ifelse(suppressWarnings(max(ppj))>xlim[2],max(ppj),xlim[2])
          if ( show.fit ) {
            if (no.facs) ppj <- pp[[style]][[ j ]][[1]] else
              ppj <- pp[[style]][i][[1]][[ j ]]
            xlim[1] <- ifelse(suppressWarnings(min(ppj))<xlim[1],min(ppj),xlim[1])
            xlim[2] <- ifelse(suppressWarnings(max(ppj))>xlim[2],max(ppj),xlim[2])
            if (!is.null(attr(pp,"dpqs"))) for ( k in 1:length(attr(pp,"dpqs")) ) {
              if (no.facs) ppk <- attr(pp,"dpqs")[[k]][[style]][[ j ]][[1]] else
                ppk <- attr(pp,"dpqs")[[k]][[style]][i][[1]][[ j ]]
              xlim[1] <- ifelse(suppressWarnings(min(ppk))<xlim[1],min(ppk),xlim[1])
              xlim[2] <- ifelse(suppressWarnings(max(ppk))>xlim[2],max(ppk),xlim[2])
            }
          }
        }
        if (!any(is.na(x.min.max))) {
          if (x.min.max[1]>xlim[1]) xlim[1] <- x.min.max[1]
          if (x.min.max[2]<xlim[2]) xlim[2] <- x.min.max[2]
        }
        for (j in inj) {
          if (j == inj[1]) plot(NA,NA,main="",type="l",ylim=ylim,xlim=xlim)
          # Draw uncertianty
          if ( show.fits && !is.null(attr(pp,"dpqs")) ) for (k in 1:length(attr(pp,"dpqs"))) {
            ppk <- attr(pp,"dpqs")[[k]]
            if (no.facs) {
              y <- as.numeric(names(ppk[[style]][[ j ]][[1]]))
              x <- ppk[[style]][[ j ]][[1]]
            } else {
              y <- as.numeric(names(ppk[[style]][i][[1]][[ j ]]))
              x <- ppk[[style]][i][[1]][[ j ]]
            }
            lines(x,y,col="grey",lty=ltys[j])
          }
          if (no.facs) ppj <- pp[[data.style]][[ j ]][[1]] else
            ppj <- pp[[data.style]][i][[1]][[ j ]]
          lines(ppj,as.numeric(names(ppj)),col="red",lty=ltys[j],lwd=lwds[j])
          if ( !any(is.na(percentiles)) & !is.null(ppj) )
            points(ppj[percentiles],as.numeric(names(ppj))[percentiles],col="red",cex=1.25)
          if ( show.fit ) {
            if (no.facs) ppj <- pp[[style]][[ j ]][[1]] else
              ppj <- pp[[style]][i][[1]][[ j ]] # ok.n[j]
            lines(ppj,as.numeric(names(ppj)),lty=ltys[j],lwd=lwds[j])
            if ( !any(is.na(percentiles)) & !is.null(ppj) )
              points(ppj[percentiles],as.numeric(names(ppj))[percentiles],pch=16,cex=.75)
          }
        }
      }
      if ( !is.na(pos) ) {
        if ( model.legend ) legend(pos,c(paste(dname,lgnd)[is.in],paste(mname,lgnd)[is.in]),
          lty=c(ltys[show.response],ltys[show.response]),
          lwd=c(lwds[show.response],lwds[show.response]),bty="n",
          col=c(rep("red", length(show.response)),rep("black",length(show.response))) ) else
        legend(pos,paste(dname,lgnd)[is.in],lty=ltys[show.response],lwd=lwds[show.response],bty="n",
            col=rep("red", length(show.response)) )

      }
    }
  }

  if (!show.fit) {
    model.legend <- FALSE
    show.fits <- FALSE
  }

  if ( !all(names(pp)[1:2]==c("pdf","cdf")) ) {
    style <- "cdf"
    dname <- paste("Average",dname)
    mname <- paste("Average",mname)
    pp <- attr(pp,"av")
  }

  if (!any(style %in% c("pdf","cdf")))
    stop("style must be either \"pdf\" or \"cdf\"")
  facs <- dimnames(pp[[1]])
  n.facs <- length(facs)
  resp <- names(pp[[1]][1][[1]])
  if (is.null(resp)) resp <- names(pp[[1]])
  n.resp <- length(resp)
  n.levs <- lapply(facs,length)
  if ( !no.layout ) {
    if ( is.null(layout) ) {
      if (class(try(pp[[data.style]][1][[1]][[2]],silent=TRUE))=="try-error")
        layout <- 1 else layout <- unlist(n.levs)
    }
    if (length(layout)==1) layout <- c(1,layout) else
      layout <- layout[1:2]
    par(mfrow=layout,mar=c(2,2,0,0))
  }
  lgnd <- NULL
  data.style <- paste("data",style,sep=".")
  if (any(is.na(show.response))) show.response <- 1:n.resp
  if (any(is.na(ltys))) ltys <- c(1:n.resp)
  if (any(is.na(lwds))) lwds <- rep(1,n.resp)
  if ( !(length(ltys)==n.resp) ) stop(paste("Must provide",n.resp,"ltys"))
  if (!(length(lwds)==n.resp)) stop(paste("Must provide",n.resp,"lwds"))

  plot.one(pp,style,data.style,pos,x.min.max,ylim,show.fits,
                  dname,mname,model.legend,lgnd,n.resp,show.response,ltys,lwds)
}


acf.dmc <- function(samples,par=NA,chain=1,lag.max=NULL,hyper=FALSE,plot=TRUE,...)
  # Plot the autocorrelation function of a chain of parameters, or of likelihood
  # if par=NA
{
  if (hyper) {
    hyper <- attr(samples,"hyper")
    if (is.na(par))
      out <- acf(window(hyper$h_log_likelihoods[,chain],...),
        main=chain,lag.max=lag.max,plot=plot) else {
        phi <- phi.as.mcmc.list(hyper)
        if (is.numeric(par)) par <- dimnames(phi[[1]])[[2]][par]
        series <- phi[[chain]][,par]
        out <- acf(window(series,...),main=par,lag.max=lag.max,plot=plot)
      }
  } else {
    if (is.na(par))
      out <- acf(window(samples$log_likelihoods[,chain],...),
        main=chain,lag.max=lag.max,plot=plot) else {
        if (is.numeric(par)) par <- dimnames(samples$theta)[[2]][par]
        series <- samples$theta[chain,par,]
        out <- acf(window(series,...),main=par,lag.max=lag.max,plot=plot)
      }
  }
  invisible(out)
}


# start=1;end=NA;hyper=FALSE;hyper1=TRUE;hyper2=TRUE;thin=NA;
# collapse.subjects=FALSE;scale.subjects=TRUE; do.plot=TRUE
pairs.dmc <- function(samples,start=1,end=NA,hyper=FALSE,hyper1=TRUE,hyper2=TRUE,thin=NA,
  collapse.subjects=FALSE,scale.subjects=TRUE,
  do.plot=TRUE,...)
  ## put histograms on the diagonal
{

  panel.hist <- function(x, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
  }

  ## put correlations on the upper panels,
  panel.cor <- function(x, y, digits = 2, prefix = "", ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y)
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    text(0.5, 0.5, txt, cex = 2)
  }

  if (hyper) {
    hyper <- attr(samples,"hyper")
    dim1 <- dim(hyper$phi[[1]])
    dim2 <- dim(hyper$phi[[2]])
    nam1 <- dimnames(hyper$phi[[1]])[[2]]
    nam2 <- dimnames(hyper$phi[[2]])[[2]]
    if (hyper1 & hyper2) {
      samples$theta <- aperm(array(c(aperm(hyper$phi[[1]],c(1,3,2)),
                       aperm(hyper$phi[[2]],c(1,3,2))),
        dim=c(dim1[c(1,3)],dim1[2]+dim2[2]),dimnames=list(NULL,NULL,
        c(paste(nam1,"h1",sep="."),paste(nam1,"h2",sep=".")))),c(1,3,2))
    } else if (hyper1) {
      samples$theta <- aperm(array(c(aperm(hyper$phi[[1]],c(1,3,2))),
        dim=c(dim1[c(1,3)],dim1[2]),dimnames=list(NULL,NULL,
        c(paste(nam1,"h1",sep=".")))),c(1,3,2))
    } else if (hyper2) {
      samples$theta <- aperm(array(c(aperm(hyper$phi[[2]],c(1,3,2))),
        dim=c(dim1[c(1,3)],dim2[2]),dimnames=list(NULL,NULL,
        c(paste(nam2,"h2",sep=".")))),c(1,3,2))
    } else stop("At least one of hyper1 and hyper2 may be true.")
  }
  if ( !any(names(samples)=="theta") ) if (collapse.subjects) {
    thetas <- lapply(samples,function(x){x$theta})
    dimall <- dim(thetas[[1]])
    i <- dimall[3]
    n <- length(samples)
    dimall[3] <- i*n
    theta <- array(dim=dimall,dimnames=dimnames(thetas[[1]]))
    for (j in 1:n) {
      if ( scale.subjects ) for (k in 1:dim(thetas[[1]])[2] )
        thetas[[j]][,k,] <- (thetas[[j]][,k,]-mean(thetas[[j]][,k,]))/sd(thetas[[j]][,k,])
      theta[,,( 1+ (j-1)*i ):(j*i)] <- thetas[[j]]
    }
    samples$theta <- theta
  } else stop("If samples not an individual subject must have hyper=TRUE or collapse.subjects=TRUE")
  if (is.na(end)) end <- dim(samples$theta)[3]
  theta <- aperm(samples$theta[,,start:end],c(1,3,2))
  if (!is.na(thin))
    theta <- theta[,c(1:dim(theta)[2])[(c(1:dim(theta)[2]) %% thin) == 0],]
  theta <- matrix(theta,ncol=dim(theta)[[3]],
                  dimnames=list(NULL,dimnames(theta)[[3]]))
  if (do.plot) pairs(theta,diag.panel = panel.hist,upper.panel = panel.cor,...)
  invisible(cor(theta))
}


plot.deviance.dmc <- function(ds=NULL,samples=NULL,digits=2,fast=TRUE,
                              main="",xlim=NA)
  # Posterior deviance histogram
{
  if (is.null(ds)) if (is.null(samples))
    stop("Must supply samples or deviance statistics") else
      ds <- Dstats.dmc(samples,TRUE,fast)
    if (any(is.na(xlim)))
      hist(ds$D,breaks="fd",xlab="Posterior Deviance",main=main) else
        hist(ds$D,breaks="fd",xlab="Posterior Deviance",main=main,xlim=xlim)
    abline(v=ds$Dmean,lty=2)
    abline(v=ds$meanD,lty=3)
    pds <- pd.dmc(ds)
    legend("topright",c(paste("NP =",ds$np),
                        paste("NPmean =",round(pds$Pmean,digits)),
                        paste("NPmin =",round(pds$Pmin,digits)),
                        paste("NPvar =",round(pds$Pvar,digits)),
                        "D(mean theta)","mean(D)"),
           bty="n",lty=c(NA,NA,NA,NA,2:3))
}


plot.score.dmc <- function(data=NULL, rnd=2, xlim=c(0,5), ymax=NA, IQR=FALSE)

  {
  if (is.null(data))
    stop("Must supply the model-data instance")

  slevs <- sort(unique(data$S))
  rlevs <- sort(unique(data$R))
  if (!all(tolower(slevs)==tolower(rlevs)))
    stop(paste("S levels (",paste(slevs,collapse=","),
                ") cant be matched with R levels (",paste(slevs,collapse=","),")"))
  correct <- tolower(data$S)==tolower(data$R)
  print(round(mean(correct),rnd))
  print(round(tapply(data$RT,list(correct),mean),rnd))
  plot.cell.density(data,C=correct,xlim=xlim)

  if (IQR)
    print(round(tapply(data$RT,list(correct),IQR),2)) # standard deviation
}


# fun=function(x){mean(x$RT)}; n.post=500; verbose=TRUE
# plot.density=TRUE;main="";bw="nrd0";report=10;SSD=Inf
ppp.dmc <- function(samples,
                    fun=function(x){mean(x$RT)},n.post=500,
                    transformp=identity,verbose=TRUE,
                    plot.density=TRUE,main="",bw="nrd0",report=10,SSD=Inf,...)
  # posterior predictive (pp) p value for function fun of data (p(observed)>pp)
{

  model <- attributes(samples$data)$model
  facs <- names(attr(model,"factors"))
  resp <- names(attr(model,"responses"))
  ns <- table(samples$data[,facs],dnn=facs)
  n.par <- dim(samples$theta)[2]
  thetas <- matrix(aperm(samples$theta,c(1,3,2)),ncol=n.par)
  colnames(thetas) <- dimnames(samples$theta)[[2]]
  posts <- thetas[sample(c(1:dim(thetas)[1]),n.post,replace=F),]
  sim <- vector(mode="list",length=n.post)
  if (verbose) {
    cat("\n")
    cat(paste("Simulating (\'.\'=",report,"): ",sep=""))
  }
  for (i in 1:n.post) {
    sim[[i]] <- simulate.dmc(transformp(posts[i,]),model,ns,SSD=SSD)
    if ( verbose && ((i %% report) == 0)) cat(".")
  }
  if (verbose) cat("\n")
  pp <- unlist(lapply(sim,fun,...))
  obs <- fun(samples$data)
  ppp <- mean(obs>pp)
  if (plot.density) {
    plot(density(pp,bw=bw),main=main)
    abline(v=obs,lty=2)
  }
  if (verbose) print(ppp)
  attr(pp,"observed") <- obs
  invisible(pp)
}

# fun=function(x){mean(x$RT)};transformp=identity
# n.post=500;save.pps=FALSE;verbose=TRUE; collapse="mean"
# plot.density=TRUE;main="";bw="nrd0";report=10;SSD=Inf
h.ppp.dmc <- function(hsamples,collapse="mean",
  fun=function(x){mean(x$RT)},transformp=identity,
  n.post=500,save.pps=FALSE,verbose=TRUE,
  plot.density=TRUE,main="",bw="nrd0",report=10,SSD=Inf,...)
  # ppp for average over subjects
{
  pps <- lapply(hsamples,ppp.dmc,fun=fun,transformp=transformp,
    n.post=n.post,plot.density=FALSE,SSD=SSD,verbose=verbose,...)
  obs <- do.call(collapse,list(unlist(lapply(pps,function(x){attr(x,"obs")}))))
  pp <- apply(do.call(rbind,pps),2,collapse)
  ppp <- mean(obs>pp)
  if (plot.density) {
    plot(density(pp,bw=bw),main=main)
    abline(v=obs,lty=2)
  }
  print(ppp)
  attr(pp,"observed") <- obs
  if (save.pps) attr(pp,"pps") <- pps
  invisible(pp)
}


### Fixed Effects

plotSpar.dmc <- function(est,p.name,len=.05)
  # plots ordered cis for p.name, est produced by summary.dmc
{
  lo <- unlist(lapply(est,function(x){x$quantiles[p.name,c("2.5%")]}))
  med <- unlist(lapply(est,function(x){x$quantiles[p.name,c("50%")]}))
  hi <- unlist(lapply(est,function(x){x$quantiles[p.name,c("97.5%")]}))
  ord <- order(med)
  n <- length(med)
  plot(1:n,med[ord],ylim=c(min(lo),max(hi)),ylab=p.name,xlab="Subject",pch=16,
    main="Medians and 95% credible intervals")
  arrows(1:n,med[ord],1:n,lo[ord],angle=90,len=len)
  arrows(1:n,med[ord],1:n,hi[ord],angle=90,len=len)
}


### Hierachical

h.profile.dmc <- function(p.name,p.num,min.p,max.p,ps,p.prior,n.point=100,
                          digits=3,ylim=NA)
  # for parameter p.name at position p.num (1 or 2) in p.prior given subject
  # parameters ps draws a likelihood profile and returns the maximum (on a grid
  # of resolution n.point)
{
  if (!(p.name %in% dimnames(ps)[[2]]))
    stop("p.name not in ps")
  if (!(p.num %in% 1:2))
    stop("p.num must be either 1 or 2")
  pps <- seq(min.p,max.p,length.out=n.point)
  ll <- numeric(n.point)
  pop <- p.prior[[p.name]]
  pop$x <- ps[,p.name]
  for (i in 1:n.point)
  {
    pop[[p.num]] <- pps[i]
    ll[i] <- sum(do.call(paste("d",attr(pop,"dist"),sep=""),pop))
  }
  names(ll) <- round(pps,digits)
  p.name <- paste(p.name,names(pop)[p.num],sep=".")
  if (any(is.na(ylim)))
    plot(pps,ll,type="l",xlab=paste(p.name,p.num,sep="."),
         ylab="log-likelihood") else
           plot(pps,ll,type="l",xlab=paste(p.name,p.num,sep="."),
                ylab="log-likelihood",ylim=ylim)
  ll[ll==max(ll)]
}


############  LUKES GGPLOTs

theme_simple <- function (base_size = 12, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.background = element_rect(fill="white"),
      panel.grid.minor.y = element_blank(),
      legend.key = element_rect(fill="white", colour= "white"),
      strip.background = element_rect(fill="white"))
}

# theme_set(theme_simple())

###NOTE: Changed post.median back to mean to accord w. old way of doing it

# sim=sim; data=data; factors=NULL; CI= c(0.025, 0.975)
# quantiles.to.get = c(0.1, 0.5, 0.9); noR = FALSE

get.fitgglist.dmc <- function (sim, data, factors=NA, noR = FALSE,  
                               quantiles.to.get = c(0.1, 0.5, 0.9), CI= c(0.025, 0.975),
                               acc.fun=function(x){as.numeric(x$S)==as.numeric(x$R)},
                               correct.only=FALSE,error.only=FALSE, custom.RTfun=NA, custom.name ="NA",
                               hPP=NA)  
  #Can also request any other function of RT (will bug with NAs unless user includes an na.rm)
  # Extracts list of data frames, pps (response proabilities) and RTs
  # from the save.simulation output of post.predict
{
  
  ##This function is to turn an array (output of tapply) into a data frame.    
  arr2df=function(arr) 
  {
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
  
  #gets the response (R) probabilities 
  #for each cell of the data frame as determined by factors
  get.ps <- function (sim, factors,R,include.na=FALSE,only.na=FALSE) 
  {
    n <- tapply(sim$RT,cbind.data.frame(sim[,factors],R=R),length)
    n[is.na(n)] <- 0 # In case some cells are empty
    nok <- tapply(sim$RT,cbind.data.frame(sim[,factors],R=R),function(x){sum(!is.na(x))})
    nok[is.na(nok)] <- 0 # In case some cells are empty
    
    # In the case that there are no factors, the apply doesn't work, and in this 
    # case we can just sum all the Rs instead.
    if(is.null(factors)) np <- sum(n) else 
      np <- rep(apply(n,1:length(factors),sum),times=length(levels(R)))
    
    if (only.na) (n-nok)/np else if (include.na) n/np else nok/np
  }  
  
  #Takes the output of tapplies to sim (array with reps), as well as vector
  #output of tapply to data, and creates a data frame with posterior median/CIs
  
  get.df <- function(simarray, dataarray, tapplyvec, len, post.mean=TRUE){
    ####Here a switch to detect when the avg is not an avg by checking sim for name 's'
    if (any (names(dimnames(simarray))=='s')){ 
      simarray <- apply(simarray,(1:length(dim(simarray)))[names(dimnames(simarray))!='s'], mean, na.rm=T)
      dataarray <- apply(dataarray,(1:length(dim(dataarray)))[names(dimnames(dataarray))!='s'], mean, na.rm=T)}
    datavec <- as.vector(dataarray)
    #Can switch this to use the posterior mean rather than median  
    if (!post.mean) summary.df <- arr2df(apply(simarray,2:len,quantile,probs=.5, na.rm=TRUE)) else{
      summary.df <- arr2df(apply(simarray,2:len,mean, na.rm=TRUE))
    }
    summary.df$lower <- as.vector(apply(simarray,2:len,quantile,probs=CI[1], na.rm=TRUE))
    summary.df$upper <- as.vector(apply(simarray,2:len,quantile,probs=CI[2], na.rm=TRUE))
    summary.df$data <- datavec
    names(summary.df)[names(summary.df)=="y"] <- "median"
    ##If there is only one factor it will end up in the row names - reverse this
    #and turn it back into a column
    if (dim(summary.df)[2]==4) {
      summary.df <- cbind.data.frame(row.names(summary.df),summary.df)
      names(summary.df)[1] <- tapplyvec[2]  
    }
    summary.df
  }
  
  ###applies a function to RT then makes a data frame with data/sim avg 
  get.df.function <- function(fun,data,sim,tapplyvec) {
    simarray <- tapply(sim$RT,sim[ ,tapplyvec],fun)
    dataarray <- tapply(data$RT,data[ ,tapplyvec[-1]],fun)
    get.df(simarray, dataarray, tapplyvec,len=len)
  }
  
  
  #new switch - in order to average in a different order
  if(!is.na(hPP)[1]){
    #gets participant index and attached to simmed data, then creates group data frame includign s names
    sim <- do.call(rbind,lapply(seq_along(hPP), function(y, n,i) { cbind(n[[i]], y[[i]]) }, y=hPP, n=names(hPP)))
    colnames(sim)[1] <- 's'
    data <- lapply(hPP, function(x) attr(x, "data"))
    data <- do.call(rbind,lapply(seq_along(data), function(y, n,i) { cbind(n[[i]], y[[i]]) }, y=data, n=names(hPP)))
    colnames(data)[1] <- 's'
  }  
  
  
  #### Use acc.fun to score df, correct.only, or error.only to filter df 
  if (correct.only & error.only)
    stop("Cant plot only correct and only error, set only one to true")
  
  if (is.function(acc.fun)) {
    try(C.sim <- acc.fun(sim),silent=TRUE)
    try(C.data <- acc.fun(data),silent=TRUE)
    if (class(C.sim)=="try-error" | class(C.data)=="try-error")
    {
      warning("Accuracy function could not score properly") 
      C.sim <- sim[,"R"]
      C.data <- data[,"R"]
      scored <- FALSE
    } else {
      C.sim <- factor(C.sim)
      C.data <- factor(C.data)
      scored <- TRUE
    }
  } else {
    C.sim <- sim[,"R"]
    C.data <- data[,"R"]
    scored <- FALSE
  }
  
  
  if ( is.null(factors) & noR ) stop ("Cannot plot when no factors and noR TRUE")
  if ( !is.null(factors) && is.na(factors[1]) ) 
    factors <-  colnames(sim) [!colnames (sim) %in% c("reps", "R", "RT","R2")] 
  if (!noR) {
    C.sim <- sim[,"R"]
    C.data <- data[,"R"]
  } 
  
  ####Here define a 'tapplyvec', which will decide the factors to use in 
  #the data summaries. 
  # Also get 'len' which is the length of factor vector + "reps" (posterior simulations)
  # + R (if there is an R factor)
  # Drop the response factor from RT calculations if noR == T
  if (noR == FALSE) { 
    len <- length(factors) +2 
    tapplyvec <- c("reps", factors,"R") 
  } else { 
    if (noR == TRUE) { 
      len <- length(factors) +1
      tapplyvec <- c("reps", factors)  
    }
  }
  
  #If subject number is part of the df change the len as s number will be subtracted out later
  if (any (names(sim)=='s')) len = len-1
  
  #Get various response probabilities
  #as noR is not relevant to response prob get another len called plen which is length ignoring noR
  if (noR) plen= len+1 else plen=len
  
  # get ps With non-responses ignored
  ps <- get.ps(sim, factors=c("reps", factors),R=C.sim)
  ps[is.nan(ps)] <- 0
  dps <- get.ps(data, factors,R=C.data)
  dps[is.nan(dps)] <- 0
  pp.df <- get.df(ps, dps , tapplyvec, len=plen)
  
  # get ps Including non-responses
  ps.na <- get.ps(sim, factors=c("reps", factors),R=C.sim,include.na=TRUE)
  ps.na[is.nan(ps.na)] <- 0
  dps.na <- get.ps(data, factors,R=C.data,include.na=TRUE)
  dps.na[is.nan(dps.na)] <- 0
  pp.df.na <- get.df(ps.na,dps.na, tapplyvec, len=plen)
  
  # get Only non-responses
  ps.pna <- get.ps(sim, factors=c("reps", factors),R=C.sim,only.na=TRUE)
  ps.pna[is.nan(ps.pna)] <- 0
  dps.pna <- get.ps(data, factors,R=C.data,only.na=TRUE)
  dps.pna[is.nan(dps.pna)] <- 0
  pp.df.pna <- get.df(ps.pna, dps.pna,  tapplyvec, len=plen)
  
  ###Filter RTs by correct/error if requested  
  if (scored) {
    if (correct.only) {
      sim <- sim[C.sim=="TRUE",]
      data <- data[C.data=="TRUE",]
    }
    if (error.only) {
      sim <- sim[C.sim!="TRUE",]
      data <- data[C.data!="TRUE",]
    }
  }
  
  #Summarise quantile RTs. To save time, calculate them all at once.
  # First calc all quantiles of the RT distribution for each rep
  all.quants <- tapply(sim$RT,sim[ ,tapplyvec],quantile, prob=quantiles.to.get,na.rm=TRUE)
  DIM <- dim(all.quants)
  DIMnames <- dimnames(all.quants)  
  all.quants <- lapply(all.quants, function(x) as.numeric(as.character(x)))
  
  # Loop through specified quantiles of the RT distribution, for each 
  # calculating posterior mean + CI then bind to data frame. 
  for (i in 1:length(quantiles.to.get)) {
    quant.array <- unlist(lapply(all.quants, function(x) x[i]))
    dim(quant.array) <-  DIM
    dimnames(quant.array) <- DIMnames 
    quant.df <- get.df(quant.array, tapply(data$RT, data[, c(tapplyvec[tapplyvec!="reps"])], 
                                           quantile, prob=quantiles.to.get[i],na.rm=TRUE), tapplyvec,
                       len=len)
    quant.df$quantile <- as.character(quantiles.to.get[i])
    if (i==1) RT.df <- quant.df  else  RT.df <- rbind(RT.df, quant.df) 
    
  }
  names(RT.df)[names(RT.df)=="y"] <- "median"
  
  
  if (!is.function(custom.RTfun)) {out <- list(pp.df, pp.df.na, pp.df.pna,RT.df,get.df.function(function(x) mean(x, na.rm=T),data,sim,tapplyvec),
                                               get.df.function(function(x) sd(x, na.rm=T),data,sim,tapplyvec))
  names(out) <- c("pps","pps.NR","NRpps","RTs","MeanRTs","SDRTs")} else {
    #If any special function of RT was requested, return that as well
    out <- list(pp.df, pp.df.na, pp.df.pna,RT.df,get.df.function(function(x) mean(x, na.rm=T),data,sim,tapplyvec),
                get.df.function(function(x) sd(x, na.rm=T),data,sim,tapplyvec), get.df.function(custom.RTfun,data,sim,tapplyvec))
    names(out) <- c("pps","pps.NR","NRpps","RTs","MeanRTs","SDRTs", custom.name)
  }
  out
}


ggplot.RP.dmc <- function (df, xaxis = 'R', panels.ppage=4, nrow=NULL, ncol=NULL)
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
      geom_point(size=3) + geom_errorbar(aes(ymax = upper, ymin = lower), width= 0.2) +
      geom_point(aes_string(x = xaxis, y= 'data'), pch=21, size=3, colour="black") +
      geom_line(aes(group = 1, y=data), linetype=2) + ylab("Response Proportion")

    if (length(grid)!=0)
      plot <- plot + facet_wrap(grid, labeller = grid_labels, scales = 'free_x',
                                nrow=nrow, ncol=ncol)
    plots[[j]] <- plot
  }

  if (length(plots)==1) plots[[1]] + theme_simple() else lapply(plots, function(x) x + theme_simple())
}

ggplot.RA.dmc <- function(df, xaxis = 'R', panels.ppage=4,ylim=NA,df.type="pps",
                          nrow=NULL, ncol=NULL)
  # Plots accuracy and returns df with only correct response rows, dropping
  # R column.
{
  if (!is.null(attr(df, "gglist"))) {
    cat ("Treating as a pp object for a single participant\n")
    df <- attr(df, "gglist")[[df.type]]
  }

  if (!is.null(attr(attr(df, "av"), "gglist"))) {
    cat ("Treating as a pp object for a group of participants\n")
    df <- attr(attr(df, "av"), "gglist")[[df.type]]
  }
  acc.df <- df[df$R=="TRUE",names(df)[names(df)!="R"]]
  if (xaxis=='R') xaxis <- names(acc.df)[1]
  if ( any(is.na(ylim)) )
    ylim <- c(min(c(acc.df$lower,acc.df$data)),max(c(acc.df$upper,acc.df$data)))
  out <-  ggplot.RP.dmc(acc.df, xaxis=xaxis,panels.ppage=panels.ppage,
                        nrow=nrow, ncol=ncol)
  #If else to account for single page (single ggplot)
  # vs multipage (list of ggplots)
  if (is.null(names(out))) {
    out <- lapply (out, function(x) x + ylab("accuracy") + ylim(ylim)  ) } else {
    out <- out+ ylab("accuracy") + ylim(ylim)
  }
  out
}

# df<-pp
# xaxis='S'
# panels.ppage=4; do.quantiles=TRUE
#                            nrow=NULLncol=NULL

ggplot.RT.dmc <- function (df, xaxis = 'R', panels.ppage=4, do.quantiles=TRUE,
                           nrow=NULL, ncol=NULL)
 
  # xaxis is the name of the factor to put on the X axis
{
 
  if (!is.null(attr(df, "gglist"))) {
    cat ("Treating as a pp object for a single participant")
    df <- attr(df, "gglist")[["RTs"]]
  }

  if (!is.null(attr(attr(df, "av"), "gglist"))) {
    cat ("Treating as a pp object for a group of participants")
    df <- attr(attr(df, "av"), "gglist")[["RTs"]]
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

    plot <- ggplot(active.df, aes_string(x = xaxis, y= 'median')) + geom_point(size=3) +
        geom_errorbar(aes(ymax = upper, ymin = lower), width= 0.2) +
        geom_point(aes_string(x = xaxis, y= 'data'), pch=21, size=3, colour="black") +
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

### Stop-signal plots: select design cell by removing irrelevant data, plots all
###                    values in data and sim (post.predict.dmc save) objects.


# minN=1;ci=.99;sim.line=FALSE; cut.on.data=TRUE
# main="Median signal-respond RTs";do.correct=FALSE;isC=NULL; only.ssrt=TRUE
# plot.it=TRUE;violin=0.05;ylab="Median SRRT (s)";xlab=NA;ylim=NA
# probs=NA;n.intervals=5;low.interval=NA; percentile.av=TRUE; xlabels=NULL
# do.correct=FALSE;Cfun=function(x){as.numeric(x$S)==(as.numeric(x$R)-1)}

plot_SS_srrt.dmc <- function(data,sim=NULL,minN=1,ci=.99,sim.line=FALSE,
  main="Median signal-respond RTs",xlab=NA,
  do.correct=FALSE,Cfun=function(x){as.numeric(x$S)==(as.numeric(x$R)-1)},
  plot.it=TRUE,violin=0.05,ylab="Median SRRT (s)",ylim=NA,only.ssrt=TRUE,
  probs=NA,n.intervals=5,low.interval=NA,percentile.av=FALSE,cut.on.data=TRUE,xlabels=NULL)
  # Plots median signal-respond RT function, optionally with violin plots and
  # invisibly returns a data frame with number of observations and Gelman p
  # values for each SSD. sim is save from post.predict.dmc. Like data, must have
  # an SSD and RT column. Recommend at least 500 reps in sim object.
  # NB1: p for n=1 is meaningless, so p=NA in such cases.
  # NB2: n.sim is number of simulated experiments where the median was defined,
  #      be cautious about p value and violins for SSDs where it is small,
  # NB3: Does not plot violins if violin=0 or FALSE, otherwise absolute value
  #      of violin gives "wex", width of density plot.
  # NB4: sim.line joins violin medians with dotted line.
  # NB5: If probs=NA then n.intervals is used to create this number of evenly
  #      spaced SSD intervals when doing an average plot. low.interval specifies
  #      a SSD = 0-low.inteval ranges to be removed first.
  # NB6: percentile.av calculates intervals for each subject based on simulation (data
  #      to unstable) then averages. Otherwise does absolute SSD based on all
  #      data (so can do for sim=NULL case). In either case can base on data SSDs
  #      (cut.on.data=TRUE, maybe unstable for percentile.av) or simulation (FALSE)

{
  fit <- !is.null(sim)
  if (is.na(xlab)) xlab <- ifelse(is.data.frame(data),"SSD (s)","SSD")
  if ( !is.data.frame(data) ) { # average over subjects
    av.plot <- TRUE
    snames <- names(data)
    if ( fit && !is.list(sim) )
      stop("If sim is for multiple subjects data must be a samples object") else
    {
      # Make sure breaks are unique
      data <- lapply(data,function(x) {
        if (only.ssrt) out <- x$data[!is.na(x$data$RT),] else out <- x$data
        if (!is.na(low.interval)) ok <- out$SSD>low.interval else
          ok <- rep(TRUE,dim(out)[1])
        out$SSD[ok] <- out$SSD[ok] + (runif(length(out$SSD[ok]))-.5)/1e6
        out
      })
      sim <- lapply(sim,function(x) {
        if (only.ssrt) x <- x[!is.na(x$RT),]
        if (!is.na(low.interval)) ok <- x$SSD>low.interval else
          ok <- rep(TRUE,dim(x)[1])
         x$SSD[ok] <- x$SSD[ok] + (runif(length(x$SSD[ok]))-.5)/1e6
         x
      })
    }
    if ( do.correct ) {
      data <- lapply(data,function(x){x[Cfun(x) | x$R=="NR",]})
      if (fit) sim <- lapply(sim,function(x){x[Cfun(x) | x$R=="NR",]})
    }
    if ( !fit & cut.on.data==FALSE )
        stop("sim must be supplied if cut.on.data=FALSE")
    if ( any(is.na(probs)) )
      probs <- seq(0,1,1/n.intervals)[-c(1,n.intervals+1)]

    if ( percentile.av ) {
      absSSD <- unlist(lapply(data,function(x){x$SSD[is.finite(x$SSD)]}))
      for ( i in snames ) {
        data[[i]] <- data[[i]][is.finite(data[[i]]$SSD),]
        if (fit) sim[[i]] <- sim[[i]][is.finite(sim[[i]]$SSD),]
        if ( cut.on.data ) {
          if ( is.na(low.interval) ) SSDok <- data[[i]]$SSD else
            SSDok <- data[[i]]$SSD[data[[i]]$SSD>low.interval]
        } else {
          if ( is.na(low.interval) ) SSDok <- sim[[i]]$SSD else
            SSDok <- sim[[i]]$SSD[sim[[i]]$SSD>low.interval]
        }
        # Make sure breaks are unique
        cuts <- c(quantile(SSDok,probs),max(SSDok))
        if ( is.na(low.interval) ) cuts <- c(0,cuts) else
          cuts <- c(0,low.interval,cuts)
        if (fit) sim[[i]]$SSD <- as.numeric(cut(sim[[i]]$SSD,cuts))
        data[[i]]$SSD <- as.numeric(cut(data[[i]]$SSD,cuts))
      }
      data <- do.call(rbind,data)
      if (fit) sim <- do.call(rbind,sim)
      SSDnames <- round(c(0,probs,1)*100)
      SSDnames <- paste("(",
        paste(SSDnames[-length(SSDnames)],SSDnames[-1],sep=","),
      "%]",sep="")
      if ( !is.na(low.interval) ) SSDnames <-
        c(paste(paste("[0",low.interval,sep=","),"]",sep=""),SSDnames)
      data$SSD <- factor(data$SSD,labels=SSDnames)
      if (fit) sim$SSD <- factor(sim$SSD,labels=SSDnames)
      absSSD <- tapply(absSSD,data$SSD,mean)
      cat("Average absolute SSD for each percentile range\n")
      print(round(absSSD,3))
      cat("\n")
    } else {

      data <- do.call(rbind,data)
      data <- data[is.finite(data$SSD),]
      if (fit) {
        sim <- do.call(rbind,sim)
        sim <- sim[is.finite(sim$SSD),]
      }

      if ( cut.on.data ) {
        if ( is.na(low.interval) ) SSDok <- data$SSD else
          SSDok <- data$SSD[data$SSD>low.interval]
      } else {
        if (is.na(low.interval)) SSDok <- sim$SSD else
          SSDok <- sim$SSD[sim$SSD>low.interval]
      }
      cuts <- c(quantile(SSDok,probs),max(SSDok))
      if (is.na(low.interval)) cuts <- c(0,cuts) else
        cuts <- c(0,low.interval,cuts)
      if (fit) sim$SSD <- cut(sim$SSD,cuts)
      data$SSD <- cut(data$SSD,cuts)
    }
  } else  {
    av.plot <- FALSE
    if (only.ssrt) {
      if (fit) sim <- sim[is.finite(sim$SSD) & !is.na(sim$RT),]
       data <- data[is.finite(data$SSD) & !is.na(data$RT),]
    } else {
      if (fit) sim <- sim[is.finite(sim$SSD),]
       data <- data[is.finite(data$SSD),]
    }
    if (do.correct) data <- data[Cfun(data) | data$R=="NR",]
  }

  # Prepare SSDs
  nrt <- tapply(!is.na(data$RT),data[,c("SS","SSD")],sum)[2,]
  SSDs <- names(nrt[nrt>=minN])

  # Remove small SSD
  tmp <- factor(data$SSD)
  levels(tmp)[!(levels(tmp) %in% SSDs)] <- NA
  data <- data[!is.na(tmp),]
  nrt <- nrt[nrt>=minN]

  # Calculate data summaries
  mrt <- tapply(data$RT,data$SSD,median,na.rm=TRUE)
  ok <- !is.na(mrt)

  # Get graph x axis
  if ( av.plot ) {
    if (is.null(xlabels)) xlabs <- levels(data$SSD) else
      if (length(xlabels) == length(levels(data$SSD))) xlabs=xlabels else
        stop(paste("xlabels not the right length (",length(levels(data$SSD))),")")
    xvals <- 1:length(xlabs)
    names(xvals) <- levels(data$SSD)
  } else {
    SSD <- as.numeric(SSDs)
    names(SSD) <- SSDs
    xvals <- SSD
    if (is.null(xlabels)) xlabs <- SSD else
      if (length(xlabels) == length(SSD)) xlabs <- xlabels else
        stop(paste("xlabels not the right length (",length(SSD),")"))
  }

  if ( plot.it & !fit ) { # Signal respond plot
    if (any(is.na(ylim))) ylim <- c(min(mrt[ok])-0.01,max(mrt[ok])+0.01)
    plot(xvals[ok],mrt[ok],ylim=ylim,pch=16,xlab=xlab,ylab=ylab,main=main,xaxt="n")
    axis(1,xvals[ok],xlabs[ok])
    points(xvals[ok],mrt[ok],pch=16)
    lines(xvals[ok],mrt[ok],lwd=2)
  }

  if ( fit ) {

    tmp <- factor(sim$SSD);levels(tmp)[!(levels(tmp) %in% SSDs)] <- NA
    sim <- sim[!is.na(tmp),]

    # Calcualte simulation summaries
    mrt.sim <- tapply(sim$RT,sim[,c("reps","SSD")],median,na.rm=TRUE)

    # Simulation ci
    mrt.lo <- apply(mrt.sim,2,quantile,probs=(1-ci)/2,na.rm=TRUE)
    mrt.hi <- apply(mrt.sim,2,quantile,probs=1-(1-ci)/2,na.rm=TRUE)

    if ( plot.it ) { # Signal respond plot
      if (any(is.na(ylim)) & fit) ylim <- c(min(c(mrt.lo,mrt[ok]),na.rm=TRUE),
        max(c(mrt.hi,mrt[ok]),na.rm=TRUE))

      plot(xvals[ok],mrt[ok],ylim=ylim,pch=16,xlab=xlab,ylab=ylab,main=main,xaxt="n")
      axis(1,xvals[ok],xlabs[ok])

      if ( as.logical(abs(violin)) ) for (i in SSDs[ok]) {
        vioplot(mrt.sim[!is.na(mrt.sim[,i]),i],at=xvals[i],add=TRUE,
                col="lightgrey",wex=abs(violin))
      }
      if (sim.line) lines(xvals[ok],apply(mrt.sim,2,median,na.rm=TRUE)[ok],lty=3)
    }
    points(xvals[ok],mrt[ok],pch=16,cex=1.5,col="red")
    lines(xvals[ok],mrt[ok],lwd=2)

    # Gelman p values
    p <- round(apply((mrt.sim-rep(mrt,each=dim(mrt.sim)[1]))>=0,2,mean,na.rm=TRUE),3)
    n.sim <- apply(mrt.sim,2,function(x){sum(!is.na(x))})
    out <- cbind.data.frame(SSD=xlabs[ok],nrt=nrt[ok],p=p[ok],n.sim=n.sim[ok])
    out$p[out$n==1] <- NA
    row.names(out) <- NULL
    cat(paste(main,"\n"))
    print(out)
  }
}



# minN=1;plot.it=TRUE;violin=.05
# ylab="p(Respond)";xlab="SSD (s)";ylim=NA;sim.line=FALSE;main="Inhibition function"
# probs=NA;n.intervals=5;low.interval=NA; cut.on.data=FALSE; percentile.av=TRUE
plot_SS_if.dmc <- function(data,sim=NULL,minN=1,plot.it=TRUE,violin=.05,
  ylab="p(Respond)",xlab=NA,ylim=NA,sim.line=FALSE,main="Inhibition function",
  probs=NA,n.intervals=5,low.interval=NA,percentile.av=TRUE,cut.on.data=TRUE,xlabels=NULL)
  # Plots inhibiiton function, optionally wit violin plots and invisibly
  # returns a data frame with number of observaitons and Gelman p values for
  # each SSD. sim is save from post.predict.dmc. Like data, must have an SSD and
  # RT column. Recommend at east 500 reps in sim object.
  # NB1: p for n=1 is meaningless, so p=NA in such cases)
  # NB2: Does not plot violins if violin=0 or FALSE, otherwise absolute value
  #      of violin gives "wex", width of density plot.)
  # NB3: sim.line joins violin medians with dotted line.
  # NB4: Other options similar to SSRT plot
{
  fit <- !is.null(sim)
  if (is.na(xlab)) xlab <- ifelse(is.data.frame(data),"SSD (s)","SSD")
  if ( !is.data.frame(data) ) { # average over subjects
    av.plot <- TRUE
    snames <- names(data)
    if ( fit && !is.list(sim) )
      stop("If sim is for multiple subjects data must be a samples object") else {
      # Make sure breaks are unique
      data <- lapply(data,function(x) {
        out <- x$data
        if (!is.na(low.interval)) ok <- out$SSD>low.interval else
          ok <- rep(TRUE,dim(out)[1])
        out$SSD[ok] <- out$SSD[ok] + (runif(length(out$SSD[ok]))-.5)/1e6
        out
      })
      sim <- lapply(sim,function(x) {
        if (!is.na(low.interval)) ok <- x$SSD>low.interval else
          ok <- rep(TRUE,dim(x)[1])
         x$SSD[ok] <- x$SSD[ok] + (runif(length(x$SSD[ok]))-.5)/1e6
         x
      })
    }
    if ( !fit & cut.on.data==FALSE )
        stop("sim must be supplied if cut.on.data=FALSE")
    if ( any(is.na(probs)) )
      probs <- seq(0,1,1/n.intervals)[-c(1,n.intervals+1)]

    if ( percentile.av ) {
      absSSD <- unlist(lapply(data,function(x){x$SSD[is.finite(x$SSD)]}))
      for ( i in snames ) {
        data[[i]] <- data[[i]][is.finite(data[[i]]$SSD),]
        if (fit) sim[[i]] <- sim[[i]][is.finite(sim[[i]]$SSD),]
        if (cut.on.data) {
          if (is.na(low.interval)) SSDok <- data[[i]]$SSD else
            SSDok <- data[[i]]$SSD[data[[i]]$SSD>low.interval]
        } else {
          if (is.na(low.interval)) SSDok <- sim[[i]]$SSD else
            SSDok <- sim[[i]]$SSD[sim[[i]]$SSD>low.interval]
        }
        cuts <- c(quantile(SSDok,probs),max(SSDok))
        if ( is.na(low.interval) ) cuts <- c(0,cuts) else
          cuts <- c(0,low.interval,cuts)
        if (fit) sim[[i]]$SSD <- as.numeric(cut(sim[[i]]$SSD,cuts))
        data[[i]]$SSD <- as.numeric(cut(data[[i]]$SSD,cuts))
      }
      data <- do.call(rbind,data)
      if (fit) sim <- do.call(rbind,sim)
      SSDnames <- round(c(0,probs,1)*100)
      SSDnames <- paste("(",
        paste(SSDnames[-length(SSDnames)],SSDnames[-1],sep=","),
      "%]",sep="")
      if ( !is.na(low.interval) ) SSDnames <-
        c(paste(paste("[0",low.interval,sep=","),"]",sep=""),SSDnames)
      data$SSD <- factor(data$SSD,labels=SSDnames)
      if (fit) sim$SSD <- factor(sim$SSD,labels=SSDnames)
      absSSD <- tapply(absSSD,data$SSD,mean)
      cat("Average absolute SSD for each percentile range\n")
      print(round(absSSD,3))
      cat("\n")
    } else {
      data <- do.call(rbind,data)
      data <- data[is.finite(data$SSD),]
      if (fit) {
        sim <- do.call(rbind,sim)
        sim <- sim[is.finite(sim$SSD),]
      }
      if ( cut.on.data ) {
        if ( is.na(low.interval) ) SSDok <- data$SSD else
          SSDok <- data$SSD[data$SSD>low.interval]
      } else {
        if (is.na(low.interval)) SSDok <- sim$SSD else
          SSDok <- sim$SSD[sim$SSD>low.interval]
      }
      cuts <- c(quantile(SSDok,probs),max(SSDok))
      if (is.na(low.interval)) cuts <- c(0,cuts) else
        cuts <- c(0,low.interval,cuts)
      if (fit) sim$SSD <- cut(sim$SSD,cuts)
      data$SSD <- cut(data$SSD,cuts)
    }
  } else  {
    av.plot <- FALSE
    if (fit) sim <- sim[is.finite(sim$SSD),]
    data <- data[is.finite(data$SSD),]
  }


  # Prepare SSDs
  n <- tapply(data$RT,data$SSD,length)
  SSDs <- names(n[n>=minN])

  # Remove small SSD
  tmp <- factor(data$SSD)
  levels(tmp)[!(levels(tmp) %in% SSDs)] <- NA
  data <- data[!is.na(tmp),]
  n <- n[n>=minN]

  # Calculate data summaries
  pR <- tapply(!is.na(data$RT),data$SSD,mean) # probability of responding

  # Get graph x axis
  if ( av.plot ) {
    if (is.null(xlabels)) xlabs <- levels(data$SSD) else
      if (length(xlabels) == length(levels(data$SSD))) xlabs=xlabels else
        stop(paste("xlabels not the right length (",length(levels(data$SSD))),")")
    xvals <- 1:length(xlabs)
    names(xvals) <- levels(data$SSD)
  } else {
    SSD <- as.numeric(SSDs)
    names(SSD) <- SSDs
    xvals <- SSD
    if (is.null(xlabels)) xlabs <- SSD else
      if (length(xlabels) == length(SSD)) xlabs=xlabels else
        stop(paste("xlabels not the right length (",length(SSD),")"))
  }

  if ( !fit & plot.it ) { # Inhibition Function
    if (any(is.na(ylim))) ylim <- c(0,1)
    plot(xvals,pR,ylim=ylim,pch=16,xlab=xlab,ylab=ylab,main=main,xaxt="n")
    axis(1,xvals,xlabs)
    points(xvals,pR,pch=16)
    lines(xvals,pR)
  }

  if ( fit ) {
    tmp <- factor(sim$SSD)
    levels(tmp)[!(levels(tmp) %in% SSDs)] <- NA
    sim <- sim[!is.na(tmp),]
    # Calcualte simulation summaries
    pR.sim <- tapply(!is.na(sim$RT),sim[,c("reps","SSD")],mean)

    if ( plot.it ) { # Inhibition Function
      if (any(is.na(ylim))) ylim <- c(0,1)
      plot(xvals,pR,ylim=ylim,pch=16,xlab=xlab,ylab=ylab,main=main,xaxt="n")
      axis(1,xvals,xlabs)
      if ( as.logical(abs(violin)) ) for (i in SSDs) {
        vioplot(pR.sim[,i],at=xvals[i],add=TRUE,col="grey",wex=abs(violin))
      }
      if (sim.line) lines(xvals,apply(pR.sim,2,median,na.rm=TRUE),lty=3)
    }
    points(xvals,pR,pch=16,cex=1.5,col="red")
    lines(xvals,pR)

    # Gelman p values
    p <- round(apply((pR.sim-rep(pR,each=dim(pR.sim)[1]))>=0,2,mean,na.rm=TRUE),3)
    out <- cbind.data.frame(SSD=xlabs,n=n,p=p)
    out$p[out$n==1] <- NA
    row.names(out) <- NULL
    cat(paste(main,"\n"))
    print(out)
  }
}


