plot.matchit <- function(x, ...)
  {
   if (identical(eval(x$call$exact),TRUE)){
    return(warning("Plot() not appropriate for exact matching.  No plots generated"))
    }
    match.matrix <- x$match.matrix
    xdata <- x$data
    treata <- model.frame(x$formula,xdata)[,1,drop=FALSE]
    pscore <- xdata[,"pscore"]
    in.sample <- x$in.sample
    covariates <- x$covariates
    treat <- as.vector(treata[,1])
    names(treat) <- row.names(treata)
   covariates <- model.frame(delete.response(terms(x$formula)),xdata)[,,drop=FALSE]
    mahvars <- eval(x$call$mahvars)
    exact <- eval(x$call$exact)

 if (!is.null(mahvars)){
      w <- mahvars%in%names(covariates)
      if(sum(w)!=length(mahvars)){
      md <- as.data.frame(as.matrix(xdata[,mahvars[!w]]))
      names(md) <- mahvars[!w]
      covariates <- cbind(covariates,md)
    }
  }
  if (!identical(TRUE, exact) & !identical(FALSE, exact)) {
    w <- exact%in%names(covariates)
    if(sum(w)!=length(exact)) {
      ed <- as.data.frame(as.matrix(xdata[,exact[!w]]))
      names(ed) <- exact[!w]
      covariates <- cbind(covariates,ed)
    }
  }
    psclass <- x$psclass
    if(is.null(psclass)){subclass <- 0}else{subclass <- qbins <- max(na.omit(x$psclass))}
    nearest <- !is.null(x$match.matrix)
    q.cut <- x$q.cut
    n <- length(treat)
    n0 <- length(treat[treat==0])
    n1 <- length(treat[treat==1])

    if(is.null(names(treat))){names(treat) <- seq(1,n)}
    labels <- names(treat)
    clabels <- labels[treat==0]
    tlabels <- labels[treat==1]
   if(is.null(in.sample)){in.sample <- rep(TRUE,n); names(in.sample) <- names(treat)}
    if(nearest==TRUE) {
      match.matrix <- match.matrix[in.sample[treat==1],,drop=F]
      num.matches <- dim(match.matrix)[2]-apply(as.matrix(match.matrix), 1, function(x) { sum(is.na(x)) })
      names(num.matches) <- tlabels[in.sample[treat==1]]
      t.units <- row.names(match.matrix)[num.matches>0]
      c.units <- na.omit(as.vector(as.matrix(match.matrix)))
      matched <-c(t.units,unique(c.units))
      weights <- rep(0,length(treat))
      names(weights) <- labels
      weights[t.units] <- 1
      
      for (cont in clabels) {
        treats <- na.omit(row.names(match.matrix)[cont==match.matrix[,1]])
        if (dim(match.matrix)[2]>1) {
          for (j in 2:dim(match.matrix)[2]) 
            treats <- c(na.omit(row.names(match.matrix)[cont==match.matrix[,j]]),treats)	
        }
	if (length(treats)==0) weights[cont] <- 0
	else for (k in 1:length(treats)) weights[cont] <- weights[cont]+1/num.matches[treats[k]] 
      }

    } else {
      	matched <- names(treat)
	weights <- rep(1,length(treat))
	names(weights) <- names(treat)
	}
    
    weights[!in.sample] <- 0 

    doverlay <- function(x,treat,weights=NULL,xlab="",main="",lines=FALSE){
      xmiss <- !is.na(x)
      x <- x[xmiss]
      treat <- treat[xmiss]
      weights <- weights[xmiss]
      minobs <- min(x)
      maxobs <- max(x)
      dx1 <- density(x[treat==1],from=minobs,to=maxobs)
      dx0 <- density(x[treat==0],from=minobs,to=maxobs)
      if(!is.null(weights))
        {
          x1 <- x[treat==1&weights!=0]
          x0 <- sample(x[treat==0], size=min(10000,(100*length(x[treat==0]))),
                       replace=TRUE,prob=weights[treat==0]/sum(weights[treat==0]))
          d1 <- density(x1,from=minobs,to=maxobs)
          bw <- d1$bw 
          d0 <- density(x0,from=minobs,to=maxobs,bw=bw)
          par(mfrow=c(2,1))
          matplot(dx0$x,cbind(dx1$y,dx0$y),type="l",ylim=range(c(dx0$y,dx1$y,d1$y,d0$y)),
                  ylab="Density",xlab=xlab,main=paste(main,": All Units",sep=""))
          legend(minobs,max(c(d1$y,d0$y,dx1$y,dx0$y)), lty=1:2, col=1:2,
                 legend=c("Treatment","Control"))
          matplot(dx0$x,cbind(d1$y,d0$y),type="l",ylim=range(c(dx0$y,dx1$y,d1$y,d0$y)),
                  ylab="Density",xlab=xlab,main=paste(main,": Matched Units",sep=""))
          legend(minobs,max(c(d1$y,d0$y,dx1$y,dx0$y)), lty=1:2, col=1:2,
                 legend=c("Treatment","Control"))
          if (lines==T) abline(v=q.cut,col="grey",lty=1)
          par(mfrow=c(1,1))
        } else{
          matplot(dx0$x,cbind(dx1$y,dx0$y),type="l",ylab="Density",xlab=xlab,main=main)
          legend(minobs,max(c(dx1$y,dx0$y)),lty=1:2,col=1:2,legend=c("Treatment","Control"))
	  if (lines==T) abline(v=q.cut, col="grey", lty=1)
        }
    }
    
    choice.menu <- function(choices,question)
      {
        k <- length(choices)-1
        Choices <- data.frame(choices)
        row.names(Choices) <- 0:k
        names(Choices) <- "Choices"
        print.data.frame(Choices,right=FALSE)
        ans <- readline(question)          
        while(!ans%in%c(0:k))
          {
            print("Not valid -- please pick one of the choices")
            print.data.frame(Choices,right=FALSE)
            ans <- readline(question)
          }
        return(ans)
      }
    if(!is.null(pscore))
      {
        densq <- ("Would you like to see density estimates of the propensity scores?")
        denschoice <- c("No","Yes")
        densplot <- choice.menu(denschoice,densq)
        if(densplot==1){
          if(nearest){
            doverlay(pscore,treat,weights,main="Propensity Score",lines=T)
            #if(length(subclass)!=1 | subclass!=0){abline(v=q.cut,col="grey",lty=1)}
	  } else {
              doverlay(pscore,treat,main="Propensity Score", lines=T) 
              #if(length(subclass)!=1 | subclass!=0){abline(v=q.cut,col="grey",lty=1)}
          }
        }
      
        jitq <- ("Would you like to see a jitterplot of the propensity scores?")
        jitchoice <- c("No","Yes")
        jitplot <- choice.menu(jitchoice,jitq)
        if(jitplot==1) {
          jitp <- jitter(rep(1,length(treat)),factor=6)-(treat==0)
          cwt <- sqrt(weights)
          plot(pscore,xlim=range(na.omit(pscore)),ylim=c(-1,2),type="n",ylab="",xlab="Propensity Score",axes=F,main="Distribution of Propensity Scores")
          if(length(subclass)!=1 | subclass!=0){abline(v=q.cut,col="grey",lty=1)}
          points(pscore[treat==1&weights!=0],jitp[treat==1&weights!=0],pch=18,cex=cwt[treat==1&weights!=0])
          points(pscore[treat==1&weights==0],jitp[treat==1&weights==0],pch=5,col="grey",cex=0.5)
          points(pscore[treat==0&weights!=0],jitp[treat==0&weights!=0],pch=18,cex=cwt[treat==0&weights!=0])
          points(pscore[treat==0&weights==0],jitp[treat==0&weights==0],pch=5,col="grey",cex=0.5)
          axis(1)
          text(sum(range(na.omit(pscore)))/2,1.5,"Treatment Units")
          text(sum(range(na.omit(pscore)))/2,-0.5,"Control Units")
          box()
          print("To identify the units, use first mouse button; to stop, use second.")
          identify(pscore,jitp,names(treat))
        }
      }

    if(nearest){
      choices <- c("No",paste("Yes : ",names(covariates)))
      question <- "Would you like to see density estimates of any other covariates?"
      ans <- -1
      while(ans!=0)
        {
          ans <- as.numeric(choice.menu(choices,question))
          if(ans!=0)
            if(nearest)
              {
                doverlay(covariates[,(ans)],treat,weights,main=names(covariates)[ans])
              } else {
                doverlay(covariates[,(ans)],treat,main=names(covariates)[ans])
              }
        }
    }
    if(length(subclass)!=1 | subclass!=0){
      choices <- c("No",paste("Yes : Subclass ", 1:qbins))
      question <- "Would you like to see density estimates of any subclass covariates?"
      ans <- -1
      while(ans!=0)
        {
          ans <- as.numeric(choice.menu(choices,question))
          if(ans!=0)
            {
              question2 <- "Which covariates?"
              choices2 <- c(paste(names(covariates)))
              ans2 <- as.numeric(choice.menu(choices2,question2))
              ans2 <- ans2+1
              if(sum(treat[psclass==ans]==1,na.rm=TRUE)<=2){
                print("Not enough treatment units in this subclass")
              } else if(sum(treat[psclass==ans]==0,na.rm=TRUE)<=2){
                print("Not enough control units in this subclass")} else{
                  doverlay(covariates[,(ans2)][psclass==ans],treat[psclass==ans],main=paste("Subclass ",ans," : ",names(covariates)[ans2]))}
            }
        }
      
    } 
    
  }

