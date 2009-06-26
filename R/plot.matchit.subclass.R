plot.matchit.subclass <- function(x, discrete.cutoff=5,
                                  type="QQ", interactive = T,
                                  subclass = NULL, which.xs=NULL,...){
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
  if(type=="QQ"){
    if(interactive){
      choices <- c("No",paste("Yes : Subclass ", 1:max(x$subclass,na.rm=T)))
      question <- "Would you like to see quantile-quantile plots of any subclasses?"
      ans <- -1
      while(ans!=0)
        {
          ans <- as.numeric(choice.menu(choices,question))
          if(ans!=0)
            {
              matchit.qqplot(x,discrete.cutoff,which.subclass=ans,
                             interactive = interactive, which.xs=which.xs,...)     
            }
        }
    } else {
      matchit.qqplot(x,discrete.cutoff,which.subclass=subclass,
                     interactive=interactive, which.xs=which.xs,...)
    }
  } else if(type=="jitter"){
    jitter.pscore(x, interactive=interactive,...)
  } else if(type=="hist"){
    hist.pscore(x,...)
  } else {
    stop("Invalid type")
  }
}
