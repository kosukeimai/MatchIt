plot.matchit.subclass <- function(x, discrete.cutoff=5,type="QQ", ...){
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
    choices <- c("No",paste("Yes : Subclass ", 1:max(x$subclass,na.rm=T)))
    question <- "Would you like to see density estimates of any subclass covariates?"
    ans <- -1
    while(ans!=0)
      {
        ans <- as.numeric(choice.menu(choices,question))
        if(ans!=0)
          {
            matchit.qqplot(x,discrete.cutoff,which.subclass=ans)     
          }
      }
  } else if(type=="jitter"){
    jitter.pscore(x)
  } else {
    stop("Invalid type")
  }
}
