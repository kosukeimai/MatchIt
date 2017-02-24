

in_frame <- function(df, col) {
  exists(col, envir= as.environment(df))
}

#----------------------------------------------------------
### Inheritance
#----------------------------------------------------------

#' @title Checks matchit Class
#' @description Function that checks if the target object is a \code{matchit} object.
#' @param object any R object
#' @return Returns \code{TRUE} if its argument has class "matchit" among its classes and
#' \code{FALSE} otherwise.
#' @export
is.matchit <- function(object) {
  inherits(object, "matchit")
}


#----------------------------------------------------------
### GET MATCHES 
#----------------------------------------------------------

#' @export
get_matches <- function(object, model_frame, id_cols= NULL, newdata= NULL) {
  UseMethod("get_pairs", object)
}

#' @title Get matches from matchit object
#' @description Get the resulting matches from a \code{matchit} model object. This function allows the 
#' user to extract the matches from the original dataset used in model building or from a new dataset
#' that has a matching set of key column(s) (\code{id_cols}).
#' @param object The \code{'matchit'} class model object
#' @param model_frame The \code{'data.frame'} class object used in creation of \code{object}.
#' @param id_cols A string indicating the ID for the datset used in the call to \code{\link{matchit}}. 
#' This can be used in combination with \code{newdata} to return the base dataset. Defaults to 
#' \code{NULL}.
#' @param newdata A new \code{data.frame} object to extract matched observations from. Used in 
#' conjunction with \code{id_cols}. Defaults to \code{NULL}.
#' @return If \code{newdata} is \code{NULL}, a subset of \code{model_frame} containing the rows 
#' corresponding to the matched treatement and control observations with weights appended. If 
#' \code{newdata} is not \code{NULL}, an equivalent subset of \code{newdata} is returned.
get_matches.matchit <- function(object, model_frame, id_cols= NULL, newdata= NULL) {
  # 00. error checking
  if (!is.matchit(object)) stop("object must be a class 'matchit' object.")
  if (!is.data.frame(model_frame)) stop("model_frame must be a data.frame")
  if (nrow(model_frame) != nrow(object$X)) 
    stop("model_frame must have the same number of rows as the data used in the call to 'matchit'.")
  if ( (!is.null(id_cols) & is.null(newdata)) | (is.null(id_cols) & !is.null(newdata)) ) {
    stop("For identity returns, both id_cols and newdata must be supplied.")
  }
  if (!is.null(id_cols) & !is.null(newdata)) {
    if (!is.character(id_cols)) stop("id_cols must be a character vector.")
    if (!all(sapply(id_cols, in_frame, df= model_frame))) 
      stop("all values in id_cols must exist in model_frame")
    if (!all(sapply(id_cols, in_frame, df= newdata))) 
      stop("all values in id_cols must exist in newdata")
  }
  
  # 01. preliminaries
  is_exact_matching <- base::grepl(x = object$call[4], pattern= "exact", ignore.case= TRUE)
  use_newdata <- ifelse(is.null(newdata), FALSE, TRUE)
  
  # 02. extract pairs, either exact matching or otherwise
  if (is_exact_matching) {
    row_idx <- as.integer(names(object$subclass)[which(!is.na(object$subclass))])
    match_wts <- object$weights[which(!is.na(object$subclass))]
    model_subset <- data.frame(model_frame[row_idx, ], weight= match_wts)
  } else {
    n_matches <- ncol(object$match.matrix)
    match_wts <- 1 / n_matches # not done via norm of distance
    
    treated_obs <- data.frame(model_frame[as.integer(rownames(object$match.matrix)), ], weight= 1)
    control_obs <- list()
    for (n in 1:n_matches) {
      control_obs[[n]] <- data.frame(model_frame[as.integer(object$match.matrix[, n]), ], weight= match_wts)
    }
    control_obs[[n + 1]] <- treated_obs 
    model_subset <- do.call("rbind", control_obs)  
  }
  
  # 03. return  
  if (!use_newdata) {
    return(model_subset)
  } else { # using newdata via id_cols
    newdata <- base::as.data.frame(newdata) # in case of data.table
    unique_ids <- base::unique(x= model_subset[, c(id_cols, "weight")])
    return(merge(newdata, unique_ids, by= id_cols, all.x=FALSE, all.y=TRUE))
  }
}

#----------------------------------------------------------
### PLOT METHODS
#----------------------------------------------------------

# Need to account for weights -- how do we do qq plots with weights
plot.matchit <- function(x, discrete.cutoff=5, type="QQ",
                         numdraws=5000, interactive = T, which.xs =
                           NULL, ...){
  if ("matchit.exact" %in% class(x)){
    stop("Not appropriate for exact matching.  No plots generated.")
  }
  if(type=="QQ"){
    matchit.qqplot(x=x,discrete.cutoff=discrete.cutoff,
                   numdraws=numdraws, interactive=interactive,
                   which.xs = which.xs, ...)
  } else if(type=="jitter"){
    if("matchit.mahalanobis" %in% class(x)){
      stop("Not appropriate for pure Mahalanobis matching.  No plots generated.")
    }
    jitter.pscore(x, interactive=interactive,...)
  } else if(type=="hist"){
    if("matchit.mahalanobis" %in% class(x)){
      stop("Not appropriate for pure Mahalanobis matching.  No plots generated.")
    }
    hist.pscore(x,...)
  } else {
    stop("Invalid type")
  }
}

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


plot.summary.matchit <- function(x, interactive = TRUE, ...) {
  if ("matchit.exact" %in% class(x)){
    stop("Not appropriate for exact matching.  No plots generated.")
  }
  
  if (!"Std. Mean Diff."%in%names(x$sum.all)){ 
    stop(paste("Not appropriate for unstandardized summary.  Run summary() with the ", 
               "standardize=TRUE option, and then plot."))
  }
  
  sd.pre <- abs(x$sum.all$"Std. Mean Diff.")
  sd.post <- abs(x$sum.matched$"Std. Mean Diff.")
  
  if (!is.null(x$q.table)) sd.post <- abs(x$sum.subclass$"Std. Mean Diff") 
  
  ases.dat <- data.frame(es.unw = sd.pre, es.w = sd.post)
  par(mfrow=c(1,1))
  plot(c(0.85, 2.15), c(0, min(3, max(unlist(ases.dat[, 
                                                      1:2]), na.rm = TRUE))), type = "n", xaxt = "n", 
       ylab = "Absolute Standardized Diff in Means", 
       xlab = "", main = "")
  abline(h = c(0.2, 0.4, 0.6, 0.8, 1.0))
  axis(side = 1, at = 1:2, labels = c("All Data", "Matched Data"))
  for (i in 1:nrow(ases.dat)) {
    points(1:2, abs(ases.dat[i, c("es.unw", "es.w")]), 
           type = "b", col = "grey", pch=19)
  }
  temp1 <- ases.dat[abs(ases.dat$es.unw) < abs(ases.dat$es.w),]
  for (i in 1:nrow(temp1)) {
    points(1:2, abs(temp1[i, c("es.unw", "es.w")]), type = "b", 
           col = "black", lwd = 2, pch=19)
  }
  if (max(ases.dat$es.w, na.rm = TRUE) > 3) 
    mtext(text = "Some standardized diffs in means > 3 after matching!", side = 3, 
          col = "red")
  
  if(interactive==TRUE) {
    print("To identify the variables, use first mouse button; to stop, use second.")
    identify(rep(1, length(sd.pre)),sd.pre,rownames(x$sum.all),atpen=T)
    identify(rep(2, length(sd.post)),sd.post,rownames(x$sum.all),atpen=T)
  }
}

#----------------------------------------------------------
### PRINT METHODS
#----------------------------------------------------------

print.matchit <- function(x, digits = getOption("digits"), ...){
  cat("\nCall: ", deparse(x$call), sep="\n")
  cat("\nSample sizes:\n")
  
  #if(any(x$weights>0)) 
  #  nn <- rbind(table(x$treat),
  #              table(x$weights>0, x$treat),
  #              c(0,0))
  #else 
  #  nn <- rbind(table(x$treat),
  #              table(x$weights>0,x$treat)[2:1,])
  
  print.table(x$nn, ...)
  invisible(x)
  cat("\n")
}

print.matchit.exact <- function(x, digits = getOption("digits"), ...){
  cat("\nCall: ", deparse(x$call), sep = "\n")
  cat("\nExact Subclasses: ", max(x$subclass, na.rm=T),"\n",sep="")
  cat("\nSample sizes:\n")
  ntab <- table(factor(!is.na(x$subclass),
                       levels=c("TRUE","FALSE")),
                x$treat)
  nn <- rbind(table(x$treat),
              ntab[c("TRUE","FALSE"),])
  dimnames(nn) <- list(c("All","Matched","Unmatched"),
                       c("Control","Treated"))
  print.table(nn, ...)
  invisible(x)
  cat("\n")
}


print.matchit.full <- function(x, digits = getOption("digits"), ...){
  cat("\nCall: ", deparse(x$call), sep = "\n")
  cat("\nSample sizes:\n")
  
  if (any(x$weights>0)) 
    nn <- rbind(table(x$treat),
                table(x$weights>0, x$treat),
                c(0,0))
  else 
    nn <- rbind(table(x$treat),
                table(x$weights>0,x$treat)[2:1,])
  dimnames(nn) <- list(c("All","Matched","Discarded"),
                       c("Control","Treated"))
  print.table(nn, ...)
  invisible(x)
  cat("\n")
}


print.matchit.subclass <- function(x, digits = getOption("digits"), ...){
  cat("\nCall: ", deparse(x$call), sep = "\n")
  cat("\nSample sizes by subclasses:\n\n")
  nsub <- table(x$subclass,x$treat)
  nn <- rbind(table(x$treat),nsub)
  dimnames(nn) <-
    list(c("All",paste("Subclass",dimnames(nsub)[[1]])),
         c("Control","Treated"))
  print.table(nn, ...)
  invisible(x)
  cat("\n")
}


print.summary.matchit.exact <- function(x, digits = max(3,
                                                        getOption("digits") - 3),
                                        ...){  
  cat("\nCall:", deparse(x$call), sep = "\n")
  cat("\nSample sizes:\n")
  print.table(x$nn,digits=digits)
  cat("\nMatched sample sizes by subclass:\n")
  print.data.frame(x$q.table, digits = digits)
  cat("\n")
  invisible(x)
}


print.summary.matchit <- function(x, digits = max(3,
                                                  getOption("digits") - 3), ...){
  
  sum.all <- x$sum.all
  sum.matched <- x$sum.matched
  q.table <- x$q.table
  xn <- x$nn
  cat("\nCall:", deparse(x$call), sep = "\n")
  cat("\nSummary of balance for all data:\n")
  print.data.frame(round(sum.all,digits))
  cat("\n")
  
  xs1 <- sum.matched
  cc <- row.names(sum.all)
  if(!is.null(x$sum.matched) | identical(eval(x$call$method),"All"))
  {
    cat("\nSummary of balance for matched data:\n")
    print.data.frame(round(xs1,digits))
    cat("\nPercent Balance Improvement:\n")
    print.data.frame(round(x$reduction,digits))
    cat("\nSample sizes:\n")
    print.table(xn, digits=digits)
    cat("\n")
  }
  invisible(x)
}


print.summary.matchit.subclass <- function(x, digits = max(3,
                                                           getOption("digits") -
                                                             3), ...){   
  sum.all <- x$sum.all
  sum.matched <- x$sum.matched
  q.table <- x$q.table
  cat("\nCall:", deparse(x$call), sep = "\n")
  cat("Summary of balance for all data:\n")
  print.data.frame(round(sum.all,digits))
  cat("\n")
  cat("\nSummary of balance by subclasses:\n")
  print.table(round(q.table, digits))
  cat("\nSample sizes by subclasses:\n")
  print.data.frame(x$qn, digits = digits)
  cat("\nSummary of balance across subclasses\n")
  print.data.frame(round(x$sum.subclass, digits))
  cat("\nPercent Balance Improvement:\n")
  print.data.frame(round(x$reduction,digits))
  cat("\n")
}


#----------------------------------------------------------
### SUMMARY METHODS
#----------------------------------------------------------

summary.matchit <- function(object, interactions = FALSE,
                            addlvariables = NULL, standardize = FALSE,
                            ...) {
  
  X <- object$X
  ## Fix X matrix so that it doesn't have any factors  
  varnames <- colnames(X)
  for(var in varnames) {
    if(is.factor(X[,var])) {
      tempX <- X[,!colnames(X)%in%c(var)]   
      form<-formula(substitute(~dummy-1,list(dummy=as.name(var))))
      X <- cbind(tempX, model.matrix(form, X))
    }
  }
  
  ## No distance output for pure Mahalanobis
  if("matchit.mahalanobis"%in%class(object)){
    XX <- X 
  } else{
    XX <- cbind(distance=object$distance,X)
  }
  if (!is.null(addlvariables)) XX <- cbind(XX, addlvariables)
  
  treat <- object$treat
  weights <- object$weights
  nam <- dimnames(XX)[[2]]
  dupnam <- duplicated(nam)
  if(sum(dupnam)>0){
    nam[dupnam] <- paste(nam[dupnam],".1",sep="")
  }
  kk <- ncol(XX)
  
  ## Summary Stats
  aa <- apply(XX,2,qoi,tt=treat,ww=weights,standardize=standardize,std=T)
  sum.all <- as.data.frame(matrix(0,kk,7))
  sum.matched <- as.data.frame(matrix(0,kk,7))
  row.names(sum.all) <- row.names(sum.matched) <- nam
  names(sum.all) <- names(sum.matched) <- names(aa[[1]])
  sum.all.int <- sum.matched.int <- NULL
  for(i in 1:kk){
    sum.all[i,] <- aa[[i]][1,]
    sum.matched[i,] <- aa[[i]][2,]
    if(interactions){
      for(j in i:kk){
        x2 <- XX[,i]*as.matrix(XX[,j])
        jqoi <- qoi(x2,tt=treat,ww=weights,standardize=standardize,std=T)
        sum.all.int <- rbind(sum.all.int,jqoi[1,])
        sum.matched.int <- rbind(sum.matched.int,jqoi[2,])
        row.names(sum.all.int)[nrow(sum.all.int)] <-
          row.names(sum.matched.int)[nrow(sum.matched.int)] <-
          paste(nam[i],nam[j],sep="x")
      }
    }
  }
  xn <- aa[[1]]$xn
  sum.all <- rbind(sum.all,sum.all.int)
  sum.matched <- rbind(sum.matched,sum.matched.int)
  
  ## Imbalance Reduction
  stat0 <- abs(cbind(sum.all[,2]-sum.all[,1],
                     sum.all[,5:7]))
  stat1 <- abs(cbind(sum.matched[,2]-sum.matched[,1],
                     sum.matched[,5:7]))
  reduction <- as.data.frame(100*(stat0-stat1)/stat0)
  if(sum(stat0==0 & stat1==0, na.rm=T)>0){
    reduction[stat0==0 & stat1==0] <- 0
  }
  if(sum(stat0==0 & stat1>0,na.rm=T)>0){
    reduction[stat0==0 & stat1>0] <- -Inf
  }
  if (standardize)
    names(reduction) <- c("Std. Mean Diff.", "eCDF Med","eCDF Mean",
                          "eCDF Max")
  else
    names(reduction) <- c("Mean Diff.", "eQQ Med","eQQ Mean", "eQQ Max")
  
  ## output
  res <- list(call=object$call, nn = object$nn, sum.all = sum.all,
              sum.matched = sum.matched, reduction = reduction)
  class(res) <- "summary.matchit"
  return(res)
}

summary.matchit.exact <- function(object, covariates = FALSE, ...) {
  XX <- object$X
  treat <- object$treat
  qbins <- max(object$subclass,na.rm=TRUE)
  if(!covariates){
    q.table <- as.data.frame(matrix(0,qbins,3))
    names(q.table) <- c("Treated","Control","Total")
    for(i in 1:qbins){
      qi <- object$subclass==i
      q.table[i,] <- c(sum(treat[qi]==1, na.rm=T), sum(treat[qi]==0, na.rm=T),
                       length(treat[qi & !is.na(qi)]))
    }
  } else {
    kk <- ncol(XX)
    q.table <- as.data.frame(matrix(0,qbins,kk+3))
    names(q.table) <- c("Treated","Control","Total",dimnames(XX)[[2]])
    for(i in 1:qbins){
      qi <- object$subclass==i & !is.na(object$subclass==i)
      q.table[i,] <- c(sum(treat[qi]==1, na.rm=T), sum(treat[qi]==0, na.rm=T),
                       length(treat[qi & !is.na(qi)]),as.numeric(XX[qi,,drop=F][1,])) 
    }
  }
  ntab <- table(factor(!is.na(object$subclass),
                       levels=c("TRUE","FALSE")), treat)
  nn <- rbind(table(treat),
              ntab[c("TRUE","FALSE"),])
  dimnames(nn) <- list(c("All","Matched","Discarded"),
                       c("Control","Treated"))
  ## output
  res <- list(q.table = q.table, nn = nn, subclass = object$subclass,
              treat = object$treat, call = object$call)
  class(res) <- c("summary.matchit.exact", "summary.matchit")
  return(res)
}

summary.matchit.full <- function(object, interactions = FALSE,
                                 addlvariables = NULL, numdraws =
                                   5000, standardize = FALSE,
                                 ...) {
  
  XX <- cbind(distance=object$distance,object$X)
  if (!is.null(addlvariables)) XX <- cbind(XX, addlvariables)
  
  treat <- object$treat
  weights <- object$weights
  nam <- dimnames(XX)[[2]]
  kk <- ncol(XX)
  
  ## Get samples of T and C units to send to qqplot
  t.plot <- sample(names(treat)[treat==1], numdraws/2, replace=TRUE, prob=weights[treat==1])
  c.plot <- sample(names(treat)[treat==0], numdraws/2, replace=TRUE, prob=weights[treat==0])
  
  ## Summary Stats
  aa <- apply(XX,2,qoi,tt=treat,ww=weights, t.plot=t.plot,
              c.plot=c.plot, standardize=standardize)
  sum.all <- as.data.frame(matrix(0,kk,6))
  sum.matched <- as.data.frame(matrix(0,kk,6))
  row.names(sum.all) <- row.names(sum.matched) <- nam
  names(sum.all) <- names(sum.matched) <- names(aa[[1]])
  sum.all.int <- sum.matched.int <- NULL
  for(i in 1:kk){
    sum.all[i,] <- aa[[i]][1,]
    sum.matched[i,] <- aa[[i]][2,]
    if(interactions){
      for(j in i:kk){
        x2 <- XX[,i]*as.matrix(XX[,j])
        names(x2) <- names(XX[,1])
        jqoi <- qoi(x2,tt=treat,ww=weights, t.plot=t.plot,
                    c.plot=c.plot, standardize=standardize)
        sum.all.int <- rbind(sum.all.int,jqoi[1,])
        sum.matched.int <- rbind(sum.matched.int,jqoi[2,])
        row.names(sum.all.int)[nrow(sum.all.int)] <-
          row.names(sum.matched.int)[nrow(sum.matched.int)] <-
          paste(nam[i],nam[j],sep="x")
      }
    }
  }
  xn <- aa[[1]]$xn
  sum.all <- rbind(sum.all,sum.all.int)
  sum.matched <- rbind(sum.matched,sum.matched.int)
  
  ## Imbalance Reduction
  stat0 <- abs(cbind(sum.all[,2]-sum.all[,1],
                     sum.all[,4:6]))
  stat1 <- abs(cbind(sum.matched[,2]-sum.matched[,1],
                     sum.matched[,4:6]))
  reduction <- as.data.frame(100*(stat0-stat1)/stat0)
  if(sum(stat0==0 & stat1==0, na.rm=T)>0){
    reduction[stat0==0 & stat1==0] <- 0
  }
  if(sum(stat0==0 & stat1>0,na.rm=T)>0){
    reduction[stat0==0 & stat1>0] <- -Inf
  }
  if (standardize)
    names(reduction) <- c("Std. Mean Diff.", "eCDF Med","eCDF Mean", "eCDF Max")
  else
    names(reduction) <- c("Mean Diff.", "eQQ Med","eQQ Mean", "eQQ Max")
  
  ## Sample sizes
  nn <- matrix(0, ncol=2, nrow=4)
  nn[1,] <- c(sum(object$treat==0), sum(object$treat==1))
  nn[2,] <- c(sum(object$treat==0 & object$weights>0), sum(object$treat==1 & object$weights>0))
  nn[3,] <- c(sum(object$treat==0 & object$weights==0 & object$discarded==0), 
              sum(object$treat==1 & object$weights==0 & object$discarded==0))
  nn[4,] <- c(sum(object$treat==0 & object$weights==0 & object$discarded==1), 
              sum(object$treat==1 & object$weights==0 & object$discarded==1))
  
  dimnames(nn) <- list(c("All","Matched","Unmatched","Discarded"),
                       c("Control","Treated"))
  
  ## output
  res <- list(call=object$call, nn = nn, sum.all = sum.all, sum.matched = sum.matched,
              reduction = reduction)
  class(res) <- c("summary.matchit.full", "summary.matchit")
  return(res)
}


summary.matchit.subclass <- function(object, interactions = FALSE,
                                     addlvariables=NULL, standardize = FALSE,
                                     ...) {
  
  X <- object$X
  ## Fix X matrix so that it doesn't have any factors
  varnames <- colnames(X)
  for(var in varnames) {
    if(is.factor(X[,var])) {
      tempX <- X[,!colnames(X)%in%c(var)]
      form<-formula(substitute(~dummy-1,list(dummy=as.name(var))))
      X <- cbind(tempX, model.matrix(form, X))
    }
  }
  
  XX <- cbind(distance=object$distance,X)
  if (!is.null(addlvariables)) XX <- cbind(XX, addlvariables)
  
  treat <- object$treat
  weights <- object$weights
  nam <- dimnames(XX)[[2]]
  kk <- ncol(XX)
  
  ## Summary Stats
  aa <- apply(XX,2,qoi,tt=treat,ww=as.numeric(weights!=0),standardize=standardize)
  sum.all <- as.data.frame(matrix(0,kk,6))
  sum.matched <- as.data.frame(matrix(0,kk,6))
  row.names(sum.all) <- row.names(sum.matched) <- nam
  names(sum.all) <- names(sum.matched) <- names(aa[[1]])
  sum.all.int <- sum.matched.int <- NULL
  for(i in 1:kk){
    sum.all[i,] <- aa[[i]][1,]
    sum.matched[i,] <- aa[[i]][2,]
    if(interactions){
      for(j in i:kk){
        x2 <- XX[,i]*as.matrix(XX[,j])
        jqoi <- qoi(x2,tt=treat,ww=as.numeric(weights!=0),standardize=standardize)
        sum.all.int <- rbind(sum.all.int,jqoi[1,])
        sum.matched.int <- rbind(sum.matched.int,jqoi[2,])
        row.names(sum.all.int)[nrow(sum.all.int)] <-
          row.names(sum.matched.int)[nrow(sum.matched.int)] <-
          paste(nam[i],nam[j],sep="x")
      }
    }
  }
  xn <- aa[[1]]$xn
  sum.all <- rbind(sum.all,sum.all.int)
  sum.matched <- rbind(sum.matched,sum.matched.int)
  
  ## By Subclass
  qbins <- max(object$subclass,na.rm=TRUE)
  if(interactions){
    q.table <- array(0,dim=c(kk+sum(1:kk),6,qbins))
    ii <- 0
    nn <- NULL
  } else {
    q.table <- array(0,dim=c(kk,6,qbins))
  }
  aa <- apply(XX,2,qoi.by.sub,tt=treat,ww=weights,
              qq=object$subclass,standardize=standardize)
  for(i in 1:kk){
    if(!interactions){
      q.table[i,,] <- as.matrix(aa[[i]]$q.table)
      nn <- names(aa)
    } else {
      ii <- ii + 1 
      q.table[ii,,] <- as.matrix(aa[[i]]$q.table)
      nn <- c(nn,names(aa)[i])
      for(j in i:kk){
        ii <- ii + 1 
        x2 <- XX[,i]*as.matrix(XX[,j])
        q.table[ii,,] <- as.matrix(qoi.by.sub(x2,tt=treat,ww=weights,qq=object$subclass,
                                              standardize=standardize)$q.table)
        nn <- c(nn,paste(nam[i],nam[j],sep="x"))
      }
    }
  }   
  qn <- aa[[1]]$qn
  dimnames(q.table) <- list(nn,row.names(aa[[i]]$q.table),paste("Subclass",1:qbins))
  
  ## Aggregate Subclass 
  if(is.null(object$call$sub.by)){
    object$call$sub.by <- "treat"
  }
  if(object$call$sub.by=="treat") {
    wsub <- qn[1,]/sum(qn[1,])
  } else if(object$call$sub.by=="control") {
    wsub <- qn[2,]/sum(qn[2,])
  } else if(object$call$sub.by=="all") {
    wsub <- qn[3,]/sum(qn[3,])
  }
  sum.subclass <- sum.all
  for(i in 1:kk){
    for(j in 1:6){
      if(j==3) {
        sum.subclass[i,j] <- sqrt(sum((wsub^2)*(q.table[i,j,]^2)))
      } else {
        sum.subclass[i,j] <- sum(wsub*q.table[i,j,])
      }
    }
  }
  
  ## Imbalance Reduction
  stat0 <- abs(cbind(sum.all[,2]-sum.all[,1],
                     sum.all[,4:6]))
  stat1 <- abs(cbind(sum.subclass[,2]-sum.subclass[,1],
                     sum.subclass[,4:6]))
  reduction <- as.data.frame(100*(stat0-stat1)/stat0)
  if(sum(stat0==0 & stat1==0, na.rm=T)>0){
    reduction[stat0==0 & stat1==0] <- 0
  }
  if(sum(stat0==0 & stat1>0,na.rm=T)>0){
    reduction[stat0==0 & stat1>0] <- -Inf
  }
  if (standardize)
    names(reduction) <- c("Std. Mean Diff.", "eCDF Med","eCDF Mean",
                          "eCDF Max")
  else
    names(reduction) <- c("Mean Diff.", "eQQ Med","eQQ Mean",
                          "eQQ Max")
  ## output
  res <- list(call=object$call, sum.all = sum.all, sum.matched = sum.matched,
              sum.subclass = sum.subclass, reduction = reduction,
              qn = qn, q.table = q.table)
  class(res) <- c("summary.matchit.subclass", "summary.matchit")
  return(res)
}

