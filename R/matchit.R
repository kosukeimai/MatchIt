matchit <- function(formula, model="logit", data, discard=0,
                    reestimate=FALSE, nearest=TRUE, replace=FALSE,
                    m.order=2, ratio=1, caliper=0, calclosest=FALSE,
                    subclass=0, sub.by="treat", mahvars=NULL, exact=FALSE,
                    counter=TRUE, full = list(full=FALSE, min.controls = 0, 
                                    max.controls = Inf,
                                    omit.fraction = NULL, tol = 0.01), ...){ 
  cl <- match.call()
  if (m.order==1 | full$full)
    require(optmatch)
  
  #Checking input format
  #data input
  is.null(data)  #there must be a better way to get a warning for data input
  if(!is.data.frame(data)){
    stop("Data ", cl$data, " must be a dataframe",call.=FALSE)}
  #check pscore / psweights names in dataframe
  if("pscore"%in%names(data)){
    stop("Dataframe contains the variable 'pscore'.  Please change this name.",call. = FALSE)
  } 
  if("psweights"%in%names(data)){
    stop("Dataframe contains the variable 'psweights'.  Please change this name.",call.=FALSE)
  }
  if("psclass"%in%names(data)){
    stop("Dataframe contains the variable 'psclass'.  Please change this name.", call.=FALSE)
  }

  #sub.by
  if(!is.vector(sub.by)|length(sub.by)!=1){
    warning(cl$sub.by," is not a valid sub.by option; sub.by=\"treat\" used instead",call.=FALSE); sub.by <- "treat"}
  if(!sub.by%in%c("treat","control","all")){
    warning(cl$sub.by, " is not a valid sub.by option; sub.by=\"treat\" used instead",call.=FALSE); sub.by <- "treat"}
  #subclass
  if(length(subclass)==1){
    if(subclass<0 | subclass==1){
      warning(cl$subclass, " is not a valid subclass; subclass=0 used instead",call.=FALSE); subclass <- 0}
  } else {
    if(!is.vector(subclass)){
      warning(cl$subclass, " is not a valid subclass; subclass=0 used instead",call.=FALSE); subclass <- 0}
    if(sum(subclass<=1 & subclass>=0)!=length(subclass)){
      warning("Subclass ", cl$subclass, " is not bounded by 0 and 1; subclass=0 used instead",
              call.=FALSE); subclass <- 0}
  }
  #nearest
  if(!(identical(nearest,TRUE) | identical(nearest,FALSE))){
    warning("nearest=",cl$nearest," is invalid; used nearest=TRUE instead",call.=FALSE);nearest=TRUE}
  #replace
  if(!(identical(replace,TRUE) | identical(replace,FALSE))){
    warning("replace=",cl$replace," is invalid; used replace=FALSE instead",call.=FALSE);replace=0}
  #m.order
  if(!(identical(m.order,2) | identical(m.order,3) |
       identical(m.order,4) | identical(m.order,1))){
    warning("m.order=",cl$m.order," is invalid; used m.order=2 instead",call.=FALSE);m.order=2}
  #ratio
  ratio <- round(ratio)
  if(!is.numeric(ratio) | ratio[1]<1 | !identical(round(length(ratio)),1)){
    warning("ratio=",cl$ratio," is invalid; used ratio=1 instead",call.=FALSE);ratio=1}
  #caliper
  if(!is.vector(caliper) | !identical(round(length(caliper)),1)){
    warning("caliper=",cl$caliper," is invalid; used caliper=0 instead",call.=FALSE);caliper=0}
  if(caliper<0){
    warning("caliper=",cl$caliper," is less than 0; used caliper=0 instead",call.=FALSE);caliper=0}
  #calclosest
  if(!(identical(calclosest,TRUE)| identical(calclosest,FALSE))){
    warning("calclosest=",cl$calclosest," is invalid; used calclosest=FALSE instead",call.=FALSE)
    calclosest=FALSE}  
  #mahvars & caliper
  if (!is.null(mahvars) & caliper[1]==0){
    warning("Must specify caliper > 0 to use Mahalanobis matching. Mahalanobis matching not done",call. = FALSE)}
  #counter
  if(!(identical(counter,TRUE)| identical(counter,FALSE))){
    warning("counter=",cl$counter," is invalid; used counter=TRUE instead",call.=FALSE);counter=TRUE}    
  #full matching
  if(!(identical(full$full,TRUE) | identical(full$full,FALSE))){
    warning("full=list(full=",cl$full,") is invalid; used full=list(full=FALSE) instead",call. = FALSE)
    full$full <- FALSE
  } else if(full & nearest){
    warning("Full matching will ignore nearest neighbor inputs",call. = FALSE)
    nearest <- F
  }
  
  # Set up for output
  tt <- terms(formula(cl))

  if (identical(TRUE, exact)) {
    if(!is.null(c(cl$discard, cl$reestimate, cl$model, cl$nearest, cl$replace, cl$m.order, cl$ratio,
                  cl$caliper, cl$calclosest, cl$subclass, cl$sub.by, cl$mahvars)))
      warning("exact=TRUE chosen; all other matching options ignored")
    data <- eval(data,parent.frame())
    mf <- match.call()
    dotsub <- eval(mf$subset,data)
    if(!is.null(dotsub)) { 
      data <- data[dotsub,]
    }

    treata <- model.frame(formula,data)[,1,drop=FALSE]
    treat <- as.vector(treata[,1])
    names(treat) <- row.names(treata)
    covariates <- model.frame(delete.response(terms(formula)),data)[,,drop=FALSE]

    b <- exactmatch(formula, data, counter=counter)
    
    dummy.pscore <- rep(0.5, dim(covariates)[1])
    dummy.in.sample <- rep(TRUE, dim(covariates)[1])
    names(dummy.pscore) <- names(treat)
    names(dummy.in.sample) <- names(treat)
    c <- diagnose(formula, match.matrix=NULL, pscore=dummy.pscore,
                  in.sample = dummy.in.sample, data=data,exact=TRUE,
                  mahvars=NULL, subclass=max(b$psclass),
                  psclass=b$psclass, nearest=FALSE, counter=counter)
    
    psweights <- c$weights
    pscore <- dummy.pscore

    if (is.null(b$psclass)) psclass <- rep(1, length(pscore))
    else psclass <- b$psclass

    data <- cbind(data,psweights,pscore, psclass)
    
    z <- list(formula=formula, psweights=c$weights,
              treat=attr(terms(formula(cl)),"variables")[[2]],
              in.sample=NULL, match.matrix=b$match.matrix,
              covariates=attr(delete.response(tt),"term.labels"), psclass=b$psclass,
              q.cut=NULL, assign.model=NULL, call=cl, matched=(c$weights!=0),
              data= data)
    class(z) <- "matchit"
    return(z)
  }
  else {
    mf <- match.call()
    dotsub <- eval(mf$subset,data)
    mf[[1]] <- as.name("distance")
    mf$nearest <- mf$replace <- mf$m.order <- mf$ratio <- mf$caliper <- 
      mf$calclosest <- mf$subclass <- mf$sub.by <- mf$mahvars <-
        mf$exact <- mf$full <- NULL 
    a <- eval(as.call(mf), sys.frame(sys.parent()))
    if(!is.null(dotsub)) { 
      data <- data[dotsub,]
    }

    b <- matchdef(formula, a$in.sample, a$pscore, nearest = nearest,
                  replace = replace, m.order = m.order, ratio = ratio,
                  caliper = caliper, calclosest = calclosest, 
                  mahvars = mahvars,data = data, exact = exact,
                  counter = counter)
    
    b1 <- subclassify(formula, data, a$in.sample, a$pscore,
                      nearest = nearest, b$match.matrix,
                      subclass = subclass, sub.by = sub.by, counter =
                      counter, full = full)
    
    c <- diagnose(formula, b$match.matrix, a$pscore, a$in.sample,
                  data=data, exact=exact, mahvars=mahvars,
                  subclass=subclass, psclass=b1$psclass,nearest=nearest,
                  q.cut=b1$q.cut, counter=counter)
 
    ##adding pscore and weights to dataframe
    pscore <- a$pscore
    psweights <- c$weights
    
    # Create psclass with subclasses
    if (is.null(b1$psclass)) psclass <- rep(1, length(pscore))
    else psclass <- b1$psclass
    
    data <- cbind(data,psweights,pscore, psclass)

    z <- list(formula=formula,psweights=psweights,
              treat=attr(terms(formula(cl)),"variables")[[2]],
              in.sample=a$in.sample, match.matrix=b$match.matrix,
              covariates=attr(delete.response(tt),"term.labels"),
              psclass=b1$psclass, q.cut=b1$q.cut, assign.model=a$assign.model,
              call=cl,matched=(c$weights!=0), data=data)
    class(z) <- "matchit"
    return(z)
  }
}
