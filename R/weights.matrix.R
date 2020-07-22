weights.matrix <- function(match.matrix, treat, discarded){

  n <- length(treat)
  labels <- names(treat)
  tlabels <- labels[treat==1]
  clabels <- labels[treat==0]

  in.sample <- !discarded
  names(in.sample) <- labels

  match.matrix <- match.matrix[tlabels,,drop=F][in.sample[tlabels],,drop=F]
  num.matches <- ncol(match.matrix) - rowSums(is.na(match.matrix))
  names(num.matches) <- tlabels[in.sample[tlabels]]

  t.units <- row.names(match.matrix)[num.matches > 0]
  c.units <- na.omit(as.vector(as.matrix(match.matrix)))

  weights <- rep(0,length(treat))
  names(weights) <- labels
  weights[t.units] <- 1

  for (cont in clabels) {
    treats <- na.omit(row.names(match.matrix)[cont==match.matrix[,1]])
    if (ncol(match.matrix)>1)
      for (j in 2:ncol(match.matrix))
        treats <- c(na.omit(row.names(match.matrix)[cont==match.matrix[,j]]),treats)
    for (k in unique(treats))
      weights[cont] <- weights[cont] + 1/num.matches[k]
  }

  if (sum(weights[clabels])==0)
    weights[clabels] <- rep(0, length(weights[clabels]))
  else
    weights[clabels] <- weights[clabels]*length(unique(c.units))/sum(weights[clabels])

  weights[!in.sample] <- 0
  if (sum(weights)==0)
    stop("No units were matched")
  else if (sum(weights[tlabels])==0)
    stop("No treated units were matched")
  else if (sum(weights[clabels])==0)
    stop("No control units were matched")

  return(weights)
}

