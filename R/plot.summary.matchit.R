plot.summary.matchit <- function(x, interactive = TRUE, ...) {
  if ("matchit.exact" %in% class(x)){
    stop("Not appropriate for exact matching.  No plots generated.")
  }

  if (!"Std. Mean Diff."%in%names(x$sum.all)){ 
	stop("Not appropriate for unstandardized summary.  Run summary() with the standardize=TRUE option, and then plot.")
	}
	
	sd.pre <- abs(x$sum.all$"Std. Mean Diff.")
	sd.post <- abs(x$sum.matched$"Std. Mean Diff.")

	if (!is.null(x$q.table)) sd.post <- abs(x$sum.subclass$"Std. Mean Diff") 

	ases.dat <- data.frame(es.unw = sd.pre, es.w = sd.post)
	par(mfrow=c(1,1))
        plot(c(0.85, 2.15), c(0, min(3, max(unlist(ases.dat[, 
            1:2]), na.rm = TRUE))), type = "n", xaxt = "n", ylab = "Absolute Standardized Diff in Means", 
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


