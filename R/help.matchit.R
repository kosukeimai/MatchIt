help.matchit <- function (object=NULL) 
{
    under.unix <- !(version$os == "Microsoft Windows" || version$os == 
        "Win32" || version$os == "mingw32")
    sys <- function(command, text = NULL) {
        cmd <- if (length(text)) 
            paste(command, text)
        else command
        if (under.unix) 
            system(cmd)
        else shell(cmd, wait = TRUE)
    }
    browser <- .Options$help.browser
    if (!length(browser)) 
        browser <- .Options$browser
    if (!length(browser)) 
        browser <- getOption("browser")
    url <- NULL
 
    if (is.null(object))
      url <- c("http://gking.harvard.edu/matchit")
    if (!is.null(object)) {
    if (object == "matchit")
        url <- c("http://gking.harvard.edu/matchit/docs/Usage.html")
    if (object == "distance")
      url <- c("http://gking.harvard.edu/matchit/docs/_TT_distance_TT.html")
    if (object == "matchdef")
      url <- c("http://gking.harvard.edu/matchit/docs/_TT_matchdef_TT.html")
    if (object == "subclassify")
      url <- c("http://gking.harvard.edu/matchit/docs/_TT_subclassify_TT.html")
    if (object == "print.matchit")
      url <- c("http://gking.harvard.edu/matchit/docs/Print.html")
    if (object == "summary.matchit")
      url <- c("http://gking.harvard.edu/matchit/docs/Summary.html")
    if (object == "plot.matchit")
      url <- c("http://gking.harvard.edu/matchit/docs/Plot.html")
    }

    if (is.null(url)) {
        cat("Error:", object, "currently not documented in help.matchit. \n Please check http://gking.harvard.edu/matchit. \n", 
            sep = " ")
	url <- c("http://gking.harvard.edu/matchit")
    }

    if (under.unix) {
        sys(paste(browser, url, "&"))
        invisible()
    }
    if (!under.unix) {
        browseURL(url, browser = browser)
        invisible("")
    }
}
