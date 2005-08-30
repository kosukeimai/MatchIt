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
        url <- c("http://gking.harvard.edu/matchit/docs/Reference_Manual.html")
    if (object == "exact")
      url <- c("http://gking.harvard.edu/matchit/docs/Exact_Matching2.html")
    if (object == "subclass")
      url <- c("http://gking.harvard.edu/matchit/docs/Subclassification2.html")
    if (object == "nearest")
      url <- c("http://gking.harvard.edu/matchit/docs/Nearest_Neighbor_Match2.html")
    if (object == "optimal")
      url <- c("http://gking.harvard.edu/matchit/docs/Optimal_Matching2.html")
    if (object == "full")
      url <- c("http://gking.harvard.edu/matchit/docs/Full_Matching2.html")
    if (object == "match.data")
      url <- c("http://gking.harvard.edu/matchit/docs/_TT_match_data_TT.html")
    if (object == "summary")
      url <- c("http://gking.harvard.edu/matchit/docs/_TT_summary_TT.html")
    if (object == "plot")
      url <- c("http://gking.harvard.edu/matchit/docs/_TT_plot_TT.html")
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
