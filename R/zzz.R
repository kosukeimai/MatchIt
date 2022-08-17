#Used to load backports functions. No need to touch, but must always be included somewhere.
.onLoad <- function(libname, pkgname) {
  backports::import(pkgname)
}