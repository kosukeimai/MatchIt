#ifndef INTERNAL_H
#define INTERNAL_H

#include <Rcpp.h>
using namespace Rcpp;

IntegerVector tabulateC_(const IntegerVector& bins,
                         const Nullable<int>& nbins = R_NilValue);


#endif