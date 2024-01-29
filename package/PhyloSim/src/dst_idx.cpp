/***************************************************************************\
 *
 * Authors: Guillaume Guénard
 * Université de Montréal - 2023
 *
 * Simulation of DNA sequences evolution as a Markov process at the scale of
 * single nucleotides, with the purpose of generating simulated data for deep
 * learning model developments.
 *
 * Calculating distance vector indices from lower triangular matrix indices.
 *
 \***************************************************************************/

#include <Rcpp.h>
using namespace Rcpp;

IntegerVector dstIdx_Cpp(const int n, const IntegerVector& from,
                         const IntegerVector& to) {
  
  int ni = from.size(), nj = to.size();
  int nn = (ni > nj) ? ni : nj;
  int i, ii, j, jj, k;
  IntegerVector idx(nn);
  
  for(i = 0, j = 0, k = 0; k < nn; i++, j++, k++) {
    if(i == ni) i = 0;
    if(j == nj) j = 0;
    ii = from[i];
    jj = to[j];
    if(jj > ii)
      idx[k] = jj + (ii - 1)*n - ii*(ii + 1)/2;
    else if(jj < ii)
      idx[k] = ii + (jj - 1)*n - jj*(jj + 1)/2;
    else
      idx[k] = NA_INTEGER;
  }
  
  return(idx);
}
