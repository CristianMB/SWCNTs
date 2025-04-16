//ADAPTED VERSION FOR VOIGT FUNCTION (all other functions are covered through the more generic Faddeeva_mex.cc, but this function has different number and type of arguments)
//#define FADDEEVA_FUNC Faddeeva::Voigt
//#define FADDEEVA_REAL 1

/* Copyright (c) 2012 Massachusetts Institute of Technology
 * 
 * [Also included are functions derived from derfc in SLATEC
 *  (netlib.org/slatec), which "is in the public domain"
 *  and hence may be redistributed under these or any terms.]
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
 */

/* This file provides generic Matlab wrappers for the different Faddeeva::
   functions.  To wrap a specific function, we create a new .cc file with

   #define FADDEEVA_FUNC ...name of function...
   #define FADDEEVA_REAL ...1 if real for real z, 0 if not...
   #include "Faddeeva_mex.cc"
*/

#include "Faddeeva.hh"

#include <mex.h>

// double FADDEEVA(Voigt)(double x, double sigma, double gamma, double relerr)

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  if (nrhs < 3 || nrhs > 4)
    mexErrMsgTxt("expecting three or four arguments");
  if (nlhs > 1)
    mexErrMsgTxt("expecting only one return value");

  if (!mxIsNumeric(prhs[0]) || !(mxIsDouble(prhs[0]) || mxIsSingle(prhs[0])))
    mexErrMsgTxt("first argument must be numeric array (double or single precision)");

  // relerr = prhs[1], if any
  double relerr;
  if (nrhs < 4)
    relerr = 0;
  else if (mxIsNumeric(prhs[3]) && mxGetM(prhs[3]) * mxGetN(prhs[3]) == 1) {
    if (mxIsDouble(prhs[3]))
      relerr = *mxGetPr(prhs[3]);
    else if (mxIsSingle(prhs[3]))
      relerr = *((float *) mxGetData(prhs[3]));
    else
      mexErrMsgTxt("fourth argument must be double or single precision");
  }
  else
    mexErrMsgTxt("fourth argument must be real scalar");

//SECOND AND THIRD ARGUMENTS: sigma, gamma:
  double sigma;
  double gamma;
  if (mxIsNumeric(prhs[1]) && mxGetM(prhs[1]) * mxGetN(prhs[1]) == 1) {
    if (mxIsDouble(prhs[1]))
      sigma = *mxGetPr(prhs[1]);
    else if (mxIsSingle(prhs[1]))
      sigma = *((float *) mxGetData(prhs[1]));
    else
      mexErrMsgTxt("second argument (sigma) must be double or single precision");
  }
  else
    mexErrMsgTxt("second argument (sigma) must be real scalar");

  if (mxIsNumeric(prhs[2]) && mxGetM(prhs[2]) * mxGetN(prhs[2]) == 1) {
    if (mxIsDouble(prhs[2]))
      gamma = *mxGetPr(prhs[2]);
    else if (mxIsSingle(prhs[2]))
      gamma = *((float *) mxGetData(prhs[2]));
    else
      mexErrMsgTxt("third argument (gamma) must be double or single precision");
  }
  else
    mexErrMsgTxt("third argument (gamma) must be real scalar");



  mwSize ndim = mxGetNumberOfDimensions(prhs[0]);
  const mwSize *dims = mxGetDimensions(prhs[0]);
  plhs[0] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
  double *wr = mxGetPr(plhs[0]);
  //double *wi = mxGetPi(plhs[0]);

  size_t N = 1;
  for (mwSize d = 0; d < ndim; ++d) N *= dims[d];  // get total size of array

  void *vzr = mxGetData(prhs[0]);
  //void *vzi = mxGetImagData(prhs[0]);
  if (mxIsDouble(prhs[0])) {
    double *zr = (double*) vzr;
      for (size_t i = 0; i < N; ++i) {
	wr[i] = Faddeeva::Voigt(zr[i],sigma,gamma,relerr);
      }
  }
  else { // single precision
    float *zr = (float*) vzr;
      for (size_t i = 0; i < N; ++i) {
	wr[i] = Faddeeva::Voigt(zr[i],sigma,gamma,relerr);
      }
  }
}
