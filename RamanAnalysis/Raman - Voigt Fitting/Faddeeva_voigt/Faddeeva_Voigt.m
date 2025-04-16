% Usage: v = Faddeeva_Voigt(x,x0,sigma,gamma [, relerr])
% 
% Compute Voigt(x,x0,sigma,gamma), the convolution of a Gaussian and a Lorentzian, 
% for an (array or matrix of) real value(s) x, central position x (single value),
% Gaussian width sigma, and Lorentzian HALF-width-at-half-maximum (HWHM) gamma.
% 
% relerr, if supplied, indicates a desired relative error tolerance in
% w; the default is 0, indicating that machine precision is requested (and
% a relative error < 1e-13 is usually achieved).  Specifying a larger
% relerr may improve performance for some z (at the expense of accuracy).
% 
% Slight modification (by Wim Wenseleers) of the code for the Faddeeva function by
% S. G. Johnson, http://ab-initio.mit.edu/Faddeeva
