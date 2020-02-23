function y = GaussianCurve(x, params)
%GAUSSMF is a curvev of a normalized Gaussian, centered about a value C,
%   and with a width of sigma.
%
%   GaussianCurve(X, PARAMS) returns a matrix which is the Gaussian
%   membership function evaluated at X. PARAMS is a 2-element vector
%   that determines the width and center of the Gaussian:
%
%   GaussianCurve(X, [SIGMA, C]) = EXP(-(X - C).^2/(2*SIGMA^2));
%  

sigma = params(1); c = params(2);
y = exp(-(x - c).^2/(2*sigma^2));
