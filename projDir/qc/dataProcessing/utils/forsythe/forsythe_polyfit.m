function [pc, ff, fc ]= forsythe_polyfit(x,y,n,varargin);
% forsythe_polyfit - perform least squares fit using forsythe polynomial method
%
%  [pc, ff, fc ]= forsythe_polyfit(x,y,n,'param1',value1,...);
%
%  x - vector of x-values
%  y - array of corresponding y-values.  Can be matrix with one of the dimensions
%      coresponding to x (must have the same size).  That dimension is assumed to be
%      the first dimension.  Set the 'dim' parameter to change this.
%  n - the order of the highest polynomial desired for the fit.
%
%  optional:
%     dim - default 1; the dimension of y corresponding to x.
%     wts - default ones(size(x)); The weights for the fit, if desired.
%
%  Outputs:
%     pc - the polynomial coefficients of the best fit.  Will have the same size as y except
%          that the size(pc,dim) will be n+1.  Basically, the data along dimension dim, will be
%          the matlab polynomial (coeffients in descending order).
%     ff - The best fit fpolynomial evaluated at x.  In other words, this is just pc evaluated at x
%          This will have the same size as y.
%     fc - these are the coefficeients of the best fit relative to the forsythe polynomials.
%          This is probably not helpful.

dim = 1;
wts = ones(size(x));

paramparse(varargin);

% capture the forsythe polynomial coeffs, the two diagnostic parameters, and the forsythe polys evaluated at x
[p,rc,ortho,px] = forsythe_poly(x,n,'wts',wts);

% now put the important dim first, and reshape to a matrix
sz = size(y);
pvec = 1:max(dim,length(sz));
pvec([1 dim]) = pvec([dim 1]);
y = permute(y,pvec);

sz2 = size(y);
y = reshape(y,sz2(1),[]);

% now y should be length(x) x ?
% px is a n x length(x) matrix of the form p_n(x)
% so *y should give all poly fits.
fc = px*diag(wts)*y;

% now fc are the coefficients in front of the forsythe polynomials
% p = n x n, and multiplying on left by p should give us the
% typical coeffes that can by used by polyfit
pc = p*fc;

% fc is a n x ?    px.' is length(x) x n  so px.'*fc is length(x) x ?
% This gives us the best fit evaluated at x.
ff = px.'*fc;

% now reshape and permute it back to what it was 
sz3 = sz2;
sz3(1) = size(fc,1);

ff = reshape(ff,sz2);
ff = permute(ff,pvec);

fc = reshape(fc,sz3);
fc = permute(fc,pvec);

pc = reshape(pc,sz3);
pc = permute(pc,pvec);

