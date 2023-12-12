function [np,rc,ortho,npx] = forsythe_poly(x,n,varargin)
% np is a matrix of poly coeffs, [Pn,Pn-1,...P0], 
%     Pn column vec here is the nth order forsythe poly for based on x
%     with decreasing order from nth down to 0th order.
%     P0, then, has only the last term in it.
% rc  reciprocal condition number of forsythe polynomials.
%     values near 1 = good, not near 1 = bad.
% ortho matrix of all inner products of the forsythe polys
%     perfect would be an identity matrix
% npx forsythe polys evaluated at x [Pn(x);Pn-1(x);...;P0(x)] ROWS!
%

wts = ones(size(x));

paramparse(varargin);


x = reshape(x,[],1);
wts = reshape(wts,[],1);

% Polynomials are just represented as a row vector of coefficients
% in decreasing order.  Length n implies a poly of order n-1.

% initialize matrix of polynomial coefficients
% the kth row corresponds the P_k-1, (indexes from 1)
% the last column is the constant coeff, the 2nd to last is the
% linear x^1 coeff, and so on.
p = zeros(n+1,n+1);

% set the first polynomial to 1;

p(1,end) = 1;

for ll = 1:n
  % get previous poly
  pp = p(ll,end-ll+1:end);
  % evaluate it at the x's
  ppxs = polyval(pp,x).^2;
  if ll==1
    % special startup case
    beta(1) = 0;
    lt = 0;
  else
    % otherwise compute beta and lt
    ppp = p(ll-1,end-ll+1:end);
    beta(ll) = sum(wts.*ppxs)/sum(wts.*polyval(ppp,x).^2);
    lt = beta(ll)*[0 ppp];
  end
  alpha(ll+1) = sum(wts.*x.*ppxs)/sum(wts.*ppxs);
  % recursion formula.  The conv multiplies 2 polynomials
  np = conv([1 -alpha(ll+1)],pp) - lt;
  % save back into the polynomial matrix
  p(ll+1,(end-ll):end) = np;
end

% Now normalize all so they have wtd L2 norm = 1
% must do this now, otherwise the recusrion formula won't work.
for ll = 1:size(p,1)
  np(ll,:) = p(ll,:)/sqrt(sum(wts.*polyval(p(ll,:),x).^2));
end

% now reorganize the polynomial matrix.  Store them as columns in
% decreasing order: [Pn,Pn-1,...P0]
np = flipud(np).';

% Now evaluate the polunomials at x as ROWS: [Pn(x);Pn-1(x);..;P0(x)]
for ll = 1:size(p,1), 
  npx(ll,:) = polyval(np(:,ll),x);
end;

% estimate how well it worked by estimating all the innerproducts
ortho = npx*diag(wts)*npx.';
% compute the reciprocal condition number
rc = rcond(ortho);
if rc<.9
  warning('The polynoials are diverging from orthogonal.  The fits probably will not be good.');
end