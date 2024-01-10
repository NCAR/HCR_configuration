function [np,rc,ortho,npx] = forsythe_poly(x,n,varargin)
% [np,rc,ortho,npx] = forsythe_poly(x,n,'param1',value1,...)
%
%  Generate the forsythe polynomials for the provided domain (x)
%  i.e. polynomials that are orthonormal wrt the inner product:
%    <f,g> = sum (f(x(i))*g(x(i)))
%             i
% Inputs:
%   x - x-values that the polynomials will be evaluated at
%   n - the highest order polynomial to compute
%  optional inputs:
%   wts - can input wts corresponding to x.  Not sure this is working properly.  default is all 1's
%   do_npx_independently - if 1 (Default) will compute npx numerically iteratively.  If 0, then
%       npx is created based on the polynomial coefficients.  For larger values of 'n', setting to
%       0 will not give stable results.
% Outputs:
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
%  

wts = ones(size(x));
do_npx_independently = 1;

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

if do_npx_independently
  % initialize npx.  the kth row corresponds the P_k-1, (indexes from 1).
  % columns correspond to x
  npx = zeros(n,length(x));
  
  % set the first polynomial to 1
  npx(1,:) = 1;
end

for ll = 1:n
  % get previous poly
  pp = p(ll,end-ll+1:end);
  % evaluate it at the x's
  
  ppxs = polyval(pp,x).^2;
  
  if do_npx_independently
    ppx_n = npx(ll,:).';
    ppxs_n = (ppx_n).^2;
  end
  
  if ll==1
    % special startup case
    beta(1) = 0;
    lt = 0;
    
    if do_npx_independently
      beta_n(1) = 0;
      ppx_n1 = 0;
      ppxs_n1 = 0;
    end
  else
    % otherwise compute beta and lt
    ppp = p(ll-1,end-ll+1:end);
    beta(ll) = sum(wts.*ppxs)/sum(wts.*polyval(ppp,x).^2);
    lt = beta(ll)*[0 ppp];
    
    if do_npx_independently
      ppx_n1 = npx(ll-1,:).';
      ppxs_n1 = (ppx_n1).^2;
      beta_n(ll) = sum(wts.*ppxs_n)/sum(wts.*ppxs_n1); %!
    end
  end
  
  alpha(ll+1) = sum(wts.*x.*ppxs)/sum(wts.*ppxs);
  
  if do_npx_independently
    alpha_n(ll+1) = sum(wts.*x.*ppxs_n)/sum(wts.*ppxs_n);
  end
  
  
  % recursion formula.  The conv multiplies 2 polynomials
  np = conv([1 -alpha(ll+1)],pp) - lt;
  % save back into the polynomial matrix
  p(ll+1,(end-ll):end) = np;

  if do_npx_independently
    npx(ll+1,:) = ( (x-alpha_n(ll+1)).*ppx_n - beta_n(ll)*ppx_n1).';
  end
  


end

% Now normalize all so they have wtd L2 norm = 1
% must do this now, otherwise the recusrion formula won't work.
for ll = 1:size(p,1)
  np(ll,:) = p(ll,:)/sqrt(sum(wts.*polyval(p(ll,:),x).^2));
end

if do_npx_independently
  npx = npx./sqrt(npx.^2 * wts);
end

% now reorganize the polynomial matrix.  Store them as columns in
% decreasing order: [Pn,Pn-1,...P0]
np = flipud(np).';

% Now evaluate the polunomials at x as ROWS: [Pn(x);Pn-1(x);..;P0(x)]
if do_npx_independently
  npx = flipud(npx);
else
  for ll = 1:size(p,1), 
    npx(ll,:) = polyval(np(:,ll),x);
  end;
end

% estimate how well it worked by estimating all the innerproducts
ortho = npx*diag(wts)*npx.';
% compute the reciprocal condition number
rc = rcond(ortho);
if rc<.9
  warning('The polynoials are diverging from orthogonal.  The fits probably will not be good.');
end