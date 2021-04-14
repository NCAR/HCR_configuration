function [Zc,PIA,kwet,b] = back_corr_W(varargin);

% Backward correction for Z
%   This code is for NCAR/EOL HCR
% Usage: [Zc,PIA,k] = back_corr_W(Z,DZ,Rng,alt,surId);
%
%   Input  : Z     => Uncorrect Z in dBZ; a mxn matrix where n is the # of 
%                     ray and m is the # of gates.  
%            DZ    => Total PIA for Nth ray in dB Using for correction; 
%                     a nx1 or 1xn vector;
%            Rng   => Range of ray in km; It has same dimension as Z.
%            altId => The index of altitude of HCR.  It is a 1xn or nx1 
%                     vector.
%            surId => An index of Rng to present the surface. 
%
%   Output : Zc  => Corrected Z in dBZ
%            PIA => Path Integrate Attenuation in dB
%            k   => specific attenuation in dB/km
%            b   => exponent using for each ray.  It is a nx1 or 1xn
%                   vector.

%  Check # of input arguments and assign them
if nargin~=5
  disp('Incorrect # of input arguments.');
  help back_corr;
  return;
end

Z = varargin{1};     DZ = varargin{2};   Rng = varargin{3};
altId = varargin{4};   surId = varargin{5}; 

%  Check the dimension of input;
[Mgate,Nray] = size(Z);
if sum(size(Z)==size(Rng))~=2
  disp('The Dimension of Z and Rng is not agree.');
  help back_corr;
  return;
end

if sum(size(altId)==size(DZ))~=2 
  disp('The length of altId and DZ is not agree.');
  help back_corr;
  return;
end

if sum(size(surId)==size(DZ))~=2 
  disp('surId and DZ should be nx1 or 1xn');
  help back_corr;
  return;
else
  if max(size(altId))~=Nray
    disp(sprintf('altId has %d rays and Z has %d rays',max(size(b)),Nray));
    help back_corr;
    return;
  end
end

%  Start Correction
PIA = zeros(size(Z));
b  = NaN*ones(size(DZ));
kwet = PIA;
Zc  = Z;
C1  = 4/(20*log10(exp(1)));
%  b1=>cloud   b2=>drizzle
b1 = 0.618;   b2 = 0.672;

for n = 1:Nray
  if DZ(n) <= 0 | isnan(DZ(n))==1
    continue;
  end
  k_wet = zeros(Mgate,1);
  pia   = k_wet;
  ZZ = squeeze(Z(:,n));
%  Compute exponent based on Z
  nb1 = length(find(ZZ<-12));   nb2 = length(find(ZZ>-20));
  b(n) = b1*(nb1/(nb1+nb2))+b2*(nb2/(nb1+nb2));
  z  = (10.^(0.1*ZZ)).^b(n);
  z(isnan(z)==1)=0;
  range = squeeze(Rng(:,n));
  dRng = abs(range(2)-range(1));
  
%  Find start and end gate
%     Start gate must be 4 gates away from HCR
  i1 = altId(n)+1;          
  i2 = surId(n)-1;
  if i2-i1<10                       % airplane just tak off
    continue;
  end
  I0 = C1*b(n)*trapz(range(i1:i2),z(i1:i2));
  CC = 10.^(0.1*b(n)*DZ(n))-1;
  for m = i1:i2
    if m == i2
      Ir = 0;
    elseif m == i2-1
      Ir = C1*b(n)*(z(m)+z(m+1))*dRng/2;
    else
      ii1 = m;   
      Ir = C1*b(n)*trapz(range(ii1:i2),z(ii1:i2));
    end
    k_wet(m) = (z(m)*CC)/(I0+CC*Ir);
    pia(m+1) = pia(m)+2*k_wet(m)*dRng;
%     disp(sprintf('m::%d   k_wet::%f   pia::%f', ...
%                  m,k_wet(m),pia(m)));
  end
  if i2<length(range)
    pia(i2+1:end) = pia(i2);
  end
  PIA(:,n) = pia;
  Zc(:,n) = Zc(:,n)+pia(end-1);
  kwet(:,n) = k_wet;
end