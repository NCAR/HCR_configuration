function [vr,vi]=getInitFallSpeed(dbz)
% Calculate initial fall speed

% V rain
% https://doi-org.cuucar.idm.oclc.org/10.1175/1520-0426(1994)011<1656:EOVIPS>2.0.CO;2
a=3.5;
b=0.084;

zlin=10.^(dbz./10);

vr=a.*zlin.^b;

% V ice
% https://doi.org/10.1175/2009JAS3132.1
vi=nan(size(dbz));

vi(dbz<30)=1.4041+0.00237.*dbz(dbz<30);
vi(dbz>=30)=-3.4+0.1864.*dbz(dbz>=30);

end