function momentsTime=calcMomentsTime(cIQ,ii,momentsTime,data)
%% V
cIQ.v=cIQ.v.*sqrt(size(cIQ.v,2));

R0v=mean(real(cIQ.v).^2+imag(cIQ.v).^2,2);
R1v=mean(cIQ.v(:,1:end-1).*conj(cIQ.v(:,2:end)),2);
R2v=mean(cIQ.v(:,1:end-2).*conj(cIQ.v(:,3:end)),2);
R3v=mean(cIQ.v(:,1:end-3).*conj(cIQ.v(:,4:end)),2);
R4v=mean(cIQ.v(:,1:end-4).*conj(cIQ.v(:,5:end)),2);

momentsTime.powerV(:,ii)=10*log10(R0v)-data.rx_gain_v;
momentsTime.vel(:,ii)=data.lambda/(4*pi*mode(data.prt))*angle(R1v);
momentsTime.width(:,ii)=data.lambda/(2*pi.*mode(data.prt)*6^.5)*abs(log(abs(R1v./R2v))).^0.5;
momentsTime.skew(:,ii)=abs(log(abs(R3v./(R2v.^3))));
momentsTime.kurt(:,ii)=abs(log(abs(R4v./(R2v.^2))));

% SNRV
noiseLinV=10.^(data.noise_v./10);
snrLinV=(R0v-noiseLinV)./noiseLinV;
snrLinV(snrLinV<0)=nan;
momentsTime.snr(:,ii)=10*log10(snrLinV);

% DBZV
data.range(data.range<0)=nan;
momentsTime.dbz(:,ii)=momentsTime.snr(:,ii)+20*log10(data.range./1000)+data.dbz1km_v;

%% H
if isfield(cIQ,'h')
    cIQ.h=cIQ.h.*sqrt(size(cIQ.h,2));
    R0h=mean(real(cIQ.h).^2+imag(cIQ.h).^2,2);

    % Power H
    momentsTime.powerH(:,ii)=10*log10(R0h)-data.rx_gain_h;

    % SNRH
    noiseLinH=10.^(data.noise_h./10);
    snrLinH=(R0h-noiseLinH)./noiseLinH;
    snrLinH(snrLinH<0)=nan;
    snrH=10*log10(snrLinH);

    % DBZH
    dbzH=snrH+20*log10(data.range./1000)+data.dbz1km_h;

    % LDR
    momentsTime.ldr(:,ii)=dbzH-momentsTime.dbz(:,ii);
end
end

