function momentsSpec=calcMomentsSpec_powerFields(specPowerLin,ii,momentsSpec,data)
%% V
% DBM
powerLinV=mean(specPowerLin.V,2);
momentsSpec.powerV(:,ii)=10*log10(powerLinV)-data.rx_gain_v;

% SNR
noiseLinV=10.^(data.noise_v./10);
snrLinV=(powerLinV-noiseLinV)./noiseLinV;
snrLinV(snrLinV<0)=nan;
momentsSpec.snr(:,ii)=10*log10(snrLinV);

% DBZ
data.range(data.range<0)=nan;
momentsSpec.dbz(:,ii)=momentsSpec.snr(:,ii)+20*log10(data.range./1000)+data.dbz1km_v;

%% H
if isfield(specPowerLin,'H')
    % DBM
    powerLinH=mean(specPowerLin.H,2);
    momentsSpec.powerH(:,ii)=10*log10(powerLinH)-data.rx_gain_h;

    % SNR
    noiseLinH=10.^(data.noise_h./10);
    snrLinH=(powerLinH-noiseLinH)./noiseLinH;
    snrLinH(snrLinH<0)=nan;
    snrH=10*log10(snrLinH);

    % DBZ
    dbzH=snrH+20*log10(data.range./1000)+data.dbz1km_h;

    % LDR
    momentsSpec.ldr(:,ii)=dbzH-momentsSpec.dbz(:,ii);
end
end