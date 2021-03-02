% Ocean scan calibration for HCR data

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

% If 1, plots for individual calibration events will be made, if 0, only
% total plots will be made

project='otrec'; %socrates, aristo, cset, otrec
quality='qc1'; %field, qc1, or qc2
dataFreq='10hz';

b_drizz = 0.52; % Z<-17 dBZ
b_rain = 0.68; % Z>-17 dBZ
alpha = 0.21;
salinity=35; % Ocean salinity for sig0model in per mille (world wide default is 35) and sensitivity to that number is low
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

directories.figdir=['/scr/sci/romatsch/liquidWaterHCR/',project,'/'];

%directories.dataDir=HCRdir(project,quality,dataFreq);
directories.dataDir=['/run/media/romatsch/RSF0006/rsf/meltingLayer/',project,'/10hz/'];

startTime=datetime(2019,9,30,16,5,0);
endTime=datetime(2019,9,30,16,17,0);

%% Get data
%[data frq]=f_load_sort_data_nadir(directories.dataDir,startTime,endTime);

% get data
fileList=makeFileList(directories.dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

data=[];

data.DBZ = [];
data.DBMVC=[];
data.DBMHX=[];
data.PRESS=[];
data.TEMP=[];
data.RH=[];
data.SST=[];
data.TOPO=[];
data.U_SURF=[];
data.V_SURF=[];
data.FLAG=[];
data.pulse_width=[];
data.rotation=[];
data.tilt=[];

dataVars=fieldnames(data);

% Load data
data=read_HCR(fileList,data,startTime,endTime);

% Check if all variables were found
for ii=1:length(dataVars)
    if ~isfield(data,dataVars{ii})
        dataVars{ii}=[];
    end
end

dataVars=dataVars(~cellfun('isempty',dataVars));
xmitPowV=ncread(fileList{1},'r_calib_xmit_power_v');
antGainV=ncread(fileList{1},'r_calib_antenna_gain_v');
xmitPowH=ncread(fileList{1},'r_calib_xmit_power_h');
antGainH=ncread(fileList{1},'r_calib_antenna_gain_h');
beamWidthH=ncread(fileList{1},'radar_beam_width_h');
beamWidthV=ncread(fileList{1},'radar_beam_width_v');
waveGuideLossV=ncread(fileList{1},'r_calib_two_way_waveguide_loss_v');
radomeLossV=ncread(fileList{1},'r_calib_two_way_radome_loss_v');
waveGuideLossH=ncread(fileList{1},'r_calib_two_way_waveguide_loss_h');
radomeLossH=ncread(fileList{1},'r_calib_two_way_radome_loss_h');
recMismatchLoss=ncread(fileList{1},'r_calib_receiver_mismatch_loss');
ksquared=ncread(fileList{1},'r_calib_k_squared_water');

data.frq=ncread(fileList{1},'frequency');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sort out bad data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reflTemp=data.DBZ;
data.reflMask=ones(size(data.time));
elev=abs(data.elevation+90);

%sort out upward pointing
outElevInd=find(elev>90);
data.reflMask(outElevInd)=0;

% sort out data from below 2500m altitude
altInd=find(data.altitude<2500);
data.reflMask(altInd)=0;

% Remove bang
reflTemp(1:17,:)=nan;

% Find ocean surface gate
[data.refl maxGate]=max(reflTemp,[],1);

% Exclude data with max gate=1
max1=find(maxGate==1);
data.reflMask(max1)=0;

%Get the linear index of the maximum reflectivity value
maxGateLin=sub2ind(size(reflTemp),maxGate,1:size(reflTemp,2));

% Calculate reflectivity sum inside and outside ocean surface
reflLin=10.^(reflTemp./10);
reflOceanLin=nan(size(data.time));
reflNoOceanLin=nan(size(data.time));

for ii=1:length(data.time)
    if ~(maxGate(ii)<10 | maxGate(ii)>size(reflLin,1)-5)
        reflRay=reflLin(:,ii);
        reflOceanLin(ii)=sum(reflRay(maxGate(ii)-5:maxGate(ii)+5),'omitnan');
        reflNoOceanLin(ii)=sum(reflRay(1:maxGate(ii)-6),'omitnan');
    end
end

% Remove data where reflectivity outside of ocean swath is more than
% 0.8
tooMuchRefl=find(reflNoOceanLin>0.8);
data.reflMask(tooMuchRefl)=0;
data.reflMask(isnan(reflOceanLin))=0;

% remove data before and after drop outs
nanInds=find(isnan(data.refl));
if ~isempty(find(nanInds==1)) | ~isempty(find(nanInds==2)) | ~isempty(find(nanInds==3))
    nanInds=nanInds(4:end);
end
if ~isempty(find(nanInds==length(data.time))) | ~isempty(find(nanInds==length(data.time)-1)) | ~isempty(find(nanInds==length(data.time)-2))
    nanInds=nanInds(1:end-3);
end
data.reflMask(nanInds+1)=0;
data.reflMask(nanInds-1)=0;
data.reflMask(nanInds+2)=0;
data.reflMask(nanInds-2)=0;
data.reflMask(nanInds+3)=0;
data.reflMask(nanInds-3)=0;

windSpd=sqrt(data.U_SURF.^2+data.V_SURF.^2);

data.reflMask(data.TOPO>0)=0;
data.SST(data.TOPO>0)=nan;
windSpd(data.TOPO>0)=nan;

%% Remove data with clouds

%dbzMeasGood=data.refl;
%dbzMeasGood(data.reflMask==0)=nan;

data.DBMVC(data.reflMask==0)=nan;
data.DBMHX(data.reflMask==0)=nan;

%% Attenuation
[ituTot,alphaTot,layer_itu,layer_ituC]= get_gas_atten(data);

%% Calculate sig0
peak_power_dBmV =xmitPowV;
peak_powerV=10^((peak_power_dBmV)*0.1);

peak_power_dBmH =xmitPowH;
peak_powerH=10^((peak_power_dBmH)*0.1);

c=299792458;
wave_len = c/data.frq/1000; %transmit wavelength in km;

tx_ant_gainV=10^(antGainV/10);
tx_ant_gainH=10^(antGainH/10);

rc1V=peak_powerV*(tx_ant_gainV*tx_ant_gainV)*wave_len*wave_len;
rc1H=peak_powerH*(tx_ant_gainH*tx_ant_gainH)*wave_len*wave_len;
rc3=beamWidthH*beamWidthV*pi*pi/(180*180);
rc4= 512*log(2)*pi*pi;
Wband_RCV=10*log10(rc1V*rc3/rc4);
Wband_RCH=10*log10(rc1H*rc3/rc4);

lossV=waveGuideLossV+radomeLossV+recMismatchLoss;
lossH=waveGuideLossH+radomeLossH+recMismatchLoss;

%calculate sig0 measured 
avg_WV_lossL= 2*alphaTot; % Two-way WV attenuation in  dB
Wband_RC_revisedLV=Wband_RCV-lossV-avg_WV_lossL;
sig0measuredV=-Wband_RC_revisedLV'+data.DBMVC(maxGateLin)+20*log10(data.range(maxGateLin)/1000)-10*log10(cosd(elev));
Wband_RC_revisedLH=Wband_RCH-lossH-avg_WV_lossL;
sig0measuredH=-Wband_RC_revisedLH'+data.DBMHX(maxGateLin)+20*log10(data.range(maxGateLin)/1000)-10*log10(cosd(elev));

%% DBM range corrected
% dbmvcLin=10.^(data.DBMVC./10);
% dbmvcLinRange=dbmvcLin.*data.range.^2;
% DBMVC=10.*log10(dbmvcLinRange);
% 
% dbmhxLin=10.^(data.DBMHX./10);
% dbmhxLinRange=dbmhxLin.*data.range.^2;
% DBMHX=10.*log10(dbmhxLinRange);

if ~max(data.reflMask)==0    
    %% Plot lines
    close all
    f1 = figure('Position',[200 500 2000 1200],'DefaultAxesFontSize',12,'renderer','painters');
    
    subplot(3,1,1)
    hold on
    l1=plot(data.time,sig0measuredV,'-b','linewidth',1);
    l2=plot(data.time,sig0measuredH+20,'-r','linewidth',1);
    ylabel('Sig0 (dB)');
    %ylim([40 50]);
    grid on
    
    yyaxis right
    l3=plot(data.time,sig0measuredV-sig0measuredH,'-g','linewidth',1);
    %ylim([22 32]);
    ylabel('Sig0 diff. (dB)');
    ax = gca;
    ax.YColor = 'k';
        
    xlim([data.time(1),data.time(end)]);
    
    legend([l1 l2 l3],{'Sig0 V','Sig0 H + 20','Sig0 V - Sig0 H'},'location','southwest');
    
    title(['Reflectivity: ',datestr(data.time(1)),' to ',datestr(data.time(end))])
    
    subplot(3,1,2)
    
    hold on
    plot(data.time,data.U_SURF,'-b','linewidth',2);
    plot(data.time,data.V_SURF,'-g','linewidth',2);
    plot(data.time,windSpd,'-r','linewidth',2);
    plot(data.time,data.altitude./1000,'--k','linewidth',2);
    ylabel('Wind (m/s), Alt (km)');
    ylim([-15 15]);
    xlim([data.time(1),data.time(end)]);
    grid on
       
    legend({'U wind','V wind','Wind speed','Altitude'},'location','northeast');
       
    subplot(3,1,3)
    hold on
    plot(data.time,data.elevation,'-k','linewidth',2);
    ylabel('Elev(deg)');
    ylim([-90.4 -89.6]);
    
    ax = gca;
    ax.YColor = 'k';
    grid on
    
    yyaxis right
    plot(data.time,data.tilt,'-b','linewidth',2);
    plot(data.time,data.rotation-180,'-r','linewidth',2);
    xlim([data.time(1),data.time(end)]);
    ylabel('Rotation-180, Tilt (deg)');
    ylim([-20 20]);
    ax = gca;
    ax.YColor = 'k';
    
    legend({'Elevation','Tilt','Rotation-180'},'location','northeast');
        
    set(gcf,'PaperPositionMode','auto')
    print(f1,[directories.figdir,project,'_lines_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
    
end