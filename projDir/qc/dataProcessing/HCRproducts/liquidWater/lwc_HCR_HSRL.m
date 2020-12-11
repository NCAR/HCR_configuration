%https://anl.box.com/s/0ow3ffpo35xbtph3clxumoz6v6rvosor

clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

indir='/run/media/romatsch/RSF0006/rsf/combined_hcr_hsrl/socrates/';
startTime=datetime(2018,2,24,2,23,0);
endTime=datetime(2018,2,24,2,31,0);

ylimits=[-0.0 2.0];

formatOut = 'yyyymmdd_HHMM';

%% Load data

fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

if isempty(fileList)
    error('No data found.')
end

%HCR data
data.HCR_DBZ=[];
data.HCR_VEL=[];
data.HCR_WIDTH=[];
data.HCR_LDR=[];
data.FLAG=[];

%HSRL data
data.HSRL_Aerosol_Backscatter_Coefficient=[];
data.HSRL_Volume_Depolarization=[];
data.HSRL_Aerosol_Extinction_Coefficient=[];

dataVars=fieldnames(data);

% Load data
data=read_HCR(fileList,data,startTime,endTime);

% Check if all variables were found
for ii=1:length(dataVars)
    if ~isfield(data,dataVars{ii})
        dataVars{ii}=[];
    end
end

data.HCR_DBZ(data.FLAG>1)=nan;
data.HCR_VEL(data.FLAG>1)=nan;
data.HCR_WIDTH(data.FLAG>1)=nan;
data.HCR_LDR(data.FLAG>1)=nan;

data.HSRL_Aerosol_Backscatter_Coefficient(data.HSRL_Aerosol_Backscatter_Coefficient<0.5e-6)=nan;

%% Calculate attenuation correction

Z_95_lin=10.^(data.HCR_DBZ*0.1);
Z_95_lin(data.HCR_DBZ < -200)=0.;

wt_coef=zeros(size(data.HCR_DBZ));
wt_exp=zeros(size(data.HCR_DBZ));
wt_coef(data.HCR_DBZ < - 20)=20.;
wt_exp(data.HCR_DBZ < - 20)=0.52;
wt_coef(-20 <data.HCR_DBZ <-15 )=1.73;
wt_exp(-20 <data.HCR_DBZ < -15 )=0.15;
wt_coef(data.HCR_DBZ > -15)=0.22;
wt_exp(data.HCR_DBZ > -15)=0.68;

att_cumul=2.*0.0192*cumsum((wt_coef.*Z_95_lin.^wt_exp),1,'omitnan');
att_cumul(data.HCR_DBZ < -200)=NaN;

dBZ_cor=data.HCR_DBZ+att_cumul;
Z_95_lin_cor=10.^(dBZ_cor*0.1);

%% Censor data

cen_ind=data.HSRL_Aerosol_Backscatter_Coefficient  <= 0;

dBZ = data.HCR_DBZ;
dBZ(cen_ind)= NaN;
beta = data.HSRL_Aerosol_Backscatter_Coefficient;
beta(cen_ind) = NaN;

%% Calculate RES and LWC

radar_lidar= dBZ_cor-10.*log10(beta);
radar_lidar_lin= 10.^(radar_lidar*0.1);

retrieved_RES=8.3478*(radar_lidar_lin).^0.3102; % in microns
ref_normRES_retrieved=Z_95_lin_cor./((1e-3*retrieved_RES).^3); % retrieved Cld.RES is in microns and it is converted to mm by the fctor 1e-3
dsdretrieved_LWC=5.238e-4*ref_normRES_retrieved+2.212e-4;

%retrieved_RLES=9.12*((radar_lidar_lin).^0.25);  % in microns
%ref_normRLES_retrieved=Z_95_lin_cor./((1e-3*retrieved_RLES).^3.5);
%dsdretrieved_LWC_a=2.33e-6*(ref_normRLES_retrieved-0.031);

%lin_beta=10.^(0.1*beta+3);

%% Plot RES and LWC

close all

figure('DefaultAxesFontSize',12,'Position',[200 500 1000 800],'DefaultAxesFontWeight','bold');
colormap('jet');

subplot(2,1,1,'align');

surf(data.time,data.asl./1000,retrieved_RES./1000,'edgecolor','none');
view(2)
hold on
plot(data.time,data.altitude,'k','LineWidth',2);
h = colorbar;
caxis([0 0.3]);
set(get(h,'title'),'string','mm');
grid on
ylabel('altitude (km)');
title({'Estimated Charactristic Diameter ';[datestr(data.time(1)) ' to ' datestr(data.time(end))]});
ylim(ylimits);

subplot(2,1,2,'align');

surf(data.time,data.asl./1000,dsdretrieved_LWC,'edgecolor','none');
view(2)
hold on
plot(data.time,data.altitude,'k','LineWidth',2);
h = colorbar;
caxis([0 3]);
set(get(h,'title'),'string','g/m^3');
grid on
ylabel('altitude (km)');
title('Liquid water content');
ylim(ylimits);

%% Plot reflectivity and backscatter

figure('DefaultAxesFontSize',12,'Position',[200 500 1000 1200],'DefaultAxesFontWeight','bold');
colormap('jet');

subplot(3,1,1,'align');

surf(data.time,data.asl./1000, data.HCR_DBZ,'edgecolor','none');
view(2)
hold on
plot(data.time,data.altitude,'k','LineWidth',2);
caxis([-40 20]);
h = colorbar;
set(get(h,'title'),'string','dBZ');
grid on
ylabel('altitude (km)');
title({'HCR Reflectivity';[datestr(data.time(1)) ' to ' datestr(data.time(end))]});
ylim(ylimits);

subplot(3,1,2,'align');

surf(data.time,data.asl./1000,dBZ_cor,'edgecolor','none');
view(2)
caxis([-40 20]);
hold on
plot(data.time,data.altitude,'k','LineWidth',2);
h = colorbar;
set(get(h,'title'),'string','dBZ');
grid on
ylabel('altitude (km)');
title('Attenuation corrected HCR reflectivity');
ylim(ylimits);

subplot(3,1,3,'align');

surf(data.time,data.asl./1000,log10(beta),'edgecolor','none');
view(2)
hold on
plot(data.time,data.altitude,'k','LineWidth',2);
h = colorbar; caxis([-7 -2]);
set(get(h,'title'),'string','log_{10}[{\beta}] m^{-1} Sr^{-1}');
grid on
ylabel('altitude (km)');
title('HSRL Backscatter{\beta}');
ylim(ylimits);


