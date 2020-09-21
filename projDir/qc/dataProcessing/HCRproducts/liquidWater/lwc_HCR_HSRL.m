%https://anl.box.com/s/0ow3ffpo35xbtph3clxumoz6v6rvosor

clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

indir='/Volumes/RSF-Vivek/SOCRATES/';
startTime=datetime(2018,2,24,2,23,0);
endTime=datetime(2018,2,24,2,31,0);

ylimits=[-0.0 2.0];

formatOut = 'yyyymmdd_HHMM';

fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

if isempty(fileList)
    error('No data found.')
end

%HCR data
data.HCR_DBZ=[];
data.HCR_VEL=[];
data.HCR_WIDTH=[];
data.HCR_LDR=[];
data.TEMP=[]; data.temp=data.TEMP+273.15;

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


[maskData antStat]=echoMask(data);
data.HCR_DBZ(maskData>1)=nan;

data.HCR_VEL(maskData>1)=nan;
data.HCR_WIDTH(maskData>1)=nan;
data.HCR_LDR(maskData>1)=nan;

data.HSRL_Aerosol_Backscatter_Coefficient(data.HSRL_Aerosol_Backscatter_Coefficient<0.5e-6)=nan;
close all
%calculate above sea level altitudes
asl=nan(size(data.range));
downInd=find(data.elevation<0);
upInd=find(data.elevation>=0);
asl(:,downInd)=-1*((data.range(:,downInd).*cosd(abs(data.elevation(downInd))-90)./1000)-data.altitude(downInd)./1000);
asl(:,upInd)=data.range(:,upInd).*cosd(abs(data.elevation(upInd))-90)./1000+data.altitude(upInd)./1000;

%% Find correct time

absolute_time=data.time;

alt_msl = data.altitude;
elevation=data.elevation;

dBZ = data.HCR_DBZ;
range = data.range;
beta = data.HSRL_Aerosol_Backscatter_Coefficient;

%% Calculate


Z_95_lin=10.^(dBZ*0.1);

wt_coef=zeros(size(dBZ));
wt_exp=zeros(size(dBZ));
Z_95_lin(dBZ < -200)=0.;
DBZ_temp=dBZ;
wt_coef(dBZ < - 20)=20.;
wt_exp(dBZ < - 20)=0.52;
wt_coef(-20 <dBZ <-15 )=1.73;
wt_exp(-20 <dBZ < -15 )=0.15;
wt_coef(dBZ > -15)=0.22;
wt_exp(dBZ > -15)=0.68;
att_cumul=2.*0.0192*cumsum((wt_coef.*Z_95_lin.^wt_exp),1,'omitnan');
att_cumul(dBZ < -200)=NaN;
dBZ_cor=dBZ+att_cumul;
Z_95_lin_cor=10.^(dBZ_cor*0.1);

cen_ind= dBZ== -9999 | beta  <= 0;
dBZ(cen_ind)= NaN;
beta(cen_ind) = NaN;

%% More variables

radar_lidar= dBZ_cor-10.*log10(beta);

radar_lidar_lin= 10.^(radar_lidar*0.1);

retrieved_RES=8.3478*(radar_lidar_lin).^0.3102; % in microns

retrieved_RLES=9.12*((radar_lidar_lin).^0.25);  % in microns


ref_normRES_retrieved=Z_95_lin_cor./((1e-3*retrieved_RES).^3); % retrieved Cld.RES is in microns and it is converted to mm by the fctor 1e-3

ref_normRLES_retrieved=Z_95_lin_cor./((1e-3*retrieved_RLES).^3.5);


dsdretrieved_LWC=5.238e-4*ref_normRES_retrieved+2.212e-4;

dsdretrieved_LWC_a=2.33e-6*(ref_normRLES_retrieved-0.031);

lin_beta=10.^(0.1*beta+3);

%% Plot retrieved quantities

f4=figure('DefaultAxesFontSize',12,'Position',[200 500 1000 800],'DefaultAxesFontWeight','bold');
colormap('jet');

subplot(2,1,1,'align');

surf(absolute_time,asl,log10(retrieved_RES),'edgecolor','none');
view(2)
hold on
plot(absolute_time,alt_msl,'k','LineWidth',2);
h = colorbar;
set(get(h,'title'),'string','log_{10}[RLED] microns');
grid on
ylabel('altitude (km)');
title({'Estimated Charactristic Diameter ';[datestr(absolute_time(1)) ' to ' datestr(absolute_time(end))]});
ylim(ylimits);

subplot(2,1,2,'align');

surf(absolute_time,asl,log10(dsdretrieved_LWC),'edgecolor','none');
view(2)
hold on
plot(absolute_time,alt_msl,'k','LineWidth',2);
h = colorbar; caxis([-3 0.5]);
set(get(h,'title'),'string','log_{10}[LWC] g/m^3');
grid on
ylabel('altitude (km)');
title('Liquid water content');
ylim(ylimits);

%% More plots

f5=figure('DefaultAxesFontSize',12,'Position',[200 500 1000 1200],'DefaultAxesFontWeight','bold');
colormap('jet');

subplot(3,1,1,'align');

surf(absolute_time,asl, dBZ,'edgecolor','none');
view(2)
hold on
plot(absolute_time,alt_msl,'k','LineWidth',2);
caxis([-40 20]);
h = colorbar;
set(get(h,'title'),'string','dBZ');
grid on
ylabel('altitude (km)');
title({'HCR Reflectivity';[datestr(absolute_time(1)) ' to ' datestr(absolute_time(end))]});
ylim(ylimits);

subplot(3,1,2,'align');

surf(absolute_time,asl,log10(beta),'edgecolor','none');
view(2)
hold on
plot(absolute_time,alt_msl,'k','LineWidth',2);
h = colorbar; caxis([-7 -2]);
set(get(h,'title'),'string','log_{10}[{\beta}] m^{-1} Sr^{-1}');
grid on
ylabel('altitude (km)');
title('HSRL Backscatter{\beta}');
ylim(ylimits);

subplot(3,1,3,'align');

surf(absolute_time,asl,dBZ_cor,'edgecolor','none');
view(2)
caxis([-40 20]);
hold on
plot(absolute_time,alt_msl,'k','LineWidth',2);
h = colorbar;
set(get(h,'title'),'string','dBZ');
grid on
ylabel('altitude (km)');
title('Attenuation corrected HCR reflectivity');
ylim(ylimits);

