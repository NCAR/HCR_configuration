% Read HCR HSRL combined data

clear all;
close all;
fs = filesep;

QF_lab = {'CLD';'SPK';'EXT';'BLE';'OOG';'BAN';'WSE';'LSE';'BLS'; ...
           'NSC';'AIT';'MIS'};
cmp_QF = nab_color_table('./particle');
% load fit_RES_LWC.mat -mat;
load fit_RES_LWC_nofilt.mat -mat;

disp('Projects :');
disp('  1. CSET    2. SOCRATES');
sel_proj = input('Select: ');
if sel_proj == 1
  proj_name = 'CSET';
  air_thres = 0.98;
elseif sel_proj == 2
  proj_name = 'SOCRATES';
  air_thres = 0.98;
else
  disp('Unknown Project')
  return;
end

% startTime=datetime(2018,1,24,3,59,00);
% endTime=datetime(2018,1,24,4,00,00);

% startTime=datetime(2018,1,23,00,00,0);
% endTime=datetime(2018,1,23,01,00,0);
sTime = input('Input Start Date and Time as [yyyy mm dd hh mm ss]: ');
eTime = input('Input End   Date and Time as [yyyy mm dd hh mm ss]: ');

startTime=datetime(sTime);
endTime=datetime(eTime);

ylimits=[-0.0 2.5];

plotlidars=0; % 1 to plot lidar data, 0 to not plot lidar
plotradars=0; % 1 to plot radar data, 0 to not plot radar
saveStatis=0; % 1 Save sea surface statistic data
filltemp=1; % 1 to fill in temperature data with workaround to masking issue

%addpath('/h/eol/romatsch/git/private/utils/');
%addpath('/h/eol/romatsch/git/private/process_HSRL/');

% indir='/scr/rain1/rsfdata/projects/socrates/hcr/qc/cfradial/hsrl_merge/2hz/';
% modeldir='/scr/sci/romatsch/data/reanalysis/ecmwf/era5interp/socrates/2hz/';
%figdir='/scr/sci/romatsch/HCR_HSRL/';

%  indir='/scr/rain1/rsfdata/projects/socrates/hcr/qc/cfradial/hsrl_merge/2hz/';
%  modeldir='/scr/sci/romatsch/data/reanalysis/ecmwf/era5interp/socrates/2hz/';
indir = sprintf('.%cData_%s%cmerge%c',fs,proj_name,fs,fs);
indir_hcr = sprintf('.%cData_%s%chcr%c',fs,proj_name,fs,fs);

 %indir='/Volumes/RSF-Vivek/ecmwf/era5interp/socrates/hcr/qc/cfradial/hsrl_merge/2hz/';
% modeldir='/Volumes/RSF-Vivek/ecmwf/era5interp/socrates/2hz/';
formatOut = 'yyyymmdd_HHMM';

fileList1=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
fileList2=makeFileList(indir_hcr,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

if ~isempty(fileList1) & ~isempty(fileList2)
    
%  HCR data
  data.HCR_DBZ=[];
  data.HCR_VEL=[];
  data.HCR_WIDTH=[];
  data.HCR_LDR=[];
    
%  HSRL data
  data.HSRL_Aerosol_Backscatter_Coefficient=[];
  data.HSRL_Volume_Depolarization=[];
  data.HSRL_Aerosol_Extinction_Coefficient=[];
  
%  HCR only
  data_hcr.DBZ=[];
  data_hcr.FLAG =[];
%   data_hcr.azimuth = [];
%   data_hcr.tilt    = [];
  
%  Others
  data.FLAG=[];
  data.TEMP = [];
  data.azimuth = [];
  data.tilt    = [];
      
  dataVars=fieldnames(data);
  dataVars_hcr=fieldnames(data_hcr);
        
% Load data
  data=read_HCR(fileList1,data,startTime,endTime);
  data_hcr=read_HCR(fileList2,data_hcr,startTime,endTime);  
    
% Check if all variables were found
  for ii=1:length(dataVars)
    if ~isfield(data,dataVars{ii})
      dataVars{ii}=[];
    end
  end
  for ii = 1:length(dataVars_hcr)
    if ~isfield(data_hcr,dataVars_hcr{ii})
      dataVars_hcr{ii}=[];
    end
  end
%   return
%   Convert HCR data only to merge time resolution
  t1 = datenum(data.time);
  t2 = datenum(data_hcr.time);
  DBZ = NaN*ones(size(data.HCR_DBZ));
  DBZ2 = data_hcr.DBZ;
  DBZ2(data_hcr.FLAG==2)=NaN;
  DBZ2 = 10.^(0.1*DBZ2);
  disp('Align HCR-only DBZ to merge DBZ');
  for n = 1:size(DBZ,2)
    try
      dt = (t1(n+1)-t1(n))/2;
    catch
      dt = (t1(n)-t1(n-1))/2;
    end
    it = find(t2>t1(n)-dt & t2<=t1(n)+dt);
    mdbz = nanmean(DBZ2(:,it),2);
    DBZ(:,n) = 10*log10(mdbz);
  end
  clear DBZ2 t1 t2;
  
%   convert datetime to second of the day
  [year mon day hh mm ss] = datevec(data.time(1));
  Nday = datenum(data.time)-datenum([year mon day]);
  Nsec = 24*60*60*Nday;
  
%  Using data.FLAG instead of echoMask
%   comment out echoMask
%   [maskData antStat]=echoMask(data);
%   data.HCR_DBZ(maskData>1)=nan;  
%   data.HCR_VEL(maskData>1)=nan;
%   data.HCR_WIDTH(maskData>1)=nan;
%   data.HCR_LDR(maskData>1)=nan;

  beta = data.HSRL_Aerosol_Backscatter_Coefficient;
  sldBZ=gradient(data.HCR_DBZ(:,1),0.2,5);
  Z_95_lin=10.^(data.HCR_DBZ*0.1);
  Z_95 = data.HCR_DBZ;
  Z_95(Z_95<=-50 | data.FLAG == 2) = NaN;
  dR = (data.range(2)-data.range(1))/1000;
  
  cen_ind= data.HCR_DBZ== -9997.5 | beta  <= 1e-9 | beta > 0.001;
%   cen_ind= data.HCR_DBZ== -9997.5 | beta > 0.001;
  data.HCR_DBZ(cen_ind)= NaN;
  beta(cen_ind) = NaN;

  wt_coef=zeros(size(data.HCR_DBZ));
  wt_exp=zeros(size(data.HCR_DBZ));
  Z_95_lin(data.HCR_DBZ < -200)=0.;
  DBZ_temp=data.HCR_DBZ;
  wt_coef(data.HCR_DBZ < - 20)=20.;
  wt_exp(data.HCR_DBZ < - 20)=0.52;
  wt_coef(-20 <data.HCR_DBZ <-15 )=1.73;
  wt_exp(-20 <data.HCR_DBZ < -15 )=0.15;
  wt_coef(data.HCR_DBZ > -15)=0.22;
  wt_exp(data.HCR_DBZ > -15)=0.68;
  
%  When elevation angles are -90 the radar is looking down
%  When elevation angles are  90 the radar is looking up
  i_el = find(data.elevation<80);
  att_cumul=2.*dR*cumsum((wt_coef.*Z_95_lin.^wt_exp),1,'omitnan');
  att_cumul(i_el)=2.*dR*cumsum((wt_coef(i_el).*Z_95_lin(i_el).^wt_exp(i_el)), ...
                             1,'reverse','omitnan');
  att_cumul(data.HCR_DBZ < -200)=NaN;
  
  dBZ_cor=data.HCR_DBZ+att_cumul;
  Z95_cor=Z_95+att_cumul;
  Z_95_lin_cor=10.^(dBZ_cor*0.1);
  att_95 = wt_coef.*Z_95_lin.^wt_exp;
  
%% More variables

  radar_lidar= dBZ_cor-10.*log10(beta);
  radar_lidar_lin= 10.^(radar_lidar*0.1);
  retrieved_RLES=9.12*((radar_lidar_lin).^0.25);  % in microns
  ref_normRLES_retrieved=Z_95_lin_cor./((0.53*1e-3*retrieved_RLES).^3.74);
  dsdretrieved_LWC_a=2.33e-6*ref_normRLES_retrieved+0.004;
  dsdretrieved_LWC_a(dsdretrieved_LWC_a<0) = NaN;

  lin_beta=10.^(0.1*beta+3);
  
%% era5 data
    
%     model.temp=[];
%     
%     model=read_model(model,modeldir,data.time(1),data.time(end));
%     
  data.temp=data.TEMP+273.15;
%% Calculate variables
  close all
  asl=nan(size(data.range));
  downInd=find(data.elevation<0);
  upInd=find(data.elevation>=0);
  asl(:,downInd)=-1*((data.range(:,downInd).*cosd(abs(data.elevation(downInd))-90)./1000)-data.altitude(downInd)./1000);
  asl(:,upInd)=data.range(:,upInd).*cosd(abs(data.elevation(upInd))-90)./1000+data.altitude(upInd)./1000;
    
  backscatLog = real(log10(data.HSRL_Aerosol_Backscatter_Coefficient));
  extLog = real(log10(data.HSRL_Aerosol_Extinction_Coefficient));
  depolLog = real(log10(data.HSRL_Volume_Depolarization));
    
  lidarRatio=10.^(extLog-backscatLog);
    
%  Particle ID
  vol_depol=data.HSRL_Volume_Depolarization./(2-data.HSRL_Volume_Depolarization);
  lin_depol=vol_depol./(2-vol_depol);
%     pid=aer_cld2c(data.HSRL_Aerosol_Backscatter_Coefficient,lin_depol,data.temp);
  pid=aer_cld2c_ver1(data.HSRL_Aerosol_Backscatter_Coefficient,lin_depol,data.temp);
    
%v[pidhcr,m]=precip_cld_a(data.HCR_DBZ,data.HCR_LDR,data.HCR_VEL_CORR,data.HCR_WIDTH_CORR,data.temp);
    
%   [pidhcr,m]=precip_cld_a(dBZ_cor,data.HCR_LDR,data.HCR_VEL,data.HCR_WIDTH,data.temp);
  [pidhcr,m]=precip_cld_a(Z95_cor,data.HCR_LDR,data.HCR_VEL,data.HCR_WIDTH,data.temp);
%vpidhcr(isnan(data.HCR_DBZ))=nan;
    
  pidhcr(isnan(Z95_cor))=nan;
    
  pid_comb=combine_pid_hcr_hsrl(pidhcr,pid,data.temp);
  
%  My SRT parts -----------------------------------------------------------
%   asl_hcr=nan(size(data_hcr.range));
%   downInd_hcr=find(data_hcr.elevation<0);
%   upInd_hcr=find(data_hcr.elevation>=0);
%   asl_hcr(:,downInd_hcr)=-1*((data_hcr.range(:,downInd_hcr).* ...
%       cosd(abs(data_hcr.elevation(downInd_hcr))-90)./1000)- ...
%       data_hcr.altitude(downInd_hcr)./1000);
%   asl_hcr(:,upInd_hcr)=data_hcr.range(:,upInd_hcr).* ...
%       cosd(abs(data_hcr.elevation(upInd_hcr))-90)./1000+ ...
%       data_hcr.altitude(upInd_hcr)./1000;

  alt_hcr   = data.altitude/1000;
  range_hcr = data.range/1000;   
  ele_hcr   = data.elevation;
  azi_hcr   = data.azimuth;
  tilt_hcr  = data.tilt;
  
%   DBZ = data_hcr.DBZ;   
  flag = data.FLAG;
%   DBZ(flag==2)=NaN;

%  Find the surface reflectivity for clear air
%   i_clear = 0 :: Cannot find sea surface reflectivity
%           = 1 :: Clear air with nadir
%           = 2 :: With precipitation and nadir
%           = 3 :: Clear air and off nadir
%           = 4 :: Precipitation and off nadir
%           = 5 :: Signal is attenuated out
  n_rays = size(DBZ,2);
  i_clear = zeros(1,n_rays);  sea_gates = i_clear;
  surInd  = i_clear;   altInd = i_clear;  
  DBZ_sea = NaN*i_clear;   
  
  rat0 = 0;
  h = waitbar(rat0,'Finding the clear air.');
  for n = 1:n_rays
    is = find(flag(:,n)==7);   ia = find(flag(:,n)==6);
    if isempty(ia)~=1 & isempty(is)~=1
      itest = find(isnan(DBZ(is,n))~=1);
      if isempty(itest) == 1
        continue;
      end
% start and end gate for attenuation correction
      m1 = is(1)+2;   m2 = ia(end)-2;  
      ray_len = (is(1))-(ia(end))+1;
      Idnan = find(isnan(DBZ(ia(end)+2:is(1)-2,n))==1);
      sea_gates(n) = length(is);
      iz = find(DBZ(is,n)==max(DBZ(is,n)));
      surInd(n) = is(1);
      altInd(n) = ia(end);
      DBZ_sea(n) = DBZ(is(iz(1)),n);
      if abs(abs(ele_hcr(n))-90)>=0.5
        if length(Idnan)/ray_len>=air_thres
          i_clear(n) = 3;
        else
          i_clear(n) = 4;
        end
        continue;
      end
      if length(Idnan)/ray_len>=air_thres
        i_clear(n) = 1;
      else
        i_clear(n) = 2;
      end
    end

%  For signal extinct
    ii = find(flag(:,n)==3);
    if isempty(ii)~=1 & isnan(altInd(n))~=1
      surInd(n) = ii(1);
      altInd(n) = ia(end);
%     DBZ_sea(n) = 1.14;                  % One STD based on cases
      dbz = DBZ(altInd(n):surInd(n),n);
      dbz = dbz(isnan(dbz)~=1);
      DBZ_sea(n) = max([dbz(end) 1.14]);
      i_clear(n) = 5;
    end
    rat = n/n_rays;
    if rat-rat0>0.001
      rat0 = rat;
      message = sprintf('Finding the clear air.  Finish %5.2f%%.',100*rat0);
      waitbar(rat0,h,message);
    end
  end
  close(h);
  
%   clear air DBZ
  clrDBZ = DBZ_sea(i_clear==1);   clrDBZ = clrDBZ(isnan(clrDBZ)~=1);
  n_clear = length(clrDBZ);
  clr_rat = 100*n_clear/n_rays;
  m_DBZ = nanmean(clrDBZ);   s_DBZ = nanstd(clrDBZ);

  if clr_rat>1
    surDBZ0 = m_DBZ;   std_DBZ0 = s_DBZ;
  else
    surDBZ0 = 47.98;   std_DBZ0 = 0.62;     % Default from Cases 
  end
  
% Output statistic data of clear sky
  disp(sprintf('# of ray::%d    %% of clear ray::%7.2f',n_rays,clr_rat));
  disp('DBZ');
  disp('==============================================================');
  disp(sprintf('Mean::%6.2f    STD::%6.2f',m_DBZ,s_DBZ));
  
%  Attenuation Correction using ZPHI
%    only the ray is looking down can be using SRT for attenuation
%    correction
%   gasDZ(isnan(gasDZ)==1)=0;   gasPia(isnan(gasPia)==1) = 0;
%   DZ = surDBZ0-DBZ_sea-gasDZ; 
  DZ = surDBZ0-DBZ_sea; 
  DZ(upInd) = 0;   DZ(DZ<=s_DBZ & i_clear==2)=2*s_DBZ;   DZ(i_clear==1)=0;
  [Zc,PIA,Aw,b] = back_corr_W(DBZ,DZ,range_hcr,altInd,surInd);
  DBZ_cor = DBZ+PIA;
  [Zc_m,PIA_m,Aw_m,b_m] = back_corr_W(Z_95,DZ,range_hcr,altInd,surInd);
  Z_95_cor = Z_95+PIA_m;
  zw = 10.^(0.1*DBZ_cor);
  Xres = (zw./Aw).^(1/3); 

% Load fitting data
zw_bnd_c  = interp1(Awbnd_c,zwbnd_c,Aw,'linear','extrap');
zw_bnd_c(zw_bnd_c<0)=0;
zw_bnd_dc = interp1(Awbnd_dc,zwbnd_dc,Aw,'linear','extrap');
zw_bnd_dc(zw_bnd_dc<0)=0;
LWC = NaN*ones(size(DBZ_cor));   RES = LWC;   

% Compute LWC, RES and classify PID
%   Cloud => PID=1;   Drizzle => PID=2;   Mixed => PID = 3;
%   Rain  => PID=4; (Not able to compute LWC and RES yet)  
%   Bang  => PID = 5;  Sea Surface => PID = 6;
%   Beyend (below sea or above radar => PID = 7;
%   Un-reasonable => PID = 8;
%   clear air (gas) => PID = NaN;
PID = NaN*ones(size(LWC));   
PD_lab = {'CLD';'DZL';'MIX';'RAN';'BAN';'WSE';'BLS';'URN'};
%   Cloud
% II = (Aw>=Aw_bnd_c)&(Aw>4)|(Aw>=Aw_bnd_c)&(Xres>0.35);
II = (zw<=zw_bnd_c)&(Aw>7.5);
PID(II)=8;
% II = (Aw>=Aw_bnd_c)&(Aw<4)&(Xres<0.35);
II = ((zw<=zw_bnd_c)&(Aw<=7.5))|(zw<min(zwbnd_dc) & Aw<=7.5);
RES(II) = P_c(1)*Xres(II)+P_c(2);
LWC(II) = Q_c(1)*Aw(II)+Q_c(2);
PID(II) = 1;

%   Mix
II = ((zw>zw_bnd_c)&(zw<zw_bnd_dc))&(Xres>Xres_thres);
PID(II)=8;
II = ((zw>zw_bnd_c)&(zw<zw_bnd_dc))&(Xres<=Xres_thres)&(Aw<=10);
%     Compute the Ratio
LWC1 = NaN*ones(size(LWC));   LWC2 = LWC1;   RES1 = LWC1;   RES2 = LWC1;
Ratio = abs(zw-zw_bnd_dc)./(abs(zw-zw_bnd_c)+abs(zw-zw_bnd_dc));
RES1(II) = P_c(1)*Xres(II)+P_c(2);
RES2(II) = P_dr1(1)*Xres(II)+P_dr1(2);
LWC1(II) = Q_c(1)*Aw(II)+Q_c(2);
LWC2(II) = Q_dr1(1)*Aw(II)+Q_dr1(2);
LWC(II)  = Ratio(II).*LWC1(II)+(1-Ratio(II)).*LWC2(II);
RES(II)  = Ratio(II).*RES1(II)+(1-Ratio(II)).*RES2(II);
PID(II)=3;
II = ((zw>zw_bnd_c)&(zw<zw_bnd_dc))&(Aw>10);
PID(II)=4;

%   Drizzle and light rain
II = (zw>=zw_bnd_dc & zw>=min(zw_bnd_dc))&(Xres>Xres_thres);
PID(II) = 4;
II = (zw>=zw_bnd_dc)&(Xres<=Xres_thres);
RES(II) = P_dr1(1)*Xres(II)+P_dr1(2);
LWC(II) = Q_dr1(1)*Aw(II)+Q_dr1(2);
PID(II)=2;

%  Aw=0
II = Aw==0 & zw>0 & zw<=zwc0_bnd;
PID(II) = 1;
LWC(II) = Q_c0(1)*zw(II).^Q_c0(2);
RES(II) = P_c0(1)*zw(II).^P_c0(2);

II = Aw==0 & zw>zwc0_bnd & zw<zwdr0_bnd;
PID(II) = 2;
LWC(II) = Q_dr0(1)*zw(II).^Q_dr0(2);
RES(II) = P_dr0(1)*zw(II).^P_dr0(2);

II = Aw==0 & zw>zwdr0_bnd;
PID(II) = 8;
LWC(II) = NaN;
RES(II) = NaN;

%  Others
PID(flag==6)=5;   PID(flag==7)=6;   PID(flag==3)=NaN;
PID(flag==9)=7;
LWC(flag==3 | flag==6 | flag==7 | flag==9)=NaN;
RES(flag==3 | flag==6 | flag==7 | flag==9)=NaN;

LWC(LWC<0) = 0;   RES(RES<0) = 0;

%  LWP
LWC_p = LWC;   LWC_mp = dsdretrieved_LWC_a;
% LWC_mp(Aw==0) = NaN;
LWC_p(isnan(LWC_mp)==1) = NaN;
LWP_merge = 1e3*dR*nansum(dsdretrieved_LWC_a,1);
LWP_hcr   = 1e3*dR*nansum(LWC_p,1);

%% Plot
  close all
  if plotlidars==1
    f1=figure('DefaultAxesFontSize',12);
    set(f1,'Position',[400 45 1000 700]);
        
    colormap(jet);
        
% Vel Raw
    subplot(3,1,1)
    hold on
    fig1=surf(data.time,asl,backscatLog);
    fig1.EdgeColor='none';
    ylim(ylimits);
    xlim([startTime,endTime]);
%         xlim('auto');
    caxis([-7 -2]);
    colorbar;
    view(2);
    ylabel('Altitude (km)');
    title(['Log10 aerosol backscatter coefficient (m^-1 sr^-1)'],'interpreter','none');
            
    subplot(3,1,2)
    hold on
    fig1=surf(data.time,asl,data.HSRL_Volume_Depolarization);
    fig1.EdgeColor='none';
    ylim(ylimits);
    caxis([0 0.8]);
    xlim([startTime,endTime]);
    colorbar;
    view(2);
    ylabel('Altitude (km)');
    title(['Aerosol Linear depolarization ratio'],'interpreter','none');
        
    subplot(3,1,3)
    hold on
    fig1=surf(data.time,asl,data.temp);
    fig1.EdgeColor='none';
    ylim(ylimits);
    xlim([startTime,endTime]);
    caxis([260 280]);
    colorbar;
    view(2);
    ylabel('Altitude (km)');
    title(['Temp, K'],'interpreter','none');
        
    cscale=[1,1,1;1,0.67,0;1,0.67,0;0,1,1;0,1,1;1,0,1;0,1,0];
    units_str={'No Signal', 'Aerosol','', 'Ice', '', 'SLD', 'Cloud droplets'};
        
    f2=figure('DefaultAxesFontSize',12);
    set(f2,'Position',[500 45 1000 700]);    
    colormap(jet);
           
    subplot(3,1,1)
    hold on
    fig1=surf(data.time,asl,backscatLog);
    fig1.EdgeColor='none';
    ylim(ylimits);
    xlim([startTime,endTime]);
    caxis([-8 -2]);
    colorbar;
    view(2);
    ylabel('Altitude (km)');
    title(['Log10 aerosol backscatter coefficient (m^-1 sr^-1)'],'interpreter','none');
    subplot(3,1,2)
    hold on
% fig1=surf(data.time,asl,data.HSRL_Volume_Depolarization);
    fig1=surf(data.time,asl,depolLog);
    fig1.EdgeColor='none';
    ylim(ylimits);
    caxis([-2 0]);
    xlim([startTime,endTime]);
    colorbar;
    view(2);
    ylabel('Altitude (km)');
    title(['Aerosol Linear depolarization ratio'],'interpreter','none');
        
    subplot(3,1,3)
    hold on
    fig1=surf(data.time,asl,pid);
    fig1.EdgeColor='none';
    ylim(ylimits);
    xlim([startTime,endTime]);
    caxis([.5 7.5]);
    cb=colorbar;
    colormap(gca,cscale);
    cb.Ticks=1:7;
    cb.TickLabels=units_str;
    view(2);
    ylabel('Altitude (km)');
    title(['Particle ID HSRL'],'interpreter','none');
  end
    
  if plotradars==1
        
    f3=figure('DefaultAxesFontSize',16);
    set(f3,'Position',[600 45 1000 700]);
    set(gca,'fontname','Timesnewroman','FontSize',14','fontweight','b');
%         colormap(jet);
        
% Vel Raw
    subplot(3,1,1)
    set(gca,'FontSize',18);
    set(gca,'FontWeight','bold');
    fig1=surf(data.time,asl,data.HCR_DBZ);
    fig1.EdgeColor='none';
    ylim(ylimits);
    xlim([startTime,endTime]); set(gca,'xticklabel',[]); %set(gca,'xtick',[])
    caxis([-30 15]);
    colorbar;
    view(2);
    ylabel('Altitude (km)','fontsize',10,'fontweight','b');
    title('Radar reflectivity (dBZ)  ', ...
          'interpreter','none');
    colormap(gca,cid_cmap2)
    set(gca,'TickLabelInterpreter','none');
    set(gca,'fontweight','bold','fontsize',16);
        
    subplot(3,1,2)
    fig1=surf(data.time,asl,data.HSRL_Volume_Depolarization);
    fig1.EdgeColor='none';
    ylim(ylimits);
    caxis([0 0.8]);
    xlim([startTime,endTime]);set(gca,'xticklabel',[]);
    colorbar;
    view(2);
    ylabel('Altitude (km)','fontsize',10,'fontweight','b');
    title(['Lidar depolarization ratio'],'interpreter','none');
    set(gca,'TickLabelInterpreter','none');
    set(gca,'fontweight','bold','fontsize',16);
        
    subplot(3,1,3)    
    cscale=[1,1,1; 0,0,1; 0,1,0.; 1,0,0; 1,0,1; 0,1,1; 1,1,0; 0.5,0,0; 1,0.67,0];
    units_str={'No signal','Droplets above freezing', 'Drizzle','Rain', 'Droplets below freezing',... 
            'Ice crystals', 'Rimed ice', 'Wet snow', 'Aerosol'};
    fig1=surf(data.time,asl,pid_comb);
    fig1.EdgeColor='none';
    ylim(ylimits);
    xlim([startTime,endTime]);
    caxis([.5 9.5]);
    cb=colorbar;
    colormap(gca,cscale);
    cb.Ticks=1:9;
    cb.TickLabels=units_str;
    view(2);
    ylabel('Altitude (km)','fontsize',10,'fontweight','b');
    xlabel('Time (UTC)','fontsize',10,'fontweight','b')
    title(['Particle ID Combined'],'interpreter','none');
    set(gca,'TickLabelInterpreter','none');
    set(gca,'fontweight','bold','fontsize',16);
        
    cscale=[1,1,1; 0,0,.5; 0,1,0.; 1,0,0; 1,0,1; 0,1,1; 1,1,0; 0.5,0,0];
        
    units_str={'No signal','Cloud droplets', 'Drizzle','Rain', 'SLD', 'Ice crystals', 'Aggregates', 'Wet snow/rimed ice'};
        
    f4=figure('DefaultAxesFontSize',16);
    set(f4,'Position',[700 45 1000 700]);
        
%     colormap(jet);
        
    subplot(3,1,1)
    hold on
    fig1=surf(data.time,asl,data.HCR_DBZ);
    fig1.EdgeColor='none';
    ylim(ylimits);
    xlim([startTime,endTime]);
    caxis([-30 15]);
    colorbar;
    view(2);
    ylabel('Altitude (km)');
    title('Reflectivity (dBZ)  ', ...
          'interpreter','none');
    colormap(gca,cid_cmap2)
        
    subplot(3,1,2)
    hold on
    fig1=surf(data.time,asl,data.HCR_VEL);
    fig1.EdgeColor='none';
    ylim(ylimits);
    xlim([startTime,endTime]);
    caxis([-2 6]);
    colorbar;
    view(2);
    ylabel('Altitude (km)');
    title(['Vel m/s)'],'interpreter','none');
    colormap(gca,cid_cmap)
            
    subplot(3,1,3)
    hold on
    fig1=surf(data.time,asl,pidhcr);
%     fig1=surf(data.time,asl,pid_comb);
    fig1.EdgeColor='none';
    ylim(ylimits);
    xlim([startTime,endTime]);
    caxis([.5 8.5]);
    cb=colorbar;
    colormap(gca,cscale);
    cb.Ticks=1:8;
    cb.TickLabels=units_str;
    view(2);
    ylabel('Altitude (km)');
    title(['Particle ID HCR'],'interpreter','none');
        
    f5=figure('DefaultAxesFontSize',12);
    subplot(2,1,1)
    hold on
    fig1=surf(data.time,asl,data.HCR_DBZ);
    fig1.EdgeColor='none';
    ylim(ylimits);
    xlim([startTime,endTime]);
    caxis([-30 15]);
    colorbar;
    view(2);
    ylabel('Altitude (km)');
    title('Reflectivity (dBZ)  ', ...
          'interpreter','none');
        
    subplot(2,1,2)
    hold on
    fig1=surf(data.time,asl,data.temp);
    fig1.EdgeColor='none';
    ylim(ylimits);
    xlim([startTime,endTime]);
    caxis([260 280]);
    colorbar;
    view(2);
    ylabel('Altitude (km)');
    title(['Temp, K'],'interpreter','none');
  end 

%   Reflectivity and beta
f1=figure('DefaultAxesFontSize',12);
set(f1,'Position',[400 45 1000 700]);        
colormap(jet);
        
% beta
subplot(3,1,1)
hold on
fig1=surf(data.time,asl,backscatLog);
fig1.EdgeColor='none';
ylim(ylimits);
xlim([startTime,endTime]);
%         xlim('auto');
caxis([-9 -3]);
hcb = colorbar;
view(2);
ylabel('Altitude (km)');
title('Log_{10} aerosol backscatter coefficient');
title(hcb,'m^{-1} sr^{-1}}');
%  HCR-HSRL HCR reflectivity
subplot(3,1,2)
hold on
fig1=surf(data.time,asl,Z_95);
fig1.EdgeColor='none';
ylim(ylimits);
xlim([startTime,endTime]);
%         xlim('auto');
caxis([-50 30]);
hcb = colorbar;
view(2);
ylabel('Altitude (km)');
title('HCR-HSRL merged HCR Reflectivity');
title(hcb,'dBZ');
%  HCR-Only HCR reflectivity
subplot(3,1,3)
hold on
fig1=surf(data.time,asl,DBZ);
fig1.EdgeColor='none';
ylim(ylimits);
xlim([startTime,endTime]);
%         xlim('auto');
caxis([-50 30]);
hcb = colorbar;
view(2);
ylabel('Altitude (km)');
title('HCR-Only HCR Reflectivity');
title(hcb,'dBZ');

%   Corrected Reflectivity
f2=figure('DefaultAxesFontSize',12);
set(f2,'Position',[400 45 1000 700]);        
colormap(jet);
        
%  HCR-HSRL HCR reflectivity
subplot(2,1,1)
hold on
fig1=surf(data.time,asl,Z95_cor);
fig1.EdgeColor='none';
ylim(ylimits);
xlim([startTime,endTime]);
%         xlim('auto');
caxis([-50 30]);
hcb = colorbar;
view(2);
ylabel('Altitude (km)');
title('HCR-HSRL merged HCR Corrected Reflectivity');
title(hcb,'dBZ');
%  HCR-Only HCR reflectivity
subplot(2,1,2)
hold on
fig1=surf(data.time,asl,DBZ_cor);
fig1.EdgeColor='none';
ylim(ylimits);
xlim([startTime,endTime]);
%         xlim('auto');
caxis([-50 30]);
hcb = colorbar;
view(2);
ylabel('Altitude (km)');
title('HCR-Only HCR Corrected Reflectivity');
title(hcb,'dBZ');
 
%   Compare PID
  figure('DefaultAxesFontSize',12);
  set(gcf,'Position',[750 45 1000 700]);
  cscale=[1,1,1; 0,0,.5; 0,1,0.; 1,0,0; 1,0,1; 0,1,1; 1,1,0; 0.5,0,0];
  units_str={'No signal','Cloud droplets', 'Drizzle','Rain', 'SLD', 'Ice crystals', 'Aggregates', 'Wet snow/rimed ice'};
  subplot(3,1,2)
  hold on
  fig1=surf(data.time,asl,pidhcr);
  fig1.EdgeColor='none';
  ylim(ylimits);
  xlim([startTime,endTime]);
  caxis([.5 8.5]);
  cb=colorbar;
  colormap(gca,cscale);
  cb.Ticks=1:8;
  cb.TickLabels=units_str;
  view(2);
  ylabel('Altitude (km)');
  title(['Particle ID HCR'],'interpreter','none');

  cscale=[1,1,1; 0,0,.5; 0,1,0.; 1,0,0; 1,0,1; 0,1,1; 1,1,0; 0.5,0,0; 1,0.67,0];
  units_str={'No signal','Cloud droplets', 'Drizzle','Rain', 'SLD', 'Ice crystals', 'Aggregates', 'Wet snow/rimed ice', 'Aerosol'};
  subplot(3,1,1)
  hold on
  fig1=surf(data.time,asl,pid_comb);
  fig1.EdgeColor='none';
  ylim(ylimits);
  xlim([startTime,endTime]);
  caxis([.5 9.5]);
  cb=colorbar;
  colormap(gca,cscale);
  cb.Ticks=1:9;
  cb.TickLabels=units_str;
  view(2);
  ylabel('Altitude (km)');
  title(['Particle ID Combined'],'interpreter','none');

  subplot(3,1,3)
  hold on
  fig1=surf(data.time,asl,PID);
  fig1.EdgeColor='none';
  ylim(ylimits);
  xlim([startTime,endTime]);
  caxis([1 9]);
  ic = [2 3 4 5 7 8 10 13];
  colormap(gca,cmp_QF(ic,:));
  hcb = colorbar;
  hcb.Ticks = 1.5:9.5;
  hcb.TickLabels = PD_lab;
  view(2);
  ylabel('Altitude (km)');
  title(['Particle ID'],'interpreter','none'); 
  
% Compare RES
  figure('DefaultAxesFontSize',12);
  set(gcf,'Position',[750 45 1000 700]);
  subplot(211)
  f1 = pcolor(data.time,asl,log10(retrieved_RLES*1e-3));
  shading 'flat';
  ylim(ylimits);
  xlim([startTime,endTime]);
  ylabel('Altitude (km)');
  title('HCR-HSRL RLES')
  cb = colorbar;
  Tick = cb.Ticks;
  for n = 1:length(Tick)
    TickLab{n} = sprintf('10^{%4.1f}',Tick(n)); 
  end
  cb.TickLabels = TickLab;
  title(cb,'mm');
  cv = caxis;
  
  subplot(212)
  f1 = pcolor(data.time,asl,log10(RES));
  shading 'flat';
  ylim(ylimits);
  xlim([startTime,endTime]);
  ylabel('Altitude (km)');
  title('HCR RES');
  caxis(cv);
  cb = colorbar;
  Tick = cb.Ticks;
  for n = 1:length(Tick)
    TickLab2{n} = sprintf('10^{%4.1f}',Tick(n)); 
  end
  cb.TickLabels = TickLab2;
  title(cb,'mm');
  
% Compare LWC
  figure('DefaultAxesFontSize',12);
  set(gcf,'Position',[750 45 1000 700]);
  subplot(211)
  f1 = pcolor(data.time,asl,log10(dsdretrieved_LWC_a));
  shading 'flat';
  ylim(ylimits);
  xlim([startTime,endTime]);
  ylabel('Altitude (km)');
  title('HCR-HSRL LWC')
  cb = colorbar;
  Tick = cb.Ticks;
  for n = 1:length(Tick)
    TickLab{n} = sprintf('10^{%4.1f}',Tick(n)); 
  end
  cb.TickLabels = TickLab;
  title(cb,'g m^{-3}');
  cv = caxis;
  
  subplot(212)
  f1 = pcolor(data.time,asl,log10(LWC));
  shading 'flat';
  ylim(ylimits);
  xlim([startTime,endTime]);
  ylabel('Altitude (km)');
  title('HCR LWC');
  caxis(cv);
  cb = colorbar;
  Tick = cb.Ticks;
  for n = 1:length(Tick)
    TickLab2{n} = sprintf('10^{%4.1f}',Tick(n)); 
  end
  cb.TickLabels = TickLab2;
  title(cb,'g m^{-3}');
  
%  LWP
  figure;
  set(gcf,'Position',[750 45 1000 400]);
  plot(data.time,LWP_merge,'b-',data.time,LWP_hcr,'r--');
  ylabel('LWP (g m^{-2})');
  legend('Merge','SRT')
end