% Calculate liquid water content from HCR ocean return

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='cset'; %socrates, aristo, cset, otrec
quality='qc2'; %field, qc1, or qc2
dataFreq='10hz';

b_drizz = 0.52; % Z<-17 dBZ
b_rain = 0.68; % Z>-17 dBZ
%alpha = 0.21;

ylimUpper=5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

%figdir=['/scr/sci/romatsch/liquidWaterHCR/'];
figdir=['/home/romatsch/plots/HCR/liquidWater/',project,'/cases/'];

%dataDir=HCRdir(project,quality,dataFreq);
dataDir=['/run/media/romatsch/RSF0006/rsf/meltingLayer/',project,'/10hz/'];

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/liquidWater_',project,'.txt'];

caseList=readtable(casefile);
caseStart=datetime(caseList.Var1,caseList.Var2,caseList.Var3, ...
    caseList.Var4,caseList.Var5,0);
caseEnd=datetime(caseList.Var6,caseList.Var7,caseList.Var8, ...
    caseList.Var9,caseList.Var10,0);

load('fit_RES_LWC_nofilt.mat')

for aa=1:length(caseStart)
    
    disp(['Case ',num2str(aa),' of ',num2str(length(caseStart))]);
    
    startTime=caseStart(aa);
    endTime=caseEnd(aa);
    
    %% Get data
    
    fileList=makeFileList(dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    data=[];
    
    data.DBZ = [];
    data.U_SURF=[];
    data.V_SURF=[];
    data.SST=[];
    data.TEMP=[];
    data.PRESS=[];
    data.RH=[];
    data.TOPO=[];
    data.FLAG=[];
    data.LDR=[];
    data.WIDTH=[];
    data.VEL_CORR=[];
    data.pitch=[];
    data.MELTING_LAYER=[];
    data.ICING_LEVEL=[];
    data.pulse_width=[];
    
    dataVars=fieldnames(data);
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
    
    frq=ncread(fileList{1},'frequency');
    
    % Check if all variables were found
    for ii=1:length(dataVars)
        if ~isfield(data,dataVars{ii})
            dataVars{ii}=[];
        end
    end
    
    dataVars=dataVars(~cellfun('isempty',dataVars));
    
    data.frq=ncread(fileList{1},'frequency');
    
    data.dbzMasked=data.DBZ;
    data.dbzMasked(data.FLAG>1)=nan;
    
    %% One way and two way gaseous attenuation
    
    [gasAttClear,gasAttCloud,gasAttClearMat,gasAttCloudMat]=get_gas_atten(frq/1e+9,data);
    gasAttCloud2=2*gasAttCloud';
    
    %% Calculate sigma0 from model and from reflectivity
    
    % Find ocean surface gate
    [linInd maxGate rangeToSurf] = hcrSurfInds(data);
    
    % Measured sig0 from surface reflectivity
    data.surfRefl=data.DBZ(linInd);
    sig0measured=calc_sig0_surfRefl(data);
    
    sig0measAtt=sig0measured+gasAttCloud2;
    
    % sig0 from models
    sig0modelAll= calc_sig0_model(data);
    %sig0modelFV=sig0modelAll(2,:);
    %sig0modelWu=sig0modelAll(5,:);
    sig0modelCM=sig0modelAll(8,:);
    
    %% Create ocean surface mask
    % 0 extinct or not usable
    % 1 cloud
    % 2 clear air
    
    surfFlag=makeSurfFlag(data,maxGate);
    
    %% Create field with reference sig0
    refSig0=makeRefSig0(sig0measAtt,sig0modelCM,surfFlag);
    
    %% 2 way path integrated attenuation from hydrometeors
    
    piaHydromet2=refSig0-sig0measAtt;
    piaHydromet2(surfFlag~=1)=nan;
    
    %% Separate warm and cold precip
    data.MELTING_LAYER(data.MELTING_LAYER<20)=10;
    data.MELTING_LAYER(data.MELTING_LAYER>19)=20;
    
    warmRefl=data.dbzMasked;
    warmRefl(data.MELTING_LAYER==20)=nan;
    coldRefl=data.dbzMasked;
    coldRefl(data.MELTING_LAYER==10)=nan;
    
    %% Make flag field that flags type of hydrometeor attenuation
    % 0 no attenuation
    % 1 liquid only
    % 2 mixed
    % 3 ice only
    
    warmFlag=any(~isnan(warmRefl),1);
    coldFlag=any(~isnan(coldRefl),1);
    
    attFlag=nan(size(data.time));
    attFlag(warmFlag & ~coldFlag)=1;
    attFlag(warmFlag & coldFlag)=2;
    attFlag(~warmFlag & coldFlag)=3;
    attFlag(surfFlag==0 | surfFlag==2)=0;
    
    %% Calculate two way ice attenuation
    
    % This equation comes from DOI: 10.1175/JTECH-D-18-0154.1 but I don't think
    % it applies to our very low reflectivities. Also, the units seem weird.
    %
    % coldReflLin=10.^(coldRefl./10);
    % iceSpecAtt=0.0325.*coldReflLin;
    %
    % iceAttAll=iceSpecAtt.*(data.range(2)-data.range(1))./1000;
    % piaIce2=sum(iceAttAll,1,'omitnan');
    
    %% Calculate liquid attenuation
    
    piaLiq2=piaHydromet2;%-piaIce2;
    
    %% Calculate specific liquid attenuation with zPhi method
    specAttLiq=nan(size(data.DBZ));
    
    C1=4/(20*log10(exp(1)));
    
    b=nan(size(data.DBZ));
    b(warmRefl<=-17)=b_drizz;
    b(warmRefl>-17)=b_rain;
    
    meanB=mode(b,1);
    
    cloudInds=find(attFlag==1);
    
    for ii=1:length(cloudInds)
        dbzRay=warmRefl(:,cloudInds(ii));
        cloudIndsRay=find(~isnan(dbzRay));
        
        if length(cloudIndsRay)>2
            % Two way specific attenuation
            % Z phi method
            dbzLinB  = (10.^(0.1.*dbzRay)).^meanB(cloudInds(ii));
            I0 = C1*meanB(cloudInds(ii))*trapz(data.range(cloudIndsRay,cloudInds(ii))./1000,dbzLinB(cloudIndsRay));
            CC = 10.^(0.1*meanB(cloudInds(ii))*piaLiq2(cloudInds(ii)))-1;
            for mm = 1:length(cloudIndsRay)
                if mm < length(cloudIndsRay)
                    Ir = C1*meanB(cloudInds(ii))*trapz(data.range(cloudIndsRay(mm:end),cloudInds(ii))./1000,dbzLinB(cloudIndsRay(mm:end)));
                else
                    Ir = 0;
                end
                specAttLiq(cloudIndsRay(mm),cloudInds(ii)) = (dbzLinB(cloudIndsRay(mm))*CC)/(I0+CC*Ir);
            end
        end
    end
    
    %% Calculate LWC and RES
    % Method taking Mie scattering into account
    piaInd=2*specAttLiq*(data.range(2)-data.range(1))./1000;
    piaGate=cumsum(piaInd,1,'omitnan');
    piaGate(isnan(specAttLiq))=nan;
      
    % Attenuation corrected reflectivity
    dbzCorr=data.dbzMasked+piaGate;
    
    % Calculate Xres
    zw=10.^(0.1*dbzCorr);
    Xres = (zw./specAttLiq).^(1/3);
    
    % Interpolate fitting data
    zw_bnd_c  = interp1(Awbnd_c,zwbnd_c,specAttLiq,'linear','extrap');
    zw_bnd_c(zw_bnd_c<0)=0;
    zw_bnd_dc = interp1(Awbnd_dc,zwbnd_dc,specAttLiq,'linear','extrap');
    zw_bnd_dc(zw_bnd_dc<0)=0;
    
    
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
%% Method without Mie scattering
    alpha=1./(4.792-3.63e-2*data.TEMP-1.897e-4*data.TEMP.^2);
    LWCorig=specAttLiq.*alpha;
    
    %% Plot
    close all
    
    %timeMat=repmat(data.time,size(data.TEMP,1),1);
    
    sig0measClear=nan(size(data.time));
    sig0measClear(surfFlag==2)=sig0measAtt(surfFlag==2);
    sig0measCloud=nan(size(data.time));
    sig0measCloud(surfFlag==1)=sig0measAtt(surfFlag==1);
    
    f1 = figure('Position',[200 500 1500 900],'DefaultAxesFontSize',12);
    
    s1=subplot(3,1,1);
    hold on
    l0=plot(data.time,sig0modelCM,'-c','linewidth',2);
    l1=plot(data.time,sig0measClear,'-b','linewidth',1);
    l2=plot(data.time,sig0measCloud,'color',[0.5 0.5 0.5],'linewidth',0.5);
    l3=plot(data.time,refSig0,'-r','linewidth',2);
    ylabel('Sig0 (dB)');
    ylim([0 20]);
    
    yyaxis right
    l4=plot(data.time,gasAttCloud2,'-k','linewidth',1);
    l5=plot(data.time,piaLiq2,'-g','linewidth',1);
    %l6=plot(data.time,piaIce2*10,'-m','linewidth',1);
    %l6=plot(data.time,piaHydromet2,'-y','linewidth',1);
    ylabel('Atten. (dB)');
    ylim([-5 15]);
    grid on
    set(gca,'YColor','k');
    
    xlim([data.time(1),data.time(end)]);
    
    legend([l0 l1 l3 l4 l5],{'sig0 model','sig0 measured','sig0 clear','2-way gaseous atten.','2-way PIA liq.'},...
        'orientation','horizontal','location','north');
    title([datestr(data.time(1)),' to ',datestr(data.time(end))])
    s1pos=s1.Position;
    
    s2=subplot(3,1,2);
    
    colormap jet
    
    hold on
    surf(data.time,data.asl./1000,data.dbzMasked,'edgecolor','none');
    view(2);
    %plot(data.time,data.ICING_LEVEL./1000,'-k','linewidth',2);
    ylabel('Altitude (km)');
    caxis([-25 25]);
    ylim([0 ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    colorbar
    grid on
    title('Reflectivity (dBZ)')
    s2pos=s2.Position;
    s2.Position=[s2pos(1),s2pos(2),s1pos(3),s2pos(4)];
    
    s3=subplot(3,1,3);
    
    colmap=jet;
    colmap=cat(1,[1 0 1],colmap);
    
    hold on
    surf(data.time,data.asl./1000,LWCorig,'edgecolor','none');
    view(2);
    colormap(s3,colmap)
    ylabel('Altitude (km)');
    caxis([0 2]);
    ylim([0 ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    colorbar
    grid on
    title('Liquid water content (g m^{-3})')
    s3pos=s3.Position;
    s3.Position=[s3pos(1),s3pos(2),s1pos(3),s3pos(4)];
    
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_lwc_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
    
end