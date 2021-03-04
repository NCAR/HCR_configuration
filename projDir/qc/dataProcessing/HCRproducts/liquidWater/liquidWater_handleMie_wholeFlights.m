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
if strcmp(project,'otrec')
    ylimRefl=15;
else
    ylimRefl=10;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

%figdir=['/scr/sci/romatsch/liquidWaterHCR/'];
figdir=['/home/romatsch/plots/HCR/liquidWater/',project,'/flights/'];

%dataDir=HCRdir(project,quality,dataFreq);
dataDir=['/run/media/romatsch/RSF0006/rsf/meltingLayer/',project,'/10hz/'];

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

%load('fit_RES_LWC_nofilt.mat')

Awbnd_c=[0.00161,0.00181,0.00582,0.00781,0.0052,0.0107,0.00985,...
    0.00723,0.00944,0.0206,0.0262,0.0364,0.0534,0.107,0.158,...
    0.254,0.405,0.655,0.893,1.62];

Awbnd_dc=[0.000633,0.00132,0.00235,0.00334,0.00486,0.00577,...
    0.00962,0.0127,0.0175,0.0284,0.0346,0.0598,0.074,0.113,...
    0.154,0.19,0.268,0.423,0.925,1.13,1.4,1.66,1.76,1.89,...
    1.97,2.48,3.41,4.09,5.56,5.64,7.29,8.78,9.71];

zwbnd_c=[1.24451461177138e-05,1.92752491319093e-05,2.98538261891796e-05,...
    4.62381021399260e-05,7.16143410212901e-05,0.000110917481526240,...
    0.000171790838715759,0.000266072505979881,0.000412097519097330,...
    0.000638263486190548,0.000988553094656938,0.00153108746168203,...
    0.00237137370566166,0.00367282300498084,0.00568852930843841,...
    0.00881048873008014,0.0136458313658892,0.0211348903983665,...
    0.0327340694878838,0.0506990708274704];

zwbnd_dc=[0.000240021419577473,0.000502631844408201,0.000727360981363162,...
    0.00105256760608252,0.00152317569097253,0.00220419493452256,...
    0.00318970118691486,0.00461583206750714,0.00667959298596078,...
    0.00966607142668268,0.0139878188719147,0.0202418405737606,...
    0.0292920657298628,0.0423886904748601,0.0613408797025088,...
    0.0887666847106126,0.128454700237875,0.185887420117083,...
    0.389269161853171,0.815173400515363,1.17964138399968,1.70706477169755,...
    2.47030171567091,3.57478560135610,5.17308959249628,7.48600305479624,...
    10.8330313508791,22.6855858864991,47.5061679731930,99.4832580824082,...
    208.329129899000,301.474094917887,436.264625838075];

P_c=[0.126782731258678,-0.000115939385913554];
P_c0=[0.161877656181692,0.240155793871566];
P_dr0=[0.331402635889504,0.229591393662394];
P_dr1=[0.189219903249597,-0.0835385698027379];

Q_c=[0.254509910846974,0.00105380577830583];
Q_c0=[0.126682917191433,0.284748888720728];
Q_dr0=[0.0136623949301112,0.581217037948342];
Q_dr1=[0.112138117429853,0.00128536965660402];

Xres_thres=3;

zwc0_bnd=0.00293494522226979;
zwdr0_bnd=2.50149707846010;

sig0ClearAll=[];

for aa=1:size(caseList,1)
    disp(['Flight ',num2str(aa)]);
    disp('Loading HCR data ...')
    disp(['Starting at ',datestr(datetime('now'),'yyyy-mm-dd HH:MM')]);
    
    clearvars -except aa Awbnd_c Awbnd_dc b_drizz b_rain caseList dataDir dataFreq figdir infile ...
        P_c P_c0 P_dr0 P_dr1 project Q_c Q_c0 Q_dr0 Q_dr1 quality sig0ClearAll Xres_thres ylimRefl ylimUpper ...
        zwbnd_c zwbnd_dc zwc0_bnd zwdr0_bnd
    
    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));
    
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
    %data.LDR=[];
    %data.WIDTH=[];
    %data.VEL_CORR=[];
    %data.pitch=[];
    data.rotation=[];
    data.MELTING_LAYER=[];
    %data.ICING_LEVEL=[];
    data.pulse_width=[];
    
    %data.DBMVC=[];
    
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
    
    data.frq=ncread(fileList{1},'frequency');
    
    %% Remove all up pointing data
    
    upInds=find(data.elevation>-85);
    upInds=cat(2,upInds,find(any(data.FLAG>9,1)),find(any(data.FLAG==3,1)));
    
    infields=fields(data);
    for bb=1:length(infields)
        if strcmp(infields{bb},'DBZ') | strcmp(infields{bb},'FLAG') | ...
                strcmp(infields{bb},'rotation') | strcmp(infields{bb},'elevation')
            currfield=data.(infields{bb});
            currfield(:,upInds)=nan;
            data.(infields{bb})=currfield;
        end
    end
    
    data.dbzMasked=data.DBZ;
    data.dbzMasked(data.FLAG>1)=nan;
    %% One way and two way gaseous attenuation
    
    disp('Calculating gaseous attenuation ...');
    [gasAttClear,gasAttCloud,gasAttClearMat,gasAttCloudMat]=get_gas_atten(data);
    gasAttCloud2=2*gasAttCloud';
    
    data=rmfield(data,'TEMP');
    data=rmfield(data,'PRESS');
    data=rmfield(data,'RH');
    
    clear gasAttClearMat
    %% Calculate sigma0 from model and from reflectivity
    
    disp('Calculating sig0 ...');
    
    % Find ocean surface gate
    [linInd maxGate rangeToSurf] = hcrSurfInds(data);
    
    % Measured sig0 from surface reflectivity
    data.surfRefl=data.DBZ(linInd);
    sig0measured=calc_sig0_surfRefl(data);
    
    sig0measAtt=sig0measured(linInd)+gasAttCloud2;
    sig0measAtt(data.elevation>-85)=nan;
    
    sig0measLin=10.^(sig0measured./10);
    sig0meas3gates=nan(size(data.time));
    for kk=1:length(data.time)
        if ~isnan(maxGate(kk))
            sig0meas3gates(kk)=sum(sig0measLin(maxGate(kk)-1:maxGate(kk)+1,kk),'omitnan');
        end
    end
    
    clear sig0measLin
    
    sig0measAtt3gates=10.*log10(sig0meas3gates)+gasAttCloud2;  
    sig0measAtt3gates(data.elevation>0)=nan;
    sig0measAtt3gates(imag(sig0measAtt3gates)~=0)=nan;
    
    % sig0 from models
    sig0modelAll= calc_sig0_model(data);
    %sig0modelFV=sig0modelAll(2,:);
    %sig0modelWu=sig0modelAll(5,:);
    sig0modelCM=sig0modelAll(8,:);
    
    %% Create ocean surface mask
    % 0 extinct or not usable
    % 1 cloud
    % 2 clear air
        
    [surfFlag1 atmFrac]=makeSurfFlag(data,gasAttCloudMat,maxGate);
    
    clear gasAttCloudMat
    %% Create field with reference sig0
    % RefFlag
    % 1 clear air
    % 2 interpolated
    % 3 model
    
    [refSig0,surfFlag,refFlag]=makeRefSig0(sig0measAtt,sig0modelCM,surfFlag1);
    
    % Find surfFlag values that have previously been clear air but are now
    % not
    surfFlag(surfFlag1~=0 & surfFlag==0 & any(data.FLAG==1,1))=1;
    
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
    
    clear coldRefl warmRefl
    %% Calculate liquid attenuation
    
    piaLiq2=piaHydromet2;%-piaIce2;
    
    %% Calculate specific liquid attenuation with zPhi method
    
    disp('Calculating specific liquid attenuation ...');
    
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
    
    specAttLiqNeg=find(specAttLiq<0);
    specAttLiq(specAttLiqNeg)=nan;
    %% Calculate LWC and RES with method taking Mie scattering into account
    
    disp('Calculating LWC, RES, and PID ...');
    
    PID=nan(size(data.DBZ));
    LWC=nan(size(data.DBZ));
    RES=nan(size(data.DBZ));
    
    piaInd=2*specAttLiq*(data.range(2)-data.range(1))./1000;
    piaGate=cumsum(piaInd,1,'omitnan');
    piaGate(isnan(specAttLiq))=nan;
      
    % Attenuation corrected reflectivity
    dbzCorr=data.dbzMasked+piaGate;
    
    % Calculate Xres
    zLin=10.^(0.1*dbzCorr);
    Xres = (zLin./specAttLiq).^(1/3);
    
    % Interpolate fitting data
    zw_bnd_c=interp1(Awbnd_c,zwbnd_c,specAttLiq,'linear','extrap');
    zw_bnd_c(zw_bnd_c<0)=0;
    zw_bnd_dc=interp1(Awbnd_dc,zwbnd_dc,specAttLiq,'linear','extrap');
    zw_bnd_dc(zw_bnd_dc<0)=0;
    
    % Cloud
    % Unreasonable cloud
    cloudIndsU=zLin<=zw_bnd_c & specAttLiq>7.5;
    PID(cloudIndsU)=5;
    
    % Cloud indices
    cloudInds=(zLin<=zw_bnd_c | zLin<min(zwbnd_dc)) & specAttLiq<=7.5;
    
    RES(cloudInds)=P_c(1)*Xres(cloudInds)+P_c(2);
    LWC(cloudInds)=Q_c(1)*specAttLiq(cloudInds)+Q_c(2);
    PID(cloudInds)=1;
    
    %   Mixed
    % Unreasonable mixed
    mixedIndsU=zLin>zw_bnd_c & zLin<zw_bnd_dc & Xres>Xres_thres;
    PID(mixedIndsU)=5;
    
    % Mixed indices
    mixedInds=zLin>zw_bnd_c & zLin<zw_bnd_dc & Xres<=Xres_thres & specAttLiq<=10;
    
    % Compute the Ratio
    LWC1=nan(size(data.DBZ));
    LWC2=nan(size(data.DBZ));
    RES1=nan(size(data.DBZ));
    RES2=nan(size(data.DBZ));
    
    Ratio = abs(zLin-zw_bnd_dc)./(abs(zLin-zw_bnd_c)+abs(zLin-zw_bnd_dc));
    
    RES1(mixedInds)=P_c(1)*Xres(mixedInds)+P_c(2);
    RES2(mixedInds)=P_dr1(1)*Xres(mixedInds)+P_dr1(2);
    LWC1(mixedInds)=Q_c(1)*specAttLiq(mixedInds)+Q_c(2);
    LWC2(mixedInds)=Q_dr1(1)*specAttLiq(mixedInds)+Q_dr1(2);
    LWC(mixedInds)=Ratio(mixedInds).*LWC1(mixedInds)+(1-Ratio(mixedInds)).*LWC2(mixedInds);
    RES(mixedInds)=Ratio(mixedInds).*RES1(mixedInds)+(1-Ratio(mixedInds)).*RES2(mixedInds);
    PID(mixedInds)=3;
    
    % Rain below rain boundary because of large attenuation
    rainInds1=zLin>zw_bnd_c & zLin<zw_bnd_dc & specAttLiq>10;
    PID(rainInds1)=4;
    
    % Rain
    rainInds2=zLin>=zw_bnd_dc & zLin>=min(zw_bnd_dc) & Xres>Xres_thres;
    PID(rainInds2) = 4;
    
    % Drizzle (above rain boundary but small Xres)
    drizzleInds=zLin>=zw_bnd_dc & Xres<=Xres_thres;
    RES(drizzleInds)=P_dr1(1)*Xres(drizzleInds)+P_dr1(2);
    LWC(drizzleInds)=Q_dr1(1)*specAttLiq(drizzleInds)+Q_dr1(2);
    PID(drizzleInds)=2;
    
    % specAttLiq=0 and cloud
    a0cInds=specAttLiq==0 & zLin>0 & zLin<=zwc0_bnd;
    PID(a0cInds)=1;
    LWC(a0cInds)=Q_c0(1)*zLin(a0cInds).^Q_c0(2);
    RES(a0cInds)=P_c0(1)*zLin(a0cInds).^P_c0(2);
    
    % specAttLiq=0 and drizzle
    a0dInds=specAttLiq==0 & zLin>zwc0_bnd & zLin<zwdr0_bnd;
    PID(a0dInds)=2;
    LWC(a0dInds)=Q_dr0(1)*zLin(a0dInds).^Q_dr0(2);
    RES(a0dInds)=P_dr0(1)*zLin(a0dInds).^P_dr0(2);
    
    % specAttLiq=0 and mixed or rain -> unreasonable
    a0mrInds=specAttLiq==0 & zLin>zwdr0_bnd;
    PID(a0mrInds)=8;
    LWC(a0mrInds)=nan;
    RES(a0mrInds)=nan;
       
    %LWC(LWC<0) = 0;   RES(RES<0) = 0;
    
    LWC(specAttLiqNeg)=-99;
    RES(specAttLiqNeg)=-99;
    
    %% Method without Mie scattering
%     alpha=1./(4.792-3.63e-2*data.TEMP-1.897e-4*data.TEMP.^2);
%     LWCorig=specAttLiq.*alpha;
    
    %% Plot reflectivity, LWC, and RES
   
    disp('Plotting ...');
    
    categories = {'Cloud';'Drizzle';'Mixed';'Rain';'Unreasonable'};
    
    sig0modelCM(upInds)=nan;
     
    sig0measAtt(upInds)=nan;
    
    sig0measClear=nan(size(data.time));
    sig0measClear(surfFlag==2)=sig0measAtt(surfFlag==2);
    sig0measCloud=nan(size(data.time));
    sig0measCloud(surfFlag==1)=sig0measAtt(surfFlag==1);
    
    refSig0(upInds)=nan;
    gasAttCloud2(upInds)=nan;
    
    sig0refMeas=nan(size(data.time));
    sig0refMeas(refFlag==1)=refSig0(refFlag==1);
    sig0refInt=nan(size(data.time));
    sig0refInt(refFlag==2)=refSig0(refFlag==2);
    sig0refMod=nan(size(data.time));
    sig0refMod(refFlag==3)=refSig0(refFlag==3);
    
    %% Plot scatter of sig0clear and altitude
    
    sig0ClearFlight=cat(2,sig0measClear(~isnan(sig0measClear)),data.altitude(~isnan(sig0measClear))./1000);
    sig0ClearAll=cat(1,sig0ClearAll,sig0ClearFlight);
    
    disp('Plotting ...');
    
    f1 = figure('Position',[200 500 700 700],'DefaultAxesFontSize',12);
    scatter(sig0ClearFlight(:,1),sig0ClearFlight(:,2),'filled');
    xlim([5,15])
    ylim([0,15])
    xlabel('sig0 clear air (db)');
    ylabel('Altitude (km)');
    title(['Flight ',num2str(aa),': sig0 clear air vs altitude'])
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_Flight',num2str(aa),'_sig0vsAlt'],'-dpng','-r0')
    
    startPlot=startTime;
    
    while startPlot<endTime
        
        close all
        
        endPlot=startPlot+minutes(15);
        timeInds=find(data.time>=startPlot & data.time<=endPlot);
        
        timePlot=data.time(timeInds);
        dbzPlot=data.DBZ(:,timeInds);
        dbzMaskedPlot=data.dbzMasked(:,timeInds);
        aslPlot=data.asl(:,timeInds);
        pidPlot=PID(:,timeInds);
        lwcPlot=LWC(:,timeInds);
        resPlot=RES(:,timeInds);
        
        f1 = figure('Position',[200 500 1500 900],'DefaultAxesFontSize',12);
        
        colormap jet
        
        s1=subplot(4,1,1);
        surf(data.time(:,timeInds),data.asl(:,timeInds)./1000,PID(:,timeInds),'edgecolor','none');
        view(2);
        ylabel('Altitude (km)');
        caxis([0.5 5.5]);
        
        ylim([0 ylimUpper]);
        xlim([data.time(timeInds(1)),data.time(timeInds(end))]);
        s1.Colormap=[0,0,1;0,1,0;1,1,0;1,0,0;0.5,0,1];
        c=colorbar('TickLabels',categories);
        grid on
        title([{[datestr(data.time(timeInds(1)),'yyyy-mm-dd HH:MM:SS'),' to ',datestr(data.time(timeInds(end)),'yyyy-mm-dd HH:MM:SS')]};{'PID'}])
        s1pos=s1.Position;
        
        s2=subplot(4,1,2);
        
        hold on
        surf(data.time(:,timeInds),data.asl(:,timeInds)./1000,data.dbzMasked(:,timeInds),'edgecolor','none');
        view(2);
        ylabel('Altitude (km)');
        caxis([-25 25]);
        ylim([0 ylimUpper]);
        xlim([data.time(timeInds(1)),data.time(timeInds(end))]);
        colorbar
        grid on
        title('Reflectivity (dBZ)')
        s2pos=s2.Position;
        s2.Position=[s2pos(1),s2pos(2),s1pos(3),s2pos(4)];
        
        s3=subplot(4,1,3);
        
        colmap=jet;
        colmap=cat(1,[1 0 1],colmap);
        
        hold on
        surf(data.time(:,timeInds),data.asl(:,timeInds)./1000,LWC(:,timeInds),'edgecolor','none');
        view(2);
        colormap(s3,colmap)
        ylabel('Altitude (km)');
        caxis([0 1]);
        ylim([0 ylimUpper]);
        xlim([data.time(timeInds(1)),data.time(timeInds(end))]);
        colorbar
        grid on
        title('Liquid water content (g m^{-3})')
        s3pos=s3.Position;
        s3.Position=[s3pos(1),s3pos(2),s1pos(3),s3pos(4)];
        
        s4=subplot(4,1,4);
        
        colmap=jet;
        colmap=cat(1,[1 0 1],colmap);
        
        hold on
        surf(data.time(:,timeInds),data.asl(:,timeInds)./1000,RES(:,timeInds),'edgecolor','none');
        view(2);
        colormap(s4,colmap)
        ylabel('Altitude (km)');
        caxis([0 0.5]);
        ylim([0 ylimUpper]);
        xlim([data.time(timeInds(1)),data.time(timeInds(end))]);
        colorbar
        grid on
        title('Radar estimated size (mm)')
        s4pos=s4.Position;
        s4.Position=[s4pos(1),s4pos(2),s1pos(3),s4pos(4)];
        
        set(gcf,'PaperPositionMode','auto')
        print(f1,[figdir,project,'_lwc_',datestr(data.time(timeInds(1)),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(timeInds(end)),'yyyymmdd_HHMMSS')],'-dpng','-r0')
        
        %% Plot debug
        
        f1 = figure('Position',[200 500 1500 900],'DefaultAxesFontSize',12);
        
        s1=subplot(4,1,1);
        hold on
        l0=plot(data.time(:,timeInds),sig0modelCM(:,timeInds),'-c','linewidth',2);
        l1=plot(data.time(:,timeInds),sig0measClear(:,timeInds),'-b','linewidth',1);
        l2=plot(data.time(:,timeInds),sig0measCloud(:,timeInds),'color',[0.5 0.5 0.5],'linewidth',0.5);
        l3=plot(data.time(:,timeInds),sig0refMeas(:,timeInds),'-r','linewidth',2);
        l4=plot(data.time(:,timeInds),sig0refInt(:,timeInds),'-','color',[0.5 0 1],'linewidth',2);
        l5=plot(data.time(:,timeInds),sig0refMod(:,timeInds),'-m','linewidth',2);
        ylabel('Sig0 (dB)');
        ylim([0 20]);
        
        yyaxis right
        l6=plot(data.time(:,timeInds),gasAttCloud2(:,timeInds),'-k','linewidth',1);
        l7=plot(data.time(:,timeInds),piaLiq2(:,timeInds),'-g','linewidth',1);
        ylabel('Atten. (dB)');
        ylim([-5 15]);
        grid on
        set(gca,'YColor','k');
        
        xlim([data.time(timeInds(1)),data.time(timeInds(end))]);
        
        legend([l0 l1 l2 l3 l4 l5 l6 l7],{'sig0 model','sig0 meas clear','sig0 meas cloud',...
            'sig0 ref meas','sig0 ref int','sig0 ref mod','2-way gas att','2-way PIA liq'},...
            'orientation','horizontal','location','north');
        title([datestr(data.time(1),'yyyy-mm-dd HH:MM:SS'),' to ',datestr(data.time(end),'yyyy-mm-dd HH:MM:SS')])
        s1pos=s1.Position;
        
        s2=subplot(4,1,2);
        
        colormap jet
        
        hold on
        surf(data.time(:,timeInds),data.asl(:,timeInds)./1000,data.DBZ(:,timeInds),'edgecolor','none');
        view(2);
        l1=plot(data.time(:,timeInds),data.altitude(:,timeInds)./1000,'-k','linewidth',2);
        ylabel('Altitude (km)');
        caxis([-25 25]);
        ylim([-0.5 ylimRefl]);
        xlim([data.time(timeInds(1)),data.time(timeInds(end))]);
        colorbar
        grid on
        legend(l1,'Altitude');
        title('Reflectivity (dBZ)')
        s2pos=s2.Position;
        s2.Position=[s2pos(1),s2pos(2),s1pos(3),s2pos(4)];
        
        s3=subplot(4,1,3);
        
        hold on
        l0=plot(data.time(:,timeInds),atmFrac(:,timeInds),'-r','linewidth',1);
        l1=plot([data.time(timeInds(1)),data.time(timeInds(end))],[1e-7,1e-7],'-k','linewidth',2);
        ylim([0 0.0000005]);
        xlim([data.time(timeInds(1)),data.time(timeInds(end))]);
        ylabel('Fraction');
        
        title('Fraction of atmosphere over ocean reflectivity')
        s3pos=s3.Position;
        s3.Position=[s3pos(1),s3pos(2),s1pos(3),s3pos(4)];
        
        s4=subplot(4,1,4);
        
        hold on
        l0=plot(data.time(:,timeInds),data.elevation(:,timeInds),'-k','linewidth',2);
        ylabel('Elevation (deg)');
        ylim([-90 -89.6]);
        
        yyaxis right
        l1=plot(data.time(:,timeInds),data.rotation(:,timeInds),'-r','linewidth',2);
        ylabel('Rotation (deg)');
        ylim([160 200]);
        grid on
        set(gca,'YColor','k');
        
        xlim([data.time(timeInds(1)),data.time(timeInds(end))]);
        
        legend([l0 l1],{'Elevation','Rotation'},'location','northeast');
        
        title('Liquid water content (g m^{-3})')
        s4pos=s4.Position;
        s4.Position=[s4pos(1),s4pos(2),s1pos(3),s4pos(4)];
        
        set(gcf,'PaperPositionMode','auto')
        print(f1,[figdir,project,'_lines_',datestr(data.time(timeInds(1)),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(timeInds(end)),'yyyymmdd_HHMMSS')],'-dpng','-r0')
        startPlot=endPlot;
    end
end