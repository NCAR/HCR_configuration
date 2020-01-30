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

addpath('/h/eol/romatsch/gitPriv/process_HCR/oceanScans/functions/');
addpath(genpath('/h/eol/romatsch/gitPriv/utils/'));

directories.figdir=['/scr/sci/romatsch/liquidWaterHCR/',project,'/'];

directories.dataDir=HCRdir(project,quality,dataFreq);

startTime=datetime(2019,8,7,13,47,20);
endTime=datetime(2019,8,7,13,53,0);

%% Get data
[data frq]=f_load_sort_data_nadir(directories.dataDir,startTime,endTime);

% get mask
fileList=makeFileList(directories.dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

data=[];

data.DBZ = [];
data.DBMVC=[];
data.DBMHX=[];
%data.PRESS=[];
%data.TEMP=[];
%data.RH=[];
data.SST=[];
data.TOPO=[];
data.U_SURF=[];
data.V_SURF=[];
data.FLAG=[];
data.pulse_width=[];

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

data.freq=ncread(fileList{1},'frequency');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sort out bad data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reflTemp=data.DBZ;
data.reflMask=ones(size(data.time));
data.elevation=abs(data.elevation+90);

%sort out upward pointing
outElevInd=find(data.elevation>90);
data.reflMask(outElevInd)=0;

% sort out data from below 2500m altitude
altInd=find(data.altitude<2500);
data.reflMask(altInd)=0;

% Remove bang
reflTemp(data.FLAG==6)=nan;

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

%% DBM range corrected
dbmvcLin=10.^(data.DBMVC./10);
dbmvcLinRange=dbmvcLin.*data.range.^2;
DBMVC=10.*log10(dbmvcLinRange);

dbmhxLin=10.^(data.DBMHX./10);
dbmhxLinRange=dbmhxLin.*data.range.^2;
DBMHX=10.*log10(dbmhxLinRange);

if ~max(data.reflMask)==0
    %% Attenuation
    %[attenuation.liebe,attenuation.liebeMat,attenuation.itu,attenuation.ituMat,attenuation.windspeed]= ...
    %    get_atten(frq/1e+9,data.time,1,data);
    
    %% Bias
    %         data.elev=data.elevation;
    %         data.pulseWidth=data.pulse_width;
    %
    %         data=calc_sig0_atten(data,frq);
    %
    %         data.sig0measured=data.sig0measured';
    %         data=f_sigma0_model(data,{windSpd},frq,data.SST,salinity);
    
    %% Remove data with clouds
    
    dbzMeasGood=data.refl;
    dbzMeasGood(data.reflMask==0)=nan;
    dbzMeasBad=data.refl;
    dbzMeasBad(data.reflMask==1)=nan;
    
    %% Clear air ocean refl
    
    meanRefl=movmedian(dbzMeasGood,1200,'omitnan');
    
    %% Calculate liquid attenuation
    
    dbzMasked=data.DBZ;
    dbzMasked(data.FLAG>1)=nan;
    
    cloudInds=find(any(~isnan(dbzMasked),1));
    
    attLiq=nan(size(meanRefl));
    clearReflUsed=nan(size(meanRefl));
    
    meanReflRay=meanRefl(cloudInds(1)-1);
    
    specAtt=nan(size(dbzMasked));
    C1=4/(20*log10(exp(1)));
    
    b=nan(size(dbzMasked));
    b(dbzMasked<=-17)=b_drizz;
    b(dbzMasked>-17)=b_rain;
    
    meanB=mean(b,1,'omitnan');
    
    for ii=1:length(cloudInds)
        dbzRay=dbzMasked(:,cloudInds(ii));
        cloudIndsRay=find(~isnan(dbzRay));
        
        if length(cloudIndsRay)>2
            
            % Total attenuation
            if ~isnan(meanRefl(cloudInds(ii)))
                meanReflRay=meanRefl(cloudInds(ii));
            end
            
            clearReflUsed(cloudInds(ii))=meanReflRay;
            
            attDiff=meanReflRay-data.refl(cloudInds(ii));
            if attDiff>=0
                attLiq(cloudInds(ii))=attDiff;
            else
                attLiq(cloudInds(ii))=0;
            end
            
            % Specific attenuation
            % Z phi method
            dbzLin  = (10.^(0.1.*dbzRay)).^meanB(cloudInds(ii));
            I0 = C1*meanB(cloudInds(ii))*trapz(data.range(cloudIndsRay,cloudInds(ii))./1000,dbzLin(cloudIndsRay));
            CC = 10.^(0.1*meanB(cloudInds(ii))*attLiq(cloudInds(ii)))-1;
            for mm = 1:length(cloudIndsRay)
                if mm < length(cloudIndsRay)
                    Ir = C1*meanB(cloudInds(ii))*trapz(data.range(cloudIndsRay(mm:end),cloudInds(ii))./1000,dbzLin(cloudIndsRay(mm:end)));
                else
                    Ir = 0;
                end
                specAtt(cloudIndsRay(mm),cloudInds(ii)) = (dbzLin(cloudIndsRay(mm))*CC)/(I0+CC*Ir);
            end
        end
    end
    
    LWC=specAtt*alpha;
    
    %% Plot lines
    close all
    f1 = figure('Position',[200 500 2000 1200],'DefaultAxesFontSize',12,'renderer','painters');
    
    subplot(3,1,1)
    hold on
    l1=plot(data.time,dbzMeasGood,'-b','linewidth',1);
    l2=plot(data.time,meanRefl,'-r','linewidth',2);
    l3=plot(data.time,dbzMeasBad,'color',[0.5 0.5 0.5],'linewidth',0.5);
    ylabel('Refl. (dBZ)');
    ylim([20 50]);
    grid on
    
    yyaxis right
    l4=plot(data.time,DBMVC(maxGateLin)-DBMHX(maxGateLin),'-g','linewidth',1);
        
    xlim([data.time(1),data.time(end)]);
    
    legend([l1 l2 l4],{'Measured','Mean Measured','DBMVC-DBMHX range corr'},'location','southwest');
    
    title(['Reflectivity: ',datestr(data.time(1)),' to ',datestr(data.time(end))])
    cb=colorbar;
    cb.Visible='off';
    
    subplot(3,1,2)
    
    hold on
    l1=plot(data.time,reflNoOceanLin,'-k','linewidth',1);
    plot([data.time(1),data.time(end)],[0.8 0.8],...
        '-r','linewidth',1);
    ylabel('Lin refl');
    ylim([0 100]);
    xlim([data.time(1),data.time(end)]);
    %yticks(0:3:15);
    ax = gca;
    ax.YColor = 'k';
    grid on
    set(gca, 'YScale', 'log')
    
    legend(l1,{'ReflAboveOcean'},'location','northeast');
    cb=colorbar;
    cb.Visible='off';
    
    subplot(3,1,3)
    yyaxis right
    hold on
    plot(data.time,data.altitude./1000,'-b','linewidth',2);
    plot(data.time,(data.range(maxGateLin)-data.altitude)./1000,'-c','linewidth',2);
    plot(data.time,data.elevation,'-k','linewidth',2);
    ylabel('Alt (km), Elev (deg)');
    if strcmp(project,'cset') | strcmp(project,'otrec')
        ylim([-0.2 16.5]);
        yticks(0:1.5:15);
    else
        ylim([-0.2 10]);
    end
    ax = gca;
    ax.YColor = 'k';
    grid on
    
    yyaxis left
    plot(data.time,windSpd,'-g','linewidth',2);
    plot(data.time,data.SST,'-r','linewidth',2);
    ylabel('Wdspd (m s^{-1}), SST (C)');
    if strcmp(project,'cset') | strcmp(project,'otrec')
        ylim([-0.4 33]);
        yticks(0:3:33);
    else
        ylim([-0.6 30]);
        yticks(0:3:30);
    end
    xlim([data.time(1),data.time(end)]);
    ax = gca;
    ax.YColor = 'k';
    
    legend({'WindSpd','SST','Altitude','SfcRange-Alt','Elevation'},'location','northeast');
    cb=colorbar;
    cb.Visible='off';
    
    set(gcf,'PaperPositionMode','auto')
    print(f1,[directories.figdir,project,'_lines_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
    
    %% Liquid attenuation
    
    f2 = figure('Position',[200 500 2000 1200],'DefaultAxesFontSize',12,'renderer','painters');
    
    subplot(3,1,1)
    hold on
    l1=plot(data.time,clearReflUsed,'-k','linewidth',2);
    l2=plot(data.time,data.refl,'-r','linewidth',2);
    l3=plot(data.time,attLiq,'-b','linewidth',2);
    ylabel('Refl. (dBZ), Atten. (dB)');
    ylim([0 50]);
    grid on
    
    xlim([data.time(1),data.time(end)]);
    
    legend([l1 l2 l3],{'Refl. clear','Refl. measured','Liquid Attenuation'},'location','east');
    
    title(['Attenuation: ',datestr(data.time(1)),' to ',datestr(data.time(end))])
    cb=colorbar;
    cb.Visible='off';
    
    subplot(3,1,2)
    
    colormap jet
    
    hold on
    surf(data.time,data.asl./1000,dbzMasked,'edgecolor','none');
    view(2);
    ylabel('Reflectivity (dBZ)');
    caxis([-25 25]);
    ylim([0 5]);
    xlim([data.time(1),data.time(end)]);
    colorbar
    grid on
    title('Reflectivity')
    
    subplot(3,1,3)
    
    colormap jet
    
    hold on
    surf(data.time,data.asl./1000,LWC,'edgecolor','none');
    view(2);
    ylabel('LWC (g m^{-3})');
    caxis([0 4]);
    ylim([0 5]);
    xlim([data.time(1),data.time(end)]);
    colorbar
    grid on
    title('Liquid water content')
    
    set(gcf,'PaperPositionMode','auto')
    print(f2,[directories.figdir,project,'_lwc_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
    
end