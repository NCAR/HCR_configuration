% Calculate liquid water content from HCR ocean return

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='cset'; %socrates, aristo, cset, otrec
quality='qc2'; %field, qc1, or qc2
dataFreq='10hz';

ylimUpper=5;
ylimRefl=15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

%figdir=['/scr/sci/romatsch/liquidWaterHCR/'];
figdir=['/home/romatsch/plots/HCR/liquidWater/',project,'/correctOffsets/'];

corrdir=['/home/romatsch/plots/HCR/liquidWater/',project,'/offsets/'];

%dataDir=HCRdir(project,quality,dataFreq);
dataDir=['/run/media/romatsch/RSF0006/rsf/gasAtt/',project,'/10hz/'];

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

sig0ClearModAll=[];
sig0ClearAll=[];

load([corrdir,project,'_corrCoeff.mat']);

for aa=4:size(caseList,1)
    disp(['Flight ',num2str(aa)]);
    disp('Loading HCR data ...')
        
    clearvars -except aa caseList dataDir dataFreq figdir infile corrections ...
        project quality ylimRefl ylimUpper sig0ClearModAll sig0ClearAll
    
    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));
    
    %% Get data
    
    fileList=makeFileList(dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    data=[];
    
    data.DBZ = [];
    data.U_SURF=[];
    data.V_SURF=[];
    data.SST=[];
    data.TOPO=[];
    data.FLAG=[];
    data.ANTFLAG=[];
    data.rotation=[];
    data.MELTING_LAYER=[];
    data.pulse_width=[];
    data.PATH_INTEGRATED_GASEOUS_ATTENUATION_2WAY=[];
    
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
    
    upInds=find(data.ANTFLAG~=0);
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
   
    %% Calculate sigma0 from model and from reflectivity
    
    disp('Calculating sig0 ...');
    
    % Find ocean surface gate
    [linInd maxGate rangeToSurf] = hcrSurfInds(data);
    
    % Measured sig0 from surface reflectivity
    data.surfRefl=data.DBZ(linInd);
    sig0measured=calc_sig0_surfRefl(data);
    
%     sig0measAtt=sig0measured(linInd)+data.PATH_INTEGRATED_GASEOUS_ATTENUATION_2WAY;
%     sig0measAtt(data.elevation>-85)=nan;
    
    sig0measLin=10.^(sig0measured./10);
    sig0meas3gates=nan(size(data.time));
    for kk=1:length(data.time)
        if ~isnan(maxGate(kk))
            sig0meas3gates(kk)=sum(sig0measLin(maxGate(kk)-1:maxGate(kk)+1,kk),'omitnan');
        end
    end
    
    clear sig0measLin sig0measured
    
    sig0measAtt3gates=10.*log10(sig0meas3gates)+data.PATH_INTEGRATED_GASEOUS_ATTENUATION_2WAY;  
    sig0measAtt3gates(data.elevation>0)=nan;
    sig0measAtt3gates(imag(sig0measAtt3gates)~=0)=nan;
    
    % sig0 from models
    sig0modelAll= calc_sig0_model(data);
    %sig0modelFV=sig0modelAll(2,:);
    %sig0modelWu=sig0modelAll(5,:);
    sig0modelCM=sig0modelAll(8,:);
    
    clear sig0modelAll
    
    %% Correct sig0
    
    sig0measCorrLow=nan(size(sig0measAtt3gates));
    lowInds=find(data.altitude./1000<corrections.maxAlt);
    sig0measCorrLow(lowInds)=sig0measAtt3gates(lowInds);
    fittedSig=polyval(corrections.fitA,data.altitude./1000);
    sig0measCorrLow=sig0measCorrLow+(corrections.maxSig0-fittedSig);
    sig0measCorrLow(isnan(sig0measCorrLow))=sig0measAtt3gates(isnan(sig0measCorrLow));
    
    % Correct for total offset
    sig0measCorr=sig0measCorrLow-corrections.bias;
    
    %% Create ocean surface mask
    % 0 extinct or not usable
    % 1 cloud
    % 2 clear air
        
    [surfFlag1 atmFrac]=makeSurfFlag(data,linInd);
    
    %% Create field with reference sig0
    % RefFlag
    % 1 clear air
    % 2 interpolated
    % 3 model
    
    [refSig0,surfFlag,refFlag]=makeRefSig0(sig0measCorr,sig0modelCM,surfFlag1);
    
    % Find surfFlag values that have previously been clear air but are now
    % not
    surfFlag(surfFlag1~=0 & surfFlag==0 & any(data.FLAG==1,1))=1;
    
    %% 2 way path integrated attenuation from hydrometeors
    
    piaHydromet2=refSig0-sig0measCorr;
    piaHydromet2(surfFlag~=1)=nan;
  
    %% Plot preparation
   
    disp('Plotting ...');
    
    %categories = {'Cloud';'Drizzle';'Mixed';'Rain';'Unreasonable'};
    
    sig0modelCM(upInds)=nan;
     
    sig0measCorr(upInds)=nan;
        
    sig0measClear=nan(size(data.time));
    sig0measClear(surfFlag==2)=sig0measCorr(surfFlag==2);
    sig0measCloud=nan(size(data.time));
    sig0measCloud(surfFlag==1)=sig0measCorr(surfFlag==1);
    
    refSig0(upInds)=nan;
    data.PATH_INTEGRATED_GASEOUS_ATTENUATION_2WAY(upInds)=nan;
    
    sig0refMeas=nan(size(data.time));
    sig0refMeas(refFlag==1)=refSig0(refFlag==1);
    sig0refInt=nan(size(data.time));
    sig0refInt(refFlag==2)=refSig0(refFlag==2);
    sig0refMod=nan(size(data.time));
    sig0refMod(refFlag==3)=refSig0(refFlag==3);
    
    %% Plot scatter of sig0clear and altitude
    
    sig0ClearFlight=cat(2,sig0measClear(~isnan(sig0measClear))',(data.altitude(~isnan(sig0measClear))./1000)');
    sig0ClearAll=cat(1,sig0ClearAll,sig0ClearFlight);
    sig0ClearMod=cat(2,sig0measClear(~isnan(sig0measClear))',(sig0modelCM(~isnan(sig0measClear)))');
    sig0ClearModAll=cat(1,sig0ClearModAll,sig0ClearMod);   
    
    if ~isempty(sig0ClearFlight) & ~isempty(sig0ClearMod)
        f1 = figure('Position',[200 500 700 700],'DefaultAxesFontSize',12);
        scatter(sig0ClearFlight(:,1)-sig0ClearMod(:,2),sig0ClearFlight(:,2),'filled');
        xlim([-10,10])
        ylim([0,15])
        grid on
        xlabel('sig0 (clear air - model) (db)');
        ylabel('Altitude (km)');
        title(['Flight ',num2str(aa),': sig0 (clear air - model) vs altitude'])
        set(gcf,'PaperPositionMode','auto')
        print(f1,[figdir,project,'_Flight',num2str(aa),'_sig0vsAlt'],'-dpng','-r0')
    end
    
    %% Plot scatter of sig0clear, and model
        
    if ~isempty(sig0ClearMod)
        f1 = figure('Position',[200 500 700 700],'DefaultAxesFontSize',12);
        hold on
        scatter(sig0ClearMod(:,2),sig0ClearMod(:,1),'filled');
        ylim([0,25])
        xlim([0,25])
        plot([0,25],[0,25],'-k','linewidth',2);
        legend(['mean=',num2str(mean(sig0ClearMod(:,1)-sig0ClearMod(:,2),'omitnan')),...
            ' dB , std=',num2str(std(sig0ClearMod(:,1)-sig0ClearMod(:,2),'omitnan')),' dB'],...
            'location','northwest');
        
        grid on
        ylabel('sig0 clear air measured (db)');
        xlabel('sig0 model (db)');
        title(['Flight ',num2str(aa),': sig0 clear air vs sig0 model'])
        set(gcf,'PaperPositionMode','auto')
        print(f1,[figdir,project,'_Flight',num2str(aa),'_sig0vsModel'],'-dpng','-r0')
    end
    
    %% Line plots
    
    startPlot=startTime;
    
    while startPlot<endTime
        
        close all
        
        endPlot=startPlot+minutes(20);
        timeInds=find(data.time>=startPlot & data.time<=endPlot);
        
        if length(timeInds)<3000
            startPlot=endPlot;
            continue
        end
        
        if max(max(~isnan(data.dbzMasked(:,timeInds))))==0
            startPlot=endPlot;
            continue
        end
        
        %% Plot debug
        
        f1 = figure('Position',[200 500 1500 900],'DefaultAxesFontSize',12);
        
        s1=subplot(4,1,1);
        hold on
        l0=plot(data.time(:,timeInds),sig0modelCM(:,timeInds),'-c','linewidth',2);
        l1=plot(data.time(:,timeInds),sig0measClear(:,timeInds),'-b','linewidth',1);
        l2=plot(data.time(:,timeInds),sig0measCloud(:,timeInds),'color',[0.5 0.5 0.5],'linewidth',0.5);
        l4=plot(data.time(:,timeInds),sig0refMeas(:,timeInds),'-r','linewidth',2);
        l5=plot(data.time(:,timeInds),sig0refInt(:,timeInds),'-','color',[0.5 0 1],'linewidth',2);
        l6=plot(data.time(:,timeInds),sig0refMod(:,timeInds),'-m','linewidth',2);
        ylabel('Sig0 (dB)');
        ylim([5 25]);
        
        yyaxis right
        l7=plot(data.time(:,timeInds),data.PATH_INTEGRATED_GASEOUS_ATTENUATION_2WAY(:,timeInds),'-k','linewidth',1);
        ylabel('Atten. (dB)');
        ylim([-5 15]);
        grid on
        set(gca,'YColor','k');
        
        xlim([data.time(timeInds(1)),data.time(timeInds(end))]);
        
        leg=legend([l0 l1 l2 l4 l5 l6 l7],{'sig0 model','sig0 meas clear','sig0 meas cloud',...
            'sig0 ref meas','sig0 ref int','sig0 ref mod','2-way gas att'},...
            'orientation','horizontal','location','north');
        leg.ItemTokenSize=[20,18];
        title([datestr(data.time(timeInds(1)),'yyyy-mm-dd HH:MM:SS'),' to ',datestr(data.time(timeInds(end)),'yyyy-mm-dd HH:MM:SS')])
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
        l1=plot([data.time(timeInds(1)),data.time(timeInds(end))],[3,3],'-k','linewidth',2);
        ylim([0 15]);
        xlim([data.time(timeInds(1)),data.time(timeInds(end))]);
        ylabel('Num cloud pix');
        
        title('Number of cloud pixels')
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
        
        title('Elevation and rotation angle')
        s4pos=s4.Position;
        s4.Position=[s4pos(1),s4pos(2),s1pos(3),s4pos(4)];
        
        set(gcf,'PaperPositionMode','auto')
        print(f1,[figdir,project,'_lines_',datestr(data.time(timeInds(1)),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(timeInds(end)),'yyyymmdd_HHMMSS')],'-dpng','-r0')
        
        startPlot=endPlot;
    end
end
%% Save

saveAll=cat(2,sig0ClearAll(:,1),sig0ClearModAll(:,2),sig0ClearAll(:,2));
save([figdir,project,'_sig0data.mat'],'saveAll');

%% Plot

edges={-30:0.1:30 -30:0.1:30};

N=hist3(cat(2,saveAll(:,2),saveAll(:,1)),'Edges',edges);

N2=hist3(cat(2,saveAll(:,3),saveAll(:,1)-saveAll(:,2)),'Edges',edges);

% Regression
xFit = -30:0.1:30;

% Linear fit
fitOrth=gmregress(saveAll(:,2),saveAll(:,1),1);
fitAll=[fitOrth(2) fitOrth(1)];
yFit = polyval(fitAll, xFit);

wi=11;
hi=5;

fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[3,100,wi,hi]);
fig1.PaperPositionMode = 'manual';
fig1.PaperUnits = 'inches';
fig1.Units = 'inches';
fig1.PaperPosition = [0, 0, wi, hi];
fig1.PaperSize = [wi, hi];
fig1.Resize = 'off';
fig1.InvertHardcopy = 'off';

set(fig1,'color','w');

colormap jet

s1=subplot(1,2,1);

hold on
surf(edges{1},edges{2},log10(N'),'edgecolor','none')
view(2)

%axis equal
ylim([0,25])
xlim([0,25])
caxis([0 3])

l1=plot([0,25],[0,25],'-g','linewidth',1.5);
l2=plot(xFit, yFit,'-k','linewidth',2);
s1.SortMethod='childorder';

text(2,22,['y=',num2str(fitAll(1)),'x+',num2str(fitAll(2))],'fontsize',12);
text(2,24.2,['meas-mod:'],'fontsize',12);
text(2,23.1,['mean=',num2str(mean(saveAll(:,1)-saveAll(:,2),'omitnan')),' dB, std=',num2str(std(saveAll(:,1)-saveAll(:,2),'omitnan')),' dB'],'fontsize',12);


grid on
ylabel('Sig0 measured (dB)');
xlabel('Sig0 model (dB)');
title(['Sig0 measured vs model'])

s2=subplot(1,2,2);

hold on
surf(edges{1},edges{2},log10(N2'),'edgecolor','none')
view(2)

%axis equal
ylim([-10,10])
xlim([0,15])
caxis([0 3])

grid on
ylabel('Sig0 (db)');
xlabel('Altitude (km)');
title(['Sig0 (meas-mod) vs altitude'])

hcb=colorbar;
hcb.Title.String='log_{10}(N)';

s1.Position=[0.1300    0.1100    0.3347    0.8150];
s2.Position=[0.5703    0.1100    0.3347    0.8150];
hcb.Position=[0.93    0.2204    0.0202    0.5939];

set(gcf,'PaperPositionMode','auto')
print([figdir,project,'_sig0vsModel'],'-dpng','-r0');