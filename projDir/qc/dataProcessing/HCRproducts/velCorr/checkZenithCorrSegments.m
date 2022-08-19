% Plot up and down pointing segments for the zenith pointing corrections

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='socrates'; % socrates, cset, aristo, otrec
quality='qc3'; % field, qc1, qc2
qcVersion='v3.1';
freqData='10hz'; % 10hz, 100hz, or 2hz
whichModel='era5';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

indir='/scr/snow2/rsfdata/projects/socrates/hcr/qc3/cfradial/v3.0_full/10hz/';

[~,indirCorr]=modelDir(project,whichModel,quality,qcVersion,dataFreq);

figdir=[indirCorr(1:end-5),'meltLayerPlots/testMat/'];

infileZ=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/velCorrZenithSegments_',project,'.txt'];
caseListZ = table2array(readtable(infileZ));

infileN=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/velCorrNadirSegments_',project,'.txt'];
caseListN = table2array(readtable(infileN));

%% Load vel
disp('Loading velCorr ...');

fileIn1=dir([indirCorr,whichModel,'.velCorr.*.Flight13.mat']);
meltLayerIn=load([indirCorr,fileIn1.name]);
meltLayerGet=meltLayerIn.meltLayer;

% Get zenith pointing
for aa=1:size(caseListZ,1)
        
    startTime=datetime(caseListZ(aa,1:6));
    endTime=datetime(caseListZ(aa,7:12));
       
    %% Load vel
    disp('Loading melt layer.');
    
    fileIn1=dir([indirCorr,whichModel,'.velCorr.*.Flight',num2str(aa),'.mat']);
    meltLayerIn=load([indirCorr,fileIn1.name]);
    meltLayerGet=meltLayerIn.meltLayer;
        
    %% Plot in hourly increments
    
    disp('Plotting ...');
    
    startPlot=startTime;
    
    while startPlot<endTime
        
        close all
        
        endPlot=startPlot+minutes(30);
        timeInds=find(data.time>=startPlot & data.time<=endPlot);
        
        dbzPlot=data.DBZ_MASKED(:,timeInds);
        
        if sum(sum(~isnan(dbzPlot)))~=0
                        
            timeIndsMelt=find(timeMELT>=data.time(timeInds(1)) & timeMELT<=data.time(timeInds(end)));
            
            meltPlot=meltLayerGet(:,timeIndsMelt);
            icePlot=iceLevelGet(timeIndsMelt);
            
            %% Get indices
            
            elevenInds=find(meltPlot==11);
            twelveInds=find(meltPlot==12);
            thirteenInds=find(meltPlot==13);
            fourteenInds=find(meltPlot==14);
            
            twentyoneInds=find(meltPlot==21);
            twentytwoInds=find(meltPlot==22);
            twentythreeInds=find(meltPlot==23);
            twentyfourInds=find(meltPlot==24);
            
            %% Plot
            
            timeMat=repmat(data.time(:,timeInds),size(data.LDR(:,timeInds),1),1);
            dbzMasked=data.DBZ_MASKED(:,timeInds);
            ldrMasked=data.LDR(:,timeInds);
            velMasked=data.VEL_UNFOLDED(:,timeInds);
            timeMasked=data.time(timeInds);
            aslMasked=data.asl(:,timeInds);
            
            close all
            
            if etime(datevec(endPlot),datevec(startPlot))<=900
                newInds=1:1:size(timeMat,2);
            elseif etime(datevec(endPlot),datevec(startPlot))<=3600
                newInds=1:10:size(timeMat,2);
            else
                newInds=1:100:size(timeMat,2);
            end
            
            % Resample for plotting
            newDBZ=dbzMasked(:,newInds);
            newLDR=ldrMasked(:,newInds);
            newVEL=velMasked(:,newInds);
            newASL=aslMasked(:,newInds);
            newFindMelt=meltPlot(:,newInds);
            newTime=timeMasked(newInds);
            
            fig1=figure('DefaultAxesFontSize',11,'position',[100,1300,1500,1200],'visible','off');
            
            ax1=subplot(4,1,1);
            hold on;
            sub1=surf(newTime,newASL./1000,newDBZ,'edgecolor','none');
            view(2);
            sub1=colMapDBZ(sub1);
            scatter(timeMat(twentyoneInds),aslMasked(twentyoneInds)./1000,10,'k','filled');
            scatter(timeMat(elevenInds),aslMasked(elevenInds)./1000,10,...
                'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7]);
            
            scatter(timeMat(twentyfourInds),aslMasked(twentyfourInds)./1000,10,...
                'MarkerEdgeColor',[0.45 0.76 0.42],'MarkerFaceColor',[0.45 0.76 0.42]);
            scatter(timeMat(twentythreeInds),aslMasked(twentythreeInds)./1000,10,...
                'MarkerEdgeColor',[0.7 0.8 0.87],'MarkerFaceColor',[0.7 0.8 0.87]);
            scatter(timeMat(twentytwoInds),aslMasked(twentytwoInds)./1000,10,...
                'MarkerEdgeColor',[0.17 0.45 0.7],'MarkerFaceColor',[0.17 0.45 0.7]);
            
            scatter(timeMat(fourteenInds),aslMasked(fourteenInds)./1000,10,'g','filled');
            scatter(timeMat(thirteenInds),aslMasked(thirteenInds)./1000,10,'c','filled');
            scatter(timeMat(twelveInds),aslMasked(twelveInds)./1000,10,'b','filled');
            
            ax = gca;
            ax.SortMethod = 'childorder';
            ylim(ylimits);
            ylabel('Altitude (km)');
            xlim([startPlot,endPlot]);
            title(['Flight ',num2str(aa),': Reflectivity and melting layer'])
            grid on
            set(gca,'xticklabel',[])
            ax1.Position=[0.06 0.765 0.87 0.21];
            
            ax2=subplot(4,1,2);
            hold on;
            sub1=surf(newTime,newASL./1000,newFindMelt,'edgecolor','none');
            ax2.Colormap=([1 0 1;1 1 0]);
            view(2);
            scatter(timeMat(twentyoneInds),aslMasked(twentyoneInds)./1000,10,'k','filled');
            scatter(timeMat(elevenInds),aslMasked(elevenInds)./1000,10,...
                'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7]);
            
            scatter(timeMat(twentyfourInds),aslMasked(twentyfourInds)./1000,10,...
                'MarkerEdgeColor',[0.45 0.76 0.42],'MarkerFaceColor',[0.45 0.76 0.42]);
            scatter(timeMat(twentythreeInds),aslMasked(twentythreeInds)./1000,10,...
                'MarkerEdgeColor',[0.7 0.8 0.87],'MarkerFaceColor',[0.7 0.8 0.87]);
            scatter(timeMat(twentytwoInds),aslMasked(twentytwoInds)./1000,10,...
                'MarkerEdgeColor',[0.17 0.45 0.7],'MarkerFaceColor',[0.17 0.45 0.7]);
            
            scatter(timeMat(fourteenInds),aslMasked(fourteenInds)./1000,10,'g','filled');
            scatter(timeMat(thirteenInds),aslMasked(thirteenInds)./1000,10,'c','filled');
            scatter(timeMat(twelveInds),aslMasked(twelveInds)./1000,10,'b','filled');
            
            plot(data.time(timeInds),icePlot./1000,'linewidth',1,'color',[0.6 0.6 0.6]);
            ax = gca;
            ax.SortMethod = 'childorder';
  
            ylim(ylimits);
            ylabel('Altitude (km)');
            xlim([startPlot,endPlot]);
            title('Reflectivity and melting layer')
            grid on
            set(gca,'xticklabel',[])
            ax2.Position=[0.06 0.525 0.87 0.21];
            
            
            %%%%%%%%%%%%%%%%%%%%%%%% LDR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ax3=subplot(4,1,3);
            hold on;
            sub3=surf(newTime,newASL./1000,newLDR,'edgecolor','none');
            view(2);
            caxis([-25 5]);
            colorbar
            ylim(ylimits);
            ylabel('Altitude (km)');
            xlim([startPlot,endPlot]);
            title('LDR')
            grid on
            set(gca,'xticklabel',[])
            ax3.Position=[0.06 0.287 0.87 0.21];
            
            ax4=subplot(4,1,4);
            ax4.Colormap=jet;
            hold on;
            sub3=surf(newTime,newASL./1000,newVEL,'edgecolor','none');
            view(2);
            caxis([-10 10]);
            colorbar
            ylim(ylimits);
            ylabel('Altitude (km)');
            xlim([startPlot,endPlot]);
            title('VEL')
            grid on
            ax4.Position=[0.06 0.05 0.87 0.21];
            
            linkaxes([ax1 ax2 ax3 ax4],'xy');
            
            formatOut = 'yyyymmdd_HHMM';
            set(gcf,'PaperPositionMode','auto')
            print([figdir,'meltRefl_',datestr(newTime(1),formatOut),'_to_',datestr(newTime(end),formatOut)],'-dpng','-r0');
                  
        end
        startPlot=endPlot;
    end
end