% Plot HCR pid from mat file in hourly plots

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='noreaster'; %socrates, aristo, cset, otrec
quality='qc2'; %field, qc1, or qc2
dataFreq='10hz';
qcVersion='v2.0';
whichModel='era5';

if strcmp(project,'otrec')
    ylimits=[0 15];
else
    ylimits=[0 10];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

indir=HCRdir(project,quality,qcVersion,dataFreq);

[~,directories.modeldir]=modelDir(project,whichModel,quality,qcVersion,dataFreq);
modeldir=directories.modeldir;

figdir=[indir(1:end-5),'meltLayerPlots/testMat/'];

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'.txt'];

caseList = table2array(readtable(infile));

for aa=1:size(caseList,1)
    disp(['Flight ',num2str(aa)]);
    disp('Loading HCR data.')
    disp(['Starting at ',datestr(datetime('now'),'yyyy-mm-dd HH:MM')]);
    
    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));
    
    %% Get data
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    data=[];
    
    data.DBZ_MASKED = [];
    data.LDR=[];
    data.VEL_UNFOLDED=[];
        
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
    
    %% Load pid
    disp('Loading melt layer.');
    
    fileIn1=dir([modeldir,whichModel,'.meltLayer.*.Flight',num2str(aa),'.mat']);
    meltLayerIn=load([modeldir,fileIn1.name]);
    meltLayerGet=meltLayerIn.meltLayer;
    
    fileIn2=dir([modeldir,whichModel,'.iceLevel.*.Flight',num2str(aa),'.mat']);
    iceLevelIn=load([modeldir,fileIn2.name]);
    iceLevelGet=iceLevelIn.iceLayer;
    
    fileIn3=dir([modeldir,whichModel,'.time.*.Flight',num2str(aa),'.mat']);
    meltTime=load([modeldir,fileIn3.name]);
    timeMELT=meltTime.timeHCR;
    
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