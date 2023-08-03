% Plot HCR pid from mat file in hourly plots

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='socrates'; %socrates, aristo, cset
quality='qc3'; %field, qc1, or qc2
qcVersion='v3.0';
freqData='combined'; % 10hz, 100hz, 2hz, or combined
whichModel='era5';

showPlot='off';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

indir=HCRdir(project,quality,qcVersion,freqData);

figdir=[indir(1:end-4),'pidPlotsComb/compare_HCR_HSRLHCR/'];

cscale_hcr=[1,0,0; 1,0.6,0.47; 0,1,0; 0,0.7,0; 0,0,1; 1,0,1; 0.5,0,0; 1,1,0; 0,1,1; 0,0,0; 0.5,0.5,0.5];
units_str_hcr={'Rain','Supercooled Rain','Drizzle','Supercooled Drizzle','Cloud Liquid','Supercooled Cloud Liquid',...
    'Mixed Phase','Large Frozen','Small Frozen','Precip','Cloud'};

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

catDiffAll=zeros(1,10);

for aa=1:15
    disp(['Flight ',num2str(aa)]);
    disp('Loading data ...')


    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));

    %% Get HCR data

    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    data=[];
    data.HCR_PID=[];
    data.PID=[];
    data.HCR_MELTING_LAYER=[];

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

    if isempty(data.time)
        startTime=endTime;
        continue
    end

    %% Same

    data.HCR_PID=round(data.HCR_PID);

    data.HCR_PID(data.HCR_MELTING_LAYER<20)=nan;
    data.PID(data.HCR_MELTING_LAYER<20)=nan;

    data.HCR_PID(data.HCR_PID>9)=nan;
    data.PID(isnan(data.HCR_PID))=nan;

    numSame=length(find(data.HCR_PID==data.PID));

    %% Different

    fillHCR=data.HCR_PID;
    fillHCR(isnan(data.HCR_PID))=-99;

    fillPID=data.PID;
    fillPID(isnan(data.PID))=-99;

    diffInds=find(fillHCR~=fillPID);
    numDiff=length(diffInds);

    catDiff=nan(1,10);

    diffOnly=nan(size(fillHCR));
    diffOnly(diffInds)=fillHCR(diffInds);

    for ii=1:9
        catDiff(ii)=length(find(diffOnly==ii));
    end

    catDiff(10)=numSame;

    catDiffAll=catDiffAll+catDiff;

    %% Plot histogram
    close all

    colBar=cat(1,cscale_hcr(1:9,:),[0.9,0.9,0.9]);
    barLabels=cat(1,units_str_hcr{1:9},{'Unchanged'});

    f1 = figure('Position',[200 500 1000 800],'DefaultAxesFontSize',12,'visible','on');

    b=bar(1:10,catDiff./sum(catDiff).*100,1);

    b.FaceColor='flat';
    for jj=1:10
        b.CData(jj,:)=colBar(jj,:);
    end

    xlim([0.5,10.5]);
    ylim([0 100]);

    xticklabels(barLabels);

    title(['Flight ',num2str(aa)]);
    set(gca, 'YGrid', 'on', 'XGrid', 'off')

    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_Flight',num2str(aa),'_catChanges.png'],'-dpng','-r0')

    %% Plot in half hourly increments
    
    disp('Plotting comparison ...');
    ylimUpper=12;
    
    startPlot=startTime;
           
    while startPlot<endTime
        
        close all
        
        endPlot=startPlot+minutes(30);
        timeInds=find(data.time>=startPlot & data.time<=endPlot);
        
        timePlot=data.time(timeInds);
        pidHCRPlot=data.HCR_PID(:,timeInds);
        pidPlot=data.PID(:,timeInds);
        
        if sum(sum(~isnan(pidPlot)))~=0
            aslPlot=data.asl(:,timeInds);
            
            timeIndsPID=find(data.time>=startPlot & data.time<=endPlot);
                            
            f1 = figure('Position',[200 500 1800 1000],'DefaultAxesFontSize',12,'visible','off');
            
            s1=subplot(2,1,1);
                        
            hold on
            surf(timePlot,aslPlot./1000,pidHCRPlot,'edgecolor','none');
            view(2);
            colormap(s1,cscale_hcr);
            cb=colorbar;
            cb.Ticks=1:11;
            cb.TickLabels=units_str_hcr;
            ylabel('Altitude (km)');
            title(['HCR particle ID']);
            caxis([.5 11.5]);
            ylim([0 ylimUpper]);
            xlim([timePlot(1),timePlot(end)]);

            grid on
            box on
            s1pos=s1.Position;
            s1.Position=[s1pos(1),s1pos(2),s1pos(3),s1pos(4)];
            
            s2=subplot(2,1,2);
            
            hold on
            surf(timePlot,aslPlot./1000,pidPlot,'edgecolor','none');
            view(2);
            colormap(s2,cscale_hcr);
            cb=colorbar;
            cb.Ticks=1:11;
            cb.TickLabels=units_str_hcr;
            ylabel('Altitude (km)');
            title(['HCR/HSRL particle ID']);
            caxis([.5 11.5]);
            ylim([0 ylimUpper]);
            xlim([timePlot(1),timePlot(end)]);
            
            grid on
            box on
            s2pos=s2.Position;
            s2.Position=[s2pos(1),s2pos(2),s1pos(3),s2pos(4)];
            
            set(gcf,'PaperPositionMode','auto')
            print(f1,[figdir,project,'_Flight',num2str(aa),'_pid_',datestr(timePlot(1),'yyyymmdd_HHMMSS'),'_to_',datestr(timePlot(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
        end
        startPlot=endPlot;
    end

end

%% Plot histogram
close all

colBar=cat(1,cscale_hcr(1:9,:),[0.9,0.9,0.9]);
barLabels=cat(1,units_str_hcr{1:9},{'Unchanged'});

f1=figure('Position',[200 500 1000 800],'DefaultAxesFontSize',12,'visible','on');

b=bar(1:10,catDiffAll./sum(catDiffAll).*100,1);

b.FaceColor='flat';
for jj=1:10
    b.CData(jj,:)=colBar(jj,:);
end

xlim([0.5,10.5]);
ylim([0 100]);

xticklabels(barLabels);

title([project,' category changes']);
set(gca, 'YGrid', 'on', 'XGrid', 'off')

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_catChanges.png'],'-dpng','-r0')

save([figdir,'catChanges.mat'],'catDiffAll');
