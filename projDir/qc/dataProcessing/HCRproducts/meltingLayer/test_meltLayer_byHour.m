% Analyze HCR clouds

clear all;
close all;

project='otrec'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz, or 2hz

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir=['/scr/snow1/rsfdata/projects/otrec/hcr/qc2/cfradial/final2/10hz/plots/testHourly/'];

if ~exist(figdir, 'dir')
    mkdir(figdir)
end

ylimits=[-0.2 8];

indir=HCRdir(project,quality,freqData);

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

for aa=1:size(caseList,1)
    disp(['Flight ',num2str(aa)]);
    
    startTime=datetime(caseList(aa,1:6));
    endTime=startTime;
    
    endTimeIn=datetime(caseList(aa,7:12));
    
    while endTime<endTimeIn
        endTime=endTime+hours(1);
        
        disp([datestr(startTime,'yyyy-mm-dd HH:MM'),' to ',datestr(endTime,'yyyy-mm-dd HH:MM')]);
        data=[];
        
        %% Load data
        
        disp('Loading data ...');
        
        data.DBZ=[];
        data.LDR=[];
        data.VEL_CORR=[];
        data.FREEZING_LEVEL=[];
        
        dataVars=fieldnames(data);
        
        % Make list of files within the specified time frame
        fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
        
        if length(fileList)==0
            disp('No data files found.');
            continue
        end
        
        % Load data
        data=read_HCR(fileList,data,startTime,endTime);
        
        % Check if all variables were found
        for ii=1:length(dataVars)
            if ~isfield(data,dataVars{ii})
                dataVars{ii}=[];
            end
        end
        
        dataVars=dataVars(~cellfun('isempty',dataVars));
        
        if isempty(data.DBZ)
            continue
        end
        
        %% Get indices
        
        zeroInds=find(data.FREEZING_LEVEL==0);
        oneInds=find(data.FREEZING_LEVEL==1);
        twoInds=find(data.FREEZING_LEVEL==2);
        threeInds=find(data.FREEZING_LEVEL==3);
        
        timeMat=repmat(data.time,size(data.DBZ,1),1);
        
        %% Plot
        close all
        
        disp('Plotting ...');
        
        % Resample for plotting
        newInds=1:20:length(data.time);
        newDBZ=data.DBZ(:,newInds);
        newLDR=data.LDR(:,newInds);
        newVEL=data.VEL_CORR(:,newInds);
        newASL=data.asl(:,newInds);
        newTime=data.time(newInds);
        
        fig3=figure('DefaultAxesFontSize',11,'position',[100,100,1500,1000]);
        
        s1=subplot(3,1,1);
        hold on;
        sub1=surf(newTime,newASL./1000,newDBZ,'edgecolor','none');
        view(2);
        sub1=colMapDBZ(sub1);
        colIndsAll=1:length(data.time);
        scatter(timeMat(zeroInds),data.asl(zeroInds)./1000,10,'k','filled');
        scatter(timeMat(oneInds),data.asl(oneInds)./1000,10,'b','filled');
        scatter(timeMat(twoInds),data.asl(twoInds)./1000,10,'c','filled');
        scatter(timeMat(threeInds),data.asl(threeInds)./1000,10,'g','filled');
        ax = gca;
        ax.SortMethod = 'childorder';
        ylim(ylimits);
        ylabel('Altitude (km)');
        xlim([data.time(1),data.time(end)]);
        title({['Flight ',num2str(aa),', ',project,', ',...
            datestr(data.time(1),'HH:MM'),' to ',datestr(data.time(end),'HH:MM')];['DBZ and freezing level']});
        grid on
        
        s2=subplot(3,1,2);
        hold on;
        sub2=surf(newTime,newASL./1000,newLDR,'edgecolor','none');
        view(2);
        colorbar
        ylim(ylimits);
        ylabel('Altitude (km)');
        xlim([data.time(1),data.time(end)]);
        title('LDR')
        grid on
        
        s3=subplot(3,1,3);
        hold on;
        sub2=surf(newTime,newASL./1000,newVEL,'edgecolor','none');
        view(2);
        s3.Colormap=jet;
        colorbar
        caxis([-5 5]);
        ylim(ylimits);
        ylabel('Altitude (km)');
        xlim([data.time(1),data.time(end)]);
        title('VEL')
        grid on
        
        formatOut = 'yyyymmdd_HHMM';
        set(gcf,'PaperPositionMode','auto')
        print([figdir,'meltRefl_',datestr(data.time(1),formatOut),'_to_',datestr(data.time(end),formatOut)],'-dpng','-r0');
        
        startTime=endTime;
    end
end