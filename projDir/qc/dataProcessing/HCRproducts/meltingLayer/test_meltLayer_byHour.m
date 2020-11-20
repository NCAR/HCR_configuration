% Analyze HCR clouds

clear all;
close all;

project='socrates'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz, 2hz, or combined

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

% figdir=['/scr/snow1/rsfdata/projects/otrec/hcr/qc2/cfradial/final2/10hz/plots/testHourly/'];
figdir='/home/romatsch/plots/HCR/meltingLayer/hourly/socrates/10hz/';

if ~exist(figdir, 'dir')
    mkdir(figdir)
end

ylimits=[-0.2 5];

%indir=HCRdir(project,quality,freqData);
indir='/run/media/romatsch/RSF0006/rsf/meltingLayer/socrates/10hz/';

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

for aa=2:size(caseList,1)
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
        
        if strcmp(freqData,'combined')
            data.HCR_DBZ=[];
            data.HCR_LDR=[];
            data.HCR_VEL=[];
        else
            data.DBZ=[];
            data.LDR=[];
            data.VEL_CORR=[];
        end
        data.MELTING_LAYER=[];
        data.ICING_LEVEL=[];
        data.FLAG=[];
        
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
        
        if strcmp(freqData,'combined')
            data.DBZ=data.HCR_DBZ;
            data=rmfield(data,'HCR_DBZ');
            data.LDR=data.HCR_LDR;
            data=rmfield(data,'HCR_LDR');
            data.VEL_CORR=data.HCR_VEL;
            data=rmfield(data,'HCR_VEL');
        end
        
        if isempty(data.DBZ)
            continue
        end
        
        %% Get indices
        
        elevenInds=find(data.MELTING_LAYER==11);
        twelveInds=find(data.MELTING_LAYER==12);
        thirteenInds=find(data.MELTING_LAYER==13);
        fourteenInds=find(data.MELTING_LAYER==14);
        
        twentyoneInds=find(data.MELTING_LAYER==21);
        twentytwoInds=find(data.MELTING_LAYER==22);
        twentythreeInds=find(data.MELTING_LAYER==23);
        twentyfourInds=find(data.MELTING_LAYER==24);
        
        %% Plot
        
        timeMat=repmat(data.time,size(data.LDR,1),1);
        dbzMasked=data.DBZ;
        dbzMasked(data.FLAG>1)=nan;
        ldrMasked=data.LDR;
        ldrMasked(data.FLAG>1)=nan;
        velMasked=data.VEL_CORR;
        velMasked(data.FLAG>1)=nan;
        
        close all
        
        if etime(datevec(endTime),datevec(startTime))<=900
            newInds=1:1:length(data.time);
        elseif etime(datevec(endTime),datevec(startTime))<=3600
            newInds=1:10:length(data.time);
        else
            newInds=1:100:length(data.time);
        end
        
        % Resample for plotting
        newDBZ=dbzMasked(:,newInds);
        newLDR=ldrMasked(:,newInds);
        newVEL=velMasked(:,newInds);
        newASL=data.asl(:,newInds);
        newFindMelt=data.MELTING_LAYER(:,newInds);
        newTime=data.time(newInds);
        
        fig1=figure('DefaultAxesFontSize',11,'position',[100,1300,1200,920]);
        
        ax1=subplot(4,1,1);
        hold on;
        sub1=surf(newTime,newASL./1000,newDBZ,'edgecolor','none');
        view(2);
        sub1=colMapDBZ(sub1);
        scatter(timeMat(twentyoneInds),data.asl(twentyoneInds)./1000,10,'k','filled');
        scatter(timeMat(elevenInds),data.asl(elevenInds)./1000,10,...
            'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7]);
        
        scatter(timeMat(twentyfourInds),data.asl(twentyfourInds)./1000,10,...
            'MarkerEdgeColor',[0.45 0.76 0.42],'MarkerFaceColor',[0.45 0.76 0.42]);
        scatter(timeMat(twentythreeInds),data.asl(twentythreeInds)./1000,10,...
            'MarkerEdgeColor',[0.7 0.8 0.87],'MarkerFaceColor',[0.7 0.8 0.87]);
        scatter(timeMat(twentytwoInds),data.asl(twentytwoInds)./1000,10,...
            'MarkerEdgeColor',[0.17 0.45 0.7],'MarkerFaceColor',[0.17 0.45 0.7]);
        
        scatter(timeMat(fourteenInds),data.asl(fourteenInds)./1000,10,'g','filled');
        scatter(timeMat(thirteenInds),data.asl(thirteenInds)./1000,10,'c','filled');
        scatter(timeMat(twelveInds),data.asl(twelveInds)./1000,10,'b','filled');
        
        ax = gca;
        ax.SortMethod = 'childorder';
        ylim(ylimits);
        ylabel('Altitude (km)');
        xlim([data.time(1),data.time(end)]);
        title(['Flight ',num2str(aa),': Reflectivity and melting layer'])
        grid on
        set(gca,'xticklabel',[])
        ax1.Position=[0.06 0.765 0.87 0.21];
        
        ax2=subplot(4,1,2);
        hold on;
        sub1=surf(newTime,newASL./1000,newFindMelt,'edgecolor','none');
        ax2.Colormap=([1 0 1;1 1 0]);
        view(2);
        scatter(timeMat(twentyoneInds),data.asl(twentyoneInds)./1000,10,'k','filled');
        scatter(timeMat(elevenInds),data.asl(elevenInds)./1000,10,...
            'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7]);
        
        scatter(timeMat(twentyfourInds),data.asl(twentyfourInds)./1000,10,...
            'MarkerEdgeColor',[0.45 0.76 0.42],'MarkerFaceColor',[0.45 0.76 0.42]);
        scatter(timeMat(twentythreeInds),data.asl(twentythreeInds)./1000,10,...
            'MarkerEdgeColor',[0.7 0.8 0.87],'MarkerFaceColor',[0.7 0.8 0.87]);
        scatter(timeMat(twentytwoInds),data.asl(twentytwoInds)./1000,10,...
            'MarkerEdgeColor',[0.17 0.45 0.7],'MarkerFaceColor',[0.17 0.45 0.7]);
        
        scatter(timeMat(fourteenInds),data.asl(fourteenInds)./1000,10,'g','filled');
        scatter(timeMat(thirteenInds),data.asl(thirteenInds)./1000,10,'c','filled');
        scatter(timeMat(twelveInds),data.asl(twelveInds)./1000,10,'b','filled');
        
        plot(data.time,data.ICING_LEVEL./1000,'linewidth',1,'color',[0.6 0.6 0.6]);
        ax = gca;
        ax.SortMethod = 'childorder';
        ylim(ylimits);
        ylabel('Altitude (km)');
        xlim([data.time(1),data.time(end)]);
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
        xlim([data.time(1),data.time(end)]);
        title('LDR')
        grid on
        set(gca,'xticklabel',[])
        ax3.Position=[0.06 0.287 0.87 0.21];
        
        ax4=subplot(4,1,4);
        ax4.Colormap=jet;
        hold on;
        sub3=surf(newTime,newASL./1000,newVEL,'edgecolor','none');
        view(2);
        caxis([-5 5]);
        colorbar
        ylim(ylimits);
        ylabel('Altitude (km)');
        xlim([data.time(1),data.time(end)]);
        title('VEL')
        grid on
        ax4.Position=[0.06 0.05 0.87 0.21];
        
        linkaxes([ax1 ax2 ax3 ax4],'xy');
        
        formatOut = 'yyyymmdd_HHMM';
        set(gcf,'PaperPositionMode','auto')
        print([figdir,'meltRefl_',datestr(data.time(1),formatOut),'_to_',datestr(data.time(end),formatOut)],'-dpng','-r0');
        
        startTime=endTime;
    end
end