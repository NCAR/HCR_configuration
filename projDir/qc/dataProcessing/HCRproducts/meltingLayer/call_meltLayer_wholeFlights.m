% Analyze HCR clouds

clear all;
close all;

project='socrates'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz
whichModel='era5';

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

%figdir=['/scr/snow1/rsfdata/projects/otrec/hcr/qc2/cfradial/final2/10hz/plots/'];
figdir='/home/romatsch/plots/HCR/meltingLayer/flights/socrates/combined/';

if ~exist(figdir, 'dir')
    mkdir(figdir)
end

ylimits=[-0.2 8];

%indir=HCRdir(project,quality,freqData);
indir='/run/media/romatsch/RSF0006/rsf/meltingLayer/socrates/combined/';

[~,directories.modeldir]=modelDir(project,whichModel,freqData);
%outdir=directories.modeldir;
outdir='/run/media/romatsch/RSF0006/rsf/meltingLayer/socratesMat/combined/';

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

zeroAdjustIn=300;
zeroAdjust=zeroAdjustIn;

for aa=1:size(caseList,1)
    disp(['Flight ',num2str(aa)]);
    disp('Loading HCR data.')
    
    clearvars -except project quality freqData whichModel figdir ...
        ylimits indir outdir caseList zeroAdjust zeroAdjustIn aa
    
    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    disp([datestr(startTime,'yyyy-mm-dd HH:MM'),' to ',datestr(endTime,'yyyy-mm-dd HH:MM')]);
    
    %% Load data
    
    disp('Loading data ...');
    
    data.DBZ=[];
    data.LDR=[];
    data.VEL_CORR=[];
    data.TEMP=[];
    data.WIDTH=[];
    data.FLAG=[];
    data.TOPO=[];
    
    dataVars=fieldnames(data);
    
    % Make list of files within the specified time frame
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if length(fileList)==0
        disp('No data files found.');
        continue
    end
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
    
    % Check if all variables were found and convert to single
    for ii=1:length(dataVars)
        if ~isfield(data,dataVars{ii})
            dataVars{ii}=[];
        end
    end
        
    dataVars=dataVars(~cellfun('isempty',dataVars));
    
    if isempty(data.DBZ)
        continue
    end
    
    %% Find melting layer
    
    for kk=1:2
        disp(['Round ',num2str(kk)]);
        
        data.dbzMasked=data.DBZ;
        data.dbzMasked(data.FLAG>1)=nan;
        
        findMelt=f_meltLayer(data,zeroAdjust);
        
        %% Calculate difference between zero degree and detected melting layer
        diffInd=[];
        
        for ii=1:size(findMelt,2)
            if data.elevation(ii)<-85
                bbRay=findMelt(:,ii);
                foundInd=min(find(bbRay==1));
                if ~isempty(foundInd)
                    zeroInd=min(find(bbRay==0));
                    diffInd=[diffInd foundInd-zeroInd];
                end
            end
        end
        
        meanDiff=mean(diffInd);
        zeroAdjust=meanDiff*(data.range(2)-data.range(1));
        
        if isnan(zeroAdjust)
            zeroAdjust=zeroAdjustIn;
        end
        
        disp(['Melting layer is on average ',num2str(meanDiff),' gates or ',num2str(zeroAdjust),' m below the zero degree level.'])
    end
    
    %% Find indices
    
    zeroInds=find(findMelt==0);
    oneInds=find(findMelt==1);
    twoInds=find(findMelt==2);
    threeInds=find(findMelt==3);
    
    timeMat=repmat(data.time,size(data.TEMP,1),1);
    
    %% Plot
    
    disp('Plotting ...');
    
    close all
    
    % Resample for plotting
    newInds=1:100:length(data.time);
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
    print([figdir,'meltRefl_Flight',num2str(aa)],'-dpng','-r0');
        
    %% Save
    disp('Saving meltLayer field ...')
    
    meltLayer=findMelt;
    save([outdir,whichModel,'.meltLayer.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'meltLayer');
    
    timeHCR=data.time;
    save([outdir,whichModel,'.time.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'timeHCR');
end