% Find maximum reflectivity in drizzle at 250 m from cloud edge and temperatures > 2C
% Only when pointing down
% In a calibrated radar these values should be ~18-20 dBZ

clear all;
close all;

project='cset'; % socrates, cset, aristo, otrec
quality='qc2'; % field, qc1, qc2
freqData='10hz'; % 10hz, 100hz, or 2hz

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

directories.figdir='/h/eol/romatsch/hcrCalib/otherCalib/drizzle/';
formatOut = 'yyyymmdd_HHMM';

indir=HCRdir(project,quality,freqData);

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/other/otherCalib/inFiles/drizzle_',project,'.txt'];

% Read file with calib events
caseList = table2array(readtable(infile));

edgesAll=[-100:1:100];
Nall=zeros(1,length(edgesAll)-1);

for ii=1:size(caseList,1)
    
    disp(['Flight ',num2str(ii)]);
    
    disp('Loading data ...');
    
    reflRangeAll=[];
    
    data=[];
    
    data.DBZ = [];  % reflectivity
    data.TEMP=[];
    data.FLAG=[];
    
    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    dataVars=fieldnames(data);
    
    if length(fileList)==0
        disp('No data files found.');
        return
    end
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
    
    % Check if all variables were found
    for kk=1:length(dataVars)
        if ~isfield(data,dataVars{kk})
            dataVars{kk}=[];
        end
    end
    
    dataVars=dataVars(~cellfun('isempty',dataVars));
    
    %% Use right data
    
    reflTemp=data.DBZ;
    tempTemp=data.TEMP;
    
    %sort out upward pointing and non cloud
    reflTemp(data.FLAG>1)=nan;
    tempTemp(data.FLAG>1)=nan;
    
    reflTemp(:,data.elevation>0)=nan;
    tempTemp(:,data.elevation>0)=nan;
    
    %% Find zero degree altitude
    tempVec=1:1:size(data.TEMP,1);
    
    tempMat=repmat(tempVec,size(data.TEMP,2),1)';
    negInds=find(tempTemp<=0 | isnan(tempTemp));
    tempMat(negInds)=nan;
    
    [mins,rowInds] = nanmin(tempMat);
    colInds=1:1:size(data.TEMP,2);
    
    linearInd = sub2ind(size(data.TEMP), rowInds, colInds);
    zTempsPos=data.asl(linearInd);
    
    zTempsPos(rowInds==1)=nan;
        
    %% Go through rays
    usedInd=[];
    edgeAltAll=[];
    
    % Go through each ray
    disp('Going through each ray...');
    for jj=1:length(data.time)
        
        % Get refl data
        reflRay=reflTemp(:,jj);
        
        if min(isnan(reflRay))==1
            continue
        end
        
        % Find cloud edge
        cloudEdge=min(find(~isnan(reflRay)));
        
         % Check if temperature is high enough
         edgeAlt=data.asl(cloudEdge,jj);
         eraAlt=zTempsPos(jj);
         
         if edgeAlt>eraAlt | isnan(eraAlt) | cloudEdge<20
             continue
         end        
        
         % Get data from cloud edge to 270 m into the cloud
         rangeRay=data.range(:,jj);
         cloudInds=find(rangeRay-rangeRay(cloudEdge)>=0 & rangeRay-rangeRay(cloudEdge)<=270);
         cloudData=reflRay(cloudInds);
         
         % Check if contiguous and if enough data
         rangeData=rangeRay(cloudInds);
         if max(isnan(cloudData))==1 | rangeData(end)-rangeData(1)<230
             continue
         end
         
         % Get data from 230 to 270 m
         coreInds=find(rangeRay-rangeRay(cloudEdge)>=230 & rangeRay-rangeRay(cloudEdge)<=270);
         reflRangeAll=[reflRangeAll;reflRay(coreInds)];
         
         usedInd=[usedInd, jj];
         edgeAltAll=[edgeAltAll,edgeAlt./1000];
    end
    
    [N,ed] = histcounts(reflRangeAll,edgesAll);
    Nall=Nall+N;
    
    close all
    
    f0=figure('DefaultAxesFontSize',12);
    set(f0,'Position',[200 500 1000 600]);
    
    hold on
    plot(data.time,data.altitude./1000,'-b','linewidth',3);
    %plot(data.time(usedInd),data.altitude(usedInd)./1000,'or');
    plot(data.time(usedInd),edgeAltAll,'og');
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    
    title([project,' flight ',num2str(ii)]);
    
    set(f0,'PaperPositionMode','auto')
    print(f0, [directories.figdir,project,'_cloudEdge_dataUsed_flight',num2str(ii)],'-dpng','-r0');
end

%% Plot histogram

close all

f1=figure('DefaultAxesFontSize',12);
set(f1,'Position',[200 500 1000 1000]);

subplot(2,1,1)

bar(edgesAll(1:end-1),Nall,1,'edgecolor','k');
title([project,' histogram of reflectivities at 250 m range'],'interpreter','none');
xlim([-50 50])
text(0,45000,['Total points: ', num2str(length(reflRangeAll))],'fontsize',14);
xlabel('Reflectivity (dBZ)');

subplot(2,1,2)

bar(edgesAll(1:end-1),Nall,1,'edgecolor','k');
title('top end');
xlim([0 60])
ylim([0 4000])
xlabel('Reflectivity (dBZ)');

set(f1,'PaperPositionMode','auto')
print(f1, [directories.figdir,project,'_cloudEdge_drizzleRefl'],'-dpng','-r0');

