% find minimum reflectivity values
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='noreaster'; % socrates, cset, aristo, otrec
quality='qc2'; % field, qc0, qc1, qc2
qcVersion='v2.1';
freqData='10hz'; % 10hz, 100hz, or 2hz
whichModel='era5'; % ecmwf or era5 or narr

era5levelFiles=1;

addTopo=0; % Set to 1 if topo data should be added and hasn't been added in separate script.
getSST=1;

formatOut = 'yyyymmdd_HHMM';

[modeldir outdir]=modelDir(project,whichModel,quality,qcVersion,freqData);

topodir=topoDir(project);

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,qcVersion,freqData);

%% Go through flights
for ii=1:size(caseList,1)
    disp(['Flight ',num2str(ii)]);
    
    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));
        
    %% HCR data
    disp('Getting HCR data ...');
    
    % Desired variables. The variable name comies after the . and must be spelled exactly
    % like in the CfRadial file
    if exist('data')
        clear data
    end
    
    data.dummy=[];
    
    dataVars=fieldnames(data);
    
    %% Load data
    % Make list of files within the specified time frame
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if length(fileList)==0
        disp('No data files found.');
        startTime=endTime;
        continue
    end
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
    
    if isempty(data.time)
        disp('No data found.');
        startTime=endTime;
        continue
    end
    
    % Check if all variables were found
    for kk=1:length(dataVars)
        if ~isfield(data,dataVars{kk})
            dataVars{kk}=[];
        end
    end

    %% We have some data where lat/lon are zero
    lonZero=find(data.longitude==0);
    if ~isempty(lonZero)
        warning([num2str(length(lonZero)),' zero longitudes replaced.']);
    end
    data.longitude(lonZero)=data.longitude(lonZero-1);

    latZero=find(data.latitude==0);
    if ~isempty(latZero)
        warning([num2str(length(latZero)),' zero latitudes replaced.']);
    end
    data.latitude(latZero)=data.latitude(latZero-1);
        
    %% Model data
    disp('Getting model data ...');
    if strcmp(whichModel,'era5')
        if era5levelFiles
            modelData=read_era5_levelFiles(modeldir,data.time(1),data.time(end),getSST);
        else
            modelData=read_era5_oneFile(modeldir,data.time(1),data.time(end),getSST);
        end
    elseif strcmp(whichModel,'ecmwf')
        modelData=read_ecmwf(modeldir,data.time(1),data.time(end),getSST);
    elseif strcmp(whichModel,'narr')
        modelData=read_narr(modeldir,data.time(1),data.time(end),getSST);
    end
    
    %% Topo data
    [modelData.topo modelData.topolon modelData.topolat]=read_gtopo30(topodir,data.longitude,data.latitude);
    
    %% Set up grid
    if strcmp(whichModel,'narr')
        lonMat=double(repmat(modelData.lon,1,1,size(modelData.z,4)));
        latMat=double(repmat(modelData.lat,1,1,size(modelData.z,4)));
    else
        lonMat=double(repmat(modelData.lon,1,size(modelData.z,2),size(modelData.z,4)));
        latMat=double(repmat(fliplr(modelData.lat'),size(modelData.z,1),1,size(modelData.z,4)));
    end
    timeMat=repmat(datenum(modelData.time),size(modelData.z,1),1,size(modelData.z,2));
    timeMat=permute(timeMat,[1,3,2]);
    
    % Remove SST data that is over land
    if strcmp(whichModel,'ecmwf')
        % Interpolate topo data to model grid
        topoModel=interpn(modelData.topolon',modelData.topolat',modelData.topo',...
            lonMat(:,:,1),latMat(:,:,1));
        
        for ll=1:size(modelData.sstSurf,3)
            tempSST=modelData.sstSurf(:,:,ll);
            tempSST(topoModel>0)=nan;
            modelData.sstSurf(:,:,ll)=tempSST;
        end
    end
    %% Interpolate
    disp('Interpolating to HCR track ...');
    
    % Make thinned out time vector
    % Get 10 second interval
    %     startMinute=datetime(year(data.time(1)),month(data.time(1)),day(data.time(1)),...
    %         hour(data.time(1)),minute(data.time(1)),0);
    %     indInt=find(data.time==startMinute+minutes(1)+seconds(10))-find(data.time==startMinute+minutes(1));
    if ~isempty(strfind(indir,'10hz'))
        indInt=100;
    elseif ~isempty(strfind(indir,'100hz'))
        indInt=1000;
    elseif ~isempty(strfind(indir,'2hz'))
        indInt=20;
    end
    timeInd=1:indInt:length(data.time);
    
    % 3D variables
    int.tempHCR=[];
    int.zHCR=[];
    int.pHCR=[];
    int.rhHCR=[];
    int.uHCR=[];
    int.vHCR=[];
    
    if strcmp(whichModel,'narr')
        
        % Make grid smaller
        minLonFlight=min(wrapTo360(data.longitude(timeInd)));
        maxLonFlight=max(wrapTo360(data.longitude(timeInd)));
        minLatFlight=min(data.latitude(timeInd));
        maxLatFlight=max(data.latitude(timeInd));
        [r c]=find(modelData.lon>=minLonFlight & modelData.lon<=maxLonFlight & ...
            modelData.lat>=minLatFlight & modelData.lat<=maxLatFlight);
        
        lonMatS=lonMat(min(r)-1:max(r)+1,min(c)-1:max(c)+1,:);
        latMatS=latMat(min(r)-1:max(r)+1,min(c)-1:max(c)+1,:);
        timeMatS=timeMat(min(r)-1:max(r)+1,min(c)-1:max(c)+1,:);
        
        xIn=cat(2,lonMatS(:),latMatS(:),timeMatS(:));
        xqIn=cat(2,wrapTo360(data.longitude(timeInd))',data.latitude(timeInd)',datenum(data.time(timeInd))');
        goodIndXQ=find(~any(isnan(xqIn),2));
        xqIn(any(isnan(xqIn),2),:)=[];        
        
        for jj=1:size(modelData.z,3)
            vIn1=squeeze(modelData.Temperature(min(r)-1:max(r)+1,min(c)-1:max(c)+1,jj,:));
            vIn=vIn1(:);
            Vq = griddatan(xIn,vIn,xqIn);
            VqOut=nan(length(timeInd),1);
            VqOut(goodIndXQ)=Vq;
            int.tempHCR=cat(1,int.tempHCR,VqOut');
            vIn1=squeeze(modelData.z(min(r)-1:max(r)+1,min(c)-1:max(c)+1,jj,:));
            vIn=vIn1(:);
            Vq = griddatan(xIn,vIn,xqIn);
            VqOut=nan(length(timeInd),1);
            VqOut(goodIndXQ)=Vq;
            int.zHCR=cat(1,int.zHCR,VqOut');
            vIn1=squeeze(modelData.p(min(r)-1:max(r)+1,min(c)-1:max(c)+1,jj,:));
            vIn=vIn1(:);
            Vq = griddatan(xIn,vIn,xqIn);
            VqOut=nan(length(timeInd),1);
            VqOut(goodIndXQ)=Vq;
            int.pHCR=cat(1,int.pHCR,VqOut');
            vIn1=squeeze(modelData.rh(min(r)-1:max(r)+1,min(c)-1:max(c)+1,jj,:));
            vIn=vIn1(:);
            Vq = griddatan(xIn,vIn,xqIn);
            VqOut=nan(length(timeInd),1);
            VqOut(goodIndXQ)=Vq;
            int.rhHCR=cat(1,int.rhHCR,VqOut');
        end
        
        % 2D variables
        xqIn=cat(2,wrapTo360(data.longitude)',data.latitude',datenum(data.time)');
        goodIndXQ=find(~any(isnan(xqIn),2));
        xqIn(any(isnan(xqIn),2),:)=[];
        
        vIn1=squeeze(modelData.pSurf(min(r)-1:max(r)+1,min(c)-1:max(c)+1,:));
        vIn=vIn1(:);
        Vq=griddatan(xIn,vIn,xqIn);
        surfData.pHCR=nan(length(data.time),1);
        surfData.pHCR(goodIndXQ)=Vq;
        vIn1=squeeze(modelData.tSurf(min(r)-1:max(r)+1,min(c)-1:max(c)+1,:));
        vIn=vIn1(:);
        Vq=griddatan(xIn,vIn,xqIn);
        surfData.tempHCR=nan(length(data.time),1);
        surfData.tempHCR(goodIndXQ)=Vq;
        vIn1=squeeze(modelData.rhSurf(min(r)-1:max(r)+1,min(c)-1:max(c)+1,:));
        vIn=vIn1(:);
        Vq=griddatan(xIn,vIn,xqIn);
        surfData.rhHCR=nan(length(data.time),1);
        surfData.rhHCR(goodIndXQ)=Vq;
        vIn1=squeeze(modelData.uSurf(min(r)-1:max(r)+1,min(c)-1:max(c)+1,:));
        vIn=vIn1(:);
        Vq=griddatan(xIn,vIn,xqIn);
        surfData.uHCR=nan(length(data.time),1);
        surfData.uHCR(goodIndXQ)=Vq;
        vIn1=squeeze(modelData.vSurf(min(r)-1:max(r)+1,min(c)-1:max(c)+1,:));
        vIn=vIn1(:);
        Vq=griddatan(xIn,vIn,xqIn);
        surfData.vHCR=nan(length(data.time),1);
        surfData.vHCR(goodIndXQ)=Vq;
        if isfield(modelData,'sstSurf')
            vIn1=squeeze(modelData.sstSurf(min(r)-1:max(r)+1,min(c)-1:max(c)+1,:));
            vIn=vIn1(:);
            Vq=griddatan(xIn,vIn,xqIn);
            surfData.sstHCR=nan(length(data.time),1);
            surfData.sstHCR(goodIndXQ)=Vq;
        end
        
    else
        for jj=1:size(modelData.z,3)
            Vq = interpn(lonMat,latMat,timeMat,squeeze(modelData.Temperature(:,:,jj,:)),...
                wrapTo360(data.longitude(timeInd)),data.latitude(timeInd),datenum(data.time(timeInd)));
            int.tempHCR=cat(1,int.tempHCR,Vq);
            Vq = interpn(lonMat,latMat,timeMat,squeeze(modelData.z(:,:,jj,:)),...
                wrapTo360(data.longitude(timeInd)),data.latitude(timeInd),datenum(data.time(timeInd)));
            int.zHCR=cat(1,int.zHCR,Vq);
            Vq = interpn(lonMat,latMat,timeMat,squeeze(modelData.p(:,:,jj,:)),...
                wrapTo360(data.longitude(timeInd)),data.latitude(timeInd),datenum(data.time(timeInd)));
            int.pHCR=cat(1,int.pHCR,Vq);
            Vq = interpn(lonMat,latMat,timeMat,squeeze(modelData.rh(:,:,jj,:)),...
                wrapTo360(data.longitude(timeInd)),data.latitude(timeInd),datenum(data.time(timeInd)));
            int.rhHCR=cat(1,int.rhHCR,Vq);
            Vq = interpn(lonMat,latMat,timeMat,squeeze(modelData.u(:,:,jj,:)),...
                wrapTo360(data.longitude(timeInd)),data.latitude(timeInd),datenum(data.time(timeInd)));
            int.uHCR=cat(1,int.uHCR,Vq);
            Vq = interpn(lonMat,latMat,timeMat,squeeze(modelData.v(:,:,jj,:)),...
                wrapTo360(data.longitude(timeInd)),data.latitude(timeInd),datenum(data.time(timeInd)));
            int.vHCR=cat(1,int.vHCR,Vq);
        end
        
        % 2D variables
        surfData.pHCR = interpn(lonMat,latMat,timeMat,modelData.pSurf,...
            wrapTo360(data.longitude),data.latitude,datenum(data.time));
        surfData.tempHCR = interpn(lonMat,latMat,timeMat,modelData.tSurf,...
            wrapTo360(data.longitude),data.latitude,datenum(data.time));
        surfData.rhHCR = interpn(lonMat,latMat,timeMat,modelData.rhSurf,...
            wrapTo360(data.longitude),data.latitude,datenum(data.time));
        surfData.uHCR = interpn(lonMat,latMat,timeMat,modelData.uSurf,...
            wrapTo360(data.longitude),data.latitude,datenum(data.time));
        surfData.vHCR = interpn(lonMat,latMat,timeMat,modelData.vSurf,...
            wrapTo360(data.longitude),data.latitude,datenum(data.time));
        if isfield(modelData,'sstSurf')
            surfData.sstHCR = interpn(lonMat,latMat,timeMat,modelData.sstSurf,...
                wrapTo360(data.longitude),data.latitude,datenum(data.time));
        end
    end
    
    % Interpolate topo
    surfData.zHCR=interpn(modelData.topolon',modelData.topolat',modelData.topo',...
        wrapTo360(data.longitude),data.latitude);
    
    if size(surfData.zHCR,1)==1
        surfData.zHCR=surfData.zHCR';
    end
    
    intFields=fields(int);
    
    % Remove nans
    for ll=1:length(intFields)
        int.(intFields{ll})= int.(intFields{ll})(any(~isnan(int.(intFields{ll})),2),:);
    end
    
    surfatInds=surfData.zHCR(timeInd);
    
    % Replace or add surface values
    for mm=1:length(data.time(timeInd))
        zLevels=int.zHCR(:,mm);
        zSurf=surfatInds(mm);
        zLevels(zLevels<=zSurf)=nan;
        nanInds=find(isnan(zLevels));
        if isempty(nanInds)
            surfInd=length(zLevels);
        else
            surfInd=(min(nanInds));
        end
        for ll=1:length(intFields);
            int.(intFields{ll})(nanInds,mm)=nan;
            int.(intFields{ll})(surfInd,mm)=surfData.(intFields{ll})(timeInd(mm));
        end
    end
    
    %% Interpolate to HCR grid
    
    disp('Interpolating to HCR grid ...');
    
    timeMatModel=repmat(datenum(data.time(timeInd)),size(int.zHCR,1),1);
    
    % Remove data that is too far below surface or too far out
    aslGood=data.asl;
    outInds=find(aslGood<-200 | aslGood>15000);
    aslGood(outInds)=nan;
    keepInds=find(~isnan(aslGood));
    
    % Output coordinates
    timeMatHCR=repmat(datenum(data.time),size(data.range,1),1);
    xq=[timeMatHCR(:) double(aslGood(:))];
    nanXq=find(any(isnan(xq),2));
    xq(nanXq,:)=[];
    
    newGrid=(0:10:15000);
    [X Y]=meshgrid(datenum(data.time(timeInd)),newGrid);
    %newTimeGrid=repmat(datenum(data.time(timeInd)),length(newGrid),1);
    
    for ll=1:length(intFields)
        if ~strcmp(intFields{ll},'zHCR')
            disp(intFields{ll});
            % First make 1d interpolation to a regular grid
            v=int.(intFields{ll});
            vq=[];
            for mm=1:size(timeMatModel,2)
                x=int.zHCR(:,mm);
                y=v(:,mm);
                xy=cat(2,x,y);
                nanXY=find(any(isnan(xy),2));
                xy(nanXY,:)=[];
                if ~isempty(xy)
                    vq1 = interp1(xy(:,1),xy(:,2),newGrid);
                    vq=cat(2,vq,vq1');
                else
                    vq=cat(2,vq,nan(length(newGrid),1));
                end
            end
            % Then grab the data points at the HCR grid
            Vq = interp2(X,Y,vq,xq(:,1),xq(:,2));
            modelvar=nan(size(data.range));
            modelvar(keepInds)=Vq;
            %                         %surf(data.time(105000:150000),data.asl(:,105000:150000),modelvar(:,105000:150000),'edgecolor','none');
            %                         surf(data.time,data.asl,modelvar,'edgecolor','none');
            %                         view(2);
            disp(['Saving ',intFields{ll},' data ...']);
            if strcmp(intFields{ll},'tempHCR')
                tempHCR=modelvar;
                save([outdir,whichModel,'.',intFields{ll},'.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
                    datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(ii),'.mat'],'tempHCR');
            elseif strcmp(intFields{ll},'pHCR')
                pHCR=modelvar;
                save([outdir,whichModel,'.',intFields{ll},'.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
                    datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(ii),'.mat'],'pHCR');
            elseif strcmp(intFields{ll},'rhHCR')
                rhHCR=modelvar;
                save([outdir,whichModel,'.',intFields{ll},'.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
                    datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(ii),'.mat'],'rhHCR');
            elseif strcmp(intFields{ll},'uHCR')
                uHCR=modelvar;
                save([outdir,whichModel,'.',intFields{ll},'.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
                    datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(ii),'.mat'],'uHCR');
            elseif strcmp(intFields{ll},'vHCR')
                vHCR=modelvar;
                save([outdir,whichModel,'.',intFields{ll},'.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
                    datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(ii),'.mat'],'vHCR');
            end
        end
    end
    
    disp(['Saving surface data ...']);
    timeHCR=data.time;
    save([outdir,whichModel,'.time.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(ii),'.mat'],'timeHCR');
    if addTopo
        topo=surfData.zHCR;
        save([outdir,whichModel,'.topo.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
            datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(ii),'.mat'],'topo');
    end
% %     aslHCR=data.asl;
% %     save([outdir,whichModel,'.asl.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
% %         datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(ii),'.mat'],'aslHCR');
%     uSurfHCR=surfData.uHCR;
%     save([outdir,whichModel,'.uSurf.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
%         datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(ii),'.mat'],'uSurfHCR');
%     vSurfHCR=surfData.vHCR;
%     save([outdir,whichModel,'.vSurf.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
%         datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(ii),'.mat'],'vSurfHCR');
    if isfield(modelData,'sstSurf')
        sstHCR=surfData.sstHCR;
        save([outdir,whichModel,'.sst.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
            datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(ii),'.mat'],'sstHCR');
    end
end

