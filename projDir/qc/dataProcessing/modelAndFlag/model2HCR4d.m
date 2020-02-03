% find minimum reflectivity values
clear all;
close all;

addpath(genpath('/h/eol/romatsch/gitPriv/utils/'));

project='socrates'; % socrates, cset, aristo, otrec
quality='qc2'; % field, qc1, qc2
freqData='10hz'; % 10hz, 100hz, or 2hz
whichModel='era5'; % ecmwf or era5

formatOut = 'yyyymmdd_HHMM';

[modeldir outdir]=modelDir(project,whichModel,freqData);

topodir=topoDir(project);

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,freqData);

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
    
    data.asl=HCRrange2asl(data.range,data.elevation,data.altitude);
    
    %% Model data
    disp('Getting model data ...');
    if strcmp(whichModel,'era5')
        modelData=read_era5(modeldir,data.time(1),data.time(end),1);
    elseif strcmp(whichModel,'ecmwf')
        modelData=read_ecmwf(modeldir,data.time(1),data.time(end),1);
    end
    
    %% Topo data
    [modelData.topo modelData.topolon modelData.topolat]=read_gtopo30(topodir,modelData.lon,modelData.lat);
    
    %% Remove sst data that is over land
    lonMat=double(repmat(modelData.lon,1,size(modelData.z,2),size(modelData.z,4)));
    latMat=double(repmat(fliplr(modelData.lat'),size(modelData.z,1),1,size(modelData.z,4)));
    timeMat=repmat(datenum(modelData.time),size(modelData.z,1),1,size(modelData.z,2));
    timeMat=permute(timeMat,[1,3,2]);
    
    % Interpolate topo data to model grid
    topoModel=interpn(modelData.topolon',modelData.topolat',modelData.topo',...
        lonMat(:,:,1),latMat(:,:,1));
    
    if strcmp(whichModel,'ecmwf')
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
    
    % Topo
    surfData.zHCR=interpn(modelData.topolon',modelData.topolat',modelData.topo',...
        wrapTo360(data.longitude),data.latitude);
    
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
    xq=[timeMatHCR(:) aslGood(:)];
    nanXq=find(any(isnan(xq),2));
    xq(nanXq,:)=[];
    
    newGrid=(0:10:15000);
    [X Y]=meshgrid(datenum(data.time(timeInd)),newGrid);
    newTimeGrid=repmat(datenum(data.time(timeInd)),length(newGrid),1);
    
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
            % Then grab the data points a the HCR grid
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
            end
        end
    end
    
    disp(['Saving surface data ...']);
    timeHCR=data.time;
    save([outdir,whichModel,'.time.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(ii),'.mat'],'timeHCR');
    topo=surfData.zHCR;
    save([outdir,whichModel,'.topo.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(ii),'.mat'],'topo');
    aslHCR=data.asl;
    save([outdir,whichModel,'.asl.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(ii),'.mat'],'aslHCR');
    uSurfHCR=surfData.uHCR;
    save([outdir,whichModel,'.uSurf.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(ii),'.mat'],'uSurfHCR');
    vSurfHCR=surfData.vHCR;
    save([outdir,whichModel,'.vSurf.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(ii),'.mat'],'vSurfHCR');
    if isfield(modelData,'sstSurf')
        sstHCR=surfData.sstHCR;
        save([outdir,whichModel,'.sst.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
            datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(ii),'.mat'],'sstHCR');
    end
end

disp(datetime('now'));
