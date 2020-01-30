% Analyze HCR clouds

clear all;
close all;

startTime=datetime(2018,1,23,3,30,0);
endTime=datetime(2018,1,23,3,40,0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

% If 1, plots for individual calibration events will be made, if 0, only
% total plots will be made
makeSingleFigs=1;
addNSCAL=0; % Noise source cal will be incorporated
addSIG0model=1; % Plot also model data from sig0 model

project='socrates'; %socrates, aristo, cset
quality='qc1'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz, or 2hz
whichModel='era5';
addName=''; % Extra name part for output files. Default is ''.

salinity=35; % Ocean salinity for sig0model in per mille (world wide default is 35) and sensitivity to that number is low
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/h/eol/romatsch/gitPriv/process_HCR/oceanScans/functions/');
addpath('/h/eol/romatsch/gitPriv/process_HCR/oceanScans/colormaps/');
addpath('/h/eol/romatsch/gitPriv/process_HCR/NSCAL/functions/');
addpath(genpath('/h/eol/romatsch/gitPriv/utils/'));

outName=[project,'_stratConv'];

outName=[outName,addName];

directories.figdir=['/h/eol/romatsch/hcrCalib/stratConv/',outName,'/2hz/'];

if ~exist(directories.figdir, 'dir')
    mkdir(directories.figdir)
end

directories.dataDir=HCRdir(project,quality,freqData);

if strcmp(whichModel,'era5')
    directories.modeldir=['/scr/sci/romatsch/data/reanalysis/ecmwf/era5interp/',project,'/',freqData,'/'];
elseif strcmp(whichModel,'ecmwf')
    directories.modeldir=['/scr/sci/romatsch/data/reanalysis/ecmwf/forecastInterp/',project,'/',freqData,'/'];
end

%% Load data

data.DBZ=[];
%data.VEL=[];
%data.VEL_RAW=[];
%data.VEL_CORR=[];
%data.WIDTH=[];
%data.WIDTH_CORR=[];
%data.DBMVC=[];
%data.SNR=[];
%data.NCP=[];
%data.LDR=[];

dataVars=fieldnames(data);

% Make list of files within the specified time frame
fileList=makeFileList(directories.dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

if length(fileList)==0
    disp('No data files found.');
    return
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

%% Reflectivity profile

refl=data.DBZ;

% Remove upward pointing data
upInd=find(data.elevation>0);
%refl(:,upInd)=nan;

% Remove bang
refl(1:20,:)=nan;

% Remove everything below surface+10 gates
[surfLinInds surfRowInds rangeToSurf]=hcrSurfInds(data);
surfIndsRM=surfRowInds-10;

maskRM=zeros(size(data.range));

for ii=1:size(maskRM,2)
    if data.elevation(ii)<0 & ~isnan(surfIndsRM(ii))
        maskRM(round(surfIndsRM(ii)):end,ii)=1;
    end
end

refl(maskRM==1)=nan;
data.DBZgood=refl;
dataVars{1}='DBZgood';

data.asl=HCRrange2asl(data.range,data.elevation,data.altitude);

%% era5 data

%model.asl=[];
%model.p=[];
model.rh=[];
model.temp=[];
%model.u=[];
%model.v=[];
%model.sst=[];
%model.time=[];

model=read_model(model,directories.modeldir,data.time(1),data.time(end));

%% Find zero degree altitude
tempVec=1:1:size(model.temp,1);
tempMat=repmat(tempVec,size(model.temp,2),1)';

tempMat(model.temp>0)=nan;
tempMat(isnan(model.temp))=nan;

[maxs,rowInds] = nanmax(tempMat);
colInds=1:1:size(model.temp,2);

linearInd = sub2ind(size(model.temp), rowInds, colInds);

[mins,rowIndsUp] = nanmin(tempMat);
colIndsUp=1:1:size(model.temp,2);

linearIndUp = sub2ind(size(model.temp), rowIndsUp, colIndsUp);

linearInd(upInd)=linearIndUp(upInd);
rowInds(upInd)=rowIndsUp(upInd);

checkTemps=model.temp(linearInd);
outInds=zeros(size(checkTemps));
outInds(checkTemps>0 | checkTemps<-2)=1;
outInds(isnan(checkTemps))=1;
linearInd(outInds==1)=[];
zTempsNeg=data.asl(linearInd);

timeMat=repmat(data.time,size(model.temp,1),1);

%% Plot

close all
figs=plot_HCR(data,dataVars,[0 9]);
figure(figs.fig1)
hold on
plot(timeMat(linearInd),zTempsNeg./1000,'k','linewidth',1.5);

ax = gca;
ax.SortMethod = 'childorder';

formatOut = 'yyyymmdd_HHMM';
set(gcf,'PaperPositionMode','auto')
print([directories.figdir,datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_zeroDegree'],'-dpng','-r0');

%% Find bright band
bbAlt=nan(size(data.time));
for ii=1:length(data.time)
    if ~outInds(ii)
        vertCol=data.DBZgood(:,ii);
        vertCol(rowInds(ii)+30:end)=nan;
        vertCol(1:rowInds(ii)-30)=nan;
        
        % Check if all nan
        if min(isnan(vertCol))==1
            continue
        end
        % Check if max refl value is at beginning or end
        vertColData=vertCol(~isnan(vertCol));
        %         if data.time(ii)>datetime(2018,2,7,21,10,0)
        %             figure
        %             plot(vertCol);
        %             test1=1;
        %         end
        if find(vertColData==max(vertColData))==1 | find(vertColData==max(vertColData))==length(vertColData)
            continue
        else
            bbAlt(ii)=find(vertCol==nanmax(vertCol));
        end
    end
end

colIndsBB=1:1:size(model.temp,2);
linearIndBB = sub2ind(size(model.temp), bbAlt, colIndsBB);
linearIndBB(isnan(bbAlt))=[];

plot(timeMat(linearIndBB),data.asl(linearIndBB)./1000,'b','linewidth',1.5);
