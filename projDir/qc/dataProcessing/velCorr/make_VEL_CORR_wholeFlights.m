% Compare data from before and after velocity correction

clear all
close all

savefig=0;

addpath(genpath('/h/eol/romatsch/gitPriv/utils/'));

project='cset'; % socrates, cset, aristo, otrec
quality='qc2'; % field, qc1, qc2
freqData='10hz'; % 10hz, 100hz, or 2hz
whichModel='era5';

plotLine=1;
plotSurf=0;

if strcmp(whichModel,'era5')
    directories.modeldir=['/scr/sci/romatsch/data/reanalysis/ecmwf/era5interp/',project,'/',freqData,'/'];
elseif strcmp(whichModel,'ecmwf')
    directories.modeldir=['/scr/sci/romatsch/data/reanalysis/ecmwf/forecastInterp/',project,'/',freqData,'/'];
end

outdir=directories.modeldir;

infile=['/h/eol/romatsch/hcrCalib/oceanScans/biasInFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

figdir=['/h/eol/romatsch/hcrCalib/velCorr/',project,'/velFigs_VEL_CORR/'];
formatOut = 'yyyymmdd_HHMM';

indir=HCRdir(project,quality,freqData);

polyTimePeriod=15; %Time period for poly fit in seconds
polyOrder=3; % Order of polynomial fit

maxEdge=14.6; % Upper edge for plotting
color_map=colormap(vel_default(29));

limits=-4.05:0.3:4.05;
limits=[-inf limits inf];

for kk=1:size(caseList,1)
    
    disp(['Flight ',num2str(kk)]);
    disp('Loading HCR data.')
    
    startTime=datetime(caseList(kk,1:6));
    endTime=datetime(caseList(kk,7:12));
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if ~isempty(fileList)
        data=[];
        
        data.vertical_velocity=[];
        data.nyquist_velocity=[];
        data.DBZ=[];
        data.VEL=[];
        data.TOPO=[];
        
        dataVars=fieldnames(data);
        
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
            
        %% Fill in extinct echo
        
        % Remove up pointing
        dbzDown=data.DBZ;
        dbzDown(:,data.elevation>0)=nan;
        
        velDown=data.VEL;
        velDown(:,data.elevation>0)=nan;
        
        % Get surface indices
        [linInd rowInd rangeToSurf] = hcrSurfInds(data);
        
        surfDBZ=dbzDown(linInd);
        surfDBZlin=10.^(surfDBZ./10);
        
        % SurfVel is the surface velocity that is used for the polinomial
        % fit
        surfVelOrig=velDown(linInd)';
        surfVel=velDown(linInd);
        
        % We remove all the data that we don't want to include in the fit
        
        % Remove all data where we don't see the surface
        surfVel(isnan(surfDBZlin))=nan;
        
        % Remove data that is out of range
        surfVel(isnan(rowInd))=nan;
        
        % Calculate standard deviation
        surfStd=movstd(surfVel,100,'includenan'); % Sets to nan when there is at least one nan
        
        % Enlarge missing data areas
        surfTemp=movmean(surfVel,100,'includenan'); % Sets to nan when there is at least one nan
        surfVel(isnan(surfTemp))=nan;
        
        % standard deviation is too high
        surfVel(surfDBZlin<10000 & surfStd>0.5)=nan;
                      
        % Fill in the missing data with moving average from before the gap
        surfMean=movmedian(surfVel,100,'omitnan'); % Ignores nans and uses the rest of the data
        surfMean(isnan(surfVel))=nan;
        
        aboveGround=data.altitude-data.TOPO;
        
        surfNan=find(isnan(surfVel) & data.elevation<=0 & aboveGround>10);
        
        for ll=1:length(surfNan)
            if ~isnan(surfMean(surfNan(ll)-1)) % At the beginning of the gap we have moving average data
                surfVel(surfNan(ll))=surfMean(surfNan(ll)-1);
            else % Once the moving average turns nan we just keep the previous one going
                surfVel(surfNan(ll))=surfVel(surfNan(ll)-1);
            end
        end
        %% Make poly fit
        velSmoothPorig=vel2vel_corr_testPoly(surfVel,data.time,polyTimePeriod,polyOrder);
        
        % Fill in data where we don't have polyfit data
        velSmoothP=movmedian(surfVelOrig,100,'omitnan');
        velSmoothP(~isnan(velSmoothPorig))=velSmoothPorig(~isnan(velSmoothPorig));
        
        % Vel corr
        velCorrP=data.VEL-velSmoothP';
        velCorrSurfP=velCorrP(linInd);
        
        %% Make output
        velCorrOut=data.VEL;
        velCorrOut(:,data.elevation<=0)=velCorrP(:,data.elevation<=0);
        
        %% Save
        disp('Saving velCorr field.')
        
        velCorr=velCorrOut;
        save([outdir,whichModel,'.velCorr.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
            datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(kk),'.mat'],'velCorr');
        
    end
end