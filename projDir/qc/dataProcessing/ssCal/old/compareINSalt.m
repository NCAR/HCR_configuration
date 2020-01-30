% Ocean scan calibration for HCR data

clear all;
close all;

project='cset';

addpath('/h/eol/romatsch/gitPriv/process_HCR/oceanScans/functions/');
addpath('/h/eol/romatsch/gitPriv/utils/');
addpath('/h/eol/romatsch/gitPriv/other/');

directories.figdir='/h/eol/romatsch/hcrCalib/oceanScans/figs/oceanAlt/';
formatOut = 'yyyymmdd_HHMM';

if strcmp(project,'socrates')
    directories.dataDir='/scr/rain1/rsfdata/projects/socrates/hcr/qc/cfradial/moments/10hz/'; %Final
    %directories.dataDir='/scr/rain1/rsfdata/projects/socrates/hcr/cfradial/moments/10hz/'; % raw data
elseif strcmp(project,'cset')
    directories.dataDir='/scr/eldora2/rsfdata/cset/hcr/qc/cfradial/moments/10hz/'; %Final
elseif strcmp(project,'aristo')
    directories.dataDir='/scr/eldora1/rsfdata/aristo-17/hcr/cfradial/moments/10hz/';
end

infile=['/h/eol/romatsch/hcrCalib/oceanScans/biasInFiles/flights_',project,'.txt'];

% Read file with calib events
caseList = table2array(readtable(infile));

altDiffAll=[];

for ii=9:size(caseList,1)
        
    PLT.time = [];
    PLT.rota = [];
    PLT.refl = [];  % reflectivity
    PLT.maxRefl=[];
    PLT.alt  = [];
    PLT.pitch = [];
    PLT.elev = [];  % elevation angle is referenced to horiz plane;
    % pitch and roll corrected.  -80 = 170 or 190 rotation
    %PLT.angdif = []; % abs(rotation - 180) - (90 + elevation) (for downward)
    PLT.roll = [];
    PLT.range = [];
    PLT.lon = [];
    PLT.lat = [];
    PLT.hdg = [];
    PLT.surfInd=[];
    
    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));
    
    oneCaseFileList=makeFileList(directories.dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    % got through all files and load data
    for kk=1:size(oneCaseFileList,2);  % step through each input file in the given set
        
        disp(['Flight ',num2str(ii),' of ',num2str(size(caseList,1)),', file ',num2str(kk),' of ',num2str(size(oneCaseFileList,2))]);
        
        inName=oneCaseFileList{kk};
        toLoc=strfind(inName,'to_');
        indataOrig=inName(1:toLoc);
        
        getInfile=dir([indataOrig '*']);
        
        if size(getInfile,1)>1
            disp('More than one file found. Taking first.');
        end
        indata=[getInfile(1).folder,'/',getInfile(1).name];
        
        % read in time and convert to datetime
        startTimeIn=ncread(indata,'time_coverage_start')';
        startTime=datetime(str2num(startTimeIn(1:4)),str2num(startTimeIn(6:7)),str2num(startTimeIn(9:10)),...
            str2num(startTimeIn(12:13)),str2num(startTimeIn(15:16)),str2num(startTimeIn(18:19)));
        timeRead=ncread(indata,'time')';
        
        timeIn=startTime+seconds(timeRead);
        
        PLT.time = [ PLT.time;  timeIn' ];
        
        %read in data that won't change
        PLT.hdg   = [ PLT.hdg;   ncread(indata,'heading') ];  % hope for no heading near North
        PLT.lat   = [ PLT.lat;   ncread(indata,'latitude') ];
        PLT.lon   = [ PLT.lon;   ncread(indata,'longitude') ];
        
        %read in temporary data
        rotaIn=ncread(indata,'rotation');
        pitchIn=ncread(indata,'pitch');
        rollIn=ncread(indata,'roll');
        elevIn=ncread(indata,'elevation');
        rangeIn=ncread(indata,'range');
        reflIn=ncread(indata,'DBZ');
        altIn=ncread(indata,'altitude');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sort out bad data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             
        %sort out scans where roll angle is too large or small
        outRollInd=find(rollIn<-2 | rollIn>2);
        reflIn(:,outRollInd)=nan;
        
        %sort out scans where rotation angle is too large or small
        outRotaInd=find(rotaIn<178 | rotaIn>182);
        reflIn(:,outRotaInd)=nan;
        
        %sort out scans where elev angle is too large or small
        outElevInd=find(elevIn<-90.5 | elevIn>-89.5);
        reflIn(:,outElevInd)=nan;
        
        % There are high reflectivities in the first gates that need to be removed.
        reflIn(1:20,:)=nan;
                
        [maxRefl maxGate]=nanmax(reflIn,[],1);
        
        %sort out scans with reflectivity values that are too high or low
        if strcmp(project,'cset')
            outRangeInd=find(maxRefl<10 | maxRefl>53);
        else
            outRangeInd=find(maxRefl<35 | maxRefl>53);
        end
        maxRefl(outRangeInd)=nan;
        
        %Get the linear index of the maximum reflectivity value
        maxGateLin=sub2ind(size(reflIn),maxGate,1:size(reflIn,2));
        
        %Check if ocean surface is at right altitude
        rangeMat=repmat(rangeIn,1,size(reflIn,2));
        sfcrng = rangeMat(maxGateLin);
        oceanSurf=sfcrng.*cosd(elevIn');
        
        wrongAltInd=find(nanmax(abs(altIn-oceanSurf'))>60);
        reflIn(:,wrongAltInd)=nan;
        maxRefl(wrongAltInd)=nan;
        
        %Get index where maximum reflectivity is empty and delete data
        emptyMax=find(isnan(maxRefl));
        if ~isempty(emptyMax)
            maxpwr(emptyMax)=nan;
            sfcrng(emptyMax)=nan;
            sfcvel(emptyMax)=nan;
            maxGate(emptyMax)=nan;
        end
        
        PLT.refl = [ PLT.refl, reflIn ];
        PLT.maxRefl = [ PLT.maxRefl; maxRefl' ];
        PLT.range = [ PLT.range; sfcrng' ];
        PLT.elev  = [ PLT.elev; elevIn ];
        PLT.rota = [ PLT.rota; rotaIn ];
        PLT.pitch = [ PLT.pitch; pitchIn ];
        PLT.roll = [ PLT.roll; rollIn ];
        PLT.alt = [ PLT.alt; altIn ];
        PLT.surfInd=[PLT.surfInd;maxGate'];
    end;
    
    %% Plot
    
    %altFromRange=PLT.range.*cosd(PLT.elev+90);
    close all
    
    f1=figure('DefaultAxesFontSize',12);
    set(f1,'Position',[200 500 2000 600]);
    
    hold on
    %plot(altFromRange);
    plot(PLT.time,PLT.alt,'-b','linewidth',1.5);
    plot(PLT.time,PLT.range,'-g','linewidth',1.5);
    ylabel('Meters');
    
    yyaxis right
    rangeMinAlt=PLT.range-PLT.alt;
    if strcmp(project,'cset')
        rangeMinAlt(rangeMinAlt<-100)=nan;
    else
        rangeMinAlt(rangeMinAlt<-10)=nan;
    end
    rangeMinAlt(rangeMinAlt>200)=nan;
    
    if strcmp(project,'socrates') & ii==12
        st=min(find(PLT.time>=datetime(2018,2,18,3,15,0)));
        en=min(find(PLT.time>=datetime(2018,2,18,3,24,0)));
        rangeMinAlt(st:en)=nan;
    end
    
    if strcmp(project,'cset') & ii==9
        st=min(find(PLT.time>=datetime(2015,7,24,23,5,0)));
        en=min(find(PLT.time>=datetime(2015,7,24,23,30,0)));
        rangeMinAlt(st:en)=nan;
    end
    
    if strcmp(project,'cset') & ii==12
        st=min(find(PLT.time>=datetime(2015,8,1,16,30,0)));
        en=min(find(PLT.time>=datetime(2015,8,1,16,40,0)));
        rangeMinAlt(st:en)=nan;
    end
    
    if strcmp(project,'cset') & ii==14
        st=min(find(PLT.time>=datetime(2015,8,7,19,0,0)));
        en=min(find(PLT.time>=datetime(2015,8,7,19,15,0)));
        rangeMinAlt(st:en)=nan;
    end
    
    if strcmp(project,'cset') & ii==16
        st=min(find(PLT.time>=datetime(2015,8,12,19,30,0)));
        en=min(find(PLT.time>=datetime(2015,8,12,19,58,0)));
        rangeMinAlt(st:en)=nan;
    end
    
    plot(PLT.time,rangeMinAlt,'-r','linewidth',1.5);
    legend('Altitude','Range to surf','Range-Alt');
    ylabel('Meters');
        
    title(['RF_',num2str(ii),', altitude of ocean surface'],'interpreter','none');
    
    xlim([PLT.time(1) PLT.time(end)]);
    
    set(f1,'PaperPositionMode','auto')
    print(f1, [directories.figdir,project,'_RF_',num2str(ii),'_oceanSurfAlt'],'-dpng','-r0');
    
    altDiffAll=[altDiffAll,nanmean(rangeMinAlt)];
end

disp(['Means: ', num2str(altDiffAll)]);
disp(['Overall mean: ', num2str(mean(altDiffAll))]);