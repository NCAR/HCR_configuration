% Find maximum reflectivity in drizzle at 250 m range and temperatures > 2C
% Only when pointing down
% In a calibrated radar these values should be ~18-20 dBZ

clear all;
close all;

project='otrec';

addpath('/h/eol/romatsch/gitPriv/process_HCR/oceanScans/functions/');
addpath('/h/eol/romatsch/gitPriv/utils/');
addpath('/h/eol/romatsch/gitPriv/other/');

directories.figdir='/h/eol/romatsch/hcrCalib/otherCalib/drizzle/';
formatOut = 'yyyymmdd_HHMM';


if strcmp(project,'socrates')
    directories.dataDir='/scr/snow2/rsfdata/projects/socrates/hcr/qc2/cfradial/moments/10hz/'; %Final
elseif strcmp(project,'cset')
    directories.dataDir='/scr/snow2/rsfdata/projects/cset/hcr/qc2/cfradial/moments/10hz/'; %Final
elseif strcmp(project,'aristo')
    directories.dataDir='/scr/eldora1/rsfdata/aristo-17/hcr/cfradial/moments/10hz/';
elseif strcmp(project,'otrec')
    directories.dataDir='/scr/snow1/rsfdata/projects/otrec/hcr/qc1/cfradial/finalERA5/10hz/';
end

infile=['/h/eol/romatsch/hcrCalib/otherCalib/inFiles/drizzle_',project,'.txt'];

% Data files with plane data
if strcmp(project,'otrec')
    directories.planeDir=['/scr/snow1/rsfdata/projects/otrec/GV/'];
    planeFiles=dir([directories.planeDir,'*rf*.nc']);
else
    directories.planeDir=['/scr/sci/romatsch/data/hiaper/',project,'/'];
    planeFiles=dir([directories.planeDir,'RF*.nc']);
end

% Read file with calib events
caseList = table2array(readtable(infile));

reflRangeAll=[];
reflRangeMedian=[];

for ii=1:size(caseList,1)
    
    PLT.time = [];
    PLT.refl = [];  % reflectivity
    PLT.elev = [];  % elevation angle is referenced to horiz plane;
    PLT.range = [];
    
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
        
        elevIn=ncread(indata,'elevation');
        rangeIn=ncread(indata,'range');
        reflIn=ncread(indata,'DBZ');
        rangeMat=repmat(rangeIn,1,size(reflIn,2));
        
        PLT.refl = [ PLT.refl, reflIn ];
        PLT.range = [ PLT.range, rangeMat ];
        PLT.elev  = [ PLT.elev; elevIn ];
    end
    
    %% Get INS data from plane
    planeFileIn=[planeFiles(ii).folder,'/',planeFiles(ii).name];
        
    temp1=ncread(planeFileIn,'ATH1');
    temp2=ncread(planeFileIn,'ATH2');
    
    if min(size(temp1))>1
        temp1=temp1(1,:)';
        temp2=temp2(1,:)';
    end
    
    planeTempOrig=nanmean(cat(2,temp1,temp2),2);
    
    planeTimeOrig=ncread(planeFileIn,'Time');
    
    info=ncinfo(planeFileIn);
    planeRefTimeIn=info.Variables(1).Attributes(3).Value;
    planeRefTime=datetime(str2num(planeRefTimeIn(15:18)),str2num(planeRefTimeIn(20:21)),str2num(planeRefTimeIn(23:24)),...
        str2num(planeRefTimeIn(26:27)),str2num(planeRefTimeIn(29:30)),str2num(planeRefTimeIn(32:33)));
    
    planeTime=planeRefTime+seconds(planeTimeOrig);
    
    % Resample temperature to hcr resolution
    
    TThcr=timetable(PLT.time,PLT.elev);
    TTplane=timetable(planeTime,planeTempOrig);
    
    % Sometimes times are not in order
    [TThcr sortIndHCR]= sortrows(TThcr);
    indDiff=sortIndHCR(2:end)-sortIndHCR(1:end-1);
    indRM=find(indDiff~=1)+1;
    
    TThcr(indRM,:)=[];
    PLT.time(indRM)=[];
    PLT.refl(:,indRM)=[];
    PLT.range(:,indRM)=[];
    PLT.elev(indRM)=[];
    
    TTplane = sortrows(TTplane);
    
    TTcombined=synchronize(TThcr,TTplane,'first','linear');
    
    planeTemp=TTcombined.planeTempOrig;
    
    %% Use right data
      
    % Sort out data with too low temperatures
    lowTempInd=find(planeTemp<5);
        
    range=PLT.range;
    reflTemp=PLT.refl;
    
    %Fill in bang
    reflTemp(1:17,:)=1;
    
    reflTemp(:,lowTempInd)=nan;
    reflTemp(range>230)=nan;
    
    for jj=1:size(reflTemp,2)
        reflRay=reflTemp(:,jj);
        nanRay=find(isnan(reflRay));
        nonNanRay=find(~isnan(reflRay));
        if min(nanRay)<max(nonNanRay)
            reflTemp(:,jj)=nan;
        end
    end
    
    refl=reflTemp;
    refl2=reflTemp;
    
    % Get right range
    refl(range<230 | range>270)=nan;
    refl2(range<200 | range>300)=nan;
    
    usedInd=any(~isnan(refl),1);
    refl2(:,usedInd==0)=nan;
    
    reflInds=find(~isnan(refl));
    reflRangeAll=[reflRangeAll;refl(reflInds)];
        
    med=nanmedian(refl2,1)';
    reflRangeMedian=[reflRangeMedian;med];
    
    close all
    
    f0=figure('DefaultAxesFontSize',12);
    set(f0,'Position',[200 500 1000 600]);
    
    hold on
    plot(PLT.time,planeTemp,'-b','linewidth',3);
    plot(PLT.time(usedInd==1),planeTemp(usedInd==1),'or');
    ylabel('Temperatrue (C)');
    xlim([PLT.time(1),PLT.time(end)]);
    
    title([project,' flight ',num2str(ii)]);
    
    set(f0,'PaperPositionMode','auto')
    print(f0, [directories.figdir,project,'_dataUsed_flight',num2str(ii)],'-dpng','-r0');
end

%% Plot histogram

edgesAll=[-70:5:40];
edgesTop=[10:1:40];

close all

f1=figure('DefaultAxesFontSize',12);
set(f1,'Position',[200 500 1000 1000]);

subplot(2,1,1)

histogram(reflRangeAll,edgesAll)
title([project,' histogram of reflectivities at 250 m range'],'interpreter','none');
xlim([-60 35])
text(0,60000,['Total points: ', num2str(length(reflRangeAll))],'fontsize',14);
xlabel('Reflectivity (dBZ)');

subplot(2,1,2)

histogram(reflRangeAll,edgesTop)
title('top end');
xlim([10 35])
ylim([0 800])
xlabel('Reflectivity (dBZ)');

set(f1,'PaperPositionMode','auto')
print(f1, [directories.figdir,project,'_drizzleRefl'],'-dpng','-r0');

f2=figure('DefaultAxesFontSize',12);
set(f2,'Position',[200 500 1000 1000]);

subplot(2,1,1)

histogram(reflRangeMedian,edgesAll)
title([project,' histogram of median of reflectivities at 250 m range'],'interpreter','none');
xlim([-60 35])
text(0,40000,['Total points: ', num2str(length(find(~isnan(reflRangeMedian))))],'fontsize',14);
xlabel('Reflectivity (dBZ)');

subplot(2,1,2)

histogram(reflRangeMedian,edgesTop)
title('top end');
xlim([10 35])
ylim([0 800])
xlabel('Reflectivity (dBZ)');

set(f2,'PaperPositionMode','auto')
print(f2, [directories.figdir,project,'_drizzleRefl_median'],'-dpng','-r0');

