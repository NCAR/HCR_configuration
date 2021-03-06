% Ocean scan calibration for HCR data

clear all;
close all;


addpath('/h/eol/romatsch/gitPriv/process_HCR/oceanScans/functions/');
addpath('/h/eol/romatsch/gitPriv/process_HCR/oceanScans/colormaps/');
addpath('/h/eol/romatsch/gitPriv/utils/');
addpath('/h/eol/romatsch/gitPriv/other/');

directories.figdir='/h/eol/romatsch/hcrCalib/oceanScans/figs/oceanDBZ/';
formatOut = 'yyyymmdd_HHMM';

%directories.dataDir='/scr/rain1/rsfdata/projects/socrates/hcr/qc/cfradial/moments/10hz/'; %Final
directories.dataDir='/scr/rain1/rsfdata/projects/socrates/hcr/cfradial/moments/10hz/'; % raw data

% Fast scans
infile='/h/eol/romatsch/hcrCalib/oceanScans/biasInFiles/cal_SOCRATES_fast.txt';
projectNameF='fastSOCRATES';

frq = 9.4406e10;

% Read file with calib events
caseList = readtable(infile,'Delimiter','space');
%convert to cell so each case has one cell
casesIn=table2array(caseList(:,1));
numCases=unique(casesIn);
uniqueCasesIn=cell(size(numCases,1),1);

for ii=1:size(numCases,1)
    caseInd=find(casesIn==ii);
    uniqueCasesIn{ii}=caseList(caseInd,:);
end

tenDegDiff=[];

for ii=1:size(uniqueCasesIn,1)
    
    uniqueCases=uniqueCasesIn{ii};
    
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
    
    % make list with files to go through
    startTimeStr=(uniqueCases.timest{:});
    endTimeStr=(uniqueCases.timend{:});
    startTime=datetime(str2num(startTimeStr(1:4)),str2num(startTimeStr(5:6)),str2num(startTimeStr(7:8)),...
        str2num(startTimeStr(10:11)),str2num(startTimeStr(12:13)),str2num(startTimeStr(14:15)));
    endTime=datetime(str2num(endTimeStr(1:4)),str2num(endTimeStr(5:6)),str2num(endTimeStr(7:8)),...
        str2num(endTimeStr(10:11)),str2num(endTimeStr(12:13)),str2num(endTimeStr(14:15)));
    
    oneCaseFileList=makeFileList(directories.dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    % got through all files and load data
    for kk=1:size(oneCaseFileList,2);  % step through each input file in the given set
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
        
        % Scans should only have reflectivity data at the very top and at the ocean surface.
        % We sort out scans with too many reflectivity values
        dataSize=sum(~isnan(reflIn),1);
        tooMuchInd=find(dataSize>50);
        
        reflIn(:,tooMuchInd)=nan;
        
        %sort out scans where roll angle is too large or small
        outRollInd=find(rollIn<-2 | rollIn>2);
        reflIn(:,outRollInd)=nan;
        
        % There are high reflectivities in the first gates that need to be removed.
        badRow=8;
        empty=false;
        
        % find first row that has less than 1% nans
        while (badRow<40 && ~empty)
            if length(find(~isnan(reflIn(badRow,:))))<size(reflIn,2)/100 | length(find(~isnan(reflIn(badRow,:))))<3;
                empty=true;
            end
            badRow=badRow+1;
        end;
        
        if badRow==40
            disp(['First row with only nans is above 40. Check data!']);
            reflIn(:,:)=nan;
        else
            reflIn(1:badRow-1,:)=nan;
        end
        
        [maxRefl maxGate]=nanmax(reflIn,[],1);
        
        %sort out scans with reflectivity values that are too high or low
        outRangeInd=find(maxRefl<0 | maxRefl>55);
        maxRefl(outRangeInd)=nan;
        
        %Get the linear index of the maximum reflectivity value
        maxGateLin=sub2ind(size(reflIn),maxGate,1:size(reflIn,2));
        
        %Check if ocean surface is at right altitude
        rangeMat=repmat(rangeIn,1,size(reflIn,2));
        sfcrng = rangeMat(maxGateLin);
        oceanSurf=sfcrng.*cosd(elevIn');
        
        wrongAltInd=find(nanmax(abs(altIn-oceanSurf'))>100);
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
    
    timeInds=[];
    for ll=1:size(uniqueCases,1)
        % take only data in right time span
        begTime=datetime(str2num(uniqueCases.timest{ll}(1:4)),str2num(uniqueCases.timest{ll}(5:6)),str2num(uniqueCases.timest{ll}(7:8)),...
            str2num(uniqueCases.timest{ll}(10:11)),str2num(uniqueCases.timest{ll}(12:13)),str2num(uniqueCases.timest{ll}(14:15)));
        endTime=datetime(str2num(uniqueCases.timend{ll}(1:4)),str2num(uniqueCases.timend{ll}(5:6)),str2num(uniqueCases.timend{ll}(7:8)),...
            str2num(uniqueCases.timend{ll}(10:11)),str2num(uniqueCases.timend{ll}(12:13)),str2num(uniqueCases.timend{ll}(14:15)));
        
        timeInds=[timeInds; find(PLT.time>begTime+1/24/60/60 & PLT.time<endTime)];
    end
    
    fieldsPLT=fields(PLT);
    
    for jj=1:length(fieldsPLT)
        if ~isempty(PLT.(fieldsPLT{jj})) & size(PLT.(fieldsPLT{jj}),2)==1
            PLT.(fieldsPLT{jj})=PLT.(fieldsPLT{jj})(timeInds);
        end
    end
    
    PLT.refl=PLT.refl(:,timeInds);
    PLT.refl=PLT.refl';
        
    %% Plot ocean surface gate vs ocean surface reflectivity
    oceanRefl=nan(length(PLT.time),25);
       
    if max(max(~isnan(PLT.refl)))>0
        
        for ii=1:length(PLT.time)
            if ~isnan(PLT.surfInd(ii));
                oceanRefl(ii,:)=PLT.refl(ii,PLT.surfInd(ii)-12:PLT.surfInd(ii)+12);               
            end
        end
        
        oceanRefl=oceanRefl';
        
        threeGates=oceanRefl(12:14,:);
        threeGatesLin=10.^(threeGates./10);
        threeSumLin=nansum(threeGatesLin,1);
        threeSum=10 * log10(threeSumLin);
        
        fiveGates=oceanRefl(11:15,:);
        fiveGatesLin=10.^(fiveGates./10);
        fiveSumLin=nansum(fiveGatesLin,1);
        fiveSum=10 * log10(fiveSumLin);
        
        sevenGates=oceanRefl(10:16,:);
        sevenGatesLin=10.^(sevenGates./10);
        sevenSumLin=nansum(sevenGatesLin,1);
        sevenSum=10 * log10(sevenSumLin);
        
        diff=threeSum-PLT.maxRefl';
        meanDiff=nanmean(diff);
        stdDiff=nanstd(threeSum-PLT.maxRefl');
        
        tenDegInd=find(abs(PLT.elev+80<0.3));
        
        tenDegDiff=[tenDegDiff,diff(tenDegInd)];
        
        close all
        
        f1=figure('DefaultAxesFontSize',12);
        set(f1,'Position',[200 500 2000 1200]);
        
        subplot(2,1,1)
        hold on
        plot(PLT.time,PLT.maxRefl,'linewidth',1.5);
        plot(PLT.time,threeSum,'linewidth',1.5);
        plot(PLT.time,fiveSum,'linewidth',1.5);
        plot(PLT.time,sevenSum,'linewidth',1.5);
        ylabel('dBZ');
        yyaxis right
        plot(PLT.time,PLT.elev);
        %plot(PLT.time,PLT.surfInd);
        xlim([PLT.time(1) PLT.time(end)]);
        legend('Refl 1 gate','Refl 3 gates','Refl 5 gates','Refl 7 gates','Elev angle');
        ylabel('Degree');
        ylimits=ylim;
        text(PLT.time(5),ylimits(1)+1,['Mean diff: ',num2str(meanDiff),' dBZ. Std: ',num2str(stdDiff),' dBZ.'],'fontsize',15);
        title('Ocean surface reflectivity');
        
        grid on
        
        subplot(2,1,2)
        fig2=surf(PLT.time,1:1:25,oceanRefl);
        fig2.EdgeColor='none';
        xlim([PLT.time(1) PLT.time(end)]);
        ylim([5 25]);
        view(2);
        colorbar('southoutside');
        title('Reflectivity (dBZ) close to ocean surface')
        ylabel('Gates')
        
        grid on
        
        set(f1,'PaperPositionMode','auto')
        print(f1, [directories.figdir,datestr(PLT.time(1),formatOut),'_to_',datestr(PLT.time(end),formatOut),'_oceanDBZ'],'-dpng','-r0');
       
    end
end

tenDegMean=nanmean(tenDegDiff);

%% Detailed analysis of one case

f2=figure('DefaultAxesFontSize',12);
set(f2,'Position',[200 500 2000 800]);

startInd=90;
endInd=120;
hold on
plot(PLT.time(startInd:endInd),PLT.maxRefl(startInd:endInd),'or','MarkerFaceColor','r');
plot(PLT.time(startInd:endInd),threeSum(startInd:endInd),'om','MarkerFaceColor','m');
plot(PLT.time(startInd:endInd),fiveSum(startInd:endInd),'og','MarkerFaceColor','g');

yyaxis right
plot(PLT.time(startInd:endInd),PLT.surfInd(startInd:endInd),'ob','MarkerFaceColor','b');

xlim([PLT.time(startInd) PLT.time(endInd)])

%ax = gca;
set(gca,'XMinorGrid','on');