% Calculate pulse pair moments from HCR time series
% Author: Ulrike Romatschke

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User input
project='socrates'; %socrates, cset, otrec, noreaster, spicule
quality='ts';
freqData=[];
qcVersion=[];

outTime=0.1; % Desired output time resolution in seconds. Must be less than or equal to one second.
sampleTime=0.1; % Length of sample in seconds.

showPlot='on';

startTime=datetime(2018,2,24,4,33,0);
endTime=datetime(2018,2,24,4,33,40);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Process
dataDirTS=HCRdir(project,quality,qcVersion,freqData);
fileListTS=makeFileList(dataDirTS,startTime+seconds(1),endTime-seconds(1),'20YYMMDDxhhmmss',1);

for bb=1:length(fileListTS)

    %% Load time series data
    if bb==1
        data=[];
        data.IVc=[];
        data.QVc=[];

        disp(['Loading time series file ',num2str(bb),' of ' num2str(length(fileListTS)),' ...']);

        data=read_TsArchive_iwrf_bulk(fileListTS{1},data);

        % Find available times
        timeTest=data.time';
        timeTest(data.time<startTime | data.time>endTime)=[];
    else
        %% Trimm first file
        data=trimFirstFile(data,endInd);

        goodTimes(1:end-1)=[];

        disp(['Loading time series file ',num2str(bb),' of ' num2str(length(fileListTS)),' ...']);

        data=read_TsArchive_iwrf_bulk(fileListTS{bb},data);

        % Find available times
        timeTest=data.time';
        timeTest(data.time<goodTimes | data.time>endTime)=[];
    end

    % Find indices in time interval
    TTdata=timetable(timeTest,ones(size(timeTest)));
    synchTT=retime(TTdata,'regular','sum','TimeStep',seconds(outTime));
    goodTimes=synchTT.timeTest(synchTT.Var1>0);

    if bb<length(fileListTS)
        beamNum=length(goodTimes)-1;
    else
        beamNum=length(goodTimes);
    end

    %% Process time series

    disp('Processing time series ...');

    momentsTimeOne.powerV=nan(size(data.range,1),beamNum);
    momentsTimeOne.velRaw=nan(size(data.range,1),beamNum);
    momentsTimeOne.width=nan(size(data.range,1),beamNum);
    momentsTimeOne.dbz=nan(size(data.range,1),beamNum);
    momentsTimeOne.snr=nan(size(data.range,1),beamNum);
    momentsTimeOne.ncp=nan(size(data.range,1),beamNum);
    momentsTimeOne.range=nan(size(data.range,1),beamNum);
    momentsTimeOne.altitude=nan(1,beamNum);
    momentsTimeOne.elevation=nan(1,beamNum);
    momentsTimeOne.azimuth_vc=nan(1,beamNum);
    momentsTimeOne.time=goodTimes(1:beamNum)';

    % Loop through beams
    for ii=1:beamNum

        % Find start and end indices for beam
        [~,startInd]=min(abs(goodTimes(ii)-seconds(sampleTime/2)-data.time));
        [~,endInd]=min(abs(goodTimes(ii)+seconds(sampleTime/2)-data.time));

        sampleNum=endInd-startInd+1;

        if sampleNum<30
            continue
        end

        % Window
        win=window(@hamming,sampleNum);  % Default window is Hamming
        winWeight=sampleNum/sum(win);
        winNorm=win*winWeight;

        % Trim data down to current beam
        dataThis=trimData(data,startInd,endInd);

        %% Create IQ
        cIQ.v=winNorm'.*(dataThis.IVc+i*dataThis.QVc)./sqrt(sampleNum);

        %% Other variables
        momentsTimeOne.range(:,ii)=dataThis.range;
        momentsTimeOne.elevation(ii)=median(dataThis.elevation);
        momentsTimeOne.azimuth_vc(ii)=median(dataThis.azimuth_vc);
        momentsTimeOne.altitude(ii)=median(dataThis.altitude);

        %% Calculate pulse pair moments
        momentsTimeOne=calcMomentsTime(cIQ,ii,momentsTimeOne,dataThis);
    end

    %% Add to output
    if bb==1
        momentsTime=momentsTimeOne;
    else
        dataFields=fields(momentsTime);

        for hh=1:length(dataFields)
            momentsTime.(dataFields{hh})=cat(2,momentsTime.(dataFields{hh}),momentsTimeOne.(dataFields{hh}));
        end

    end

end

%% Plot

close all

disp('Plotting ...');

momentsTime.asl=HCRrange2asl(momentsTime.range,momentsTime.elevation,momentsTime.altitude);

ylims=[0,2.5];

climsPow=[-105,-35];
climsDbz=[-40,30];
climsVel=[-10,10];
climsWidth=[0,3];
colMap=jet;

% Figure
f1 = figure('Position',[200 500 1600 830],'DefaultAxesFontSize',12,'visible',showPlot);

t = tiledlayout(2,2,'TileSpacing','tight','Padding','tight');

s1=nexttile(1);

hold on
surf(momentsTime.time,momentsTime.asl./1000,momentsTime.powerV,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsPow);
s1.Colormap=colMap;
colorbar
grid on
box on
title('Time domain DBMV (dBm)')
ylim(ylims);
xlim([momentsTime.time(1),momentsTime.time(end)]);

s2=nexttile(2);

hold on
surf(momentsTime.time,momentsTime.asl./1000,momentsTime.dbz,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsDbz);
s2.Colormap=colMap;
colorbar
grid on
box on
title('Time domain DBZ (dB)')
ylim(ylims);
xlim([momentsTime.time(1),momentsTime.time(end)]);

s3=nexttile(3);

hold on
surf(momentsTime.time,momentsTime.asl./1000,momentsTime.velRaw,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsVel);
s3.Colormap=colMap;
colorbar
grid on
box on
title('Time domain velocity (m s^{-1})')
ylim(ylims);
xlim([momentsTime.time(1),momentsTime.time(end)]);

s4=nexttile(4);

hold on
surf(momentsTime.time,momentsTime.asl./1000,momentsTime.width,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsWidth);
s4.Colormap=colMap;
colorbar
grid on
box on
title('Spectrum width (m s^{-1})')
ylim(ylims);
xlim([momentsTime.time(1),momentsTime.time(end)]);