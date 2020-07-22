% Plots flight and HCR info to be used in determination of surface
% reflectivity.  Determines bias according to Rilling technique.
% Uses a calibration input file, with attenuations pre-computed
% by Romatschke.
%
% Creates multi-panel plots over a selected, contiguous, number of
% data files.  Panels include:
%
%   -- scan rotation angle (rotation is in a plane perpendicular
%                           to flight direction)
%   -- max reflectivity (when pointing down, only)
%   -- range to max reflectivity (when pointing down)
%   -- Altitude
%   -- Pitch
%
% all parameters are plotted vs time, where time is in seconds from
% the start of the first file, or seconds from the start of the day.
%
% Note that for sea surface scans, we operate in vertical receive, only.
function [outTable PLTall]= f_processOceanScans_addNSCAL_dBZ(directories,infile,makeFigs,projName,nsCalFile,refPodTemp)

frq = 9.4406e10;

reflXlims=[0 25];
reflYlims=[15,55];
sigXlims=[0 25];
sigYlims=[-20 15];

% Read file with calib events
caseList = readtable(infile,'Delimiter','space');
%convert to cell so each case has one cell
casesIn=table2array(caseList(:,1));
numCases=unique(casesIn);
uniqueCases=cell(size(numCases,1),1);

for ii=1:size(numCases,1)
    caseInd=find(casesIn==ii);
    uniqueCases{ii}=caseList(caseInd,:);
end

%% Attenuation

attenuation.liebe=nan(size(uniqueCases,1),1);
attenuation.itu=nan(size(uniqueCases,1),1);
attenuation.windspeed=nan(size(uniqueCases,1),1);
attenuation.winddir=nan(size(uniqueCases,1),1);

for ii=1:size(uniqueCases,1); %Loop through all cases
    
    if strcmp(uniqueCases{ii,1}.soundg{1},'nan')
        uniqueCases{ii,1}.soundg{1}=str2num(uniqueCases{ii,1}.soundg{1});
    end
    
    if ~isnan(uniqueCases{ii,1}.soundg{1})
        fileType=uniqueCases{ii,1}.soundg{1}(end-2:end);
        if strcmp(fileType,'eol')
            soundFile=1;
        elseif strcmp(fileType,'.nc')
            soundFile=2;
        else
            disp('Wrong sounding or model file extension.');
        end
        
        if soundFile==1
            % calculate attenuation from sounding
            [attenuation.liebe(ii),attenuation.itu(ii),attenuation.windspeed(ii),attenuation.winddir(ii)] ...
                =f_atten_layers_sfcWnd([directories.sondedir,uniqueCases{ii,1}.soundg{1}],frq/1e+9);
        elseif soundFile==2
            % calculate attenuation from ecmwf model data
            %get lat lon alt
            % make list with files to go through
            startTimeStr=(uniqueCases{ii,1}.timest{:});
            endTimeStr=(uniqueCases{ii,1}.timend{:});
            startTime=datetime(str2num(startTimeStr(1:4)),str2num(startTimeStr(5:6)),str2num(startTimeStr(7:8)),...
                str2num(startTimeStr(10:11)),str2num(startTimeStr(12:13)),str2num(startTimeStr(14:15)));
            endTime=datetime(str2num(endTimeStr(1:4)),str2num(endTimeStr(5:6)),str2num(endTimeStr(7:8)),...
                str2num(endTimeStr(10:11)),str2num(endTimeStr(12:13)),str2num(endTimeStr(14:15)));
            
            oneCaseFileList=makeFileList(directories.dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
            
            indata=oneCaseFileList{1};
            
            lonS=ncread(indata,'longitude');
            latS=ncread(indata,'latitude');
            altS=ncread(indata,'altitude');
            
            %model
            [attenuation.liebe(ii),attenuation.itu(ii),attenuation.windspeed(ii),attenuation.winddir(ii)]= ...
                f_atten_layers_sfcWnd_ecmwf([directories.modeldir,uniqueCases{ii,1}.soundg{1}],...
                frq/1e+9,nanmean(lonS),nanmean(latS),nanmean(altS));
        end
    end
end

%% Load data and calculate bias
PLTall=cell(size(uniqueCases,1),1);
bias.liebe=nan(size(uniqueCases,1),1);
bias.itu=nan(size(uniqueCases,1),1);
bias.stdLiebe=nan(size(uniqueCases,1),1);
bias.N=nan(size(uniqueCases,1),1);
meanHeading=nan(size(uniqueCases,1),1);

for ii=1:size(uniqueCases,1); %Loop through all cases
    disp(['Case ',num2str(ii),' of ',num2str(size(uniqueCases,1))]);
    % Load data
    PLTall{ii}=f_load_sort_data(uniqueCases{ii,1},directories.dataDir);
    % correct reflectivity data with nscal results
    PLTall{ii}=f_reflNSCAL(PLTall{ii},nsCalFile,directories.highResTempDir,refPodTemp);
    % calculate bias
    [PLTall{ii},bias.liebe(ii),bias.itu(ii),bias.stdLiebe(ii),bias.N(ii)]=f_determine_bias_dBZ(PLTall{ii},attenuation.liebe(ii),attenuation.itu(ii),frq);
    
    % calculate mean heading
    headingSmooth=rad2deg(unwrap(deg2rad(PLTall{ii}.hdg)));
    if max(headingSmooth-PLTall{ii}.hdg)>0.1
        disp('Heading corrected');
    end
    meanHeading(ii)=mean(headingSmooth);
end

%% Create table

outTable=cell2table(cell(size(uniqueCases,1)+2,12), ...
    'VariableNames', {'StartDate','EndDate','BiasLiebe','BiasITU','Std','N','OneWayAttLiebe','OneWayAttITU', ...
    'MeanHdg','Sounding','SfcWindSpd','SfcWindDir'});

for ii=1:size(outTable,1)-2
    outTable(ii+2,1)=uniqueCases{ii,1}.timest(1);
    outTable(ii+2,2)=uniqueCases{ii,1}.timend(end);
    soundingTemp=uniqueCases{ii,1}.soundg(1);
    if isnan(soundingTemp{:})
        outTable(ii+2,10)={'nan'};
    else
        outTable(ii+2,10)=soundingTemp;
    end
end

outTable.OneWayAttLiebe=[nan;nan;attenuation.liebe];
outTable.OneWayAttITU=[nan;nan;attenuation.itu];
outTable.BiasLiebe=[nan;nan;bias.liebe];
outTable.BiasITU=[nan;nan;bias.itu];
outTable.Std=[nan;nan;bias.stdLiebe];
outTable.N=[nan;nan;bias.N];
outTable.MeanHdg=[nan;nan;meanHeading];
outTable.SfcWindSpd=[nan;nan;attenuation.windspeed];
outTable.SfcWindDir=[nan;nan;attenuation.winddir];

% add means to table
outTable(1,1)={'Mean'};
outTable(2,1)={'Std'};
% mean
outTable(1,3)={nanmean(bias.liebe)};
outTable(1,4)={nanmean(bias.itu)};
outTable(1,7)={nanmean(attenuation.liebe)};
outTable(1,8)={nanmean(attenuation.itu)};
%std
outTable(2,3)={nanstd(bias.liebe)};
outTable(2,4)={nanstd(bias.itu)};
outTable(2,7)={nanstd(attenuation.liebe)};
outTable(2,8)={nanstd(attenuation.itu)};

outTable(1:2,2)={'nan'};
outTable(1:2,10)={'nan'};

writetable(outTable,[directories.figdir,projName,'_dataTable.dat'],'WriteRowNames',true, 'Delimiter','\t')

%% Plot individual event figures
if makeFigs
    
    % check if directories exist
    if ~exist([directories.figdir 'timeSeries'], 'dir')
        mkdir([directories.figdir 'timeSeries'])
    end
    if ~exist([directories.figdir 'sig0measured'], 'dir')
        mkdir([directories.figdir 'sig0measured'])
    end
    if ~exist([directories.figdir 'reflectivity'], 'dir')
        mkdir([directories.figdir 'reflectivity'])
    end
    
    close all
    for ii=1:size(PLTall,1)
        disp(['Case ',num2str(ii),' of ',num2str(size(uniqueCases,1))]);
        close all
        % time series
        f_plot_series(PLTall{ii},[directories.figdir 'timeSeries/',projName,'_timeSeries_' uniqueCases{ii}.timest{1}]);
        
        % sig0 measured and reflectivity
        PLT=PLTall{ii};
        titleIn=[projName,' measured radar cross section vs elevation angle ' uniqueCases{ii}.timest{1}];
        titleIn2=[projName,' reflectivity vs elevation angle ' uniqueCases{ii}.timest{1}];
        if min(isnan(PLT.sig0measured))==0
            f_plot_sig0noBias(PLT,[directories.figdir 'sig0measured/',projName,'_sig0measured_' uniqueCases{ii}.timest{1}],...
                sigXlims,sigYlims,titleIn);
            f_plot_refl(PLT,[directories.figdir 'reflectivity/',projName,'_refl_' uniqueCases{ii}.timest{1}],reflXlims,reflYlims,titleIn2);
        end
    end
end

% Plot all data combined
plot_all_combined(PLTall,sigXlims,sigYlims,reflXlims,reflYlims,projName,directories,uniqueCases)
end