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
function PLTall= f_processOceanScans_simple(directories,infile,makeFigs,projName,...
    quality,nsCalFile,attenuationFrom,...
    reflXlims,reflYlims,sigXlims,sigYlims,addName,sig0model,oceanTemp,salinity)

% Read file with calib events
uniqueCases = table2array(readtable(infile,'Delimiter','space'));
startCase=1;
uniqueCases=uniqueCases(startCase:end,:);
numCases=size(uniqueCases,1);
%numCases=1; % This can be set for testing if we don't want to run all cases

% Add last row if not in file
if size(uniqueCases,2)<14
    uniqueCases=cat(2,uniqueCases,ones(size(uniqueCases,1),1));
end

%% Load data
PLTall=cell(numCases,1);
for ii=1:numCases; %Loop through all cases
    disp(['Loading case ',num2str(ii),' of ',num2str(numCases)]);
    % Load data
    [PLTall{ii},frq]=f_load_sort_data(uniqueCases(ii,:),directories.dataDir);
end

%% Attenuation

attenuation.liebe={};
attenuation.itu={};
attenuation.windspeed={};
attenuation.u={};
attenuation.v={};
%eraInds=cell(numCases,1);
SST=cell(numCases,1);

for ii=1:numCases; %Loop through all cases
    disp(['Attenuation for case ',num2str(ii),' of ',num2str(numCases)]);
    soundFile=uniqueCases(ii,13);
    
    if soundFile==0
        attenuation.liebe{end+1}=nan;
        attenuation.itu{end+1}=nan;
        attenuation.windspeed{end+1}=nan;
        attenuation.u{end+1}=nan;
        attenuation.v{end+1}=nan;
    elseif soundFile==1
        % calculate attenuation from sounding
        [attenuation.liebe{end+1},attenuation.itu{end+1},attenuation.windspeed{end+1},attenuation.u{end+1},attenuation.v{end+1}] ...
            =f_atten_layers_sfcWnd_dropsonde(directories.sondedir,frq/1e+9,PLTall{ii,1}.time,PLTall{ii,1}.alt);
    elseif soundFile==2
        if ~isempty(strfind(attenuationFrom,'era5')) | ~isempty(strfind(attenuationFrom,'ecmwf'))
            % calculate attenuation from model data
            [attenuation.liebe{end+1},attenuation.itu{end+1},attenuation.windspeed{end+1},attenuation.u{end+1},attenuation.v{end+1},SST{ii}]= ...
                f_atten_layers_modelInterp(directories.modeldir,frq/1e+9,PLTall{ii,1}.time,1);
        end
    elseif soundFile==3
        if ~isempty(strfind(attenuationFrom,'file'))
            % calculate attenuation from ecmwf model data
            [attenuation.liebe{end+1},attenuation.itu{end+1},attenuation.windspeed{end+1},attenuation.u{end+1},attenuation.v{end+1},SST{ii}]= ...
                f_atten_layers_file(frq/1e+9,PLTall{ii,1}.time,1,PLTall{ii});
        end
    end
end

%% Calculate bias

bias.liebe=nan(numCases,1);
bias.itu=nan(numCases,1);
bias.stdLiebe=nan(numCases,1);
bias.N=nan(numCases,1);
meanHeading=nan(numCases,1);

for ii=1:numCases; %Loop through all cases
    disp(['Calculating bias for case ',num2str(ii),' of ',num2str(numCases)]);
    if ~isnan(nsCalFile)
        % correct reflectivity data with nscal results
        disp('Correcting for NSCAL results.');
        PLTall{ii}=f_reflNSCAL_simple(PLTall{ii},nsCalFile,directories.highResTempDir);
    end
    % calculate bias
    [PLTall{ii},bias.liebe(ii),bias.itu(ii),bias.stdLiebe(ii),bias.N(ii)]=f_determine_bias_simple(PLTall{ii},attenuation.liebe{1,ii},attenuation.itu{1,ii},frq);
    
    % calculate mean heading
    headingSmooth=rad2deg(unwrap(deg2rad(PLTall{ii}.hdg)));
    if max(headingSmooth-PLTall{ii}.hdg)>0.1
        disp('Heading corrected');
    end
    meanHeading(ii)=mean(headingSmooth);
    
    %calculate sigma0 model
    if sig0model
        if isempty(SST{ii})
            disp(['No sea surface temperature data found. Proceeding with default ',num2str(oceanTemp),' C.']);
            SST{ii}=nan(size(PLTall{ii}.time));
            SST{ii}(:,:)=oceanTemp;
        elseif min(isnan(SST{ii}))==1
            disp(['No sea surface temperature data found. Proceeding with default ',num2str(oceanTemp),' C.']);
            SST{ii}(:,:)=oceanTemp;
        end
        PLTall{ii}=f_sigma0_model(PLTall{ii},attenuation.windspeed(ii),frq,SST{ii},salinity);
    end
end

% %% Create table
% outTable=cell2table(cell(numCases+2,13), ...
%     'VariableNames', {'StartDate','EndDate','BiasLiebe','BiasITU','Std','N','OneWayAttLiebe','OneWayAttITU', ...
%     'MeanHdg','Sounding','SfcWindSpd','SfcWindDir','SST'});
% attLiebe=nan(numCases,1);
% attITU=nan(numCases,1);
meanWindDir=nan(numCases,1);
% 
% for ii=1:size(outTable,1)-2
%     outTable(ii+2,1)={datestr(uniqueCases(ii,1:6),'yyyymmdd_HHMMSS')};
%     outTable(ii+2,2)={datestr(uniqueCases(ii,7:12),'yyyymmdd_HHMMSS')};
%     if uniqueCases(ii,13)==0
%         outTable(ii+2,10)={'nan'};
%     elseif uniqueCases(ii,13)==1
%         outTable(ii+2,10)={'dropSonde'};
%     else
%         outTable(ii+2,10)={'model'};
%     end
%     attLiebe(ii)=mean(attenuation.liebe{1,ii},'omitnan');
%     attITU(ii)=mean(attenuation.itu{1,ii},'omitnan');
%     outTable(ii+2,11)={num2str(mean(attenuation.windspeed{1,ii},'omitnan'))};
%     try
%         meanWindDir(ii)=wrapTo360(270-atan2(mean(attenuation.v{1,ii},'omitnan'),mean(attenuation.u{1,ii},'omitnan'))*180/pi);
%     catch
%         meanWindDir(ii)=nan;
%     end
%     outTable(ii+2,12)={num2str(meanWindDir(ii))};
%     outTable(ii+2,13)={num2str(mean(SST{ii},'omitnan'))};
% end
% 
% outTable.BiasLiebe=[nan;nan;bias.liebe];
% outTable.BiasITU=[nan;nan;bias.itu];
% outTable.OneWayAttLiebe=[nan;nan;attLiebe];
% outTable.OneWayAttITU=[nan;nan;attITU];
% outTable.Std=[nan;nan;bias.stdLiebe];
% outTable.N=[nan;nan;bias.N];
% outTable.MeanHdg=[nan;nan;meanHeading];
% 
% % add means to table
% outTable(1,1)={'Mean'};
% outTable(2,1)={'Std'};
% % mean
% outTable(1,3)={mean(bias.liebe,'omitnan')};
% outTable(1,4)={mean(bias.itu,'omitnan')};
% outTable(1,7)={mean(attLiebe,'omitnan')};
% outTable(1,8)={mean(attITU,'omitnan')};
% %std
% outTable(2,3)={std(bias.liebe,'omitnan')};
% outTable(2,4)={std(bias.itu,'omitnan')};
% outTable(2,7)={std(attLiebe,'omitnan')};
% outTable(2,8)={std(attITU,'omitnan')};
% 
% outTable(1:2,2)={'nan'};
% outTable(1:2,10:12)={'nan'};
% 
% if ~isnan(nsCalFile)
%     outName=[directories.figdir,projName,'_',quality,'_',attenuationFrom,'_addNScal_dataTable',addName,'.dat'];
% else
%     outName=[directories.figdir,projName,'_',quality,'_',attenuationFrom,'_dataTable',addName,'.dat'];
% end
% 
% writetable(outTable,outName,'WriteRowNames',true, 'Delimiter','\t')

%% Make poly fits to data and model
if sig0model
    for ii=1:size(PLTall,1)
        disp(['Poly fit for case ',num2str(ii),' of ',num2str(numCases)]);
        PLT=PLTall{ii};
        if min(isnan(PLT.sig0measured))==0 & min(min(isnan(PLT.sig0model)))==0
            PLTall{ii}=f_polySig0(PLTall{ii});
        end
    end
end
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
        disp(['Plotting case ',num2str(ii),' of ',num2str(numCases)]);
        close all
        
        outName=[projName,'_',quality];
        
        if ~isnan(nsCalFile)
            outName=[outName,'_addNScal'];
        end
        if isfield(PLTall{ii},'sig0model') & max(~isnan(PLTall{ii}.sig0model))>0
            outName=[outName,'_sig0model'];
        end 
        
        outName=[outName,addName];
        
        % time series
        f_plot_series_simple(PLTall{ii},[directories.figdir 'timeSeries/',outName,'_timeSeries_',datestr(uniqueCases(ii,1:6),'yyyymmdd_HHMMSS')]);
        
        % sig0 measured and reflectivity
        PLT=PLTall{ii};        
        titleIn=[{outName};{['Measured radar cross section vs elevation angle ',datestr(uniqueCases(ii,1:6),'yyyymmdd_HHMMSS')]}];
        titleIn2=[{outName};{['Reflectivity vs elevation angle ',datestr(uniqueCases(ii,1:6),'yyyymmdd_HHMMSS')]}];
        if min(isnan(PLT.sig0measured))==0
            try
                headingFrom=wrapTo360(PLT.hdg-180);
                planeU=mean(-sin(headingFrom * pi/180),'omitnan');
                planeV=mean(-cos(headingFrom * pi/180),'omitnan');
                hdgMinus180=wrapTo360(270-atan2(planeV,planeU)*180/pi);
            catch
                hdgMinus180=nan;
            end
            try
                meanWindDir=wrapTo360(270-atan2(mean(attenuation.v{1,ii},'omitnan'),mean(attenuation.u{1,ii},'omitnan'))*180/pi);
            catch
                meanWindDir=nan;
            end
            f_plot_sig0noBias(PLT,[directories.figdir 'sig0measured/',outName,'_sig0measured_',datestr(uniqueCases(ii,1:6),'yyyymmdd_HHMMSS')],...
                sigXlims,sigYlims,titleIn,mean(attenuation.windspeed{1,ii},'omitnan'),mean(SST{ii},'omitnan'),meanWindDir,hdgMinus180,...
                mean(attenuation.liebe{1,ii},'omitnan'),uniqueCases(ii,14));
            f_plot_refl(PLT,[directories.figdir 'reflectivity/',outName,'_refl_',datestr(uniqueCases(ii,1:6),'yyyymmdd_HHMMSS')],...
                reflXlims,reflYlims,titleIn2);
        end
    end
end

end