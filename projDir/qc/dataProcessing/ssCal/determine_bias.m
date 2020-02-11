% Ocean scan calibration for HCR data

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

% If 1, plots for individual calibration events will be made, if 0, only
% total plots will be made
makeSingleFigs=1;
addNSCAL=0; % Noise source cal will be incorporated
addSIG0model=1; % Plot also model data from sig0 model

project='otrec'; %socrates, aristo, cset, otrec
quality='qc1'; %field, qc1, or qc2
dataFreq='10hz';
addName=''; % Extra name part for output files. Default is ''.

attenuationFrom='era5file'; %Name part of the file that infile that states which attenuation data should be used: 
% ecmwf, ecmwfsonde, era5, era5sonde, sonde, era5file, ecmwffile

oceanTemp=27;
salinity=35; % Ocean salinity for sig0model in per mille (world wide default is 35) and sensitivity to that number is low
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

outName=[project,'_',quality,'_',attenuationFrom];

if addNSCAL
    outName=[outName,'_addNScal'];
end
if addSIG0model
    outName=[outName,'_sig0model'];
end

outName=[outName,addName];

directories.figdir=['/h/eol/romatsch/hcrCalib/oceanScans/figsCompare/',outName,'/'];

if ~exist(directories.figdir, 'dir')
    mkdir(directories.figdir)
end

directories.dataDir=HCRdir(project,quality,dataFreq);

directories.sondedir=['/scr/sci/romatsch/data/dropSonde/',project,'/'];

if strcmp(attenuationFrom,'era5') | strcmp(attenuationFrom,'era5sonde')
    [~,directories.modeldir]=modelDir(project,'era5',freqData);
elseif strcmp(attenuationFrom,'ecmwf') | strcmp(attenuationFrom,'ecmwfsonde')
    [~,directories.modeldir]=modelDir(project,'ecmwf',freqData);
end
%directories.wavedir=['/scr/sci/romatsch/data/waves/',project,'/'];
%directories.SST='/scr/sci/romatsch/data/SST/';

%File with start date, end date, and number for attenuation: 3=from the cfradial file, 2=model,
%1=sounding, 0=no sounding (data will not be used)
% Column 14: 1, if case is good, 0, if case is bad
infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/ssCal/biasInFiles/cal_',project,'_',attenuationFrom,'.txt'];

reflXlims=[0 25];
reflYlims=[15,55];
sigXlims=[0 25];
sigYlims=[-20 15];

%% If NSCAL results should be used

if addNSCAL
    if ~strcmp(quality,'field')
        disp('NS cal evaluation should be done on raw data.')
        return
    end
    if strcmp(project,'otrec')
        directories.highResTempDir='/scr/snow2/rsfdata/projects/otrec/hcr/txt/OTREC.temperatures.txt';
    elseif strcmp(project,'socrates')
        directories.highResTempDir='/scr/snow2/rsfdata/projects/socrates/hcr/qc/temperatures1s/';
    elseif strcmp(project,'cset')
        directories.highResTempDir='/h/eol/romatsch/data/hcrCalib/temps/cset_temps.txt';
    elseif strcmp(project,'aristo')
        directories.highResTempDir='/h/eol/romatsch/data/hcrCalib/temps/aristo_temps.txt';
    end
    
    nscalFile=['/h/eol/romatsch/hcrCalib/nsCal/inFiles/tempsTable_',project,'.txt'];
else
    nscalFile=nan;
    refPodTemp=nan;
end

%% Run processing

PLT=f_processOceanScans_simple(directories,infile,makeSingleFigs,project,...
    quality,nscalFile,attenuationFrom,...
    reflXlims,reflYlims,sigXlims,sigYlims,addName,addSIG0model,oceanTemp,salinity);

%% Plot results

uniqueCases = table2array(readtable(infile,'Delimiter','space'));

% Add last row if not in file
if size(uniqueCases,2)<14
    uniqueCases=cat(2,uniqueCases,ones(size(uniqueCases,1),1));
end

f_plot_all(PLT,sigXlims,sigYlims,reflXlims,reflYlims,outName,directories,uniqueCases)
