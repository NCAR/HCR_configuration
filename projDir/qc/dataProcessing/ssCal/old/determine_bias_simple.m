% Ocean scan calibration for HCR data

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

% If 1, plots for individual calibration events will be made, if 0, only
% total plots will be made
makeSingleFigs=0;
addNSCAL=1;

project='socrates'; %socrates
quality='raw'; %raw or final

attenuationFrom='model'; %model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/h/eol/romatsch/gitPriv/process_HCR/oceanScans/functions/');
addpath('/h/eol/romatsch/gitPriv/process_HCR/oceanScans/colormaps/');
addpath('/h/eol/romatsch/gitPriv/process_HCR/NSCAL/functions/');
addpath('/h/eol/romatsch/gitPriv/utils/');

directories.figdir=['/h/eol/romatsch/hcrCalib/oceanScans/figs/',project,'Compare/'];

if strcmp(project,'socrates')
    if strcmp(quality,'raw')
        directories.dataDir='/scr/rain1/rsfdata/projects/socrates/hcr/cfradial/moments/10hz/'; % raw data
    elseif strcmp(quality,'final')
        directories.dataDir='/scr/rain1/rsfdata/projects/socrates/hcr/qc/cfradial/moments/10hz/'; %Final
    end
end

directories.sondedir=['/scr/sci/romatsch/data/dropSonde/',project,'/'];
directories.modeldir=['/scr/sci/romatsch/data/reanalysis/ecmwf/era5/',project,'/'];
%directories.wavedir=['/scr/sci/romatsch/data/waves/',project,'/'];

%File with start date, end date, and number for attenuation: 1=era5,
%2=sounding
infile=['/h/eol/romatsch/hcrCalib/oceanScans/biasInFiles/cal_',project,'_',attenuationFrom,'.txt'];

if addNSCAL
    if ~strcmp(quality,'raw')
        disp('NS cal evaluation should be done on raw data.')
        return
    end
    if strcmp(project,'socrates')
        directories.highResTempDir='/scr/rain1/rsfdata/projects/socrates/hcr/qc/temperatures1s/';
        [refPodTemp refPodTempStd]=f_getRefTemp([directories.dataDir '20171113/'],datetime(2017,11,13,23,11,30),datetime(2017,11,13,23,32,17),'podTemp');
    elseif strcmp(project,'cset')
        directories.highResTempDir='/h/eol/romatsch/data/hcrCalib/temps/';
        [refPodTemp refPodTempStd]=f_getRefTemp(['/scr/eldora2/rsfdata/cset/hcr/cfradial/moments/10hz/20140825/'],datetime(2014,8,25,22,30,0),datetime(2014,8,25,23,0,0),'podTemp');
    elseif strcmp(project,'aristo')
        directories.highResTempDir='/h/eol/romatsch/data/hcrCalib/temps/';
        [refPodTemp refPodTempStd]=f_getRefTemp(['/scr/eldora2/rsfdata/cset/hcr/cfradial/moments/10hz/20140825/'],datetime(2014,8,25,22,30,0),datetime(2014,8,25,23,0,0),'podTemp');
    end
    
    nscalFile=['/h/eol/romatsch/hcrCalib/nsCal/inFiles/meansTable_',project,'.txt'];
else
    nscalFile=nan;
    refPodTemp=nan;
end

[outTable PLT]=f_processOceanScans(directories,infile,makeSingleFigs,project,nscalFile,refPodTemp);

uniqueCases = table2array(readtable(infile,'Delimiter','space'));

reflXlims=[0 25];
reflYlims=[15,55];
sigXlims=[0 25];
sigYlims=[-20 15];

plot_all(PLT,sigXlims,sigYlims,reflXlims,reflYlims,project,directories,uniqueCases)

