% Test utils files

clear all
close all

%% User input

yLimits=[-0.2,2]; % Altitude limits for plots in km

project='socrates'; % socrates, cset, aristo
quality='qc2'; % field, qc1, qc2

startTime=datetime(2018,2,18,4,55,0); % Start time as year,month,day,hour,minute,second
endTime=datetime(2018,2,18,5,20,0); % End time as year,month,day,hour,minute,second

saveFig=1; % If 1, figures will be saved, if 0, figures will not be saved
outDir='/h/eol/romatsch/testOrder/';
figNames={'socratesQC1'};

% Desired variables. The variable name comies after the . and must be spelled exactly
% like in the CfRadial file

data.DBZ=[];
data.VEL=[];
%data.VEL_RAW=[];
%data.VEL_CORR=[];
%data.WIDTH=[];
%data.WIDTH_CORR=[];
%data.DBMVC=[];
%data.SNR=[];
%data.NCP=[];
%data.LDR=[];

%% Data directory

dataVars=fieldnames(data);

if strcmp(project,'socrates')
    if strcmp(quality,'field')
        indir='/scr/snow2/rsfdata/projects/socrates/hcr/cfradial/moments/10hz/';
    elseif strcmp(quality,'qc1')
        indir='/scr/snow2/rsfdata/projects/socrates/hcr/qc/cfradial/velcorr/10hz/';
    elseif strcmp(quality,'qc2')
        indir='/scr/snow2/rsfdata/projects/socrates/hcr/qc2/cfradial/moments/10hz/';
    end
elseif strcmp(project,'cset')
    if strcmp(quality,'field')
        indir='/scr/snow2/rsfdata/projects/cset/hcr/cfradial/moments/10hz/';
    elseif strcmp(quality,'qc1')
        indir='/scr/snow2/rsfdata/projects/cset/hcr/qc/cfradial/velcorr/10hz/';
    elseif strcmp(quality,'qc2')
        indir='/scr/snow2/rsfdata/projects/cset/hcr/qc2/cfradial/velcorr/10hz/';
    end
elseif strcmp(project,'aristo')
    if strcmp(quality,'field')
        indir='/scr/snow2/rsfdata/projects/aristo-17/hcr/cfradial/moments/10hz/';
    elseif strcmp(quality,'qc1')
        disp('There are no qc1 data for aristo.')
        return
    elseif strcmp(quality,'qc2')
        indir='/scr/snow2/rsfdata/projects/aristo-17/hcr/qc2/cfradial/velcorr/10hz/';
    end
end

%% Process
% Make list of files within the specified time frame
fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

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

%% Plot data
close all
figs=plot_HCR(data,dataVars,yLimits);

if saveFig
    for ii=1:length(figs)
        print(figs.(['fig',num2str(ii)]), [outDir,figNames{ii}],'-dpng','-r0');
    end
end