% Test utils files

clear all
close all

%% User input

yLimits=[-0.2,2]; % Altitude limits for plots in km

startTime=datetime(2018,2,18,4,55,0); % Start time as year,month,day,hour,minute,second
endTime=datetime(2018,2,18,5,20,0); % End time as year,month,day,hour,minute,second

saveFig=1; % If 1, figures will be saved, if 0, figures will not be saved
outDir='./'; % Output directories for figures
figNames={'socrates1','socrates2','socrates3','socrates4'}; % Figure names

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

% Data directory
indir='/scr/snow2/rsfdata/projects/socrates/hcr/qc2/cfradial/moments/10hz/';
 
%% Process

dataVars=fieldnames(data);

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