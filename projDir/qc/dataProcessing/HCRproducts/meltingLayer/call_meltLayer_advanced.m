% Analyze HCR clouds

clear all;
close all;

project='socrates'; %socrates, aristo, cset
quality='qc3'; %field, qc1, or qc2
qcVersion='v3.2';
freqData='10hz'; % 10hz, 100hz, or 2hz

%unfoldVelocity=1;
offsetIn=-300;

thresholds.LDRlimits=[-16,-7]; % SOCRATES, OTREC, CSET default: [-16,-7]
thresholds.LDRspeclePix=50; % SOCRATES, OTREC, CSET default: [] (not used)
thresholds.LDRsolidity=0.5; % SOCRATES, OTREC, CSET default: [] (not used)
thresholds.LDRsearchPix=25; % SOCRATES, OTREC, CSET default: 18
thresholds.LDRstd=150; % SOCRATES, OTREC, CSET default: 100
thresholds.LDRaltDiff=100; % SOCRATES, OTREC, CSET default: 50
thresholds.LDRareaPix=[];

thresholds.VELsearchPix=80; % SOCRATES, OTREC, CSET default: 50
thresholds.VELstd=70; % SOCRATES, OTREC, CSET default: 35
thresholds.VELaltDiff=100; % SOCRATES, OTREC, CSET default: 100
thresholds.VELudDiff=-0.7; % SOCRATES, OTREC, CSET default: -0.7
thresholds.VEL_LDRdiff=600; % SOCRATES, OTREC, CSET default: 200

thresholds.outlier=350; % SOCRATES, OTREC, CSET default: 50
thresholds.length=10; % SOCRATES, OTREC, CSET default: 20

% Determines plot zoom.
if strcmp(project,'otrec')
    ylimits=[-0.2 15];
elseif strcmp(project,'socrates')
    ylimits=[-0.2 5];
elseif strcmp(project,'spicule')
    ylimits=[-0.2 13];
elseif strcmp(project,'cset')
    ylimits=[-0.2 9];
end

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/meltLayer_',project,'.txt'];

indir=HCRdir(project,quality,qcVersion,freqData);

figdir=[indir(1:end-5),'meltLayer/cases/'];

if ~exist(figdir, 'dir')
    mkdir(figdir)
end

% Loop through cases

caseList=readtable(casefile);
caseStart=datetime(caseList.Var1,caseList.Var2,caseList.Var3, ...
    caseList.Var4,caseList.Var5,0);
caseEnd=datetime(caseList.Var6,caseList.Var7,caseList.Var8, ...
    caseList.Var9,caseList.Var10,0);

for aa=1:length(caseStart)
    
    disp(['Case ',num2str(aa),' of ',num2str(length(caseStart))]);
    
    startTime=caseStart(aa);
    endTime=caseEnd(aa);
    
    %% Load data
    
    disp('Loading data ...');
    
    data=[];
    
    data.DBZ_MASKED=[];
    data.VEL_MASKED=[];
    data.TEMP=[];
    data.WIDTH=[];
    data.FLAG=[];
    data.TOPO=[];
        
    % Make list of files within the specified time frame
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
        
    % Check if LDR_MASKED is available
    try
        velTest=ncread(fileList{1},'LDR_MASKED');
        data.LDR_MASKED=[];
    catch
        data.LDR=[];
    end
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
    
    if isfield(data,'LDR_MASKED')
        data.LDR=data.LDR_MASKED;
        data=rmfield(data,'LDR_MASKED');
    end

    data.LDR(data.FLAG~=1)=nan;
    
    %% Find melting layer

    f_meltLayer_advanced(data,offsetIn,thresholds,figdir);
       
end