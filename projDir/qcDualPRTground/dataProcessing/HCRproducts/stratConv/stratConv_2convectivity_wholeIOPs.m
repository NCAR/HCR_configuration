% Calculate convectivity and convective/stratiform partitioning
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/'));

project='meow';
quality='qc1';
freqData='10hz_combined';
qcVersion='v1.0';
whichModel='hrrr';

infile=['~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/scriptsFiles/iops_',project,'.txt'];

saveTime=0;

if strcmp(project,'meow')
    upperLimDBZ=14; %12
    upperLimVEL=5;
    pixRadDBZ=300; % Radius over which texture is calculated in pixels. Default is 50.
    pixRadVEL=300;
else
    error('Set upperLimDBZ and upperLimVEL')
end

indir=HCRdir(project,quality,qcVersion,freqData);

[~,outdir]=modelDir(project,whichModel,quality,qcVersion,freqData);

caseList = table2array(readtable(infile));

for aa=1:size(caseList,1)
    disp(['IOP ',num2str(aa)]);
    disp(['Starting at ',datestr(datetime('now'),'yyyy-mm-dd HH:MM')]);
    disp('Loading data ...');
    
    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));
    
    %% Get data
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    data=[];
    
    data.DBZ=[];
    data.VEL=[];
    data.TEMP=[];
    data.MELTING_LAYER=[];
    data.FLAG_short=[];
                
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);

    data.TOPO=data.altitude;
    data.FLAG=data.FLAG_short;
    data=rmfield(data,'FLAG_short');

    %% Truncate to non missing
    gapSecs=10;
    nonMissingInds=findNonMissingInds(data,gapSecs);

    dataInVars=fields(data);

    dataShort=[];
    for ii=1:length(dataInVars)
        dataShort.(dataInVars{ii})=data.(dataInVars{ii})(:,nonMissingInds==1);
    end

    % Take care of up pointing VEL
    dataShort.VEL(:,dataShort.elevation<0)=-dataShort.VEL(:,dataShort.elevation<0);

    %% Texture from reflectivity and velocity

    disp('Calculating reflectivity texture ...');

    dbzBase=-10; % Reflectivity base value which is subtracted from DBZ.

    dbzText=f_reflTexture(dataShort.DBZ,pixRadDBZ,dbzBase);

    disp('Calculating velocity texture ...');

    velBase=-20; % VEL base value which is subtracted from DBZ.

    velText=f_velTexture(dataShort.VEL,pixRadVEL,velBase);

    %% Convectivity

    % Convectivity
    convDBZ=1/upperLimDBZ.*dbzText;

    convVEL=1/upperLimVEL.*velText;

    convectivityShort=convDBZ.*convVEL;
    convectivityShort(convectivityShort>1)=1;
    convectivityShort(isnan(convectivityShort))=convDBZ(isnan(convectivityShort));

    %% Basic classification

    disp('Basic classification ...');

    stratMixed=0.4; % Convectivity boundary between strat and mixed.
    mixedConv=0.5; % Convectivity boundary between mixed and conv.

    classBasic=f_classBasicBoth(convectivityShort,stratMixed,mixedConv,dataShort.MELTING_LAYER,dataShort.elevation);

    %% Sub classification

    disp('Sub classification ...');

    classSub=f_classSubBoth(classBasic,dataShort.asl,dataShort.TOPO,dataShort.MELTING_LAYER,dataShort.TEMP,dataShort.elevation);
          
    %% Save
    disp('Saving stratConv field ...')

    convectivity=nan(size(data.DBZ));
    convectivity(:,nonMissingInds==1)=convectivityShort;
    
    convStrat=nan(size(data.DBZ));
    convStrat(:,nonMissingInds==1)=classSub;
    convStrat1D=max(convStrat,[],1);
    
    save([outdir,whichModel,'.convectivity.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.IOP',num2str(aa),'.mat'],'convectivity','-v7.3');
    
    save([outdir,whichModel,'.convStrat.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.IOP',num2str(aa),'.mat'],'convStrat','-v7.3');
    
    save([outdir,whichModel,'.convStrat1D.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.IOP',num2str(aa),'.mat'],'convStrat1D','-v7.3');
    
    if saveTime
        timeHCR=data.time;
        save([outdir,whichModel,'.time.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
            datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.IOP',num2str(aa),'.mat'],'timeHCR','-v7.3');
    end
end