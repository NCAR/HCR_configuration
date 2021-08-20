% Calculate liquid water content from HCR ocean return

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='otrec'; %socrates, aristo, cset
% quality='qc2'; %field, qc1, or qc2
% qcVersion='v2.2';
% freqData='10hz'; % 10hz, 100hz, 2hz, or combined
whichModel='era5';

saveTime=1;

if strcmp(project,'otrec')
    indir='/scr/sleet2/rsfdata/projects/otrec/hcr/qc2/cfradial/development/convStrat/10hz/';
elseif strcmp(project,'socrates')
    indir='/scr/snow2/rsfdata/projects/socrates/hcr/qc2/cfradial/development/convStrat/10hz/';
end

outdir=[indir(1:end-47),'mat/convStrat/10hz/'];

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

for aa=7:size(caseList,1)
    disp(['Flight ',num2str(aa)]);
    disp(['Starting at ',datestr(datetime('now'),'yyyy-mm-dd HH:MM')]);
    disp('Loading data ...');
    
    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
      
    %% Get data
       
    data=[];
    
    data.DBZ = [];
    data.FLAG=[];
    data.TOPO=[];
    data.CLOUD_PUZZLE=[];
    data.TEMP=[];
    data.MELTING_LAYER=[];
    
    dataVars=fieldnames(data);
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
    
    %frq=ncread(fileList{1},'frequency');
    
    % Check if all variables were found
    for ii=1:length(dataVars)
        if ~isfield(data,dataVars{ii})
            dataVars{ii}=[];
        end
    end
    
    dataVars=dataVars(~cellfun('isempty',dataVars));
    
    %% Cloud puzzle
    %disp('Making cloud puzzle');
    
    %cloudPuzzle=f_cloudPuzzle_radial(data);
    cloudPuzzle=data.CLOUD_PUZZLE;
    uClouds=unique(cloudPuzzle(~isnan(cloudPuzzle)));
    uClouds(uClouds==0)=[];
    %cloudCount=length(uClouds);
    
    %% Calculate reflectivity texture and convectivity
    %stratConvNearSurf=nan(size(data.DBZ));
    
    dbzText=nan(size(data.DBZ));
    convectivity=nan(size(data.DBZ));
    classBasic=nan(size(data.DBZ));
    classSub=nan(size(data.DBZ));
        
    pixRad=50; % Radius over which texture is calculated in pixels. Default is 50.
    dbzBase=-10; % Reflectivity base value which is subtracted from DBZ.
    
    upperLim=14; % Upper limit for convectivity mapping. Texture above that will be set to 1.
    stratMixed=0.4; % Convectivity boundary between strat and mixed.
    mixedConv=0.5; % Convectivity boundary between mixed and conv.
    
    for jj=1:length(uClouds)
        disp(['Calculating texture for cloud ',num2str(jj),' of ',num2str(length(uClouds))]);
        dbzIn=data.DBZ;
        dbzIn(data.FLAG>1)=nan;
        dbzIn(cloudPuzzle~=uClouds(jj))=nan;
        
        % Shrink to good data area
        nonNanCols=find(any(~isnan(dbzIn),1));
        dbzIn=dbzIn(:,nonNanCols);
        
        dbzTextOne=f_reflTexture(dbzIn,pixRad,dbzBase);
                
        % Convectivity        
        convOne=1/upperLim.*dbzTextOne;
        
        % Basic classification       
        classBasicOne=f_classBasic(convOne,stratMixed,mixedConv);
                        
        % Fill into large matrix
        dbzTextLarge=nan(size(dbzText));
        dbzTextLarge(:,nonNanCols)=dbzTextOne;
        dbzText(~isnan(dbzTextLarge))=dbzTextLarge(~isnan(dbzTextLarge));
        
        convLarge=nan(size(convectivity));
        convLarge(:,nonNanCols)=convOne;
        convectivity(~isnan(convLarge))=convLarge(~isnan(convLarge));
        
        classBasicLarge=nan(size(classBasic));
        classBasicLarge(:,nonNanCols)=classBasicOne;
        classBasic(~isnan(classBasicLarge))=classBasicLarge(~isnan(classBasicLarge));
        
        % Prepare asl and topo and temperature
        topoIn=data.TOPO(nonNanCols);
        aslIn=data.asl(:,nonNanCols);
        aslIn(isnan(dbzIn))=nan;
        tempIn=data.TEMP(:,nonNanCols);
        meltIn=data.MELTING_LAYER(:,nonNanCols);
        elevIn=data.elevation(nonNanCols);
        
        classSubOne=f_classSub(classBasicOne,aslIn,topoIn,meltIn,tempIn,elevIn);
        
        classSubLarge=nan(size(classSub));
        classSubLarge(:,nonNanCols)=classSubOne;
        classSub(~isnan(classSubLarge))=classSubLarge(~isnan(classSubLarge));
    end
    
    convStrat=classSub;
    convStrat1D=max(convStrat,[],1);
    
    %% Save
    disp('Saving stratConv field ...')
    
    save([outdir,whichModel,'.convectivity.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'convectivity');
    
    save([outdir,whichModel,'.convStrat.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'convStrat');
    
    save([outdir,whichModel,'.convStrat1D.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'convStrat1D');
    
    if saveTime
        timeHCR=data.time;
        save([outdir,whichModel,'.time.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
            datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'timeHCR');
    end
end