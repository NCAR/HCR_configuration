% Compare different data sets to make sure the data is complete

clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='cset'; % socrates, cset, aristo, otrec
qualityGood='qc2'; % field, qc0, qc1, qc2
qualityTest='qc2'; % field, qc0, qc1, qc2
freqGood='10hz'; % 10hz, 2hz, 2hzMerged
freqTest='2hzMerged'; % 10hz, 2hz, 2hzMerged

figdir=['/h/eol/romatsch/hcrCalib/checkData/'];

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'.txt'];

caseList = table2array(readtable(infile));

indirGood=HCRdir(project,qualityGood,freqGood);
indirTest=HCRdir(project,qualityTest,freqTest);

%% Run processing

compareVars={'DBZ','HCR_DBZ',1;
    'WIDTH','HCR_WIDTH',1;
    'SNR','HCR_SNR',1;
    'DBMVC','HCR_DBMVC',0;
    'NCP','HCR_NCP',1;
    'LDR','HCR_LDR',1;
    'PRESS','PRESS',0;
    'TEMP','TEMP',0;
    'RH','RH',0;
    'FLAG','FLAG',0;
    'VEL_CORR','HCR_VEL',1;
    'SST','SST',0;
    'TOPO','TOPO',0;
    'U_SURF','U_SURF',0;
    'V_SURF','V_SURF',0};

% Go through flights
for ii=1:size(caseList,1)
    disp(['Flight ',num2str(ii),' of ',num2str(size(caseList,1))]);
    
    startTime=datetime(caseList(ii,1:6));
    %endTime=datetime(caseList(ii,7:12));
    endTime=startTime+hours(1);
    
    % Desired variables. The variable name comes after the . and must be spelled exactly
    % like in the CfRadial file
    dataGood=[];
    dataTest=[];
    
    for jj=1:size(compareVars,1)
        dataGood.(compareVars{jj,1})=[];
    end
    for jj=1:size(compareVars,1)
        dataTest.(compareVars{jj,2})=[];
    end
       
    dataVarsGood=fieldnames(dataGood);
    dataVarsTest=fieldnames(dataTest);
    
    %% Load good data
    disp('Loading good data.');
    
    fileListGood=makeFileList(indirGood,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if length(fileListGood)==0
        disp('No good data files found.');
        startTime=endTime;
        continue
    end
    
    % Load data
    dataGood=read_HCR(fileListGood,dataGood,startTime,endTime);
    
    if isempty(dataGood.time)
        disp('No good data found.');
        startTime=endTime;
        continue
    end
    
    % Check if all variables were found
    for kk=1:length(dataVarsGood)
        if ~isfield(dataGood,dataVarsGood{kk})
            dataVarsGood{kk}=[];
        end
    end
       
    %% Load test data
    disp('Loading test data.');
    
    fileListTest=makeFileList(indirTest,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if length(fileListTest)==0
        disp('No test data files found.');
        startTime=endTime;
        continue
    end
    
    % Load data
    dataTest=read_HCR(fileListTest,dataTest,startTime,endTime);
    
    if isempty(dataTest.time)
        disp('No test data found.');
        startTime=endTime;
        continue
    end
    
    % Check if all variables were found
    for kk=1:length(dataVarsTest)
        if ~isfield(dataTest,dataVarsTest{kk})
            dataVarsTest{kk}=[];
        end
    end
    
    % Apply flag field
    for jj=1:size(compareVars,1)
        if compareVars{jj,3}==1
            varChange=dataGood.(compareVars{jj,1});
            varChange(dataGood.FLAG>1)=nan;
            dataGood.(compareVars{jj,1})=varChange;
        end
    end
    %% Count seconds
    goodNans=array2table(nan(length(dataGood.time),size(dataVarsGood,1)),'VariableNames',dataVarsGood);
    goodNans=cat(2,table(dataGood.time'),goodNans);
    goodNans.Properties.VariableNames{'Var1'}='time';
    
    testNans=array2table(nan(length(dataTest.time),size(dataVarsTest,1)),'VariableNames',dataVarsTest);
    testNans=cat(2,table(dataTest.time'),testNans);
    testNans.Properties.VariableNames{'Var1'}='time';
    
    TTgood=table2timetable(goodNans);
    TTtest=table2timetable(testNans);
    
    sumGood=nan(size(dataVarsGood));
    sumTest=nan(size(dataVarsTest));
    
    for jj=1:size(dataVarsGood,1)
        TTgood.(dataVarsGood{jj})=double(any(~isnan(dataGood.(dataVarsGood{jj})),1))';
        sumGood(jj)=sum(TTgood.(dataVarsGood{jj}));
    end
    
    for jj=1:size(dataVarsTest,1)
        TTtest.(dataVarsTest{jj})=double(any(~isnan(dataTest.(dataVarsTest{jj})),1))';
        sumTest(jj)=sum(TTtest.(dataVarsTest{jj}));
    end
    
    sumGood=cat(1,size(TTgood,1),sumGood);
    sumTest=cat(1,size(TTtest,1),sumTest);
    
    if strcmp(freqGood,'10hz')
        sumGood=sumGood./10;
    elseif strcmp(freqGood,'2hz') | strcmp(freqGood,'2hzMerged')
        sumGood=sumGood./2;
    end
    
    if strcmp(freqTest,'10hz')
        sumTest=sumTest./10;
    elseif strcmp(freqTest,'2hz') | strcmp(freqTest,'2hzMerged')
        sumTest=sumTest./2;
    end
    
    
    %% Find where differences are
    if strcmp(freqTest,'2hz') | strcmp(freqTest,'2hzMerged')
        dt = seconds(0.5);
        TTgood2 = retime(TTgood,'regular','fillwithmissing','TimeStep',dt);
    else
        TTgood2=TTgood;
    end
    
    TTall=synchronize(TTgood2,TTtest,'first','fillwithmissing');
    tableAll=timetable2table(TTall);
    dataAll=table2array(tableAll(:,2:end));
    
    missingAll=nan(size(TTall,1),size(dataVarsGood,1));
    
    for jj=1:size(dataVarsGood,1)
       missingAll(find(dataAll(:,jj)==1 & dataAll(:,jj+size(dataVarsGood,1))==0),jj)=jj+2;
    end
    
    missingAll=cat(2,nan(size(TTall,1),2),missingAll);
    missingAll(find(isnan(dataAll(:,1))),1)=1;
    
    missingAll(find(isnan(dataAll(:,size(dataVarsGood,1)+1)) & missingAll(:,1)~=1),2)=2;
    missingSingle=nan(size(missingAll));
    missingDouble=nan(size(missingAll));
    
    %% Remove single rays
    maskMiss=zeros(size(missingAll));
    maskMiss(~isnan(missingAll))=1;
    
    for jj=1:size(missingAll,2)
        column=maskMiss(:,jj);
        diffCol=diff(column);
        
        startInds=find(diffCol==1);
        startInds=startInds+1;
        endInds=find(diffCol==-1);
        
        if ~isempty(startInds) & ~isempty(endInds)
            if endInds(1)<startInds(1)
                startInds=[1;startInds];
            end
            if length(startInds)~=length(endInds);
                endInds=[endInds;length(column)];
            end
        end
        
        for kk=1:length(startInds)
            if endInds(kk)-startInds(kk)==0
                missingAll(startInds(kk),jj)=nan;
                missingSingle(startInds(kk),jj)=jj;
            elseif endInds(kk)-startInds(kk)==1
                missingAll(startInds(kk):endInds(kk),jj)=nan;
                missingDouble(startInds(kk):endInds(kk),jj)=jj;
            end
        end
    end
    
    %% Plot
    
    makeLabels=cat(1,{'time'},dataVarsGood);
    makeLabels2=cat(1,{'allMissing'},makeLabels);
    
    close all
    
    f1=figure('DefaultAxesFontSize',12,'renderer','painters','defaultAxesTickLabelInterpreter','none');
    set(f1,'Position',[5 5 1700 1000]);
    
    subplot(2,1,1)
    hold on
    bar(cat(2,sumGood,sumTest));
    ylabel('Seconds');
    xticks(1:size(makeLabels,1));
    xticklabels(makeLabels);
    xtickangle(45);
    
    legend(freqGood,freqTest,'location','best');
    title([project,' RF ',num2str(ii),' ',qualityTest,' ',freqTest,' vs ',qualityGood,' ',freqGood]);
    
    subplot(2,1,2)
    hold on
    
    for jj=1:size(missingAll,2)
        hold on
        s3=scatter(TTall.time,missingSingle(:,jj),'g+');
        s2=scatter(TTall.time,missingDouble(:,jj),'r+');
        s1=scatter(TTall.time,missingAll(:,jj),'b+');       
    end
    
    legend([s1 s2 s3],{'Multi','Double','Single'});
    
    xlim([TTall.time(1) TTall.time(end)]);
    grid on
    
    yticks(1:size(missingAll,2));
    yticklabels(makeLabels2);
    
    title('Missing times');
    grid on
    
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_RF',num2str(ii),'_dataCheck_',qualityTest,'_',freqTest,'_vs_',qualityGood,'_',freqGood,'.png'],'-dpng','-r0')
    
end
