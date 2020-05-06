% Compare different data sets to make sure the data is complete

clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='socrates'; % socrates, cset, aristo, otrec
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
    
    dataVarsGood={};
    dataVarsTest={};
    
    for jj=1:size(compareVars,1)
        dataVarsGood{end+1}=compareVars{jj,1};
        dataVarsTest{end+1}=compareVars{jj,2};
    end
    
    dataVarsGood=dataVarsGood';
    dataVarsTest=dataVarsTest';

% Go through flights
for ii=1:size(caseList,1)
    disp(['Flight ',num2str(ii),' of ',num2str(size(caseList,1))]);
    
    startTime=datetime(caseList(ii,1:6));
    %endTime=datetime(caseList(ii,7:12));
    endTime=startTime+hours(1);
    
    %% Load good data
    disp('Loading good data.');
    
    fileListGood=makeFileList(indirGood,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if length(fileListGood)==0
        disp('No good data files found.');
        startTime=endTime;
        continue
    end
    
    % Go through each variable
    for ll=1:size(compareVars,1)
        
        dataGood=[];
        
        if ll==1 & any([compareVars{:,3}])==1
            dataGood.FLAG=[];
        end
        dataGood.(compareVars{ll,1})=[];
        
        % Load data
        dataGood=read_HCR(fileListGood,dataGood,startTime,endTime);
        
        if ll==1 & any([compareVars{:,3}])==1
            FLAG=dataGood.FLAG;
        end
        
        % Apply flag field
        if compareVars{ll,3}==1
            varChange=dataGood.(compareVars{ll,1});
            varChange(FLAG>1)=nan;
            dataGood.(compareVars{ll,1})=varChange;
        end
        
        % Set up matrix
        if ll==1            
            goodNans=array2table(nan(length(dataGood.time),size(dataVarsGood,1)),'VariableNames',dataVarsGood);
            goodNans=cat(2,table(dataGood.time'),goodNans);
            goodNans.Properties.VariableNames{'Var1'}='time';
            
            TTgood=table2timetable(goodNans);
            sumGood=nan(size(dataVarsGood));
        end
        
        nanMask=zeros(size(dataGood.(dataVarsGood{ll})));
        nanMask(~isnan(dataGood.(dataVarsGood{ll})))=1;
        
        TTgood.(dataVarsGood{ll})=sum(nanMask,1)';
        sumGood(ll)=length(find(TTgood.(dataVarsGood{ll})>0));
    
    end
      
     %% Load test data
    disp('Loading test data.');
    
    fileListTest=makeFileList(indirTest,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if length(fileListTest)==0
        disp('No test data files found.');
        startTime=endTime;
        continue
    end
    
    % Go through each variable
    for ll=1:size(compareVars,1)
        
        dataTest=[];
        
        if ll==1 & any([compareVars{:,3}])==1
            dataTest.FLAG=[];
        end
        
        dataTest.(compareVars{ll,2})=[];
        
        % Load data
        dataTest=read_HCR(fileListTest,dataTest,startTime,endTime);
        
        % Set up matrix
        if ll==1            
            testNans=array2table(nan(length(dataTest.time),size(dataVarsTest,1)),'VariableNames',dataVarsTest);
            testNans=cat(2,table(dataTest.time'),testNans);
            testNans.Properties.VariableNames{'Var1'}='time';
            
            TTtest=table2timetable(testNans);
            sumTest=nan(size(dataVarsTest));
        end
        
        nanMask=zeros(size(dataTest.(dataVarsTest{ll})));
        nanMask(~isnan(dataTest.(dataVarsTest{ll})))=1;
        
        TTtest.(dataVarsTest{ll})=sum(nanMask,1)';
        sumTest(ll)=length(find(TTtest.(dataVarsTest{ll})>0));
    
    end
    
    %% Process
    
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
    dataAllNum=table2array(tableAll(:,2:end));
    
    dataAll=nan(size(dataAllNum));
    dataAll(dataAllNum==0)=0;
    dataAll(dataAllNum>0)=1;
    
    %% Missing good data
    missingGood=nan(size(TTall,1),size(dataVarsGood,1));
    
    for jj=1:size(dataVarsGood,1)
       missingGood(find(dataAll(:,jj)==0 & dataAll(:,jj+size(dataVarsGood,1))==1),jj)=jj+2;
    end
    
    missingGood=cat(2,nan(size(TTall,1),2),missingGood);
    missingGood(find(isnan(dataAll(:,1))),1)=1;
    
    missingGood(find(isnan(dataAll(:,size(dataVarsGood,1)+1)) & missingGood(:,1)~=1),2)=2;
    missingGoodSinglePoint=nan(size(missingGood));
    missingGoodSingle=nan(size(missingGood));
    missingGoodDouble=nan(size(missingGood));
    
    % Remove single rays
    maskMiss=zeros(size(missingGood));
    maskMiss(~isnan(missingGood))=1;
    
    for jj=1:size(missingGood,2)
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
            
            for kk=1:length(startInds)
                if endInds(kk)-startInds(kk)==0
                    missingGood(startInds(kk),jj)=nan;
                    if jj>2 & (dataAll(startInds(kk),jj-2)-dataAll(startInds(kk),jj+size(dataVarsGood,1)-2)<0 & ...
                            dataAll(startInds(kk),jj-2)-dataAll(startInds(kk),jj+size(dataVarsGood,1)-2)>=-2)
                        missingGoodSinglePoint(startInds(kk),jj)=jj;
                    else
                        missingGoodSingle(startInds(kk),jj)=jj;
                    end
                elseif endInds(kk)-startInds(kk)==1
                    missingGood(startInds(kk):endInds(kk),jj)=nan;
                    if jj>2 & (dataAllNum(startInds(kk),jj-2)-dataAllNum(startInds(kk),jj+size(dataVarsGood,1)-2)<0 & ...
                            dataAllNum(startInds(kk),jj-2)-dataAllNum(startInds(kk),jj+size(dataVarsGood,1)-2)>=-2 & ...
                            dataAllNum(endInds(kk),jj-2)-dataAllNum(endInds(kk),jj+size(dataVarsGood,1)-2)<0 & ...
                            dataAllNum(endInds(kk),jj-2)-dataAllNum(endInds(kk),jj+size(dataVarsGood,1)-2)>=-2)
                        missingGoodSinglePoint(startInds(kk):endInds(kk),jj)=jj;
                    else
                        missingGoodDouble(startInds(kk):endInds(kk),jj)=jj;
                    end
                end
            end
        end
    end
    
     %% Missing test data
    missingTest=nan(size(TTall,1),size(dataVarsGood,1));
    
    for jj=1:size(dataVarsGood,1)
       missingTest(find(dataAll(:,jj)==1 & dataAll(:,jj+size(dataVarsGood,1))==0),jj)=jj+2;
    end
    
    missingTest=cat(2,nan(size(TTall,1),2),missingTest);
    missingTest(find(isnan(dataAll(:,1))),1)=1;
    
    missingTest(find(isnan(dataAll(:,size(dataVarsGood,1)+1)) & missingTest(:,1)~=1),2)=2;
    missingTestSinglePoint=nan(size(missingTest));
    missingTestSingle=nan(size(missingTest));
    missingTestDouble=nan(size(missingTest));
    
    % Remove single rays
    maskMiss=zeros(size(missingTest));
    maskMiss(~isnan(missingTest))=1;
    
    for jj=1:size(missingTest,2)
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
            
            for kk=1:length(startInds)
                if endInds(kk)-startInds(kk)==0
                    missingTest(startInds(kk),jj)=nan;
                    if jj>2 & (dataAllNum(startInds(kk),jj-2)-dataAllNum(startInds(kk),jj+size(dataVarsGood,1)-2)>0 & ...
                            dataAllNum(startInds(kk),jj-2)-dataAllNum(startInds(kk),jj+size(dataVarsGood,1)-2)<=2)
                        missingTestSinglePoint(startInds(kk),jj)=jj;
                    else
                        missingTestSingle(startInds(kk),jj)=jj;
                    end
                elseif endInds(kk)-startInds(kk)==1
                    missingTest(startInds(kk):endInds(kk),jj)=nan;
                    if jj>2 & (dataAllNum(startInds(kk),jj-2)-dataAllNum(startInds(kk),jj+size(dataVarsGood,1)-2)>0 & ...
                            dataAllNum(startInds(kk),jj-2)-dataAllNum(startInds(kk),jj+size(dataVarsGood,1)-2)<=2 & ...
                            dataAllNum(endInds(kk),jj-2)-dataAllNum(endInds(kk),jj+size(dataVarsGood,1)-2)>0 & ...
                            dataAllNum(endInds(kk),jj-2)-dataAllNum(endInds(kk),jj+size(dataVarsGood,1)-2)<=2)
                        missingTestSinglePoint(startInds(kk):endInds(kk),jj)=jj;
                    else
                        missingTestDouble(startInds(kk):endInds(kk),jj)=jj;
                    end
                end
            end
        end
    end
    
    %% Plot
    
    makeLabels=cat(1,{'time'},dataVarsGood);
    makeLabels2=cat(1,{'allMissing'},makeLabels);
    
    close all
    
    f1=figure('DefaultAxesFontSize',12,'renderer','painters','defaultAxesTickLabelInterpreter','none');
    set(f1,'Position',[5 5 1700 1000]);
    
    ax1=subplot(3,1,1);
    ax1.Position = [0.1300    0.7093    0.7750    0.25];

    hold on
    bar(cat(2,sumGood,sumTest));
    ylabel('Seconds');
    xticks(1:size(makeLabels,1));
    xticklabels(makeLabels);
    xtickangle(45);
    
    text(1:size(makeLabels,1),zeros(size(makeLabels)),num2str(sumGood-sumTest),...
        'HorizontalAlignment','center','VerticalAlignment','bottom','BackgroundColor','w',...
        'Margin',0.2,'EdgeColor','k');
    
    legend(freqGood,freqTest,'location','best');
    title([project,' RF ',num2str(ii),' ',qualityTest,' ',freqTest,' vs ',qualityGood,' ',freqGood]);
    
    ax2=subplot(3,1,2);
    ax2.Position = [0.1300    0.36    0.7750    0.25];
    hold on
    
    for jj=1:size(missingTest,2)
        hold on
        s4=scatter(TTall.time,missingTestSinglePoint(:,jj),'g+');
        s3=scatter(TTall.time,missingTestSingle(:,jj),'b+');
        s2=scatter(TTall.time,missingTestDouble(:,jj),'k+');
        s1=scatter(TTall.time,missingTest(:,jj),'r+');       
    end
        
    xlim([TTall.time(1) TTall.time(end)]);
    ylim([0 size(missingGood,2)+1]);
    grid on
    
    yticks(1:size(missingTest,2));
    yticklabels(makeLabels2);
    
    title(['Missing times ',qualityTest,' ',freqTest]);
    grid on
        
    ax3=subplot(3,1,3);
    ax3.Position = [0.1300    0.05    0.7750    0.25];
    hold on
    
    for jj=1:size(missingGood,2)
        hold on
        s4=scatter(TTall.time,missingGoodSinglePoint(:,jj),'g+');
        s3=scatter(TTall.time,missingGoodSingle(:,jj),'b+');
        s2=scatter(TTall.time,missingGoodDouble(:,jj),'k+');
        s1=scatter(TTall.time,missingGood(:,jj),'r+');       
    end
    
    legend([s1 s2 s3 s4],{'> 2 Rays','2 Rays, > 2 Gates','1 Ray, > 2 Gates','1 or 2 Rays, 1 Gate'},'location','best');
    
    xlim([TTall.time(1) TTall.time(end)]);
    ylim([0 size(missingGood,2)+1]);
    grid on
    
    yticks(1:size(missingGood,2));
    yticklabels(makeLabels2);
    
    title(['Missing times ',qualityGood,' ',freqGood]);
    grid on
    
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_RF',num2str(ii),'_dataCheck_',qualityTest,'_',freqTest,'_vs_',qualityGood,'_',freqGood,'.png'],'-dpng','-r0')
    
end
