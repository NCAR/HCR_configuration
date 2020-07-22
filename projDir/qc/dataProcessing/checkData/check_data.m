% Compare different data sets to make sure the data is complete

clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='cset'; % socrates, cset, aristo, otrec
qualityGood='qc2'; % field, qc0, qc1, qc2
qualityTest='qc2'; % field, qc0, qc1, qc2
freqGood='10hz'; % 10hz, 2hz, 2hzMerged
freqTest='2hzMerged'; % 10hz, 2hz, 2hzMerged

thresh12=[10 50];

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
    endTime=datetime(caseList(ii,7:12));
    %endTime=startTime+hours(1);
        
    fileListGood=makeFileList(indirGood,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    fileListTest=makeFileList(indirTest,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    %% Load good data
        
    if length(fileListGood)==0
        disp('No good data files found.');
        startTime=endTime;
        continue
    end
    if length(fileListTest)==0
        disp('No test data files found.');
        startTime=endTime;
        continue
    end
    
    % Go through each variable
    for ll=1:size(compareVars,1)
        disp(['Loading good data ',compareVars{ll,1}]);
        
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
        
        % Resample
        if ll==1
            startTimeReal=dateshift(dataGood.time(1),'start','second');
            endTimeReal=dateshift(dataGood.time(end),'end','second');
            if strcmp(freqGood,'10hz') & strcmp(freqTest,'10hz')
                newTime=startTimeReal:seconds(0.1):endTimeReal;
            elseif strcmp(freqGood,'2hz') | strcmp(freqTest,'2hz') | strcmp(freqGood,'2hzMerged') | strcmp(freqTest,'2hzMerged')
                newTime=startTimeReal:seconds(0.5):endTimeReal;
            else
                didsp('Wrong data resolution')
                return
            end
            
            [~,iaG,ibG]=intersect(dataGood.time,newTime);
        end
        
        newDataG=nan(size(dataGood.(compareVars{ll,1}),1),length(newTime));
        newDataG(:,ibG)=dataGood.(compareVars{ll,1})(:,iaG);
        
        %% Load test data
        disp(['Loading test data ',compareVars{ll,2}]);
        
        dataTest=[];
        
        dataTest.(compareVars{ll,2})=[];
        
        % Load data
        dataTest=read_HCR(fileListTest,dataTest,startTime,endTime);
        
        % Resample
        if ll==1
            [~,iaT,ibT]=intersect(dataTest.time,newTime);
        end
        
        newDataT=nan(size(dataTest.(compareVars{ll,2}),1),length(newTime));
        newDataT(:,ibT)=dataTest.(compareVars{ll,2})(:,iaT);
        
        % Mask with data
        maskG=zeros(size(newDataG));
        maskT=zeros(size(newDataT));
        
        maskG(isnan(newDataT) & ~isnan(newDataG))=1;
        maskT(isnan(newDataG) & ~isnan(newDataT))=1;
        
        % Find speckle
        countG=zeros(size(maskG));
        countT=zeros(size(maskT));
        
        CCG = bwconncomp(maskG);
        CCT = bwconncomp(maskT);
        
        for kk=1:CCG.NumObjects
            area=CCG.PixelIdxList{kk};
            countG(area)=length(area);
        end
        
        for kk=1:CCT.NumObjects
            area=CCT.PixelIdxList{kk};
            countT(area)=length(area);
        end
        
        maxCountG=max(countG,[],1);
        maxCountT=max(countT,[],1);
        
        % Set up matrix
        if ll==1
            plotMatG=nan(length(newTime),size(dataVarsGood,1),3);
            plotMatT=nan(length(newTime),size(dataVarsTest,1),3);
            
            sumGood=nan(size(dataVarsGood));
            sumTest=nan(size(dataVarsTest));
        end
        
        plotMatG((maxCountG>0 & maxCountG<=thresh12(1)),ll,1)=ll;
        plotMatG((maxCountG>thresh12(1) & maxCountG<=thresh12(2)),ll,2)=ll;
        plotMatG(maxCountG>thresh12(2),ll,3)=ll;
        
        plotMatT((maxCountT>0 & maxCountT<=thresh12(1)),ll,1)=ll;
        plotMatT((maxCountT>thresh12(1) & maxCountT<=thresh12(2)),ll,2)=ll;
        plotMatT(maxCountT>thresh12(2),ll,3)=ll;
        
        % Sum of good data
        allDataG=any(~isnan(newDataG),1);
        sumGood(ll)=sum(double(allDataG));
        
        allDataT=any(~isnan(newDataT),1);
        sumTest(ll)=sum(double(allDataT));
    end
    
    %% Add times
    plotMatG=plotMatG+2;
    plotMatT=plotMatT+2;
    plotMatG=cat(2,nan(size(plotMatG,1),2,size(plotMatG,3)),plotMatG);
    plotMatT=cat(2,nan(size(plotMatT,1),2,size(plotMatT,3)),plotMatT);
    
    timeDummyG=ones(size(dataGood.time));
    timeDummyT=ones(size(dataTest.time));
    
    newDataG=nan(size(dataGood.time,1),length(newTime));
    newDataG(:,ibG)=timeDummyG(:,iaG);
    newDataT=nan(size(dataTest.time,1),length(newTime));
    newDataT(:,ibT)=timeDummyT(:,iaT);
    
    % Mask with data
    maskG=zeros(size(newDataG));
    maskT=zeros(size(newDataT));
    
    maskG(isnan(newDataT) & ~isnan(newDataG))=1;
    maskT(isnan(newDataG) & ~isnan(newDataT))=1;
    
    % Find speckle
    countG=zeros(size(maskG));
    countT=zeros(size(maskT));
    
    CCG = bwconncomp(maskG);
    CCT = bwconncomp(maskT);
    
    for kk=1:CCG.NumObjects
        area=CCG.PixelIdxList{kk};
        countG(area)=length(area);
    end
    
    for kk=1:CCT.NumObjects
        area=CCT.PixelIdxList{kk};
        countT(area)=length(area);
    end
    
    maxCountG=max(countG,[],1);
    maxCountT=max(countT,[],1);
    
    plotMatG((maxCountG>0 & maxCountG<=thresh12(1)),2,1)=2;
    plotMatG((maxCountG>thresh12(1) & maxCountG<=thresh12(2)),2,2)=2;
    plotMatG(maxCountG>thresh12(2),2,3)=2;
    
    plotMatT((maxCountT>0 & maxCountT<=thresh12(1)),2,1)=2;
    plotMatT((maxCountT>thresh12(1) & maxCountT<=thresh12(2)),2,2)=2;
    plotMatT(maxCountT>thresh12(2),2,3)=2;
    
    timeSmallG=dataGood.time(:,iaG);
    [~,iaMissG] = setdiff(newTime,timeSmallG);
    
    plotMatG(iaMissG,1,3)=1;
    
    timeSmallT=dataTest.time(:,iaT);
    [~,iaMissT] = setdiff(newTime,timeSmallT);
    
    plotMatT(iaMissT,1,3)=1;
    
    allDataG=any(~isnan(newDataG),1);
    sumGood=cat(1,sum(double(allDataG)),sumGood);
    
    allDataT=any(~isnan(newDataT),1);
    sumTest=cat(1,sum(double(allDataT)),sumTest);
    
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
    ylabel('Number of rays');
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
    
    for jj=1:size(plotMatT,2)
        hold on
        s1=scatter(newTime,plotMatT(:,jj,1),'g+');
        s2=scatter(newTime,plotMatT(:,jj,2),'b+');
        if jj~=1
            s3=scatter(newTime,plotMatT(:,jj,3),'r+');
        else
            s4=scatter(newTime,plotMatT(:,jj,3),'m+');
        end
    end
    
    xlim([newTime(1) newTime(end)]);
    ylim([0 size(plotMatT,2)+1]);
    grid on
    
    yticks(1:size(plotMatT,2));
    yticklabels(makeLabels2);
    
    title(['Missing contiguous points ',qualityTest,' ',freqTest]);
    grid on
    
    ax3=subplot(3,1,3);
    ax3.Position = [0.1300    0.05    0.7750    0.25];
    hold on
    
    for jj=1:size(plotMatG,2)
        hold on
        s1=scatter(newTime,plotMatG(:,jj,1),'g+');
        s2=scatter(newTime,plotMatG(:,jj,2),'b+');
        if jj~=1
            s3=scatter(newTime,plotMatG(:,jj,3),'r+');
        else
            s4=scatter(newTime,plotMatG(:,jj,3),'m+');
        end
    end
    
    xlim([newTime(1) newTime(end)]);
    ylim([0 size(plotMatG,2)+1]);
    grid on
    
    yticks(1:size(plotMatG,2));
    yticklabels(makeLabels2);
    
    legend([s1 s2 s3 s4],{['<=',num2str(thresh12(1))],['>',num2str(thresh12(1)),', <=',num2str(thresh12(2))],...
        ['>',num2str(thresh12(2))],'Dropouts'},'location','best');
    
    title(['Missing contiguous points ',qualityGood,' ',freqGood]);
    grid on
    
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_RF',num2str(ii),'_dataCheck_',qualityTest,'_',freqTest,'_vs_',qualityGood,'_',freqGood,'.png'],'-dpng','-r0')
    
end
