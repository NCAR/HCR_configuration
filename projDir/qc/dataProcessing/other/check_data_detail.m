% Compare different data sets to make sure the data is complete

clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='socrates'; % socrates, cset, aristo, otrec
qualityGood='qc2'; % field, qc0, qc1, qc2
qualityTest='qc2'; % field, qc0, qc1, qc2
freqGood='10hz'; % 10hz, 2hz, 2hzMerged
freqTest='2hzMerged'; % 10hz, 2hz, 2hzMerged

startTime=datetime(2018,1,16,3,25,0);
endTime=datetime(2018,1,16,3,40,0);

figdir=['/h/eol/romatsch/hcrCalib/checkData/check_detail/'];

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'.txt'];

caseList = table2array(readtable(infile));

indirGood=HCRdir(project,qualityGood,freqGood);
indirTest=HCRdir(project,qualityTest,freqTest);

%% Run processing

compareVars={'DBMVC','HCR_DBMVC',0};

dataVarsGood={};
dataVarsTest={};

for jj=1:size(compareVars,1)
    dataVarsGood{end+1}=compareVars{jj,1};
    dataVarsTest{end+1}=compareVars{jj,2};
end

dataVarsGood=dataVarsGood';
dataVarsTest=dataVarsTest';

%% Load good data
disp('Loading good data.');

fileListGood=makeFileList(indirGood,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

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
    
    TTgood.(dataVarsGood{ll})=double(any(~isnan(dataGood.(dataVarsGood{ll})),1))';
    sumGood(ll)=sum(TTgood.(dataVarsGood{ll}));
    
end

%% Load test data
disp('Loading test data.');

fileListTest=makeFileList(indirTest,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

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
    
    TTtest.(dataVarsTest{ll})=double(any(~isnan(dataTest.(dataVarsTest{ll})),1))';
    sumTest(ll)=sum(TTtest.(dataVarsTest{ll}));
    
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
dataAll=table2array(tableAll(:,2:end));

%% Missing good data
missingGood=nan(size(TTall,1),size(dataVarsGood,1));

for jj=1:size(dataVarsGood,1)
    missingGood(find(dataAll(:,jj)==0 & dataAll(:,jj+size(dataVarsGood,1))==1),jj)=jj+2;
end

missingGood=cat(2,nan(size(TTall,1),2),missingGood);
missingGood(find(isnan(dataAll(:,1))),1)=1;

missingGood(find(isnan(dataAll(:,size(dataVarsGood,1)+1)) & missingGood(:,1)~=1),2)=2;
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
                missingGoodSingle(startInds(kk),jj)=jj;
            elseif endInds(kk)-startInds(kk)==1
                missingGood(startInds(kk):endInds(kk),jj)=nan;
                missingGoodDouble(startInds(kk):endInds(kk),jj)=jj;
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
                missingTestSingle(startInds(kk),jj)=jj;
            elseif endInds(kk)-startInds(kk)==1
                missingTest(startInds(kk):endInds(kk),jj)=nan;
                missingTestDouble(startInds(kk):endInds(kk),jj)=jj;
            end
        end
    end
end

%% Plot

makeLabels=cat(1,{'time'},dataVarsGood);
makeLabels2=cat(1,{'allMissing'},makeLabels);

close all

f1=figure('DefaultAxesFontSize',12,'defaultAxesTickLabelInterpreter','none');
set(f1,'Position',[5 5 1700 1000]);

colormap jet;

ax1=subplot(3,1,1);
ax1.Position = [0.1300    0.7093    0.7750    0.25];

hold on
surf(dataTest.time,dataTest.asl./1000,dataTest.(compareVars{1,2}),'edgecolor','none');
        view(2);
ylabel('Altitude (km)');
c1=colorbar('Location','east');
c1.Position=[0.94    0.7197    0.0127    0.2285];

xlim([TTall.time(1) TTall.time(end)]);
grid on

ylim1=ylim;
caxis1=caxis;

title([compareVars{1,2},' ',qualityTest,' ',freqTest],'interpreter','none');

ax2=subplot(3,1,2);
ax2.Position = [0.1300    0.38    0.7750    0.25];

hold on
surf(dataGood.time,dataGood.asl./1000,dataGood.(compareVars{1,1}),'edgecolor','none');
        view(2);
ylabel('Altitude (km)');
c2=colorbar('Location','east');
c2.Position=[0.94    0.3907    0.0127    0.2285];

xlim([TTall.time(1) TTall.time(end)]);
ylim(ylim1);

caxis(caxis1);
grid on

title([compareVars{1,1},' ',qualityGood,' ',freqGood],'interpreter','none');

ax3=subplot(6,1,5);
ax3.Position = [0.1300    0.2    0.7750    0.1026];
hold on

for jj=1:size(missingTest,2)
    hold on
    s3=scatter(TTall.time,missingTestSingle(:,jj),'g+');
    s2=scatter(TTall.time,missingTestDouble(:,jj),'b+');
    s1=scatter(TTall.time,missingTest(:,jj),'r+');
end

xlim([TTall.time(1) TTall.time(end)]);
ylim([0 size(missingGood,2)+1]);
grid on

yticks(1:size(missingTest,2));
yticklabels(makeLabels2);

title(['Missing times ',qualityTest,' ',freqTest]);
grid on

ax4=subplot(7,1,7);
ax4.Position = [0.1300    0.05    0.7750    0.1026];
hold on

for jj=1:size(missingGood,2)
    hold on
    s3=scatter(TTall.time,missingGoodSingle(:,jj),'g+');
    s2=scatter(TTall.time,missingGoodDouble(:,jj),'b+');
    s1=scatter(TTall.time,missingGood(:,jj),'r+');
end

legend([s1 s2 s3],{'Multi','Double','Single'},'location','best');

xlim([TTall.time(1) TTall.time(end)]);
ylim([0 size(missingGood,2)+1]);
grid on

yticks(1:size(missingGood,2));
yticklabels(makeLabels2);

title(['Missing times ',qualityGood,' ',freqGood]);
grid on

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_',datestr(startTime,'yyyymmdd_HHMMSS'),'_to_',...
    datestr(endTime,'yyyymmdd_HHMMSS'),'_dataCheckDetail_',...
    qualityTest,'_',freqTest,'_vs_',qualityGood,'_',freqGood,'_',compareVars{1,1},'.png'],'-dpng','-r0')

