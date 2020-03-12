% Calculate matrix output of Lee et al. 1994
% i.e. elevation and azimuth angles in earth coordinates

clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='cset'; % socrates, cset, aristo, otrec
quality='qc2'; % field, qc0, qc1, qc2
freq='10hz';

figdir=['/h/eol/romatsch/hcrCalib/extinctionStats/'];

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'.txt'];

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,freq);

%% Run processing

cloudSum=[];
extinctSum=[];
clearSum=[];
inCloudSum=[];
totalSum=[];

% Go through flights
for ii=1:size(caseList,1)
    disp(['Flight ',num2str(ii),' of ',num2str(size(caseList,1))]);
    
    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));
    
    % Desired variables. The variable name comies after the . and must be spelled exactly
    % like in the CfRadial file
    if exist('data')
        clear data
    end
   
    data.FLAG=[];
    
    dataVars=fieldnames(data);
    
    %% Load data
    % Make list of files within the specified time frame
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if length(fileList)==0
        disp('No data files found.');
        startTime=endTime;
        continue
    end
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
    
    if isempty(data.time)
        disp('No data found.');
        startTime=endTime;
        continue
    end
    
    % Check if all variables were found
    for kk=1:length(dataVars)
        if ~isfield(data,dataVars{kk})
            dataVars{kk}=[];
        end
    end
       
    %% Calculate extinct percent
    
    cloudDataIn=any(data.FLAG==1,1);
    cloudData=double(cloudDataIn);
    cloudData(cloudData==1)=2;
    extinctData=any(data.FLAG==3,1);
    
    % Calculate in cloud sum
    beyondBang=data.FLAG(18:22,:);
    inCloud3=zeros(size(beyondBang));
    inCloud3(beyondBang==1)=1;
    inCloudAll=sum(inCloud3,1);
    inCloud=zeros(size(data.time));
    inCloud(inCloudAll==5)=1;
    
    cloudSum=[cloudSum,sum(cloudDataIn)];
    extinctSum=[extinctSum,sum(extinctData)];
    clearSum=[clearSum,length(find(cloudData==0))];
    totalSum=[totalSum,length(data.time)];
    inCloudSum=[inCloudSum,sum(inCloud)];
    %% Plot
   
    close all
    
    f1=figure('DefaultAxesFontSize',12);
    set(f1,'Position',[200 500 2000 600]);
    
    subplot(2,1,1)
    hold on
    bar(data.time,cloudData,1,'facecolor',[0.4,0.8,1]);
    bar(data.time,extinctData,1,'facecolor',[0.5,0,0.5]);
       
    ylim([0 2.8]);
    yticks([-1,5]);
    yticklabels('');           
    xlim([data.time(1) data.time(end)]);

    legend('Cloud','Extinct');
        
    text(data.time(1)+minutes(15),2.3,['Extinct/cloud ratio: ',num2str(extinctSum(ii)/cloudSum(ii),2)],'fontsize',12,'fontweight','bold');
    title([project,' RF ',num2str(ii)]);
    
    subplot(2,1,2)
    hold on
    bar(data.time,inCloud,1,'facecolor','g');
           
    ylim([0 1.3]);
    yticks([-1,5]);
    yticklabels('');           
    xlim([data.time(1) data.time(end)]);

    legend('Flying in cloud');
        
    text(data.time(1)+minutes(15),1.1,['In cloud/total cloud ratio: ',num2str(inCloudSum(ii)/cloudSum(ii),2)],'fontsize',12,'fontweight','bold');
        
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_RF',num2str(ii),'_cloudExtinct'],'-dpng','-r0')
   
end
%% Summary plot
close all
f1=figure('DefaultAxesFontSize',12,'renderer','painters');
set(f1,'Position',[200 500 1000 600]);

subplot(2,1,1)
hold on

yyaxis left
plot(extinctSum./cloudSum*100,'color',[0.5,0,0.5],'linewidth',2);

ax = gca;
ax.YColor = [0.5,0,0.5];

ylim([0 20]);
xlim([1 length(cloudSum)]);
xticks(1:length(cloudSum));

xlabel('Flight');
ylabel('(%)');

title([project,': Extinct echo and total cloud percentage'])
text(2,18,['Total extinct/cloud ratio: ',num2str(sum(extinctSum)/sum(cloudSum),2)],'fontsize',12);

yyaxis right
plot(cloudSum./totalSum*100,'color',[0.4,0.8,1],'linewidth',2);
ylim([0 100]);
ax = gca;
ax.YColor = [0.4,0.8,1];
yticks(0:25:100)
grid on

legend('Extinct','Cloud','location','southwest');

subplot(2,1,2)
plot(inCloudSum./cloudSum*100,'color','g','linewidth',2);

ylim([0 100]);
xlim([1 length(cloudSum)]);
xticks(1:length(cloudSum));

xlabel('Flight');
ylabel('(%)');
grid on

title(['Flying in cloud percentage'])
text(2,90,['Total in cloud/cloud ratio: ',num2str(sum(inCloudSum)/sum(cloudSum),2)],'fontsize',12);

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_cloudExtinct'],'-dpng','-r0')