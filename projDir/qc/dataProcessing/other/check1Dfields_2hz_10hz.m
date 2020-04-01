% Calculate matrix output of Lee et al. 1994
% i.e. elevation and azimuth angles in earth coordinates

clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='socrates'; % socrates, cset, aristo, otrec
quality='qc2'; % field, qc0, qc1, qc2
freq10='10hz';
freq2='2hz';

figdir=['/h/eol/romatsch/hcrCalib/check1Dfields/'];

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'.txt'];

caseList = table2array(readtable(infile));

indir10=HCRdir(project,quality,freq10);
indir2=HCRdir(project,quality,freq2);

%% Run processing

% Go through flights
for ii=1:size(caseList,1)
    disp(['Flight ',num2str(ii),' of ',num2str(size(caseList,1))]);
    
    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));
    
    % Desired variables. The variable name comies after the . and must be spelled exactly
    % like in the CfRadial file
    data10=[];
    data2=[];
       
    data10.SST=[];
    data10.TOPO=[];
    data10.U_SURF=[];
    data10.V_SURF=[];
    data10.ANTFLAG=[];
    
    dataVars10=fieldnames(data10);
    
    data2.SST=[];
    data2.TOPO=[];
    data2.U_SURF=[];
    data2.V_SURF=[];
    data2.ANTFLAG=[];
    
    dataVars2=fieldnames(data2);
    
    %% Load 10hz data
    disp('Loading 10hz data.');
    
    fileList10=makeFileList(indir10,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if length(fileList10)==0
        disp('No data 10hz files found.');
        startTime=endTime;
        continue
    end
    
    % Load data
    data10=read_HCR(fileList10,data10,startTime,endTime);
    
    if isempty(data10.time)
        disp('No 10hz data found.');
        startTime=endTime;
        continue
    end
    
    % Check if all variables were found
    for kk=1:length(dataVars10)
        if ~isfield(data10,dataVars10{kk})
            dataVars10{kk}=[];
        end
    end
       
    %% 2hz data
    disp('Loading 2hz data.');
    
    fileList2=makeFileList(indir2,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if length(fileList2)==0
        disp('No data 2hz files found.');
        startTime=endTime;
        continue
    end
    
    % Load data
    data2=read_HCR(fileList2,data2,startTime,endTime);
    
    if isempty(data2.time)
        disp('No 2hz data found.');
        startTime=endTime;
        continue
    end
    
    % Check if all variables were found
    for kk=1:length(dataVars2)
        if ~isfield(data2,dataVars2{kk})
            dataVars2{kk}=[];
        end
    end
    %% Plot
   
    close all
    
    f1=figure('DefaultAxesFontSize',12,'renderer','painters');
    set(f1,'Position',[200 500 1000 1000]);
    
    subplot(2,1,1)
    hold on
    plot(data10.time,data10.TOPO,'-g','linewidth',2);
    plot(data2.time,data2.TOPO,'-k','linewidth',1);
    ylabel('TOPO (m)');
       
    xlim([data10.time(1) data10.time(end)]);
    
    yyaxis right
    plot(data10.time,data10.SST,'-b','linewidth',2);
    plot(data2.time,data2.SST,'-r','linewidth',1);
    plot(data10.time,sqrt(data10.U_SURF.^2+data10.V_SURF.^2),'-m','linewidth',2);
    plot(data2.time,sqrt(data2.U_SURF.^2+data2.V_SURF.^2),'-c','linewidth',1);
    ylabel('SST (C), Wind speed (m/s)');
    ax = gca;
    ax.YAxis(2).Color = 'k';

    legend('TOPO 10hz','TOPO 2hz','SST 10hz','SST 2hz','Wind 10hz','Wind 2hz');        
    title([project,' RF ',num2str(ii),' TOPO, SST and wind speed']);
    
    subplot(2,1,2)
    hold on
    
    plot(data10.time,data10.ANTFLAG,'-c','linewidth',2);
    plot(data2.time,data2.ANTFLAG,'-m','linewidth',1);
    ylabel('Antenna flag');
   
    xlim([data10.time(1) data10.time(end)]);
    ylim([-1 5])
    grid on

    legend('ANTFLAG 10hz','ANTFLAG 2hz');        
    title('ANTFLAG');
            
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_RF',num2str(ii),'_1Dfields'],'-dpng','-r0')
   
end
