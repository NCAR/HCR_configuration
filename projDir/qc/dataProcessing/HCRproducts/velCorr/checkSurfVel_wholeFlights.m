% Calculate matrix output of Lee et al. 1994
% i.e. elevation and azimuth angles in earth coordinates

clear all;
close all;

addpath(genpath('/h/eol/romatsch/gitPriv/utils/'));

% upwind_limit=2;
% crosswind_limit=2;

% If attitude correction was sometimes turned off we need to filter that
% data
filterAttitudeCorrOff=1;

project='otrec'; % socrates, cset, aristo, otrec
quality='qc1'; % field, qc0, qc1, qc2
freq='10hz';

figdir=['/h/eol/romatsch/hcrCalib/velCorr/',project,'/velFigsCheckSurf/'];
formatOut = 'yyyymmdd_HHMM';

infile=['/h/eol/romatsch/hcrCalib/oceanScans/biasInFiles/flights_',project,'.txt'];

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,freq);

lineCols=lines;

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
    
    data.roll=[];
    data.pitch=[];
    data.drift=[];
    
    data.VEL=[];
    data.VEL_CORR=[];
    data.DBZ=[];
    
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
    
    % Remove zenith pointing
    zenith=find(data.elevation>0);
    
    for kk=1:length(dataVars)
        dataTemp=data.(dataVars{kk});
        dataTemp(:,zenith)=nan;
        data.(dataVars{kk})=dataTemp;
    end
    
     [linInd rowInd rangeToSurf] = hcrSurfInds(data);
     
     velSurf=data.VEL(linInd);
     velCorrSurf=data.VEL_CORR(linInd);

    
    %%
    close all
    
    f1=figure('DefaultAxesFontSize',12);
    set(f1,'Position',[200 500 1200 1300]);
    
    subplot(2,1,1)
    hold on
    yyaxis right
    plot(data.time,data.roll,'-','linewidth',2,'color',lineCols(1,:));
    plot(data.time,data.pitch,'-','linewidth',2,'color',lineCols(2,:));
    plot(data.time,data.drift,'-','linewidth',2,'color',lineCols(3,:));
           
    ylabel('Angles (deg)');
    
    ylim([-14 14]);
    yticks(-20:4:20);
    ax = gca;
    ax.SortMethod = 'childorder';
    ax.YColor=[0 0 0];
    grid on
    
    yyaxis left
    plot(data.time,movmean(data.elevation+90,100),'-','linewidth',2,'color',lineCols(4,:));
    ylim([-0.175 0.175]);
    yticks(-5:0.05:5);
    ylabel('Elev+90 (deg)');
    ax = gca;
    ax.SortMethod = 'childorder';
    ax.YColor=lineCols(4,:);
    
    xlim([data.time(1) data.time(end)]);

    legend('Elev+90','Roll','Pitch','Drift');
        
    title([project,' RF ',num2str(ii)]);
        
    subplot(2,1,2)
    hold on
    plot(data.time,movmedian(velSurf,100),'linewidth',2);
    plot(data.time,movmedian(velCorrSurf,100),'linewidth',2);
    
    legend('Surf Vel','Surf Vel Corr');
    
    ylabel('Velocity (m/s)');
    
    xlim([data.time(1) data.time(end)]);
    ylim([-1 1]);
    yticks(-10:0.2:10);
    
    grid on
       
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_RF',num2str(ii),'_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
    
end