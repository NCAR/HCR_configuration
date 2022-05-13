% Analyze HCR clouds

clear all;
close all;

startTime=datetime(2019,8,25,17,40,0);
endTime=datetime(2019,8,25,17,45,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='otrec'; %socrates, aristo, cset
quality='qc3'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz, or 2hz
qcVersion='v3.0';

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir=['/scr/sleet2/rsfdata/projects/otrec/hcr/qc3/cfradial/v3.0_full/flagPlots/'];

if ~exist(figdir, 'dir')
    mkdir(figdir)
end

directories.dataDir=HCRdir(project,quality,qcVersion,freqData);

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/flag_',project,'.txt'];

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

    data=[];

    data.DBZ=[];
    data.WIDTH=[];
    data.DBMVC=[];
    data.TOPO=[];

    dataVars=fieldnames(data);

    % Make list of files within the specified time frame
    fileList=makeFileList(directories.dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    if length(fileList)==0
        disp('No data files found.');
        return
    end

    % Load data
    data=read_HCR(fileList,data,startTime,endTime);

    % Check if all variables were found
    for ii=1:length(dataVars)
        if ~isfield(data,dataVars{ii})
            dataVars{ii}=[];
        end
    end

    dataVars=dataVars(~cellfun('isempty',dataVars));

    %% Mask data

    [maskData antStat]=echoMask(data);

    refl=data.DBZ;
    refl(maskData>1)=nan;

    maskPlot=maskData;
    maskPlot(maskData==0)=nan;

    data.FLAG=maskData;

    %% Plot
    ylimits=[-0.5 15];

    ytickLabels={'Cloud (1)';'Speckle (2)';'Extinct (3)';'Backlobe (4)';'Out of range (5)';...
        'Bang (6)';'Water (7)';'Land (8)';'Below surf. (9)';...
        'NS cal (10)';'Missing (11)'};

    colMask=[0.4,0.8,1;
        0,0,0;
        0.5,0,0.5;
        0,1,0;
        0.2,0.6,0.2;
        1,0,0;
        0,0,0.6;
        0.7065,0.4415,0.2812;
        0.5,0.5,0.5;
        0.9290,0.8940,0.1250;
        1,0.6,0];

    close all

    figure('DefaultAxesFontSize',11,'position',[1,100,1800,1200],'renderer','Zbuffer');

    ax1=subplot(3,1,1);
    fig1=surf(data.time,data.asl./1000,data.DBZ,'edgecolor','none');
    view(2);
    fig1=colMapDBZ(fig1);
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('Reflectivity (all data)')

    ax3=subplot(3,1,3);
    fig3=surf(data.time,data.asl./1000,refl,'edgecolor','none');
    view(2);
    fig3=colMapDBZ(fig3);
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('Reflectivity (cloud only, i.e. DBZ(FLAG>1)=NAN)')

    ax2=subplot(3,1,2);
    fig2=surf(data.time,data.asl./1000,maskPlot,'edgecolor','none');
    view(2);
    caxis([1 11]);
    colormap(ax2,colMask);
    hcb=colorbar;
    set(hcb,'ytick',[1.5:0.9:11.5]);
    set(hcb,'YTickLabel',ytickLabels);
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('Flag field')

    linkaxes([ax1,ax2,ax3],'xy');

    formatOut = 'yyyymmdd_HHMM'; set(gcf,'PaperPositionMode','auto')
    print([figdir,datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_echoMask'],...
        '-dpng','-r0');

    %% Plot antenna status
    figure('DefaultAxesFontSize',11,'position',[1,100,1800,300],'renderer','painters');

    plot(data.time,antStat,'linewidth',2)
    xlim([data.time(1),data.time(end)]);
    title('Antenna status')
    yticks(1:6)
    yticklabels({'Down (1)','Up (2)','Pointing (3)','Scanning (4)','Transition (5)','Failure(6)'})
    ylim([0 7])

    set(gcf,'PaperPositionMode','auto')
    print([figdir,datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_antStat'],...
        '-dpng','-r0');
end