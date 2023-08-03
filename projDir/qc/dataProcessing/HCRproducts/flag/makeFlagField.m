% Analyze HCR clouds

clear all;
close all;

% startTime=datetime(2015,7,9,18,45,0);
% endTime=datetime(2015,7,9,18,51,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='spicule'; %socrates, aristo, cset
quality='qc1'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz, or 2hz
qcVersion='v1.1';

showPlot='off';

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

indir=HCRdir(project,quality,qcVersion,freqData);

figdir=[indir(1:end-5),'flagPlots/cases/'];

mkdir(figdir);

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
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    % Load data
    data=read_HCR(fileList,data,startTime,endTime);

    %% Mask data

    [maskData,antStat]=echoMask(data);

    refl=data.DBZ;
    refl(maskData>1)=nan;

    maskPlot=maskData;
    maskPlot(maskData==0)=nan;

    data.FLAG=maskData;

    %% Plot
    ylimits=[-0.5 7];

    timeInds=1:round(length(data.time)/2000):length(data.time);

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

    disp('Plotting ...');

    figure('DefaultAxesFontSize',11,'position',[1,100,1200,1200],'renderer','Zbuffer','Visible',showPlot);

    ax1=subplot(3,1,1);
    fig1=surf(data.time(timeInds),data.asl(:,timeInds)./1000,data.DBZ(:,timeInds),'edgecolor','none');
    view(2);
    fig1=colMapDBZ(fig1);
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('Reflectivity (all data)')

    ax3=subplot(3,1,3);
    fig3=surf(data.time(timeInds),data.asl(:,timeInds)./1000,refl(:,timeInds),'edgecolor','none');
    view(2);
    fig3=colMapDBZ(fig3);
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('Reflectivity (cloud only, i.e. DBZ(FLAG>1)=NAN)')

    ax2=subplot(3,1,2);
    fig2=surf(data.time(timeInds),data.asl(:,timeInds)./1000,maskPlot(:,timeInds),'edgecolor','none');
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

    ax3pos=ax3.Position;

    formatOut = 'yyyymmdd_HHMM'; set(gcf,'PaperPositionMode','auto')
    print([figdir,datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_echoMask'],...
        '-dpng','-r0');

    %% Plot antenna status
    figure('DefaultAxesFontSize',11,'position',[1,100,1200,700],'renderer','painters','Visible',showPlot);

    s1=subplot(2,1,1);
    plot(data.time,antStat,'linewidth',2)
    xlim([data.time(1),data.time(end)]);
    title('Antenna status')
    yticks(1:6)
    yticklabels({'Down (1)','Up (2)','Pointing (3)','Scanning (4)','Transition (5)','Failure(6)'})
    ylim([0 7])

    s1pos=s1.Position;

    s1.Position=[s1pos(1),s1pos(2),ax3pos(3),s1pos(4)];
    grid on
    box on

    s2=subplot(2,1,2);
    plot(data.time,data.elevation,'linewidth',1.5)
    xlim([data.time(1),data.time(end)]);
    title('Elevation angle')
    ylabel('Elevation (deg)')
    ylims=s2.YLim;
    ylim([ylims(1)-2,ylims(2)+2]);
    
    s2pos=s2.Position;

    s2.Position=[s2pos(1),s2pos(2),ax3pos(3),s2pos(4)];
    grid on
    box on

    set(gcf,'PaperPositionMode','auto')
    print([figdir,datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_antStat'],...
        '-dpng','-r0');
end