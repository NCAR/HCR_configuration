
% Calculate liquid water content from HCR ocean return

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

showPlot='on';

meltAlt=4.7; % Estimated altitude of melting layer in km
divAlt=8; % Estimated altitude of divergence level in km

dataDir='/scr/virga1/rsfdata/projects/nasa-apr3/data/3D/';

figdir='/scr/virga1/rsfdata/projects/nasa-apr3/ecco-v-Figs/';

casefile='eccoCases_apr3.txt';

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

    %% Read data

    disp("Getting data ...");

    fileList=makeFileList_apr3(dataDir,startTime,endTime,'xxxxxxxxxxxxxxxxxxxxxxxxxxxxx20YYMMDDxhhmmss',1);

    data=[];

    data.lores_zhh14=[];
    data.lores_vel14c=[];
    data.lores_Topo_Hm=[];

    % Load data
    data=read_apr3_3D(fileList,data,startTime,endTime);

    %% Prepare data

    % Remove surface echo
    data.DBZ(data.asl<100)=nan;

    % Create melting layer
    data.MELTING_LAYER=nan(size(data.DBZ));
    data.MELTING_LAYER(data.asl>=meltAlt.*1000)=20;
    data.MELTING_LAYER(data.asl<meltAlt.*1000)=10;

    % Create fake temperature profile
    data.TEMP=nan(size(data.DBZ));
    data.TEMP(data.asl>=divAlt.*1000)=-30;
    data.TEMP(data.asl<divAlt.*1000)=10;

    %% Texture from reflectivity and velocity

    disp('Calculating reflectivity texture ...');

    pixRadDBZ=5; % Radius over which texture is calculated in pixels. Default is 50.
    dbzBase=0; % Reflectivity base value which is subtracted from DBZ.

    dbzText=f_reflTexture(data.DBZ,pixRadDBZ,dbzBase);

    %% Convectivity

    % Convectivity
    upperLimDBZ=30; % 20
    convDBZ=1/upperLimDBZ.*dbzText;

    %% Basic classification

    disp('Basic classification ...');

    stratMixed=0.4; % Convectivity boundary between strat and mixed.
    mixedConv=0.5; % Convectivity boundary between mixed and conv.

    classBasic=f_classBasic(convDBZ,stratMixed,mixedConv,data.MELTING_LAYER);

    %% Sub classification

    disp('Sub classification ...');

    classSub=f_classSub(classBasic,data.asl,data.TOPO,data.MELTING_LAYER,data.TEMP);

    %% Plot strat conv

    disp('Plotting ...');

    classSubPlot=classSub;
    classSubPlot(classSub==14)=1;
    classSubPlot(classSub==16)=2;
    classSubPlot(classSub==18)=3;
    classSubPlot(classSub==25)=4;
    classSubPlot(classSub==30)=5;
    classSubPlot(classSub==32)=6;
    classSubPlot(classSub==34)=7;
    classSubPlot(classSub==36)=8;
    classSubPlot(classSub==38)=9;

    % 1D
    stratConv1D=max(classSubPlot,[],1);
    time1D=data.time(~isnan(stratConv1D));
    stratConv1D=stratConv1D(~isnan(stratConv1D));

    colmapSC=[0,0.1,0.6;
        0.38,0.42,0.96;
        0.65,0.74,0.86;
        0.32,0.78,0.59;
        1,0,0;
        1,0,1;
        1,1,0;
        0.99,0.77,0.22;
        0.7,0,0];

    col1D=colmapSC(stratConv1D,:);

    ylimUpper=(max(data.asl(~isnan(data.DBZ)))./1000)+0.5;
    textAlt=ylimUpper-1;
    textDate=data.time(1)+seconds(5);

    close all

    f1 = figure('Position',[200 500 1600 1100],'DefaultAxesFontSize',12,'visible',showPlot);

    colormap('jet');

    s1=subplot(4,1,1);

    hold on
    surf(data.time,data.asl./1000,data.DBZ,'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    clim([-10 60]);
    ylim([0 ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    set(gca,'XTickLabel',[]);
    cb1=colorbar;
    grid on
    box on

    text(textDate,textAlt,'(a) Reflectivity (dBZ)','FontSize',11,'FontWeight','bold');

    s2=subplot(4,1,2);

    hold on
    surf(data.time,data.asl./1000,convDBZ,'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    clim([0 1]);
    ylim([0 ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    cb2=colorbar;
    set(gca,'XTickLabel',[]);
    grid on
    box on
    text(textDate,textAlt,'(b) Convectivity','FontSize',11,'FontWeight','bold');

    s5=subplot(30,1,30);

    hold on
    scat1=scatter(time1D,ones(size(time1D)),10,col1D,'filled');
    set(gca,'clim',[0,1]);
    set(gca,'YTickLabel',[]);
    s5.Colormap=colmapSC;
    xlim([data.time(1),data.time(end)]);
    grid on
    box on

    s4=subplot(4,1,3);

    hold on
    surf(data.time,data.asl./1000,classSubPlot,'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    clim([0 10]);
    ylim([0 ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    s4.Colormap=colmapSC;
    clim([0.5 9.5]);
    cb4=colorbar;
    cb4.Ticks=1:9;
    cb4.TickLabels={'Strat Low','Strat Mid','Strat High','Mixed',...
        'Conv','Conv Elev','Conv Shallow','Conv Mid','Conv Deep'};
    set(gca,'XTickLabel',[]);
    grid on
    box on
    text(textDate,textAlt,'(c) Echo type','FontSize',11,'FontWeight','bold');

    s1.Position=[0.049 0.69 0.82 0.29];
    s2.Position=[0.049 0.38 0.82 0.29];
    s4.Position=[0.049 0.07 0.82 0.29];
    s5.Position=[0.049 0.035 0.82 0.02];

    cb1.Position=[0.875,0.69,0.02,0.29];
    cb2.Position=[0.875,0.38,0.02,0.29];
    cb4.Position=[0.875,0.07,0.02,0.29];

    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,'apr3_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS'),'.png'],'-dpng','-r0')
end