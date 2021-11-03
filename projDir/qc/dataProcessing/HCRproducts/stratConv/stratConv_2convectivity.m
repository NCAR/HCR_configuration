% Calculate liquid water content from HCR ocean return

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='socrates'; %socrates, aristo, cset, otrec
quality='qc3'; %field, qc1, or qc2
freqData='10hz';
qcVersion='v3.0';

showPlot='off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

dataDir=HCRdir(project,quality,qcVersion,freqData);

% if strcmp(project,'otrec')
%     dataDir='/scr/sleet2/rsfdata/projects/otrec/hcr/qc2/cfradial/development/convStrat/10hz/';
% elseif strcmp(project,'socrates')
%     dataDir='/scr/snow2/rsfdata/projects/socrates/hcr/qc2/cfradial/development/convStrat/10hz/';
% end

figdir=[dataDir(1:end-5),'convStratPlots/cases/'];

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/stratConv_',project,'.txt'];

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
    %% Get data
    
    disp("Getting data ...");

    fileList=makeFileList(dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    data=[];
    
    data.DBZ_MASKED=[];
    data.VEL_MASKED=[];
    data.TOPO=[];
    data.TEMP=[];
    data.MELTING_LAYER=[];
    data.FLAG=[];
        
    dataVars=fieldnames(data);
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
    
    % Check if all variables were found
    for ii=1:length(dataVars)
        if ~isfield(data,dataVars{ii})
            dataVars{ii}=[];
        end
    end
    
    dataVars=dataVars(~cellfun('isempty',dataVars));
    
    ylimUpper=(max(data.asl(~isnan(data.DBZ_MASKED)))./1000)+0.5;

    % Take care of up pointing VEL
    data.VEL_MASKED(:,data.elevation>0)=-data.VEL_MASKED(:,data.elevation>0);

    %% Texture from reflectivity and velocity

    disp('Calculating reflectivity texture ...');

    pixRadDBZ=50; % Radius over which texture is calculated in pixels. Default is 50.
    dbzBase=-10; % Reflectivity base value which is subtracted from DBZ.

    dbzText=f_reflTexture(data.DBZ_MASKED,pixRadDBZ,dbzBase);

    disp('Calculating velocity texture ...');

    pixRadVEL=50;
    velBase=-20; % VEL base value which is subtracted from DBZ.

    velText=f_velTexture(data.VEL_MASKED,data.elevation,pixRadVEL,velBase);

    %% Convectivity

    % Convectivity
    upperLimDBZ=12;%14
    convDBZ=1/upperLimDBZ.*dbzText;

    upperLimVEL=5;%6
    convVEL=1/upperLimVEL.*velText;

    convBoth=convDBZ.*convVEL;
    convBoth(isnan(convBoth))=convDBZ(isnan(convBoth));

    %% Basic classification

    disp('Basic classification ...');

    stratMixed=0.4; % Convectivity boundary between strat and mixed.
    mixedConv=0.5; % Convectivity boundary between mixed and conv.

    classBasic=f_classBasicBoth(convBoth,stratMixed,mixedConv,data.MELTING_LAYER,data.elevation);

    %% Sub classification

    disp('Sub classification ...');

    classSub=f_classSubBoth(classBasic,data.asl,data.TOPO,data.MELTING_LAYER,data.TEMP,data.elevation,data.FLAG);
       
    %% Plot strat conv
    
    disp('Plotting conv/strat ...');
    
    close all
    
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
        
    close all
    
    f1 = figure('Position',[200 500 1500 1100],'DefaultAxesFontSize',12,'visible',showPlot);
    
    s1=subplot(5,1,1);
    
    colormap jet
    
    hold on
    surf(data.time,data.asl./1000,data.DBZ_MASKED,'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    caxis([-35 25]);
    ylim([0 ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    colorbar
    grid on
    title('Reflectivity (dBZ)')
    s1pos=s1.Position;
    s1.Position=[s1pos(1),s1pos(2),s1pos(3),s1pos(4)];
    
    s2=subplot(5,1,2);
    hold on
    surf(data.time,data.asl./1000,data.VEL_MASKED,'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    caxis([-5 5]);
    ylim([0 ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    colorbar
    grid on
    title('Radial velocity')
    s2pos=s2.Position;
    s2.Position=[s2pos(1),s2pos(2),s1pos(3),s2pos(4)];
    
    s3=subplot(5,1,3);
    
    hold on
    surf(data.time,data.asl./1000,convBoth,'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    caxis([0 1]);
    ylim([0 ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    colorbar
    grid on
    title('Convectivity')
    s3pos=s3.Position;
    s3.Position=[s3pos(1),s3pos(2),s1pos(3),s3pos(4)];
    
    s4=subplot(5,1,4);
    
    hold on
    surf(data.time,data.asl./1000,classBasic,'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    caxis([0.5 3.5]);
    colormap(s4,[0,0,1;0.32,0.78,0.59;1 0 0]);
    cb=colorbar;
    cb.Ticks=[1,2,3];
    cb.TickLabels=cat(2,{'Stratiform','Mixed','Convective'});
    ylim([0 ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    grid on
    s5pos=s4.Position;
    s4.Position=[s5pos(1),s5pos(2),s1pos(3),s5pos(4)];
        
    s6=subplot(30,1,30);
    
    hold on
    scat1=scatter(time1D,ones(size(time1D)),10,col1D,'filled');
    set(gca,'clim',[0,1]);
    set(gca,'YTickLabel',[]);
    s6.Colormap=colmapSC;
    xlim([data.time(1),data.time(end)]);
    s6pos=s6.Position;
    s6.Position=[s6pos(1),s6pos(2)-0.023,s1pos(3),s6pos(4)];
    
    s5=subplot(5,1,5);
    
    hold on
    surf(data.time,data.asl./1000,classSubPlot,'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    caxis([0 10]);
    ylim([0 ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    s5.Colormap=colmapSC;
    caxis([0.5 9.5]);
    cb=colorbar;
    cb.Ticks=1:9;
    cb.TickLabels={'Strat Low','Strat Mid','Strat High','Mixed',...
        'Conv','Conv Elev','Conv Shallow','Conv Mid','Conv Deep'};
    set(gca,'XTickLabel',[]);
    grid on
    title('Stratiform/convective partitioning')
    s5pos=s5.Position;
    s5.Position=[s5pos(1),s5pos(2),s1pos(3),s5pos(4)];
    
    linkaxes([s1 s2 s3 s4 s5],'xy');
    
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_convStrat_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
    
    %% Plot convectivities
    
    disp('Plotting convectivities ...');
    
    f1 = figure('Position',[200 500 1500 1200],'DefaultAxesFontSize',12,'visible',showPlot);
    
    s1=subplot(5,1,1);
    
    colormap jet
    
    hold on
    surf(data.time,data.asl./1000,data.DBZ_MASKED,'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    caxis([-35 25]);
    ylim([0 ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    colorbar
    grid on
    title('Reflectivity (dBZ)')
    s1pos=s1.Position;
    s1.Position=[s1pos(1),s1pos(2),s1pos(3),s1pos(4)];
    
    s2=subplot(5,1,2);
    hold on
    surf(data.time,data.asl./1000,convDBZ,'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    caxis([0 1]);
    ylim([0 ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    colorbar
    grid on
    title('Reflectivity convectivity')
    s2pos=s2.Position;
    s2.Position=[s2pos(1),s2pos(2),s1pos(3),s2pos(4)];

    s3=subplot(5,1,3);
    
    hold on
    surf(data.time,data.asl./1000,convBoth,'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    caxis([0 0.7]);
    ylim([0 ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    colorbar
    grid on
    title('Reflectivity convectivity * velocity convectivity');
    s3pos=s3.Position;
    s3.Position=[s3pos(1),s3pos(2),s1pos(3),s3pos(4)];

    s4=subplot(5,1,4);
    
    hold on
    surf(data.time,data.asl./1000,convVEL,'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    caxis([0 1]);
    ylim([0 ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    colorbar
    grid on
    title('Velocity convectivity')
    s4pos=s4.Position;
    s4.Position=[s4pos(1),s4pos(2),s1pos(3),s4pos(4)];

    s5=subplot(5,1,5);
    hold on
    surf(data.time,data.asl./1000,data.VEL_MASKED,'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    caxis([-5 5]);
    ylim([0 ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    colorbar
    grid on
    title('Radial velocity')
    s5pos=s5.Position;
    s5.Position=[s5pos(1),s5pos(2),s1pos(3),s5pos(4)];
    
    linkaxes([s1 s2 s3 s4 s5],'xy');
    
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_convectivity_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
    
end