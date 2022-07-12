% Calculate air velocity from VEL and DBZ
% https://doi.org/10.1175/2009JAS3132.1

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='spicule'; %socrates, aristo, cset, otrec
quality='qc1'; %field, qc1, or qc2
freqData='10hz';
qcVersion='v1.1';

showPlot='on';
ylimUpper=7.5;
ylimLower=-0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

dataDir=HCRdir(project,quality,qcVersion,freqData);

figdir=[dataDir(1:end-5),'figsAirVel/dbz/'];

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/airVel_',project,'.txt'];

% Loop through cases

caseList=readtable(casefile);
caseStart=datetime(caseList.Var1,caseList.Var2,caseList.Var3, ...
    caseList.Var4,caseList.Var5,caseList.Var6);
caseEnd=datetime(caseList.Var7,caseList.Var8,caseList.Var9, ...
    caseList.Var10,caseList.Var11,caseList.Var12);

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
    data.TEMP=[];
    data.PRESS=[];
    data.RH=[];
    data.MELTING_LAYER=[];
    data.ICING_LEVEL=[];
    data.ECHO_TYPE_2D=[];
                
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

    %ylimUpper=(max(data.asl(~isnan(data.DBZ_MASKED)))./1000)+0.5;

    % Take care of up pointing VEL
    data.VEL_MASKED(:,data.elevation>0)=-data.VEL_MASKED(:,data.elevation>0);

    %% Get unweighted fall speed

    [vr,vi]=getInitFallSpeed(data);

    %% Get weights

    [wr,wi]=getFallSpeedWeights(data.ECHO_TYPE_2D,data.ICING_LEVEL,data.MELTING_LAYER,data.TEMP,data.asl);

    %% Fall speed
    vf=wr.*vr+wi.*vi;

    %% Air vel

    airVel=data.VEL_MASKED-vf;

    %% Plot vel
    close all

    %data.VEL_MASKED(:,data.elevation>0)=-data.VEL_MASKED(:,data.elevation>0);

    f1 = figure('Position',[200 500 2300 900],'DefaultAxesFontSize',12,'visible',showPlot);
    colM=colormap(velCols);
    colormap(colM);

    s1=subplot(2,2,1);
    surf(data.time,data.asl./1000,data.DBZ_MASKED,'edgecolor','none');
    view(2);
    ylabel('Range (km)');
    caxis([-40 20]);
    ylim([ylimLower ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    colormap(s1,jet);
    colorbar
    grid on
    box on
    title('Reflectivity (dBZ)');

    s2=subplot(2,2,2);
    surf(data.time,data.asl./1000,data.VEL_MASKED,'edgecolor','none');
    view(2);
    ylabel('Range (km)');
    caxis([-16 16]);
    ylim([ylimLower ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    colorbar
    grid on
    box on
    title('Velocity (m s^{-1})')

    s3=subplot(2,2,3);
    surf(data.time,data.asl./1000,vf,'edgecolor','none');
    view(2);
    ylabel('Range (km)');
    caxis([-16 16]);
    ylim([ylimLower ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    colorbar
    grid on
    box on
    title('Particle fall speed (m s^{-1})')

    s4=subplot(2,2,4);
    surf(data.time,data.asl./1000,airVel,'edgecolor','none');
    view(2);
    ylabel('Range (km)');
    caxis([-16 16]);
    ylim([ylimLower ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    colorbar
    grid on
    box on
    title('Air vel (m s^{-1})')

    linkaxes([s1 s2 s3 s4],'xy')

    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_fallSpeed_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');

end