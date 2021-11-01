% Analyze HCR clouds

clear all;
close all;

project='otrec'; %socrates, otrec, cset
quality='qc3'; %field, qc1, or qc2
qcVersion='v3.0';
freqData='10hz'; % 10hz, 100hz
whichModel='era5';

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

ylimits=[-0.2 10];

indir=HCRdir(project,quality,qcVersion,freqData);

[~,directories.modeldir]=modelDir(project,whichModel,quality,qcVersion,freqData);
outdir=directories.modeldir;

figdir=[indir(1:end-5),'meltLayerPlots/process/'];

saveTime=0;
saveOffset=1;

% Default params
% !!! DO NOT CHANGE !!!
params.fixedOffsetFlight=0; % If fixed offset in should be changed for each flight, set to 1. Will need input file with offsets.
params.offsetIn=-200; % Not needed if fixedOffsetFlight=1

% Fixed offset
params.adjustOffset=0; % Set to zero if fixed offset should be used and not from previous flight.

% If no data is found within one flight, take mean of previous flight (0)
% or fixed offset
params.prevFixed=0; % Only needed if adjustOffset=1

params.LDRlimits=[-16,-7]; % SOCRATES, OTREC, CSET default: [-16,-7]
params.LDRareaPix=[];
params.LDRspeclePix=[]; % SOCRATES, OTREC, CSET default: [] (not used)
params.LDRsolidity=[]; % SOCRATES, OTREC, CSET default: [] (not used)
params.LDRsearchPix=18; % SOCRATES, OTREC, CSET default: 18
params.LDRstd=100; % SOCRATES, OTREC, CSET default: 100
params.LDRaltDiff=50; % SOCRATES, OTREC, CSET default: 50

params.VELsearchPix=50; % SOCRATES, OTREC, CSET default: 50
params.VELstd=35; % SOCRATES, OTREC, CSET default: 35
params.VELaltDiff=100; % SOCRATES, OTREC, CSET default: 100
params.VELudDiff=-0.7; % SOCRATES, OTREC, CSET default: -0.7
params.VEL_LDRdiff=200; % SOCRATES, OTREC, CSET default: 200

params.outlier=50; % SOCRATES, OTREC, CSET default: 50
params.length=20; % SOCRATES, OTREC, CSET default: 20

% Read params file
try
    run([project,'_params.m']);
    paramFields=fields(paramsIn);
    
    for ii=1:length(paramFields)
        params.(paramFields{ii})=paramsIn.(paramFields{ii});
    end
catch
    warning('No parameter file found. Using defaults.');
end

if ~exist(figdir, 'dir')
    mkdir(figdir)
end

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

% Get fixed offset
if params.fixedOffsetFlight
    offsetFixed=readtable(['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/meltingLayer/paramsFiles/',project,'_offsetsIn.txt']);
else
    offsetFixed=params.offsetIn;
end

for aa=1:size(caseList,1)
    disp(['Flight ',num2str(aa)]);
    disp('Loading HCR data.')
    disp(['Starting at ',datestr(datetime('now'),'yyyy-mm-dd HH:MM')]);
    
    clearvars -except project quality freqData whichModel figdir ...
        ylimits indir outdir caseList offsetFixed zeroAdjust aa ...
        saveOffset saveTime params
    
    % Take care of fixed offset
    if size(offsetFixed,1)>1
        thisFixed=offsetFixed.Var2(aa);
    else
        thisFixed=offsetFixed;
    end
        
    if aa==1
        OffsetM=nan(size(caseList,1),1);
        FlightNames=table('Size',[size(caseList,1) 1],'VariableTypes',"string");
        for ii=1:size(caseList,1)
            FlightNames.Var1(ii)=['Flight ',num2str(ii)];
        end
        offsetIn1=array2table(OffsetM);
        offset=cat(2,FlightNames,offsetIn1);
        offset.Properties.VariableNames{'Var1'}='Flight';
        zeroAdjust=thisFixed;
    else
        offset=readtable([outdir,whichModel,'.offset.',project,'.txt']);
        if params.adjustOffset==0
            zeroAdjust=thisFixed;
        elseif isnan(offset.OffsetM(aa-1)) & params.prevFixed==1
            zeroAdjust=thisFixed;
        end
    end
        
    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    disp([datestr(startTime,'yyyy-mm-dd HH:MM'),' to ',datestr(endTime,'yyyy-mm-dd HH:MM')]);
    
    %% Load data
        
    data.DBZ=[];
    data.TEMP=[];
    data.WIDTH=[];
    data.FLAG=[];
    data.TOPO=[];
    
    dataVars=fieldnames(data);
    
    % Make list of files within the specified time frame
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    % Check if VEL_UNFOLDED is available
    try
        velTest=ncread(fileList{1},'VEL_UNFOLDED');
        data.VEL_UNFOLDED=[];
    catch
        data.VEL_CORR=[];
    end
    
    % Check if LDR_MASKED is available
    try
        velTest=ncread(fileList{1},'LDR_MASKED');
        data.LDR_MASKED=[];
    catch
        data.LDR=[];
    end
    
    if length(fileList)==0
        disp('No data files found.');
        continue
    end
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
    
    % Check if all variables were found and convert to single
    for ii=1:length(dataVars)
        if ~isfield(data,dataVars{ii})
            dataVars{ii}=[];
        end
    end
    
    dataVars=dataVars(~cellfun('isempty',dataVars));
    
    if isempty(data.DBZ)
        continue
    end
    
    if isfield(data,'VEL_UNFOLDED')
        data.VEL_CORR=data.VEL_UNFOLDED;
        data=rmfield(data,'VEL_UNFOLDED');
    end
    
    if isfield(data,'LDR_MASKED')
        data.LDR=data.LDR_MASKED;
        data=rmfield(data,'LDR_MASKED');
    end
    
    %% Find melting layer
    
    [meltLayer iceLayer offset.OffsetM(aa)]=f_meltLayer(data,zeroAdjust,params);
    
    if ~isnan(offset.OffsetM(aa))
        zeroAdjust=offset.OffsetM(aa);
    end        
    
    elevenInds=find(meltLayer==11);
    twelveInds=find(meltLayer==12);
    thirteenInds=find(meltLayer==13);
    fourteenInds=find(meltLayer==14);
    
    twentyoneInds=find(meltLayer==21);
    twentytwoInds=find(meltLayer==22);
    twentythreeInds=find(meltLayer==23);
    twentyfourInds=find(meltLayer==24);
    
    %% Plot
    
    timeMat=repmat(data.time,size(data.TEMP,1),1);
    ldrMasked=data.LDR;
    ldrMasked(data.FLAG>1)=nan;
    velMasked=data.VEL_CORR;
    velMasked(data.FLAG>1)=nan;
    
    close all
    
    if etime(datevec(endTime),datevec(startTime))<=900
        newInds=1:1:length(data.time);
    elseif etime(datevec(endTime),datevec(startTime))<=3600
        newInds=1:10:length(data.time);
    else
        newInds=1:100:length(data.time);
    end
    
    % Resample for plotting
    newDBZ=data.DBZ(:,newInds);
    newLDR=data.LDR(:,newInds);
    newVEL=data.VEL_CORR(:,newInds);
    newASL=data.asl(:,newInds);
    newTEMP=data.TEMP(:,newInds);
    newFindMelt=meltLayer(:,newInds);
    newTime=data.time(newInds);
    
    fig1=figure('DefaultAxesFontSize',11,'position',[100,1300,1200,920]);
    
    ax1=subplot(4,1,1);
    hold on;
    sub1=surf(newTime,newASL./1000,newDBZ,'edgecolor','none');
    view(2);
    sub1=colMapDBZ(sub1);
    scatter(timeMat(twentyoneInds),data.asl(twentyoneInds)./1000,10,...
        'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7]);    
    scatter(timeMat(elevenInds),data.asl(elevenInds)./1000,10,'k','filled');
           
    scatter(timeMat(twentyfourInds),data.asl(twentyfourInds)./1000,10,...
        'MarkerEdgeColor',[0.45 0.76 0.42],'MarkerFaceColor',[0.45 0.76 0.42]);
    scatter(timeMat(twentythreeInds),data.asl(twentythreeInds)./1000,10,...
        'MarkerEdgeColor',[0.7 0.8 0.87],'MarkerFaceColor',[0.7 0.8 0.87]);
    scatter(timeMat(twentytwoInds),data.asl(twentytwoInds)./1000,10,...
        'MarkerEdgeColor',[0.17 0.45 0.7],'MarkerFaceColor',[0.17 0.45 0.7]);
        
    scatter(timeMat(fourteenInds),data.asl(fourteenInds)./1000,10,'g','filled');
    scatter(timeMat(thirteenInds),data.asl(thirteenInds)./1000,10,'c','filled');
    scatter(timeMat(twelveInds),data.asl(twelveInds)./1000,10,'b','filled');
    
    ax = gca;
    ax.SortMethod = 'childorder';
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title(['Flight ',num2str(aa),': Reflectivity and melting layer'])
    grid on
    set(gca,'xticklabel',[])
    ax1.Position=[0.06 0.765 0.87 0.21];
    
    ax2=subplot(4,1,2);
    hold on;
    sub1=surf(newTime,newASL./1000,newFindMelt,'edgecolor','none');
    ax2.Colormap=([1 0 1;1 1 0]);
    view(2);
    scatter(timeMat(twentyoneInds),data.asl(twentyoneInds)./1000,10,...
        'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7]);    
    scatter(timeMat(elevenInds),data.asl(elevenInds)./1000,10,'k','filled');
           
    scatter(timeMat(twentyfourInds),data.asl(twentyfourInds)./1000,10,...
        'MarkerEdgeColor',[0.45 0.76 0.42],'MarkerFaceColor',[0.45 0.76 0.42]);
    scatter(timeMat(twentythreeInds),data.asl(twentythreeInds)./1000,10,...
        'MarkerEdgeColor',[0.7 0.8 0.87],'MarkerFaceColor',[0.7 0.8 0.87]);
    scatter(timeMat(twentytwoInds),data.asl(twentytwoInds)./1000,10,...
        'MarkerEdgeColor',[0.17 0.45 0.7],'MarkerFaceColor',[0.17 0.45 0.7]);
        
    scatter(timeMat(fourteenInds),data.asl(fourteenInds)./1000,10,'g','filled');
    scatter(timeMat(thirteenInds),data.asl(thirteenInds)./1000,10,'c','filled');
    scatter(timeMat(twelveInds),data.asl(twelveInds)./1000,10,'b','filled');
    
    plot(data.time,iceLayer./1000,'linewidth',1,'color',[0.6 0.6 0.6]);
    ax = gca;
    ax.SortMethod = 'childorder';
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('Reflectivity and melting layer')
    grid on
    set(gca,'xticklabel',[])
    ax2.Position=[0.06 0.525 0.87 0.21];
    
    %%%%%%%%%%%%%%%%%%%%%%%% LDR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax3=subplot(4,1,3);
    hold on;
    sub3=surf(newTime,newASL./1000,newLDR,'edgecolor','none');
    view(2);
    caxis([-25 5]);
    colorbar
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('LDR')
    grid on
    set(gca,'xticklabel',[])
    ax3.Position=[0.06 0.287 0.87 0.21];
    
    ax4=subplot(4,1,4);
    ax4.Colormap=jet;
    hold on;
    sub3=surf(newTime,newASL./1000,newVEL,'edgecolor','none');
    view(2);
    caxis([-5 5]);
    colorbar
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('VEL')
    grid on
    ax4.Position=[0.06 0.05 0.87 0.21];
    
    linkaxes([ax1 ax2 ax3 ax4],'xy');
    
    set(gcf,'PaperPositionMode','auto')
    print([figdir,'meltRefl_Flight',num2str(aa)],'-dpng','-r0');
    
    disp(['Melting layer is on average ',num2str(offset.OffsetM(aa)),' m from the zero degree isotherm.'])
    
    %% Save
    disp('Saving meltLayer field ...')
    
    save([outdir,whichModel,'.meltLayer.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'meltLayer');
    
    if saveTime
        timeHCR=data.time;
        save([outdir,whichModel,'.time.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
            datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'timeHCR');
    end
    
    save([outdir,whichModel,'.iceLevel.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'iceLayer');
    
    if saveOffset
        writetable(offset,[outdir,whichModel,'.offset.',project,'.txt'],'delimiter',' ');
    end
end

%% Calculate and save mean offset

meanOffset=mean(offset.OffsetM,'omitnan');
stdOffset=std(offset.OffsetM,'omitnan');
disp(['Offset over all flights: ',num2str(meanOffset),' +/- ',num2str(stdOffset),' m']);

offset=cat(1,offset,offset(1,:));
offset(end,1)={'Mean'};
offset.OffsetM(end)=meanOffset;

offset=cat(1,offset,offset(1,:));
offset(end,1)={'StDev'};
offset.OffsetM(end)=stdOffset;

if saveOffset
    writetable(offset,[outdir,whichModel,'.offset.',project,'.txt'],'delimiter',' ');
end