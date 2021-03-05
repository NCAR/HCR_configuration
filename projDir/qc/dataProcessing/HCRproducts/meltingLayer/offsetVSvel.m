% Analyze HCR clouds

clear all;
close all;

project='socrates'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz, 2hz, or combined

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

% figdir=['/scr/snow1/rsfdata/projects/otrec/hcr/qc2/cfradial/final2/10hz/plots/testHourly/'];
figdir=['/home/romatsch/plots/HCR/meltingLayer/offsetVel/'];

if ~exist(figdir, 'dir')
    mkdir(figdir)
end

ylimits=[-0.2 7];

%indir=HCRdir(project,quality,freqData);
indir=['/run/media/romatsch/RSF0006/rsf/meltingLayer/',project,'/',freqData,'/'];

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

velAbove=[];
velBelow=[];
offset=[];

for aa=1:size(caseList,1)
    disp(['Flight ',num2str(aa)]);
    
    startTime=datetime(caseList(aa,1:6));
    endTime=startTime;
    
    endTimeIn=datetime(caseList(aa,7:12));
    
    while endTime<endTimeIn
        endTime=endTime+hours(1);
        
        disp([datestr(startTime,'yyyy-mm-dd HH:MM'),' to ',datestr(endTime,'yyyy-mm-dd HH:MM')]);
        data=[];
        
        %% Load data
        
        disp('Loading data ...');
        
        if strcmp(freqData,'combined')
            data.HCR_VEL=[];
        else
            data.VEL_CORR=[];
        end
        data.MELTING_LAYER=[];
        %data.ICING_LEVEL=[];
        data.FLAG=[];
        
        dataVars=fieldnames(data);
        
        % Make list of files within the specified time frame
        fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
        
        if length(fileList)==0
            disp('No data files found.');
            continue
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
        
        if strcmp(freqData,'combined')
            data.VEL_CORR=data.HCR_VEL;
            data=rmfield(data,'HCR_VEL');
        end
                
        %% Get indices
        
%         elevenInds=find(data.MELTING_LAYER==11);
%         twelveInds=find(data.MELTING_LAYER==12);
%         thirteenInds=find(data.MELTING_LAYER==13);
%         fourteenInds=find(data.MELTING_LAYER==14);
%         
%         twentyoneInds=find(data.MELTING_LAYER==21);
%         twentytwoInds=find(data.MELTING_LAYER==22);
%         twentythreeInds=find(data.MELTING_LAYER==23);
%         twentyfourInds=find(data.MELTING_LAYER==24);
        
        %% Calculate offset
        
        % Detected
        altDet=data.asl;
        altDet(data.MELTING_LAYER~=12)=nan;
        
        [minDet,rowDet]=min(altDet,[],1);
        rowDet(isnan(minDet))=[];
        minDet(isnan(minDet))=[];
        colDet=find(sum(~isnan(altDet),1)>0);
        
        if isempty(minDet)
            startTime=endTime;
            continue
        end
        
        % 0 deg
        
        alt0=data.asl;
        alt0(find(data.MELTING_LAYER~=11 & data.MELTING_LAYER~=21))=nan;
        
        alt0det=alt0(:,colDet);
        minAlt0=min(alt0det,[],1,'omitnan');
        
        offsetTime=minAlt0-minDet;
        
        rowDet(isnan(offsetTime))=[];
        colDet(isnan(offsetTime))=[];
        offsetTime(isnan(offsetTime))=[];
        
        %% Velocity
        
        vel=data.VEL_CORR;
        vel(data.elevation>0)=-vel(data.elevation>0);
        vel(data.FLAG>1)=nan;
        
        vel=vel(:,colDet);
        
        velAboveTime=nan(size(colDet));
        velBelowTime=nan(size(colDet));
        
        elev=data.elevation(colDet);
        
        for jj=1:length(colDet)
            if elev(jj)>0 % up
                vel10b=vel(rowDet(jj)-10:rowDet(jj)-1,jj);
                vel10a=vel(rowDet(jj)+1:rowDet(jj)+10,jj);
            elseif elev(jj)<0 % down
                vel10a=vel(rowDet(jj)-10:rowDet(jj)-1,jj);
                vel10b=vel(rowDet(jj)+1:rowDet(jj)+10,jj);
            end
            velAboveTime(jj)=mean(vel10a,'omitnan');
            velBelowTime(jj)=mean(vel10b,'omitnan');
        end
        
        velAbove=cat(1,velAbove,velAboveTime');
        velBelow=cat(1,velBelow,velBelowTime');
        offset=cat(1,offset,offsetTime');
        
        startTime=endTime;
    end
end
%% Plot

close all
edges={-10:0.5:10 -100:20:1000};

N=hist3(cat(2,velAbove,offset),'Edges',edges);
%N(N==0)=nan;

N2=hist3(cat(2,velBelow,offset),'Edges',edges);
%N2(N2==0)=nan;

f1 = figure('Position',[200 500 1300 600],'DefaultAxesFontSize',12);
colormap jet

s1=subplot(1,2,1);

hold on
%surf(edges{1},edges{2},log10(N'),'edgecolor','none')
surf(edges{1},edges{2},N','edgecolor','none')
view(2)

%axis equal
xlim([0,8])
ylim([-100,700])
caxis([0 3000])
%xticks(-40:20:60);
%yticks(-40:20:60);

grid on
xlabel('Velocity above melting layer (m s^{-1})');
ylabel('Offset (m)');
title(['Velocity vs offset'])
s1pos=s1.Position;

% Regression
% fitOrth=gmregress(compTablePlot.DBZsur,compTablePlot.DBZrhi,1);
% fitAll=[fitOrth(2) fitOrth(1)];
% xFit = -40:0.1:60;
% yFit = polyval(fitAll, xFit);
% 
% plot(xFit, yFit,'-b','linewidth',2);
% 
% ax2.SortMethod='childorder';

s2=subplot(1,2,2);

hold on
%surf(edges{1},edges{2},log10(N'),'edgecolor','none')
surf(edges{1},edges{2},N2','edgecolor','none')
view(2)

xlim([0,8])
ylim([-100,700])
caxis([0 3000])
colorbar

grid on
xlabel('Velocity above melting layer (m s^{-1})');
ylabel('Offset (m)');
title(['Velocity vs offset'])
s2pos=s2.Position;
s2.Position=[s2pos(1) s1pos(2) s1pos(3) s1pos(4)];

formatOut = 'yyyymmdd_HHMM';
set(gcf,'PaperPositionMode','auto')
print([figdir,'offsetVSvel_',project],'-dpng','-r0');

save([figdir,project,'_velOffset.mat'],'velAbove','velBelow','offset');