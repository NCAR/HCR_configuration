% Analyze HCR clouds

clear all;
close all;

% startTime=datetime(2018,2,7,18,0,0);
% endTime=datetime(2018,2,8,12,0,0);

startTime=datetime(2019,10,2,15,0,0);
endTime=datetime(2019,10,2,15,59,0);

plotFields=1;
plotWholeFlight=0;

project='otrec'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz, or 2hz

% Expected bright band altitude. Determines plot zoom.
if strcmp(project,'otrec')
    expBBalt=5;
    ylimits=[-0.2 15];
elseif strcmp(project,'socrates')
    expBBalt=2;
    ylimits=[-0.2 9];
elseif strcmp(project,'otrec')
    expBBalt=5;
    ylimits=[-0.2 9];
end

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir=['/h/eol/romatsch/hcrCalib/clouds/cloudID/',project,'/'];

if ~exist(figdir, 'dir')
    mkdir(figdir)
end

indir=HCRdir(project,quality,freqData);

%% Load data

disp('Loading data ...');

data.DBZ=[];
data.TEMP=[];
data.FLAG=[];

dataVars=fieldnames(data);

% Make list of files within the specified time frame
fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

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

cloudIDout=nan(size(data.DBZ));

%% Mask non cloud

refl=data.DBZ;
refl(data.FLAG>1)=nan;

%% Add extinct back in

flagTemp=data.FLAG;
flagTemp(data.FLAG==8)=7;
surfMask=ones(1,size(flagTemp,2));
surfMask(find(any(data.FLAG==7 | data.FLAG==8,1)))=0;
surfDiff=diff(surfMask);

% Add surface echo back in
holeStart=find(surfDiff==1);
holeStart=holeStart+1;
holeEnd=find(surfDiff==-1);

if holeStart(1)>holeEnd(1)
    holeStart=[1 holeStart];
end
if length(holeStart)~=length(holeEnd)
    holeEnd=[holeEnd,size(flagTemp,2)-1];
end

for ii=1:length(holeStart)
    flagColStart=flagTemp(:,holeStart(ii)-1);
    flagColEnd=flagTemp(:,holeEnd(ii)+1);
    
    minStart=min(find(flagColStart==7));
    minEnd=min(find(flagColEnd==7));
    
    if isempty(minStart)
        minStart=minEnd;
    end
    if isempty(minEnd)
        minEnd=minStart;
    end
    
    newVec=holeStart(ii):holeEnd(ii);
    holeLength=holeEnd(ii)-holeStart(ii);
    
    if holeLength==0
        dataFill=holeStart(ii);
    else        
        x=[holeStart(ii),holeEnd(ii)];
        v=[minStart,minEnd];        
        dataFill=round(interp1(x,v,newVec));
    end
    for jj=1:length(dataFill)
        flagTemp(dataFill(jj):end,newVec(jj))=7;
    end
end

% Fill in extinct echo with median of column above
extInd=find(any(data.FLAG==3,1));

for ii=1:length(extInd)
    reflCol=refl(:,extInd(ii));
    reflMed=median(reflCol,'omitnan');
    
    extCol=flagTemp(:,extInd(ii));
    extData=find(extCol==3);
    refl(extData,extInd(ii))=reflMed;
end

%% Smooth with convolution

% Create mask for convolution
radius=5;
numPix=radius*2+1;
[rr cc] = meshgrid(1:numPix);
cirMask = sqrt((rr-(radius+1)).^2+(cc-(radius+1)).^2)<=radius;
cirMask=double(cirMask);

% Normalize
cirMask=cirMask./(sum(reshape(cirMask,1,[])));

% Convolution
reflConv=nanconv(refl,cirMask);

%% Remove contiguous clouds that are too small

reflLarge=reflConv;

reflMask=zeros(size(reflConv));
reflMask(~isnan(reflConv))=1;

pixCut=5000;
CC = bwconncomp(reflMask);

cloudNum=nan(size(data.DBZ));
countCloud=1;

for ii=1:CC.NumObjects
    area=CC.PixelIdxList{ii};
    if length(area)<=pixCut
        reflLarge(area)=nan;
        cloudIDout(area)=1;
    else
        cloudNum(area)=countCloud;
        countCloud=countCloud+1;
    end
end

%% Split up individual clouds
numMax=max(reshape(cloudNum,1,[]),[],'omitnan');

thresh12=[-25 0];

for ii=1:numMax
    
    cloudInds=find(cloudNum==ii);
    
    if length(cloudInds)>100000
        close all;
        cloudRefl=reflLarge(cloudInds);
        
        reflMapBig=nan(size(cloudNum));
        reflMapBig(cloudInds)=cloudRefl;
        
        [clR clC]=ind2sub(size(cloudNum),cloudInds);
        
        reflMap=reflMapBig(min(clR):max(clR),min(clC):max(clC));
        aslMap=data.asl(min(clR):max(clR),min(clC):max(clC));
        timeMap=data.time(min(clC):max(clC));
        
        BW=zeros(size(reflMap));
        BW(reflMap>=thresh12(1) & reflMap<thresh12(2))=1;
                
        [B,L,N,A] = bwboundaries(BW);
        
        cutBound=100;
        goodB={};
        figure; imshow(BW); hold on;
        % Loop through object boundaries
        for k = 1:N
            % Boundary k is the parent of a hole if the k-th column
            % of the adjacency matrix A contains a non-zero element
            if (nnz(A(:,k)) > 0)
                %                 boundary = B{k};
                %                 plot(boundary(:,2),...
                %                     boundary(:,1),'r','LineWidth',2);
                % Loop through the children of boundary k
                for l = find(A(:,k))'
                    boundary = B{l};
                    if size(boundary,1)>cutBound
                        goodB{end+1}=boundary;
                        
                        plot(boundary(:,2),...
                            boundary(:,1),'g','LineWidth',2);
                    else
                        boundary = B{k};
                        % Check if it already exists
                        goodTemp=goodB;
                        goodTemp{end+1}=boundary;
                        M = arrayfun(@(l) sum(arrayfun(@(m) isequal(goodTemp{l}, goodTemp{m}), 1:numel(goodTemp))), 1:numel(goodTemp));
                        if size(boundary,1)>cutBound & max(M)==1
                            plot(boundary(:,2),...
                                boundary(:,1),'r','LineWidth',2);
                            goodB{end+1}=boundary;
                        end
                    end
                end
            else
                boundary = B{k};
                if size(boundary,1)>cutBound
                    plot(boundary(:,2),...
                        boundary(:,1),'r','LineWidth',2);
                    goodB{end+1}=boundary;
                end
            end
        end
        
       

%        I = mat2gray(reflMap);
        
        %         % Adjust data to span data range.
%         Iadj = imadjust(I);
%         
%         % Threshold image - adaptive threshold
%         BW = imbinarize(Iadj, 'adaptive', 'Sensitivity', 0.500000, 'ForegroundPolarity', 'dark');
%         BW(isnan(reflMap))=0;
%         
%         BW = imfill(BW,'holes');
%         BW = bwareaopen(BW,1000);
%         
%         threshMask=double(BW);
%         threshMask(threshMask==0)=nan;
%         
%         properties = regionprops(BW, {'Area', 'Eccentricity', 'EquivDiameter',...
%             'Extent', 'FilledArea', 'MajorAxisLength', 'MinorAxisLength',...
%             'Orientation', 'Perimeter', 'Solidity'});
%         
%         imageRegionAnalyzer(BW);

%         
        %imageMorphology(I)
        
        % Plot
%         timeMat=repmat(timeMap,size(reflMap,1),1);
%         
%         fig1=figure('DefaultAxesFontSize',11,'position',[100,100,1300,900]);
%         
%         s1=subplot(2,1,1);
%         hold on
%         sub1=surf(timeMap,aslMap./1000,reflMap,'edgecolor','none');
%         view(2);
%         sub1=colMapDBZ(sub1);
%         ylabel('Altitude (km)');
%         xlim([timeMap(1),timeMap(end)]);
%         title('Reflectivity')
%         grid on
%         pos1=s1.Position;
%         
%         s2=subplot(2,1,2);
%         hold on
%         sub2=surf(timeMap,aslMap./1000,threshMask,'edgecolor','none');
%         view(2);
%         ylabel('Altitude (km)');
%         xlim([timeMap(1),timeMap(end)]);
%         grid on
%         pos2=s2.Position;
%         s2.Position=[pos1(1),pos2(2),pos1(3),pos1(4)];
         show1=1;
    end
end
%% Plot

close all

if plotFields
    fig1=figure('DefaultAxesFontSize',11,'position',[100,100,1300,900]);
    
    %%%%%%%%%%%%%%%%%%%%%%%% DBZ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax1=subplot(2,1,1);
    hold on;
    
    sub1=surf(data.time,data.asl./1000,refl,'edgecolor','none');
    view(2);
    sub1=colMapDBZ(sub1);
    scatter(timeMat(zeroInds),data.asl(zeroInds)./1000,10,'k','filled');
    ax = gca;
    ax.SortMethod = 'childorder';
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('Reflectivity and ERA5 freezing level')
    grid on
    
    ax2=subplot(2,1,2);
    hold on;
    sub2=surf(data.time,data.asl./1000,cloudNum,'edgecolor','none');
    view(2);
    ax = gca;
    ax.SortMethod = 'childorder';
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('Reflectivity and melting layer')
    grid on
    
    formatOut = 'yyyymmdd_HHMM';
    set(gcf,'PaperPositionMode','auto')
    % print([figdir,'refl_',datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_zeroDegree'],'-dpng','-r0');
end

