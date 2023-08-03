function velText=f_velStdTexture(VELorig,elev,pixRad,velBase)
% Adjust for up pointing
VEL=VELorig;

% Shrink velocity areas to remove outliers at the edges
VELmask=zeros(size(VEL));
VELmask(~isnan(VELorig))=1;

VELsmall=imerode(VELmask,strel('disk',5));

VEL(VELmask==1 & VELsmall==0)=nan;

% Calculate velocity texture
velText=nan(size(VEL));

% Pad data at start and end
velPadded=padarray(VEL,[0 pixRad],nan);

% Fill in areas with no data
velPadded=fillmissing(velPadded,'linear',2,'EndValues','nearest');

% Adjust velocity with base value
velPadded=velPadded-velBase;

% Loop through data points in time direction and pull out right window
for ii=1:size(velPadded,2)-pixRad*2-1
    velBlock=velPadded(:,ii:ii+pixRad*2);
        
    % Calculate and remove slope of reflectivity
    % Calculate fit
    x1=1:size(velBlock,2);
    X=repmat(x1,size(velBlock,1),1);
    
    sumX=sum(X,2,'omitnan');
    sumY=sum(velBlock,2,'omitnan');
    sumXY=sum((velBlock.*X),2,'omitnan');
    sumX2=sum(X.^2,2,'omitnan');
    sumY2=sum(velBlock.^2,2,'omitnan');
    
    N=size(velBlock,2);
    
    a=(sumY.*sumX2-sumX.*sumXY)./(N.*sumX2-sumX.^2);
    b=(N.*sumXY-sumX.*sumY)./(N.*sumX2-sumX.^2);
    
    newY=a+b.*X;
    
    % Remove slope
    velCorr=velBlock-newY+mean(velBlock,2,'omitnan');
    velCorr(velCorr<1)=1;
    
    % Calculate texture
    %tvel=sqrt(std(velCorr,[],2,'omitnan'));
    tvel=std(velCorr,[],2,'omitnan');
    
%     % Remove data points with not enough data
%     countNonNan=sum(~isnan(velBlock),2);
%     dataFrac=countNonNan/size(velBlock,2);
%     tvel(dataFrac<goodDataFrac)=nan;
    
    velText(:,ii)=tvel;
end
velText(isnan(VEL))=nan;


% f1 = figure('Position',[200 500 1500 900],'DefaultAxesFontSize',12);
% s1=subplot(3,1,1);
% surf(flipud(velText),'edgecolor','none');
% view(2);
% s1.Colormap=lines(12);
% caxis([0 12]);
% colorbar;
% 
% s2=subplot(3,1,2);
% surf(flipud(VEL),'edgecolor','none');
% view(2);
% s2.Colormap=jet;
% caxis([-20 20]);
% colorbar;

end

