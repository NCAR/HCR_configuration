function dbzText=f_reflTexture(data,pixRad,goodDataFrac)
% Calculate reflectivity texture
dbzText=nan(size(data.dbzMasked));

% Pad data at start and end
dbzPadded=padarray(data.dbzMasked,[0 pixRad],nan);

% figure
% surf(dbzPadded,'edgecolor','none');
% view(2)
% colormap(jet)
% caxis([-35 25])
% 
% % Analyze edges
% dbzBW=zeros(size(dbzPadded));
% dbzBW(~isnan(dbzPadded))=1;
% 
% % Fill small holes so they don't get bigger
% dbzBW=~bwareaopen(~dbzBW, 20);
% 
% % Shrink
% dbzErode=imerode(dbzBW, strel('disk', 50));
% 
% edgeArea=dbzBW+dbzErode;
% edgeArea(edgeArea>1)=0;
% 
% edges=dbzPadded;
% edges(edgeArea==0)=nan;
% 
% figure
% surf(edges,'edgecolor','none');
% view(2)
% colormap(jet)
% caxis([-35 25])

% Fill in areas with no data
dbzPadded=fillmissing(dbzPadded,'linear',2,'EndValues','nearest');

% Loop through data points in time direction and pull out right window
for ii=1:size(dbzPadded,2)-pixRad*2-1
    dbzBlock=dbzPadded(:,ii:ii+pixRad*2);
        
    % Calculate and remove slope of reflectivity
    % Calculate fit
    x1=1:size(dbzBlock,2);
    X=repmat(x1,size(dbzBlock,1),1);
    
    sumX=sum(X,2,'omitnan');
    sumY=sum(dbzBlock,2,'omitnan');
    sumXY=sum((dbzBlock.*X),2,'omitnan');
    sumX2=sum(X.^2,2,'omitnan');
    sumY2=sum(dbzBlock.^2,2,'omitnan');
    
    N=size(dbzBlock,2);
    
    a=(sumY.*sumX2-sumX.*sumXY)./(N.*sumX2-sumX.^2);
    b=(N.*sumXY-sumX.*sumY)./(N.*sumX2-sumX.^2);
    
    newY=a+b.*X;
    
    % Remove slope
    dbzCorr=dbzBlock-newY;
    
    % Calculate texture
    tdbz=sqrt(std(dbzCorr.^2,[],2,'omitnan'));
    
%     % Remove data points with not enough data
%     countNonNan=sum(~isnan(dbzBlock),2);
%     dataFrac=countNonNan/size(dbzBlock,2);
%     tdbz(dataFrac<goodDataFrac)=nan;
    
    dbzText(:,ii)=tdbz;
end
dbzText(isnan(data.dbzMasked))=nan;
end

