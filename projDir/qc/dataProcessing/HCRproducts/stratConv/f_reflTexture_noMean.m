function dbzText=f_reflTexture_noMean(DBZ,pixRad,dbzThresh)
% Calculate reflectivity texture
dbzText=nan(size(DBZ));

DBZ(DBZ<dbzThresh)=nan;

% Pad data at start and end
dbzPadded=padarray(DBZ,[0 pixRad],nan);

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
dbzText(isnan(DBZ))=nan;
end

