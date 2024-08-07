function dbzText=f_reflTexture(DBZ,pixRad,dbzBase)
% Calculate reflectivity texture
dbzText=nan(size(DBZ));

% Pad data at start and end
dbzPadded=padarray(DBZ,[0 pixRad],nan);

% Fill in areas with no data
dbzPadded=fillmissing(dbzPadded,'linear',2,'EndValues','nearest');

% Adjust reflectivity with base value
dbzPadded=dbzPadded-dbzBase;

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
    dbzCorr=dbzBlock-newY+mean(dbzBlock,2,'omitnan');
    dbzCorr(dbzCorr<1)=1;

    % Calculate texture
    tdbz=sqrt(std(dbzCorr.^2,[],2,'omitnan'));

    dbzText(:,ii)=tdbz;
end
dbzText(isnan(DBZ))=nan;
end