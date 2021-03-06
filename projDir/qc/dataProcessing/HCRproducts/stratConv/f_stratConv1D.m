function stratConv1D=f_stratConv1D(stratConv,meltingLayer,icingL,asl,topo)
% Convert 2D stratiform convective field to 1D
stratConv1D=nan(1,size(stratConv,2));
% Initialize all as stratiform
stratConv1D(:)=2;

% If vertical extent of warm cloud is less than 50% and the lowest pixel
% is higher than half the distance between earth surface and melting layer
% make whole cloud cold for 1D classification purposes

% If whole atmosphere is cold, make low convective convective

convMask=zeros(size(stratConv));
convMask(stratConv==11 | stratConv==12)=1;

convAreas=bwconncomp(convMask);
stratConvTemp=zeros(size(stratConv));

for ii=1:convAreas.NumObjects
    stratConvCloud=nan(size(stratConv));
    stratConvCloud(convAreas.PixelIdxList{ii})=stratConv(convAreas.PixelIdxList{ii});
    
    % Calculate fraction
    aslA=asl(stratConvCloud==12);
    aslB=asl(stratConvCloud==11);
        
    extA=max(aslA)-min(aslA);
    extTot=max(aslA)-min(aslB);
    
    aboveFrac=extA/extTot;
    
    % Make sure lowest point is not intersecting plane altitude
    [rI cI]=ind2sub(size(stratConv),convAreas.PixelIdxList{ii});
    aslAll=asl(convAreas.PixelIdxList{ii});
    [minAsl idxMin]=min(aslAll);
    minRow=rI(idxMin);
    
    % Half point between cloud and surface
    distIceTopo=median(icingL(cI)-topo(cI))/2;
    halfPoint=median(topo(cI))+distIceTopo;
    
    if aboveFrac>0.5 & min(aslB)>halfPoint & minRow>19
        stratConv(convAreas.PixelIdxList{ii})=12;
    end
    
    if min(aslA)<median(topo(cI))+220
        stratConvTemp(convAreas.PixelIdxList{ii})=1;
    end
end

% Check that areas of isolated convection are not all cold
convMask=zeros(size(stratConv));
convMask(stratConv==10)=1;

convAreas=bwconncomp(convMask);

for ii=1:convAreas.NumObjects
    stratConvCloud=nan(size(stratConv));
    stratConvCloud(convAreas.PixelIdxList{ii})=stratConv(convAreas.PixelIdxList{ii});
    
    % Lowest Asl
    aslTemp=asl;
    aslTemp(isnan(stratConvCloud))=nan;
    aslTemp(meltingLayer~=20)=nan;
    aslA=asl(~isnan(aslTemp));
        
%     extA=max(aslA)-min(aslA);
%     extTot=max(aslA)-min(aslB);
%     
%     aboveFrac=extA/extTot;
    
    % Make sure lowest point is not intersecting plane altitude
    [rI cI]=ind2sub(size(stratConv),convAreas.PixelIdxList{ii});
%     aslAll=asl(convAreas.PixelIdxList{ii});
%     [minAsl idxMin]=min(aslAll);
%     minRow=rI(idxMin);
%     
%     % Half point between cloud and surface
%     distIceTopo=median(icingL(cI)-topo(cI))/2;
%     halfPoint=median(topo(cI))+distIceTopo;
%     
%     if aboveFrac>0.5 & min(aslB)>halfPoint & minRow>19
%         stratConv(convAreas.PixelIdxList{ii})=12;
%     end
    
    if min(aslA)<median(topo(cI))+220
        stratConvTemp(convAreas.PixelIdxList{ii})=1;
    end
end

stratConvBelow=stratConv;
stratConvBelow(meltingLayer~=10)=nan;
stratConvBelow(stratConvTemp==1)=10;

stratConv1D(any(stratConvBelow<12,1))=1;

%stratConv1D(any(stratConv==10,1))=1;

% If there is only convective echo, make it convective
stratCount=sum(stratConv>19,1);
stratConv1D(stratConv1D==2 & stratCount==0)=1;

% Remove rays with no data
stratConv1D(all(isnan(stratConv),1))=nan;
end