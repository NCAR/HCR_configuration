%% Find top and bottom altitudes
meltMaskTurned=meltMaskInit';
mFirst=size(meltMaskTurned,2);
[valFirst,locFirst]=max(meltMaskTurned,[],2);
kFirst=mFirst+1-locFirst;
kFirst(valFirst==0)=nan;
kFirst=abs(size(meltMaskInit,1)-kFirst);
linFirst=sub2ind(size(meltMaskInit),kFirst',1:size(meltMaskInit,2));
linFirst(isnan(linFirst))=[];
altFirstNoNan=dataShort.asl(linFirst);
altFirst=nan(size(dataShort.time));
altFirst(~isnan(kFirst))=altFirstNoNan;

mLast=size(meltMaskTurned,2);
[valLast,locLast]=max(fliplr(meltMaskTurned),[],2);
kLast=mLast+1-locLast;
kLast(valLast==0)=nan;
linLast=sub2ind(size(meltMaskInit),kLast',1:size(meltMaskInit,2));
linLast(isnan(linLast))=[];
altLastNoNan=dataShort.asl(linLast);
altLast=nan(size(dataShort.time));
altLast(~isnan(kLast))=altLastNoNan;

% Smooth slightly
smoothVal2=15;
altFirstS=movmedian(altFirst,smoothVal2,'omitnan');
altLastS=movmedian(altLast,smoothVal2,'omitnan');

% Remove edges
flMask=~isnan(altFirst);
flMask=imclose(flMask,strel('line',5,0));
altFirstS(flMask==0)=nan;
altLastS(flMask==0)=nan;

%% Final mask
maskInds=find(flMask==1);
meltMask=zeros(size(meltMaskInit));
for jj=1:length(maskInds)
    altCol=dataShort.asl(:,maskInds(jj));
    if dataShort.elevation(maskInds(jj))>=0
        rightInds=find(altCol<=altLastS(maskInds(jj)) & altCol>=altFirstS(maskInds(jj)));
    else
        rightInds=find(altCol>=altLastS(maskInds(jj)) & altCol<=altFirstS(maskInds(jj)));
    end
    meltMask(rightInds,maskInds(jj))=1;
end