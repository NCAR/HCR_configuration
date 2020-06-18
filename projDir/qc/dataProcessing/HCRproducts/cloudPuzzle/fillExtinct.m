function refl = fillExtinct(data)
% Fill in extinct echo

refl=data.DBZ;
refl(data.FLAG>1)=nan;

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
end

