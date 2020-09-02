function [cloud refl] = fillExtinct(data,cloud,reflIn)
% Fill in extinct echo
refl=reflIn;
%refl(data.FLAG>1)=nan;

flagTemp=data.FLAG;
flagTemp(data.FLAG==8)=7;
flagTemp(flagTemp==3)=99;
surfMask=ones(1,size(flagTemp,2));
surfMask(find(any(data.FLAG==7 | data.FLAG==8,1)))=0;

% Add in up pointing
surfMask(data.elevation>0)=0;
surfDiff=diff(surfMask);

% Add surface echo back in
holeStart=find(surfDiff==1);
holeStart=holeStart+1;
holeEnd=find(surfDiff==-1);

if ~isempty(holeStart) | ~isempty(holeEnd)
    if holeStart(1)>holeEnd(1)
        holeStart=[1 holeStart];
    end
    if length(holeStart)~=length(holeEnd)
        holeEnd=[holeEnd,size(flagTemp,2)-1];
    end
    
    for ii=1:length(holeStart)
        
        flagColEnd=flagTemp(:,holeEnd(ii)+1);
        minEnd=min(find(flagColEnd==7));
        
        if holeStart(ii)>1
            flagColStart=flagTemp(:,holeStart(ii)-1);
            minStart=min(find(flagColStart==7));
        else
            minStart=minEnd;
        end
        
        if isempty(minStart)
            minStart=minEnd;
        end
        if isempty(minEnd)
            minEnd=minStart;
        end
        
        newVec=holeStart(ii):holeEnd(ii);
        holeLength=holeEnd(ii)-holeStart(ii);
        
        if holeLength==0
            dataFill=minStart;
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
    extInd=find(any(flagTemp==99,1));
    
    for ii=1:length(extInd)
        cloudCol=cloud(:,extInd(ii));
        
        extCol=flagTemp(:,extInd(ii));
        extData=find(extCol==99);
        cloudNumOut=cloudCol(min(extData-1));
        cloud(extData,extInd(ii))=cloudNumOut;
        
        reflCol=refl(:,extInd(ii));
        reflCol(cloudCol~=cloudNumOut)=nan;
        reflMed=median(cloudCol,'omitnan');
        refl(extData,extInd(ii))=reflMed;
    end
end
end

