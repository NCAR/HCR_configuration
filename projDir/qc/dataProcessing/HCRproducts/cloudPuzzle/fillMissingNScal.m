function refl = fillMissingNScal(refl,data)
% Fill in missing and NScal
flagTemp=data.FLAG;

surfMask=zeros(1,size(flagTemp,2));
surfMask(find(any(data.FLAG==10 | data.FLAG==12,1)))=1;
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
        if etime(datevec(data.time(holeEnd(ii))),datevec(data.time(holeStart(ii))))<10
            reflStart=refl(:,holeStart(ii)-1);
            refl(:,holeStart(ii):holeEnd(ii))=repmat(reflStart,1,holeEnd(ii)-holeStart(ii)+1);
        end
    end
end
end

