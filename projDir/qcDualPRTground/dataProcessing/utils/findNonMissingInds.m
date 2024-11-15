function nonMissingInds = findNonMissingInds(data,gapSecs)
% Fill in missing and NScal
nonMissingInds=ones(size(data.time));

flagTemp=data.FLAG_short;

surfMask=zeros(1,size(flagTemp,2));
surfMask(find(any(data.FLAG_short==10 | data.FLAG_short==11,1)))=1;
surfDiff=diff(surfMask);

% Add surface echo back in
holeStart=find(surfDiff==1);
holeStart=holeStart+1;
holeEnd=find(surfDiff==-1);

if ~isempty(holeStart) & ~isempty(holeEnd)
    if holeStart(1)>holeEnd(1)
        holeStart=[1 holeStart];
    end
    if length(holeStart)~=length(holeEnd)
        holeEnd=[holeEnd,size(flagTemp,2)-1];
    end
    
    for ii=1:length(holeStart)
        if etime(datevec(data.time(holeEnd(ii))),datevec(data.time(holeStart(ii))))<gapSecs &...
                sign(data.elevation(holeStart(ii)))==sign(data.elevation(holeEnd(ii)))
            nonMissingInds(holeStart(ii):holeEnd(ii))=0;
        end
    end
end
end

