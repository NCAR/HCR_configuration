function [classOut] = postProcessPID(PID,data)

classOut=PID;

% No frozen in warm region
replaceMat=zeros(size(classOut));
replaceMat(data.HCR_MELTING_LAYER==10 & (classOut==8 | classOut==9))=1;

% No mixed in warm region where no LDR
replaceMat(data.HCR_MELTING_LAYER==10 & classOut==7 & isnan(data.HCR_LDR))=1;

% Replace with closest warm pixel
[oldR oldC]=find(~isnan(classOut) & replaceMat==0 & data.HCR_MELTING_LAYER==10);
[addR addC]=find(replaceMat==1);
idx = knnsearch([oldR oldC], [addR addC]);
nearest_OldValue = classOut(sub2ind(size(classOut), oldR(idx), oldC(idx)));
classOut(sub2ind(size(classOut), addR, addC))=nearest_OldValue;

% Check if first HSRL gates agree with HCR gates. If not, change HCR to no
% phase
replaceVec=ones(1,size(PID,2));

hcrFirst=PID(18:20,:);
hsrlFirst=PID(21:23,:);

replaceVec(any(isnan(hsrlFirst),1))=0; % Do not replace if hsrl has nan

goodInds=find(replaceVec==1);
for ii=1:length(goodInds)
    uhsrl=unique(hsrlFirst(:,goodInds(ii)));
    if length(uhsrl)==1
        hcrPix=hcrFirst(:,goodInds(ii));
        % Precip
        precPix=find((hcrPix==1 | hcrPix==2 | hcrPix==4 | hcrPix==5) & hcrPix~=uhsrl);
        if ~isempty(precPix)
            classOut(17+precPix,goodInds(ii))=10;
        end
        % Cloud
        cloudPix=find((hcrPix==3 | hcrPix==6) & hcrPix~=uhsrl);
        if ~isempty(precPix)
            classOut(17+cloudPix,goodInds(ii))=11;
        end
    end
end
end