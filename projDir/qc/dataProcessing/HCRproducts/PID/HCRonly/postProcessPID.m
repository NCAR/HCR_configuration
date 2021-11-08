function [classOut] = postProcessPID(PID,data)

classOut=PID;

% No small frozen in warm region
replaceMat=zeros(size(classOut));
replaceMat(data.MELTING_LAYER==10 & classOut==9)=1;

% No mixed in warm region where no LDR
replaceMat(data.MELTING_LAYER==10 & classOut==7 & isnan(data.LDR))=1;

% Replace with closest warm pixel
[oldR oldC]=find(~isnan(classOut) & replaceMat==0 & data.MELTING_LAYER==10);
[addR addC]=find(replaceMat==1);
idx = knnsearch([oldR oldC], [addR addC]);
nearest_OldValue = classOut(sub2ind(size(classOut), oldR(idx), oldC(idx)));
classOut(sub2ind(size(classOut), addR, addC))=nearest_OldValue;

% No large frozen below mixed below melting layer
belowMelt=nan(size(PID));
belowMelt(data.MELTING_LAYER==10)=classOut(data.MELTING_LAYER==10);

replaceMat2=zeros(size(classOut));

for ii=1:size(belowMelt,2)
    ray=belowMelt(:,ii);
    if sum(~isnan(ray))==0 | isempty(find(ray==7))
        continue
    end
    if data.elevation(ii)>0 % Up
        lastMelt=min(find(ray==7));
        lf=find(ray(1:lastMelt)==8);
        if ~isempty(lf)
            replaceMat2(lf,ii)=1;
        end
    else
        lastMelt=max(find(ray==7));
        ray(1:lastMelt)=nan;
        lf=find(ray==8);
        if ~isempty(lf)
            replaceMat2(lf,ii)=1;
        end
    end
end

% Replace with closest warm pixel
[oldR2 oldC2]=find(~isnan(classOut) & replaceMat2==0 & data.MELTING_LAYER==10 & classOut~=7 & classOut~=8 & classOut~=9);
[addR2 addC2]=find(replaceMat2==1);
idx = knnsearch([oldR2 oldC2], [addR2 addC2]);
nearest_OldValue2 = classOut(sub2ind(size(classOut), oldR2(idx), oldC2(idx)));
classOut(sub2ind(size(classOut), addR2, addC2))=nearest_OldValue2;

% Remove frozen from clouds that are completely below the melting layer
pidMask=~isnan(classOut);
maskAreas=bwconncomp(pidMask);

replaceMat3=zeros(size(classOut));

for ii=1:maskAreas.NumObjects
    pixArea=maskAreas.PixelIdxList{ii};
    belowFrac=length(find(data.MELTING_LAYER(pixArea)==10))./length(pixArea);
    if belowFrac>0.9
        pidArea=nan(size(classOut));
        pidArea(pixArea)=classOut(pixArea);
        replaceMat3(pidArea>=7 & pidArea<=9)=1;
    end
end

% Replace with closest warm pixel
[oldR3 oldC3]=find(~isnan(classOut) & replaceMat3==0 & data.MELTING_LAYER==10 & classOut~=7 & classOut~=8 & classOut~=9);
[addR3 addC3]=find(replaceMat3==1);
idx = knnsearch([oldR3 oldC3], [addR3 addC3]);
nearest_OldValue3 = classOut(sub2ind(size(classOut), oldR3(idx), oldC3(idx)));
classOut(sub2ind(size(classOut), addR3, addC3))=nearest_OldValue3;

end