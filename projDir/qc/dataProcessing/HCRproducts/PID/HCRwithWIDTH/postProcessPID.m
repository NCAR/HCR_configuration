function [classOut] = postProcessPID(PID,data)

classOut=PID;

% No small frozen in warm region
replaceMat=zeros(size(classOut));
replaceMat(data.MELTING_LAYER==10 & classOut==9)=1;

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
[oldR2 oldC2]=find(~isnan(classOut) & replaceMat==0 & data.MELTING_LAYER==10 & classOut~=7 & classOut~=8 & classOut~=9);
[addR2 addC2]=find(replaceMat2==1);
idx = knnsearch([oldR2 oldC2], [addR2 addC2]);
nearest_OldValue2 = classOut(sub2ind(size(classOut), oldR2(idx), oldC2(idx)));
classOut(sub2ind(size(classOut), addR2, addC2))=nearest_OldValue2;
end