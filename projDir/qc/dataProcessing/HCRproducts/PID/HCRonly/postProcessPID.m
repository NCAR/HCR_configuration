function [PIDpost] = postProcessPID(PID,data)

% No small frozen in warm region
replaceMat=zeros(size(classOut));
replaceMat(data.MELTING_LAYER==10 & classOut==9)=1;

% Replace with closest warm pixel
[oldR oldC]=find(~isnan(classOut) & replaceMat==0 & data.MELTING_LAYER==10);
[addR addC]=find(replaceMat==1);
idx = knnsearch([oldR oldC], [addR addC]);
nearest_OldValue = classOut(sub2ind(size(classOut), oldR(idx), oldC(idx)));
classOut(sub2ind(size(classOut), addR, addC))=nearest_OldValue;

end