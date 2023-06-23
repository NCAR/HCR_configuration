function [classOut] = postProcessPID(PID,data)

classOut=PID;

% No frozen in warm region
replaceMat=zeros(size(classOut));
replaceMat(data.MELTING_LAYER==9 & (classOut==8 | classOut==9))=1;

% Replace with closest warm pixel
[oldR oldC]=find(~isnan(classOut) & replaceMat==0 & data.MELTING_LAYER==9);
[addR addC]=find(replaceMat==1);
idx = knnsearch([oldR oldC], [addR addC]);
nearest_OldValue = classOut(sub2ind(size(classOut), oldR(idx), oldC(idx)));
classOut(sub2ind(size(classOut), addR, addC))=nearest_OldValue;

% Below -40 C temps
% SC cloud liquid to cloud
classOut(data.TEMP<-40 & classOut==6)=11;
% SC rain and SC drizzle to precip
classOut(data.TEMP<-40 & (classOut==4 | classOut==2))=10;

end