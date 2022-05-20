function [classOut] = postProcessPID(PID,data)

classOut=PID;

% No frozen in warm region
replaceMat=zeros(size(classOut));
replaceMat(data.MELTING_LAYER==10 & (classOut==8 | classOut==9))=1;

% No mixed in warm region where no LDR
replaceMat(data.MELTING_LAYER==10 & classOut==7 & isnan(data.LDR))=1;

% Replace with closest warm pixel
[oldR oldC]=find(~isnan(classOut) & replaceMat==0 & data.MELTING_LAYER==10);
[addR addC]=find(replaceMat==1);
idx = knnsearch([oldR oldC], [addR addC]);
nearest_OldValue = classOut(sub2ind(size(classOut), oldR(idx), oldC(idx)));
classOut(sub2ind(size(classOut), addR, addC))=nearest_OldValue;

% Below -40 C temps
% SC drizzle and and SC cloud liquid to small ice
lowLiq=find(data.TEMP<-40 & (classOut==)

end