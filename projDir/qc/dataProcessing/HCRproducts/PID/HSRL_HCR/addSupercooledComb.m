function classOut = addSupercooledComb(pid,data)
%% Add supercooled
classOut=pid;

classOut(classOut==6)=9;
classOut(classOut==5)=8;
classOut(classOut==4)=7;
classOut(classOut==3)=5;
classOut(classOut==2)=3;

% Supercooled rain
classOut(classOut==1 & data.TEMP<0)=2;

% Supercooled drizzle
classOut(classOut==3 & data.TEMP<0)=4;

% Supercooled cloud liquid
classOut(classOut==5 & data.TEMP<0)=6;
end