function [outtable]=txtTable2matTable(filename, del1)
% Reads Mikes txt tables and convert to matlab table

fid=fopen(filename,'r');
slurp=fscanf(fid,'%c');
fclose(fid);

M=strread(slurp,'%s','delimiter','\n');

commentInds = find(contains(M,'#'));

header=strread(M{1},'%s','delimiter',del1)';

goodInds=1:length(M);
goodInds(commentInds)=[];

for i=1:length(goodInds)
    try
        temp=strread(M{goodInds(i)},'%f','delimiter',del1);
    catch
        temp=strread(M{goodInds(i)},'%s','delimiter',del1);
    end
    for j=1:length(temp)
        MM(i,j)=temp(j);
    end;
end;

try
    outtable=array2table(MM,'VariableNames',header(2:end));
catch
    outtable=array2table(MM);
end