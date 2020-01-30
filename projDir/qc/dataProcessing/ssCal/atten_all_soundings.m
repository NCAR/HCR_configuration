%Read in list of sounding files and calculate zenith one way attenuation 
%Calculates attenuation according to the Liebe and the ITU-R method

clear all;
close all;

addpath('/h/eol/romatsch/codes/matlab/hcrCalib/functions/');

indir='/scr/sci/romatsch/data/dropSonde/otrec/';
outFile='/h/eol/romatsch/hcrCalib/notes/oneWayAtt_all_otrec.dat';

inlistIn=table2array(readtable(['/h/eol/romatsch/data/hcrCalib/soundings/list_otrec.dat'],...
    'Format','%s','delimiter','\n','ReadVariableNames',false));

%inlist=cell2mat(inlistIn);

f = 94.406172672; %Frequency in GHz
%cutoff=6000; %height up to which the attenuation is calculated (optional)

listLenght=size(inlistIn);

outMat=NaN(listLenght(1),4);

for i=1:listLenght(1);
    soundFile=inlistIn{i,1:end};    
    disp(soundFile);
    passFile=[indir soundFile];
    
    %[outMat(i,1),outMat(i,2),outMat(i,3),outMat(i,4)]=f_atten_layers(passFile,f,cutoff);
    [outMat(i,1),outMat(i,2),outMat(i,3),outMat(i,4)]=f_atten_layers(passFile,f);
end

inlistIn(any(isnan(outMat),2),:)=[];
outMat(any(isnan(outMat),2),:)=[];

outTable = array2table(round(outMat,2), 'RowNames', inlistIn);
outTable.Properties.VariableNames={'Att0Liebe' 'Att0ITU' 'dryAtt0ITU' 'moistAtt0ITU'};

writetable(outTable,outFile,'WriteRowNames',true, 'Delimiter','\t')

