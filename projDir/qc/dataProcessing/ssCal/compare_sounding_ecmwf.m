%Compare attenuation calculated from sounding data vs ecmwf fake sounding

clear all;
close all;

addpath('/h/eol/romatsch/codes/matlab/hcrCalib/functions/');

indir='/scr/rain1/rsfdata/projects/socrates/hcr/cfradial/moments/10hz/';
sondedir='/h/eol/romatsch/data/hcrCalib/soundings/socrates2018.preliminary_qc_drops/';
modeldir='/h/eol/romatsch/data/reanalysis/ecmwf/socrates/soundIn/';

filedir='/h/eol/romatsch/hcrCalib/biasInFiles/';
infile='sounding_model.txt';

frq = 9.4406e10;

% read file with cases
caseListIn = readtable([filedir,infile],'Delimiter', ' ');
caseList=table2array(caseListIn);

outMat=nan(size(caseList,1),8);

for ii=1:size(caseList,1)
    % sounding
    [outMat(ii,1),outMat(ii,3),outMat(ii,5),outMat(ii,7)]=f_atten_layers_sfcWnd([sondedir,caseList{ii,1}],frq/1e+9);
    
    %get lat lon alt
    lon=ncread([indir,caseList{ii,3}],'longitude');
    lat=ncread([indir,caseList{ii,3}],'latitude');
    alt=ncread([indir,caseList{ii,3}],'altitude');
    
    %model
    [outMat(ii,2),outMat(ii,4),outMat(ii,6),outMat(ii,8)]=f_atten_layers_sfcWnd_ecmwf([modeldir,caseList{ii,2}],...
        frq/1e+9,nanmean(lon),nanmean(lat),nanmean(alt));
end