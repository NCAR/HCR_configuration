% Create composite plot for spicule clouds

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='spicule'; %socrates, aristo, cset
quality='qc0'; %field, qc1, or qc2
qcVersion='v0.1';
freqData='10hz'; % 10hz, 100hz, 2hz, or combined

ylimits=[0 10];

figdir=['/scr/sleet2/rsfdata/projects/spicule/hcr/qc0/cfradial/v0.1/cloudPlots/'];

indir=HCRdir(project,quality,qcVersion,freqData);

inHour=[18,18,18,18,18,19,19,19,19,19,19];
inMin=[33,40,45,52,58,03,08,15,21,27,33];

inYear=repmat(2021,1,length(inHour));
inMonth=repmat(5,1,length(inHour));
inDay=repmat(29,1,length(inHour));
inSec=repmat(0,1,length(inHour));

inTimes=datetime(inYear,inMonth,inDay,inHour,inMin,inSec);