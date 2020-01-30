function [PLT] = f_reflNSCAL(PLT,infile,highResTempDir,refPodTemp)

startTime=PLT.time(1);
endTime=PLT.time(end);

% Read table with results from noise source calibration
meanTable=readtable(infile);

phaseShift=round(meanTable.lnaTempLagSecs(1));

% make list with correct temperature file
tempFile=highResTempFiles_socrates(startTime,endTime,highResTempDir);

% read temperature data
indata=txtTable2matTable(tempFile,',');
EikTemp=indata.EikTemp;
PolSwitchTemp=indata.PolarizationSwitchTemp;
RfDetTemp=indata.RfDetectorTemp;
NoisSourceTemp=indata.NoiseSourceTemp;
VLnaTemp=indata.VLnaTemp;

allTemps=cat(2,EikTemp,PolSwitchTemp,RfDetTemp,NoisSourceTemp); % Calculate mean over all relevant temperatures
Temp=mean(allTemps,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correct for gain change due to VLNA temperature change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Get temperature times
timeTemp=datetime(indata.year,indata.month,indata.day,indata.hour,indata.min,indata.sec);
%adjust for phase shift
timeTempShift=timeTemp+seconds(phaseShift);
%remove data with wrong times
outOfTimeInds=find(timeTempShift<startTime | timeTempShift>(endTime));
Temp(outOfTimeInds)=[];
VLnaTemp(outOfTimeInds)=[];
timeTempShift(outOfTimeInds)=[];

% Smooth LNA and pod temperature data and synchronize with reflectivity
% data
meanVLNAt=movmean(VLnaTemp,20);
meanPodT=movmean(Temp,80);

reflTable=timetable(PLT.time,PLT.refl);
tempsTable=timetable(timeTempShift,meanVLNAt,meanPodT);

syncData=synchronize(reflTable,tempsTable,'first','linear');
syncData.Properties.VariableNames{'Var1'} = 'refl';

lnaDeltaGain=(syncData.meanVLNAt-meanTable.lnaTempReference(1))*meanTable.lnaGainChangePerC(1);
rxDeltaGain=(syncData.meanPodT-refPodTemp)*meanTable.rxGainChangePerC(1);
deltaGain=lnaDeltaGain+rxDeltaGain;

PLT.refl=PLT.refl-deltaGain;
end

