% Find ocean scan data
%!!!!!!!!!!!!!!!!! set break point at line 95 !!!!!!!!!!!!!!!!!!!!!!

clear all
close all

addpath('/h/eol/romatsch/gitPriv/process_HCR/oceanScans/functions/');
addpath('/h/eol/romatsch/gitPriv/process_HCR/oceanScans/colormaps/');
addpath('/h/eol/romatsch/gitPriv/process_HCR/NSCAL/functions/');
addpath(genpath('/h/eol/romatsch/gitPriv/utils/'));

% aristo
%direc = '/scr/snow2/rsfdata/projects/aristo-17/hcr/cfradial/moments/10hz/20170302/';

% cset
%direc = '/scr/snow2/rsfdata/projects/cset/hcr/cfradial/moments/10hz/20150812/';

%socrates
%'/scr/snow2/rsfdata/projects/socrates/hcr/cfradial/moments/10hz/20180116/';

%otrec
direc = '/scr/snow1/rsfdata/projects/otrec/hcr/qc1/cfradial/velcorr/10hz/20191002/';

%fileList=dir([direc, '/cfrad*SUR.nc']);
fileList=dir([direc, '/cfrad*.nc']);

%% Load data

timeData=[];
rotaData=[];

klim = size(fileList,1);
fprintf(' %d files have been selected.\n',klim);

for kk=1:klim;
    %disp(['File ',num2str(kk),' from ',num2str(klim)]);
    
    infile=[fileList(kk).folder,'/',fileList(kk).name];
    
    startTimeIn=ncread(infile,'time_coverage_start')';
    startTime=datetime(str2num(startTimeIn(1:4)),str2num(startTimeIn(6:7)),str2num(startTimeIn(9:10)),...
        str2num(startTimeIn(12:13)),str2num(startTimeIn(15:16)),str2num(startTimeIn(18:19)));
    timeRead=ncread(infile,'time')';
    
    timeIn=startTime+seconds(timeRead);
    
    try
        rotaIn=ncread(infile,'rotation');
        
        if length(timeIn)~=length(rotaIn)
            disp('Sizes do not match.');
        end
        
        % Get only data from downward pointing radar
        downInds=find(rotaIn<215 & rotaIn>145);
        
        timeData=cat(1,timeData,timeIn(downInds)');
        rotaData=cat(1,rotaData,rotaIn(downInds));
    end
end

%% Figure

close all

fig1=figure;
fig1.Position=[200 900 3000 400];
%plot(rotaData);
plot(timeData,rotaData)
xData=(1:length(rotaData));
    
endLoop=1;
while endLoop==1
   
    x = input('Enter start index:');
    y = input('Enter end index:');
   
    plot(xData(x:y),rotaData(x:y));
    
    endLoop=input('End loop (0) or continue (1)?');
end
%% 

rotaSmall=rotaData(x:y);
timeSmall=timeData(x:y);
%openvar('rotaSmall');
%refreshdata
fig2=figure;
fig2.Position=[200 400 1000 400];
plot(rotaSmall);
xData2=(1:length(rotaSmall));
    
endLoop=1;
while endLoop==1
   
    x2 = input('Enter start index:');
    y2 = input('Enter end index:');
   
    plot(xData2(x2:y2),rotaSmall(x2:y2));
    
    disp([timeSmall(x2),' ',timeSmall(y2)]);
    %endLoop=input('End loop (0) or continue (1)?');
end
