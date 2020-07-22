%Read in sounding file and create layers fo compare to liebe code
%calculate attenuation for each layer and sum up

clear all;
close all;

%soundFile='/h/eol/romatsch/data/radar/soundings/D20150724_165608_P.QC.eol';
soundFile='/h/eol/romatsch/data/radar/soundings/D20150724_165608ed_P.QC.eol';

f = 94.4; %94.4; %Frequency in GHz

soundDataIn=importdata(soundFile,' ',13);
soundData=soundDataIn.data;

clear soundDataIn;

ptrhalt=cat(2,soundData(:,5),soundData(:,6),soundData(:,8),soundData(:,14));

%remove missing data
ptrhalt(any(ptrhalt==-999,2),:) = [];

%thin out layers
spacing=1;
oneSeq=(1:1:size(ptrhalt,1));
indVect=find(mod(oneSeq,spacing)==0);
if indVect(1)~=1
    indVect=indVect+1;
    indVect=cat(2,1,indVect);
end
flip_ptrhalt=flip(ptrhalt);
flip_ptrhalt=flip_ptrhalt(indVect,:);

%create output file which will be read in by the liebe code
liebeFile=strcat('/h/eol/romatsch/fortran/radar/atmosnTest',num2str(spacing),'.dat')
firstLine=cat(2,f,flip_ptrhalt(1,1)/10,size(flip_ptrhalt,1)-1);
outMat=cat(2,flip_ptrhalt(:,4)./1000,flip_ptrhalt(:,2)+273.5,flip_ptrhalt(:,3)./100,zeros(size(flip_ptrhalt,1),1));

dlmwrite(liebeFile,firstLine,'delimiter',' ');
dlmwrite(liebeFile,outMat,'-append','delimiter',' ');

ptrhalt=flip(flip_ptrhalt);

layer_p=mean(cat(2,ptrhalt(1:size(ptrhalt,1)-1,1),ptrhalt(2:size(ptrhalt,1),1)),2);
layer_t=mean(cat(2,ptrhalt(1:size(ptrhalt,1)-1,2),ptrhalt(2:size(ptrhalt,1),2)),2);
layer_rh=mean(cat(2,ptrhalt(1:size(ptrhalt,1)-1,3),ptrhalt(2:size(ptrhalt,1),3)),2);
layer_depth=abs(ptrhalt(1:size(ptrhalt,1)-1,4)-ptrhalt(2:size(ptrhalt,1),4));

[spGamma, spGamma0, spGammaW]=f_atten_ITUR(f,layer_p,layer_t,layer_rh);

layer_gamma=spGamma.*layer_depth/1000;

gammaTot=sum(layer_gamma)