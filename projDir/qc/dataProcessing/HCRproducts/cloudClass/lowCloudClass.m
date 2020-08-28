function cloudFlag=lowCloudClass(cloudParams)
% Classify clouds into Low clouds: Deep (1), Ns (2), Cu (3), Sc (4), St (5), As (6), Ac (7)
cloudFlag=nan;

if cloudParams.meanLat<-65 & cloudParams.meanBaseTemp<-30 & cloudParams.meanTopTemp<-35
    cloudFlag=6;
    % There is a check for cloud fraction from MODIS in the CloudSat
    % algorithm
elseif cloudParams.inhomo<0.5 & cloudParams.maxAgl<3 & cloudParams.meanMinAgl<1.8 & cloudParams.intPrecip==0 % CloudSat checks also cloud fracion
    cloudFlag=5;
elseif cloudParams.inhomo>0.3 & cloudParams.maxAgl>3 & cloudParams.minTopTemp>-35 & ...
        cloudParams.meanMaxRefl<2 & cloudParams.meanThickness<8 & ...
        (cloudParams.meanMinAgl>1.8 | (cloudParams.meanMinAgl>1 & cloudParams.maxAgl>3.5) | ...
        (cloudParams.meanMaxRefl<-5 & cloudParams.maxAgl>3.5 & cloudParams.meanThickness>2))
    cloudFlag=7;
elseif cloudParams.meanThickness>8 & cloudParams.meanMaxRefl<0
    cloudFlag=6;
elseif (cloudParams.meanThickness>2 | cloudParams.intPrecip==0 | cloudParams.maxAgl>=3) & ...
        cloudParams.meanThickness<7 & cloudParams.meanMaxRefl>0 & cloudParams.lengthKM<100
    cloudFlag=3;
elseif (cloudParams.meanThickness>2 | cloudParams.maxAgl>=4) & ...
        cloudParams.meanThickness<7 & cloudParams.meanMaxRefl>-5 & cloudParams.lengthKM>56
    cloudFlag=2;
elseif cloudParams.meanThickness>8
    cloudFlag=1;
else
    cloudFlag=4;
end
end

