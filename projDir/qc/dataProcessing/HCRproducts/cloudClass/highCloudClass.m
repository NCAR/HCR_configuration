function cloudFlag=highCloudClass(cloudParams)
% Classify clouds into High clouds: Deep (1), Cu (3), As (6), High (8)
cloudFlag=nan;

if cloudParams.meanMaxRefl<0.05 & cloudParams.meanMaxReflTemp<-22.5 & ...
        cloudParams.minAgl>5 & cloudParams.meanThickness < 6.1 & ...
        cloudParams.meanMinAgl>5.5
    cloudFlag=8;
elseif cloudParams.meanMinAgl>2
    cloudFlag=6;
elseif cloudParams.meanThickness<6
    cloudFlag=3;
else
    cloudFlag=1;
end
end