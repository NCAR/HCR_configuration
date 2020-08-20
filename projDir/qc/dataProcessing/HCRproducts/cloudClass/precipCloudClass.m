function cloudFlag=precipCloudClass(cloudParams)
% Classify precipitating clouds into Deep (1), Ns (2), Cu (3), Sc (4), St (5), Ac (7)
cloudFlag=nan;

%% Deep flag
deepFlag=0;

if abs(cloudParams.meanLat)<23.5
    if (cloudParams.maxAgl>12 & cloudParams.max10dbzAgl>8.2) | ...
            (cloudParams.maxAgl>14 & cloudParams.meanThickness>12) | ...
            (cloudParams.meanMaxAgl>8.5 & cloudParams.max10dbzAgl>8.4)
        deepFlag=1;
    end
elseif (cloudParams.maxAgl>10 & cloudParams.max10dbzAgl>7.2) | ...
        cloudParams.max10dbzAgl>7.8 | ...
        ((cloudParams.maxAgl-cloudParams.max10dbzAgl)<1 & cloudParams.max10dbzAgl>7.5)
    deepFlag=1;
end

%% Conv flag

convFlag=0;

if deepFlag==1
    convFlag=1;
elseif (cloudParams.maxMaxRefl>14 | cloudParams.meanMaxRefl>4) & ... % CloudSat also checks for intense_prep_flag
        (cloudParams.lengthKM<80 | cloudParams.stdMaxAgl>0.5 | cloudParams.meanMaxRefl>4) & ...
        (cloudParams.meanMaxAgl-cloudParams.max10dbzAgl)<0.34 & ...
        cloudParams.max10dbzAgl>3 & cloudParams.meanThickness<5
    convFlag=1;
end
end

