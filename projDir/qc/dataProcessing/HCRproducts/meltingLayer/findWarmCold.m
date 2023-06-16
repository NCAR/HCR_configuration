function [warmCold,zeroAltsAdj,meltOnly]=findWarmCold(data)

%% Zero deg alts
[layerAltsOut,~,layerVelsOut,meltOnly,tempOut]=zeroDegIso(data);

%% Adjust freezing level alts
a=168;
b=40;

% No zero vel found
zeroVelNan=b/a;

zeroAltsAdj=nan(size(layerAltsOut));
for ii=1:size(layerAltsOut,1)
    zeroAlt=layerAltsOut(ii,:);
    zeroVel=layerVelsOut(ii,:);

    % Find offset
    zeroVelMask=~isnan(zeroVel);
    zeroVelMask=bwareaopen(zeroVelMask,15);
    zeroVel(zeroVelMask==0)=nan;
    zeroVel(isnan(zeroVel))=zeroVelNan;
    zeroVel(isnan(zeroAlt))=nan;

    % Offset
    offset=a.*zeroVel-b;
    %offset(offset<0)=0;

    offsetSmooth=movmedian(offset,25,'omitnan');
    offsetSmooth(isnan(offset))=nan;
    zeroAltsAdj(ii,:)=zeroAlt-offsetSmooth;
end

%% Mat with warm and cold regions
warmCold=zeros(size(data.TEMP));
warmCold(tempOut<=0)=0;
warmCold(tempOut>0)=2;
warmCold(isnan(data.TEMP))=nan;

%% Plot
plotYes=0;
if plotYes
    newInds=1:round(length(data.time)/2000):length(data.time);
    newASL=data.asl(:,newInds);
    newTime=data.time(newInds);
    newWarmCold=warmCold(:,newInds);

    close all
    fig1=figure('DefaultAxesFontSize',11,'position',[100,1300,1500,300]);
    s1=subplot(1,1,1);
    surf(newTime,newASL./1000,newWarmCold,'edgecolor','none');
    view(2);
    caxis([-0.5 2.5]);
    colormap=jet(3);
    ylim([0,6]);
    hold on
    for ii=1:size(layerAltsOut,1)
        plot(data.time,layerAltsOut./1000,'-r','LineWidth',2);
        plot(data.time,zeroAltsAdj./1000,'-g','LineWidth',2);
    end
    s1.SortMethod='childorder';
end

end