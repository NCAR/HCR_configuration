function [highLayerOutVel,lowLayerOutVel,highLayerOutPow,lowLayerOutPow]=processMultiRegs(majorVel,majorPow,minorVel,minorPow,multLayers,inMask)
majorVel(repmat(inMask,1,1,size(majorVel,3))==0)=nan;
majorPow(repmat(inMask,1,1,size(majorPow,3))==0)=nan;
minorVel(repmat(inMask,1,1,size(minorVel,3))==0)=nan;
minorPow(repmat(inMask,1,1,size(minorPow,3))==0)=nan;

highLayerVel=nan(size(majorVel,1),size(majorVel,2));
lowLayerVel=highLayerVel;
highLayerPow=highLayerVel;
lowLayerPow=highLayerVel;

multL3D=repmat(multLayers,1,1,size(majorVel,3));
majorVel(isnan(multL3D))=nan;
majorPow(isnan(multL3D))=nan;

% Figue out if dual distribution
divider=[];

velDistrib=majorVel(~isnan(majorVel));
maxPerc=prctile(velDistrib,95);
minPerc=prctile(velDistrib,5);

[N1,e1]=histcounts(velDistrib,floor(minPerc):0.5:ceil(maxPerc));
[N2,e2]=histcounts(velDistrib,floor(minPerc)-0.25:0.5:ceil(maxPerc)+0.25);

plotY=0;
if plotY
    close all
    figure
    s1=subplot(4,1,1);
    bar(e1(1:end-1)+0.25,N1,1)
    xlim([e2(1),e2(end)])
    s2=subplot(4,1,2);

    bar(e2(1:end-1)+0.25,N2,1)
    xlim([e2(1),e2(end)])
    ylims=s2.YLim;
end

[min1,p1]=islocalmin(N1);
[min2,p2]=islocalmin(N2);

if sum(min1)~=0 & sum(min2)~=0
    min1=find(p1==max(p1));
    min1=min1(1);
    min2=find(p2==max(p2));
    min2=min2(1);
    if abs(e1(min1)-e2(min2))<=1
        if plotY
            hold on
            plot([e1(min1)+0.25,e2(min2)+0.25],ylims,'-r','LineWidth',2);
            hold off
        end
        divider=mean([e1(min1)+0.25,e2(min2)+0.25]);
        max1=islocalmax(N1);
        max2=islocalmax(N2);
        max1=e1(find(max1==1))+0.25;
        max2=e2(find(max2==1))+0.25;
        maxes=[max1,max2];
        maxLow=floor(mean(maxes(find(maxes<divider))));
        maxHigh=ceil(mean(maxes(find(maxes>divider))));
    end

    if ~isnan(divider)
        % Find high layer
        majorHigh=majorVel;
        majorHigh(majorVel<divider)=nan;
        majorHighP=majorPow;
        majorHighP(majorVel<divider)=nan;

        [highLayerVel,highLayerPow,highRM,highRMP]=buildLayer(majorHigh,highLayerVel,majorHighP,highLayerPow,maxHigh,ceil(maxPerc),divider);

        highMinor=minorVel;
        highMinor(minorVel<divider)=nan;

        highMinor=cat(3,highRM,highMinor);

        highMinorP=minorPow;
        highMinorP(minorVel<divider)=nan;

        highMinorP=cat(3,highRMP,highMinorP);

        [highLayerVel,highLayerPow,~]=buildLayer(highMinor,highLayerVel,highMinorP,highLayerPow,maxHigh,ceil(max(highMinor(:),[],'omitmissing')),divider);
       
        % Find low layer
        majorLow=majorVel;
        majorLow(majorVel>=divider)=nan;
        majorLowP=majorPow;
        majorLowP(majorVel>=divider)=nan;

        [lowLayerVel,lowLayerPow,lowRM,lowRMP]=buildLayer(majorLow,lowLayerVel,majorLowP,lowLayerPow,maxLow,divider,floor(minPerc));

        lowMinor=minorVel;
        lowMinor(minorVel>=divider)=nan;

        lowMinor=cat(3,lowRM,lowMinor);

        lowMinorP=minorPow;
        lowMinorP(minorVel>=divider)=nan;

        lowMinorP=cat(3,lowRMP,lowMinorP);

        [lowLayerVel,lowLayerPow,~]=buildLayer(lowMinor,lowLayerVel,lowMinorP,lowLayerPow,maxLow,divider,floor(min(lowMinor(:),[],'omitmissing')));

        if plotY
            figure
            subplot(2,1,1)
            surf(highLayerVel,'edgecolor','none')
            view(2)
            caxis([-5,15])
            colormap(jet(20))
            subplot(2,1,2)
            surf(lowLayerVel,'edgecolor','none')
            view(2)
            caxis([-5,15])
            colormap(jet(20))
        end
    end
end
highLayerOutVel=highLayerVel(inMask==1);
lowLayerOutVel=lowLayerVel(inMask==1);
highLayerOutPow=highLayerPow(inMask==1);
lowLayerOutPow=lowLayerPow(inMask==1);
end