function dualParticles=findDualParticles_test(powerIn,specVelIn)
% Find maxima and minima in spectra
close all
dualParticles=[];
sampleNum=size(powerIn,2);

dataInds=find(any(~isnan(powerIn),2));

% Remove transmitter pulse
dataInds(dataInds<=12)=[];

for jj=1:length(dataInds)
    ii=dataInds(jj);

    powerOrig=powerIn(ii,:);

    % Find start and end of line segments
    lineMask=~isnan(powerOrig);
    diffLine=diff(lineMask);

    startLine=find(diffLine==1)+1;
    endLine=find(diffLine==-1);

    if isempty(startLine) & ~isempty(endLine)
        startLine=1;
    end
    if ~isempty(startLine) & isempty(endLine)
        endLine=length(powerOrig);
    end
    if startLine(1)>endLine(1)
        startLine=[1,startLine];
    end
    if startLine(end)>endLine(end)
        endLine=[endLine,length(powerOrig)];
    end

    % Find maxima and minima
    [locsMax,prom]=islocalmax(powerOrig);
    locsMax=find(locsMax==1);
    prom=prom(locsMax);
    locsMin=islocalmin(powerOrig);
    locsMin=find(locsMin==1);

    locsMin=cat(2,startLine,endLine,locsMin);
    locsMin=unique(locsMin);

    maxEnd=ismember(locsMax,[startLine,endLine]);
    locsMax(maxEnd==1)=[];

    % Find change points
    diffCurve=cat(2,nan,diff(powerOrig));
    [diffMax,promD]=islocalmax(diffCurve);
    diffMax=find(diffMax==1);
    promD=promD(diffMax);
    diffMin=islocalmin(diffCurve);
    diffMin=find(diffMin==1);

    plot(powerOrig,'-b','linewidth',2);
    hold on
    scatter(locsMin,powerOrig(locsMin),'filled','MarkerFaceColor','r')
    scatter(locsMax,powerOrig(locsMax),'filled','MarkerFaceColor','g')
    scatter(diffMin,powerOrig(diffMin),'filled','MarkerFaceColor','y')
    scatter(diffMax,powerOrig(diffMax),'filled','MarkerFaceColor','k')
    
    yyaxis right
    plot(diffCurve,'-c','linewidth',1.5);
    hold on
    scatter(diffMin,diffCurve(diffMin),'filled','MarkerFaceColor','y')
    scatter(diffMax,diffCurve(diffMax),'filled','MarkerFaceColor','k')
    ylim([-0.5,0.5]);
    grid on

    cla
    yyaxis left
    cla
end
end