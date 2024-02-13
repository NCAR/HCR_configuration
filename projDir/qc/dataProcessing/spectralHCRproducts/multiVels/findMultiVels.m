function [majorVels,majorPows,minorVels,minorPows,leftSvel,rightSvel,leftPvel,rightPvel]= ...
    findMultiVels(powerIn,specVelIn,powerRaw,majorVels,majorPows,minorVels,minorPows,leftSvel,rightSvel,leftPvel,rightPvel,nn)
% Find maxima and minima in spectra
velMajor=nan(size(powerIn,1),35);
powMajor=nan(size(powerIn,1),35);
velMinor=nan(size(powerIn,1),35);
powMinor=nan(size(powerIn,1),35);

% Remove lower part
maxIn=max(powerIn,[],2,'omitmissing');
minIn=min(powerIn,[],2,'omitmissing');

powSpread=(maxIn-minIn);
powerIn(powerIn<minIn+powSpread/3)=nan;

dataInds=find(any(~isnan(powerIn),2));

% Remove transmitter pulse
dataInds(dataInds<=12)=[];

[locsMaxAll,maxPromAll]=islocalmax(powerIn,2);

diff1Curve=cat(2,nan(size(powerIn,1),1),diff(powerIn,1,2));
diff1MaxAll=islocalmax(diff1Curve,2);
diff1MinAll=islocalmin(diff1Curve,2);

% For shoulder points
diff2Curve=cat(2,diff(diff1Curve,1,2),nan(size(powerIn,1),1));
maxIn2=max(powerIn,[],2,'omitmissing');
minIn2=min(powerIn,[],2,'omitmissing');

powSpread2=(maxIn2-minIn2);
diff2Curve(powSpread2<10,:)=nan;

diff2MinAll=islocalmin(diff2Curve,2);

for jj=1:length(dataInds)
    ii=dataInds(jj);

    powerOrig=powerIn(ii,:);

    % Find maxima % and minima
    majorMax=find(locsMaxAll(ii,:)==1);
    maxProm=maxPromAll(ii,majorMax);

    % Find change points
    diff1Max=find(diff1MaxAll(ii,:)==1);
    diff1Min=find(diff1MinAll(ii,:)==1);

    % Find steep points
    diff2Min=find(diff2MinAll(ii,:)==1);
    leftShoulder=[];
    rightShoulder=[];
    if ~isempty(diff2Min)
        leftShoulder=diff2Min(1);
        rightShoulder=diff2Min(end);
        if leftShoulder==rightShoulder
            leftShoulder=[];
            rightShoulder=[];
        end
    end
    if ~isempty(leftShoulder)
        maxProm(majorMax<=leftShoulder | majorMax>=rightShoulder)=[];
        majorMax(majorMax<=leftShoulder | majorMax>=rightShoulder)=[];
      
        % Output
        leftSvel(ii,nn)=specVelIn(ii,leftShoulder);
        rightSvel(ii,nn)=specVelIn(ii,rightShoulder);
        leftPvel(ii,nn)=powerIn(ii,leftShoulder);
        rightPvel(ii,nn)=powerIn(ii,rightShoulder);
    end

    % Find minor peaks
    minorMax=[];

    for kk=1:length(diff1Max)-1
        if powerOrig(diff1Max(kk+1))>powerOrig(diff1Max(kk))
            findMin=find(diff1Min>diff1Max(kk) & diff1Min<diff1Max(kk+1));
            if length(findMin)==1
                hasPeak=sum(majorMax>diff1Max(kk) & majorMax<diff1Max(kk+1));
                if hasPeak==0
                    minorMax=cat(1,minorMax,diff1Min(findMin));
                end
            end
        end
    end
    for kk=1:length(diff1Min)-1
        if powerOrig(diff1Min(kk))>powerOrig(diff1Min(kk+1))
            findMax=find(diff1Max>diff1Min(kk) & diff1Max<diff1Min(kk+1));
            if length(findMax)==1
                hasPeak=sum(majorMax>diff1Min(kk) & majorMax<diff1Min(kk+1));
                if hasPeak==0
                    minorMax=cat(1,minorMax,diff1Max(findMax));
                end
            end
        end
    end

    % Move peaks with low prominence to minor
    minorMax=cat(1,minorMax,(majorMax(maxProm<3))');
    majorMax(maxProm<3)=[];
    minorMax=sort(minorMax);

    %% Add to output

    velMajor(ii,1:length(majorMax))=specVelIn(ii,majorMax);
    powMajor(ii,1:length(majorMax))=powerIn(ii,majorMax);
    velMinor(ii,1:length(minorMax))=specVelIn(ii,minorMax);
    powMinor(ii,1:length(minorMax))=powerIn(ii,minorMax);

    if nn==1
        scatter(specVelIn(ii,:),powerRaw(ii,:),'MarkerEdgeColor','k');
        xlim([-8,8]);
        hold on
        plot(specVelIn(ii,:),powerOrig,'-b','linewidth',2);
        scatter(specVelIn(ii,majorMax(~isnan(majorMax))),powerOrig(majorMax(~isnan(majorMax))),80,'filled','MarkerFaceColor','r')
        scatter(specVelIn(ii,leftShoulder),powerOrig(leftShoulder),80,'filled','MarkerFaceColor','m')
        scatter(specVelIn(ii,rightShoulder),powerOrig(rightShoulder),80,'filled','MarkerFaceColor','m')
        scatter(specVelIn(ii,diff1Min),powerOrig(diff1Min),'filled','MarkerFaceColor','y')
        scatter(specVelIn(ii,diff1Max),powerOrig(diff1Max),'filled','MarkerFaceColor','k')
        if ~isempty(minorMax)
            scatter(specVelIn(ii,minorMax),powerOrig(minorMax),80,'filled','MarkerFaceColor','g')
        end

        % yyaxis right
        % plot(specVelIn(ii,:),diff2Curve(ii,:),'-c','linewidth',1.5);
        % hold on
        % scatter(diffMin,diffCurve(diffMin),'filled','MarkerFaceColor','y')
        % scatter(diffMax,diffCurve(diffMax),'filled','MarkerFaceColor','k')
        %ylim([-0.5,0.5]);
        grid on

        cla
        % yyaxis left
        % cla
    end
end

% Major
findEmpty=all(isnan(velMajor),1);
velMajor(:,findEmpty)=[];
powMajor(:,findEmpty)=[];

checkDims=size(majorVels,3)-size(velMajor,2);
if checkDims<0
    majorVels=padarray(majorVels,[0,0,abs(checkDims)],nan,'post');
    majorPows=padarray(majorPows,[0,0,abs(checkDims)],nan,'post');
end

dim3=size(velMajor,2);
majorVels(:,nn,1:dim3)=single(velMajor);
majorPows(:,nn,1:dim3)=single(powMajor);

% Minor
findEmpty=all(isnan(velMinor),1);
velMinor(:,findEmpty)=[];
powMinor(:,findEmpty)=[];

checkDims=size(minorVels,3)-size(velMinor,2);
if checkDims<0
    minorVels=padarray(minorVels,[0,0,abs(checkDims)],nan,'post');
    minorPows=padarray(minorPows,[0,0,abs(checkDims)],nan,'post');
end

dim3=size(velMinor,2);
minorVels(:,nn,1:dim3)=single(velMinor);
minorPows(:,nn,1:dim3)=single(powMinor);
end