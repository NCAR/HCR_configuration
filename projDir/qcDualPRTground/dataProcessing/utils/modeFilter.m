function fieldFilt=modeFilter(fieldIn,winSize,pixPerc)
% 2D median filter
% winSize: size of 7d box
% pixPerc: percentage of required non-nan pixels

fieldIn=single(fieldIn);

% Pad with nans
pidPadded=padarray(fieldIn,[floor(winSize/2) floor(winSize/2)],nan);

% Only short time periods can be processed in one junk
junkSize=2000;
if size(pidPadded,2)<junkSize
    % Each sliding window becomes a column
    pidCol=im2col(pidPadded,[winSize winSize],'sliding');
    
    % Count non nans
    sumNonNans=sum(~isnan(pidCol),1);
    
    % Set columns that have less than pixPer valid pixels to nan
    nonNanPerc=sumNonNans./size(pidCol,1);
    pidCol(:,(nonNanPerc<(1-pixPerc)))=nan;
    
    % Mode
    modeVec=mode(pidCol,1);
    
    % Reshape
    fieldFilt=reshape(modeVec,size(fieldIn));
else
    % Process data in junks
    junks=1+floor(winSize/2):junkSize:size(pidPadded,2);
    startJunks=junks-floor(winSize/2);
    endJunks=junks(2:end)+floor(winSize/2)-1;
    if endJunks(end)~=size(pidPadded,2)
        endJunks=[endJunks,size(pidPadded,2)];
    else
        startJunks(end)=[];
    end
    
    fieldFilt=nan(size(fieldIn,1),size(pidPadded,2));
    
    for ii=1:length(startJunks)
        % Each sliding window becomes a column
        pidJunk=pidPadded(:,startJunks(ii):endJunks(ii));
        
        pidCol=im2col(pidJunk,[winSize winSize],'sliding');
        
        % Count non nans
        sumNonNans=sum(~isnan(pidCol),1);
        
        % Set columns that have less than pixPer valid pixels to nan
        nonNanPerc=sumNonNans./size(pidCol,1);
        pidCol(:,(nonNanPerc<(1-pixPerc)))=nan;
        
        % Mode
        modeJunk=mode(pidCol,1);
        
        % Reshape
        fieldJunk=reshape(modeJunk,size(pidJunk,1)-2*floor(winSize/2),size(pidJunk,2)-2*floor(winSize/2));
        fieldFilt(:,startJunks(ii)+floor(winSize/2):endJunks(ii)-floor(winSize/2))=fieldJunk;            
    end
    fieldFilt=fieldFilt(:,startJunks(1)+floor(winSize/2):endJunks(end)-floor(winSize/2));
end

% Set all pixels that were nan in input to nan
fieldFilt(isnan(fieldIn))=nan;
end