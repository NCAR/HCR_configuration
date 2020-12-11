function fieldFilt=modeFilter(fieldIn,winSize,pixPerc)
% 2D median filter
% winSize: size of 7d box
% pixPerc: percentage of required non-nan pixels
 
 fieldIn=single(fieldIn);
   
 % Pad with nans
 pidPadded=padarray(fieldIn,[floor(winSize/2) floor(winSize/2)],nan);
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
 
 % Set all pixels that were nan in input to nan
 fieldFilt(isnan(fieldIn))=nan;
end