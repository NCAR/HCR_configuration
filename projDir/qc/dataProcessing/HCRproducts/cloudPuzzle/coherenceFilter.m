function pidFilt=coherenceFilter(pidIn,winSize,pixPerc);
% (a) If more than 35 of 49 pixels are classified as clear, then the central pixel is set to clear. 
 %(b) If the central pixel is not set to clear and there are more than 7 of 49 pixels 
 % with the same type as the central pixel, it is left unchanged; 
 % otherwise, the central pixel is set to the classification type 
 % that is most plentiful in the 7 X 7 box.
 
 pidIn=single(pidIn);
 pidInVec=reshape(pidIn,[],1);
 
 %finalVec=nan(size(pidInVec));
 
 % Pad with nans
 pidPadded=padarray(pidIn,[floor(winSize/2) floor(winSize/2)],nan);
 % Each sliding window becomes a column
 pidCol=im2col(pidPadded,[winSize winSize],'sliding');
 
 % Count non nans
 sumNonNans=sum(~isnan(pidCol),1);
 
 % Set columns that have less than pixPer valid pixels to nan
 nonNanPerc=sumNonNans./size(pidCol,1); 
 pidCol(:,(nonNanPerc<(1-pixPerc)))=nan;
 pidInVec(nonNanPerc<(1-pixPerc))=nan;
 
 % Shrink data to only non-nan values
 goodInds=find(~isnan(pidInVec)); 
 pidVecShrink=pidInVec(goodInds);
 pidCol=pidCol(:,goodInds);
 
 % Find number of pixels with same value as center pixel
 for ii=1:length(pidVecShrink)
     numCenterPix=sum(pidCol(:,ii)==pidVecShrink(ii));
     if numCenterPix<winSize
         pidVecShrink(ii)=mode(pidCol(:,ii));
     end
 end
 
 % Put it back in big vector
 pidInVec(goodInds)=pidVecShrink;
 
 % Reshape
 pidFilt=reshape(pidInVec,size(pidIn));
end