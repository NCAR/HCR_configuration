function pidFilt=coherenceFilter(pidIn,winSize,pixPerc);
% (a) If more than 35 of 49 pixels are classified as clear, then the central pixel is set to clear. 
 %(b) If the central pixel is not set to clear and there are more than 7 of 49 pixels 
 % with the same type as the central pixel, it is left unchanged; 
 % otherwise, the central pixel is set to the classification type 
 % that is most plentiful in the 7 X 7 box.
end

