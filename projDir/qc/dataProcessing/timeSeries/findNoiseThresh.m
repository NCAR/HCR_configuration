function [noiseThresh,meanNoiseOut,R2]=findNoiseThresh(powIn,avNum)
% Find noise threshold and mean noise following
% Hildebrand and Sekhon, 1974 https://doi.org/10.1175/1520-0450(1974)013%3C0808:ODOTNL%3E2.0.CO;2

% low10=prctile(powIn,10);
% noiseThresh=10.^((low10+20)./10);
% powLin=10.^(powIn./10);
% powLin(isnan(powLin))=[];
% 
% R2=0;
% R2prev=0;
% stepDown=0.00001;
% 
% while R2<1
% 
%     noiseThresh=noiseThresh-stepDown;
% 
%     powLin(powLin>noiseThresh)=[];
%     sampleNum=length(powLin);
%     if isempty(powLin)
%         noiseThresh=noiseThresh+stepDown;
%         if stepDown>0.0000001
%             stepDown=max([stepDown/10,0.0000001]);
%             powLin=10.^(powIn./10);
%             powLin(isnan(powLin))=[];
%         else
%             break
%         end
%     else
%         % R2
%         meanNoise=sum(powLin)/sampleNum;
%         Q=sum(powLin.^2/sampleNum)-meanNoise.^2;
%         R2=meanNoise.^2/(Q*avNum);
% 
%         if R2>=1 & stepDown>0.0000001
%             noiseThresh=noiseThresh+stepDown;
%             stepDown=max([stepDown/10,0.0000001]);
%             powLin=10.^(powIn./10);
%             powLin(isnan(powLin))=[];
%             R2=0;
%             R2prev=0;
%         end
% 
%         if R2>R2prev
%             noiseThreshMax=noiseThresh;
%             meanNoiseMax=meanNoise;
%         end
% 
%         R2prev=R2;
% 
%         % plot(powLin)
%         % hold on
%         % plot([1,sampleNum],[noiseThresh,noiseThresh],'-c','LineWidth',1.5);
%         % hold off
% 
%     end
% end
% 
% if stepDown>0.0000001
%     warning('Missed noise threshold')
% end
% noiseThresh=real(10*log10(noiseThreshMax));
% meanNoiseOut=10*log10(meanNoiseMax);

%% Third version
% powLin2=10.^(powIn./10);
% powLin2(isnan(powLin2))=[];
% 
% testThresh=sort(powLin2);
% testThresh=repmat(testThresh',1,size(testThresh,2));
% 
% testPow=testThresh';
% 
% testPow(testPow>testThresh)=nan;
% sampleNum2=sum(~isnan(testPow),2);
% meanNoiseMat=sum(testPow,2,'omitmissing')./sampleNum2;
% Qmat=sum(testPow.^2./sampleNum2,2,'omitmissing')-meanNoiseMat.^2;
% R2mat=meanNoiseMat.^2./(Qmat*avNum);
% 
% R2ind=find(R2mat>=1,1,'last');
% noiseThreshMax2=testThresh(R2ind,1);
% meanNoiseMax2=meanNoiseMat(R2ind);
% 
% noiseThresh3=real(10*log10(noiseThreshMax2));
% meanNoiseOut3=10*log10(meanNoiseMax2);
% R23=R2mat(R2ind);

%% Third version loop

powLin2=10.^(powIn./10);

testThresh=sort(powLin2,'descend');

testPow=testThresh;

ii=1;
R2=0;

while R2<1
testPow(testPow>testThresh(ii))=[];
sampleNum=length(testPow);
meanNoise=sum(testPow)./sampleNum;
Q=sum(testPow.^2./sampleNum)-meanNoise.^2;
R2=meanNoise.^2./(Q*avNum);
ii=ii+1;
end

noisePow=powLin2;
noisePow(noisePow>testThresh(ii-1))=[];
meanNoisePow=sum(noisePow)./sampleNum;
noiseThresh=real(10*log10(testThresh(ii-1)));
meanNoiseOut=10*log10(meanNoisePow);
end