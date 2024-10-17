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

%% Second version

low10=prctile(powIn,10);
noiseTop=10.^((low10+20)./10);
powLin2=10.^(powIn./10);
powLin2(isnan(powLin2))=[];

noiseBottom=min(powLin2);
spacing=0.00000001;

testThresh=noiseBottom:spacing:noiseTop;
testThresh=repmat(testThresh',1,size(powIn,2));
testPow=repmat(powLin2,size(testThresh,1),1);

testPow(testPow>testThresh)=nan;
sampleNum2=sum(~isnan(testPow),2);
%sampleNum=repmat(sampleNum,1,size(powIn,2));
meanNoiseMat=sum(testPow,2,'omitmissing')./sampleNum2;
Qmat=sum(testPow.^2./sampleNum2,2,'omitmissing')-meanNoiseMat.^2;
R2mat=meanNoiseMat.^2./(Qmat*avNum);

R2ind=find(R2mat>=1,1,'last');
noiseThreshMax2=testThresh(R2ind,1);
meanNoiseMax2=meanNoiseMat(R2ind);

noiseThresh=real(10*log10(noiseThreshMax2));
meanNoiseOut=10*log10(meanNoiseMax2);
R2=R2mat(R2ind);
end