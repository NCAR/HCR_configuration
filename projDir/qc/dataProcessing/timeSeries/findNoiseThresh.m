function [noiseThresh,meanNoise,R2]=findNoiseThresh(powIn,avNum)
% Find noise threshold and mean noise following
% Hildebrand and Sekhon, 1974 https://doi.org/10.1175/1520-0450(1974)013%3C0808:ODOTNL%3E2.0.CO;2

low10=prctile(powIn,10);
noiseThresh=10.^((low10+20)./10);
powLin=10.^(powIn./10);
powLin(isnan(powLin))=[];

R2=0;
R2prev=0;
stepDown=0.00001;

while R2<1
    
    noiseThresh=noiseThresh-stepDown;

    powLin(powLin>noiseThresh)=[];
    sampleNum=length(powLin);
    if isempty(powLin)
        noiseThresh=noiseThresh+stepDown;
        if stepDown>0.0000001
            stepDown=max([stepDown/10,0.0000001]);
            powLin=10.^(powIn./10);
            powLin(isnan(powLin))=[];
        else
            break
        end
    else
        % R2
        meanNoise=sum(powLin)/sampleNum;
        Q=sum(powLin.^2/sampleNum)-meanNoise.^2;
        R2=meanNoise.^2/(Q*avNum);

        if R2>=1 & stepDown>0.0000001
            noiseThresh=noiseThresh+stepDown;
            stepDown=max([stepDown/10,0.0000001]);
            powLin=10.^(powIn./10);
            powLin(isnan(powLin))=[];
            R2=0;
            R2prev=0;
        end

        if R2>R2prev
            noiseThreshMax=noiseThresh;
            meanNoiseMax=meanNoise;
        end

        R2prev=R2;

        % plot(powLin)
        % hold on
        % plot([1,sampleNum],[noiseThresh,noiseThresh],'-c','LineWidth',1.5);
        % hold off

    end
end

if stepDown>0.0000001
    warning('Missed noise threshold')
end
noiseThresh=real(10*log10(noiseThreshMax));
meanNoise=10*log10(meanNoiseMax);
end