function [dBmean, dBstdUp, dBstdDown] = dB_meanStd(dBin)
%calculate mean and standard deviation in dB space
inNonDB = 10.^(dBin/10.);
meanNonDB=nanmean(inNonDB);
stdNonDB=nanstd(inNonDB);
dBmean=10*log10(meanNonDB);
dBstdUpTot=10*log10(meanNonDB+stdNonDB);
dBstdDownTot=10*log10(meanNonDB-stdNonDB);
dBstdUp=dBstdUpTot-dBmean;
dBstdDown=dBmean-dBstdDownTot;
end

