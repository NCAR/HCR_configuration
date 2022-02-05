function [powerSpecFilt] = filterPowerSpecPerc(powerSpecLarge,sampleNum)
% Remove lowest 10% of median

powerSpecMed=movmedian(powerSpecLarge,sampleNum/10,2);

powerSpecMin=min(powerSpecMed,[],2,'omitnan');

spread=max(powerSpecMed,[],2,'omitnan')-powerSpecMin;

powerSpecFilt=powerSpecMed;
powerSpecFilt(powerSpecMed<powerSpecMin+spread./2)=nan;
powerSpecFilt(spread<5,:)=nan;
% 
% figure
% 
% for ii=1:length(spread)
%     plot(powerSpecLarge(ii,:));
%     hold on
%     %plot([1,size(powerSpecFilt,2)],[]);
%     plot(powerSpecFilt(ii,:));
%     hold off
% end
end