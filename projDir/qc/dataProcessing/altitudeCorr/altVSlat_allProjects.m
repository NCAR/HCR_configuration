% Ocean scan calibration for HCR data

clear all;
close all;

addpath('/h/eol/romatsch/gitPriv/process_HCR/oceanScans/functions/');
addpath('/h/eol/romatsch/gitPriv/process_HCR/oceanScans/colormaps/');
addpath('/h/eol/romatsch/gitPriv/process_HCR/NSCAL/functions/');
addpath(genpath('/h/eol/romatsch/gitPriv/utils/'));

load('/h/eol/romatsch/hcrCalib/latVSalt/socrates_qc2/socrates_table.mat');
socrates=allData;
clear allData

load('/h/eol/romatsch/hcrCalib/latVSalt/cset_qc1/cset_table.mat');
cset=allData;
clear allData

load('/h/eol/romatsch/hcrCalib/latVSalt/otrec_raw/otrec_table.mat');
otrec=allData;
clear allData

allData=cat(1,cset,socrates,otrec);
%% Scatter plot
close all
% % Remove outliers
% TF=isoutlier(allData.altDiff);
% allData(TF,:)=[];
% 
% if strcmp(project,'otrec')
%     TF=isoutlier(allData.altDiff);
%     allData(TF,:)=[];
% end

f1 = figure('Position',[200 500 1200 800],'DefaultAxesFontSize',12);

hold on
% scatter(abs(cset.lat),cset.altDiff);
% scatter(abs(socrates.lat),socrates.altDiff);
% scatter(abs(otrec.lat),otrec.altDiff);

scatter(cset.lat,cset.altDiff);
scatter(socrates.lat,socrates.altDiff);
scatter(otrec.lat,otrec.altDiff);


% hold on
% xlimits=xlim;
% ylimits=ylim;
% 
% fitOrth=gmregress(allData.lat,allData.altDiff,1);
% fitAll=[fitOrth(2) fitOrth(1)];
% xFit = xlimits(1):0.1:xlimits(2);
% yFit = polyval(fitAll, xFit);
% 
% plot(xFit, yFit,'-k','linewidth',3);
% 
% xlim(xlimits);
% ylim(ylimits);

xlabel('abs(Latitude) [deg]');
ylabel('GGALT-HCRins [m]');
legend('cset','socrates','otrec');
title(['Latitude vs altitude difference']);

%textbp(['y = ',num2str(fitAll(1)),' x + ',num2str(fitAll(2))],'FontSize',14);
%ylim([0 55])

set(f1,'PaperPositionMode','auto')
print(['/h/eol/romatsch/hcrCalib/latVSalt/latVSaltdiff2'],'-dpng','-r0')
