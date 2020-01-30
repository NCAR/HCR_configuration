% Compare data from before and after velocity correction

clear all
close all

savefig=0;

startTime=datetime(2018,1,15,22,0,0);
endTime=datetime(2018,1,16,6,0,0);

addpath('/h/eol/romatsch/git/private/utils/');

%indir='/scr/rain1/rsfdata/projects/socrates/hcr/qc/cfradial/velcorr/10hz/';
indir='/scr/rain1/rsfdata/projects/socrates/hcr/qc/cfradial/moments/100hz/';

figdir='/h/eol/romatsch/hcrCalib/otherCalib/negativeSNR/';
formatOut = 'yyyymmdd_HHMM';

%% Load HCR data

edges=-25:0.5:0;

fileListHCR=makeFileList(indir,startTime,endTime,1);

if ~isempty(fileListHCR)
    
    histSnrNeg=zeros(1,length(edges)-1);
    
    % Get uncorrected data
    for ii=1:size(fileListHCR,2)
        if mod(ii,10)==0
            disp([num2str(ii),' of ',num2str(size(fileListHCR,2))]);
        end
        infile=fileListHCR{ii};
        
        snr=ncread(infile,'SNR');
        snrNegInd=find(snr<0);
        histIn=histcounts(snr(snrNegInd),edges);
        histSnrNeg=histSnrNeg+histIn;
    end
end

%% Plot histogram
close all

f1=figure('DefaultAxesFontSize',14);
set(f1,'Position',[200 500 1000 1000]);

bar(histSnrNeg,1);
xticks(0.5:10:50.5);
xticklabels(-25:5:0)
xlim([0.5 50.5]);
xlabel('SNR (dB)');
ylabel('Number of samples');
title(['Negative SNR samples ',datestr(startTime,formatOut),' to ',datestr(endTime,formatOut)],'interpreter','none');

set(gcf,'PaperPositionMode','auto')
print(f1, [figdir,datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_negSNR_100Hz'],'-dpng','-r0');
