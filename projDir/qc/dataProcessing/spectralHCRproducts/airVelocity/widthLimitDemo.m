% Test the limits of the standard aircraft broadening width correction

clear all
close all

delta=0.1:0.1:1;

widthRaw=0:0.005:5;

figdir='/scr/virga1/rsfdata/projects/spicule/hcr/qc1/cfradial/v1.2_full/airMotion/';

f1=figure('Position',[200 500 600 600],'DefaultAxesFontSize',12,'Renderer','painters');
hold on

leg={};

cols=jet(length(delta));

for ii=1:length(delta)
    squares=widthRaw.^2-delta(ii).^2;
    squares(squares<0.1)=0.01;

    widthCorr=sqrt(squares);

    plot(widthRaw,widthCorr,'-','LineWidth',1.5,'Color',cols(ii,:));

    leg{end+1}=['Delta ',num2str(delta(ii))];
end

xlabel('Width raw (m s^{-1})')
ylabel('Width corrected (m s^{-1})')
legend(leg,'Location','northwest')

title('Spectrum width aircraft motion correction')

xticks(0:5)
yticks(0:5)

grid on
box on

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,'widthCorrModel'],'-dpng','-r0');