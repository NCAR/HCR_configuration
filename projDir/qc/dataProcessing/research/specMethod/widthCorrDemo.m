% Test the limits of the standard aircraft broadening width correction

clear all
close all

delta=0.1:0.1:1;

widthRaw=0:0.005:4;

figdir='/scr/virga1/rsfdata/projects/spicule/hcr/qc1/cfradial/v1.2_full/specParams/specPaperFigs/';

f1=figure('Position',[200 500 400 400],'DefaultAxesFontSize',12,'Renderer','painters');
t = tiledlayout(1,1,'TileSpacing','tight','Padding','tight');

s1=nexttile(1);

hold on

leg={};

cols=jet(length(delta));

for ii=1:length(delta)
    squares=widthRaw.^2-delta(ii).^2;
    squares(squares<0.01)=0.01;

    widthCorr=sqrt(squares);

    plot(widthRaw,widthCorr,'-','LineWidth',1.5,'Color',cols(ii,:));

    leg{end+1}=['BBS ',num2str(delta(ii))];
end

xlabel('WIDTH_{raw} (m s^{-1})')
ylabel('WIDTH (m s^{-1})')
legend(leg,'Location','northwest')

xticks(0:5)
yticks(0:5)

grid on
box on

set(gcf,'PaperPositionMode','auto')
%print(f1,[figdir,'widthCorrModel'],'-dpng','-r0');
exportgraphics(f1,[figdir,'widthCorrModel.png'],'Resolution','300');