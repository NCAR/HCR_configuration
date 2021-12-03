% Plot hit miss table

clear all
close all

figdir=['/scr/snow2/rsfdata/projects/socrates/hcr/qc3/cfradial/hcr_hsrl_merge/v3.0/pidPlotsComb/paperFigs/'];

indir='/scr/snow2/rsfdata/projects/socrates/hcr/qc3/cfradial/hcr_hsrl_merge/v3.0/pidPlotsComb/comparePID_UW_largest_overlap/';

outTableIn=load([indir,'compTable.mat']);
outTableAll=outTableIn.outTableAll;

units_str_hcr={'Rain','SC Rain','Drizzle','SC Drizzle','Cloud Liquid','SC Cloud Liquid',...
    'Melting','Large Frozen','Small Frozen','Precip','Cloud'};

%centers={0.05:0.1:0.95 0.05:0.1:0.95};
centers={0.125:0.25:0.875 0.125:0.25:0.875};
cm=turbo(12);

wi=10;
hi=4.6;

fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[3,100,wi,hi]);
fig1.PaperPositionMode = 'manual';
fig1.PaperUnits = 'inches';
fig1.Units = 'inches';
fig1.PaperPosition = [0, 0, wi, hi];
fig1.PaperSize = [wi, hi];
fig1.Resize = 'off';
fig1.InvertHardcopy = 'off';

set(fig1,'color','w');

colormap(cm(2:11,:));
s={};

tickLoc=0.5:1:4.5;
tickLabel={'0','0.25','0.5','0.75','1'};

panel={'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'};
panelColor={'k','k','k','k','k','k','k','w','w'};

for ii=1:9
    pidTI=find(outTableAll.pidHCR==ii);

    lfU=outTableAll.numLiqHCR(pidTI)./outTableAll.numAllHCR(pidTI);
    lfH=outTableAll.numLiqLargestP(pidTI)./outTableAll.numAllLargestP(pidTI);

    lfUH=cat(2,lfU,lfH);
    lfUH(any(isnan(lfUH),2),:)=[];

    if ~isempty(lfUH)
        s{end+1}=subplot(2,5,ii);
        if ~isnan(lfUH)
            N=hist3(lfUH,'Ctrs',centers);
            N(N==0)=nan;
            Nnorm=N'./size(lfUH,1).*100;
            imagesc(Nnorm,'AlphaData',~isnan(Nnorm))
            set(gca,'YDir','normal');
            set(gca,'Xtick',tickLoc);
            set(gca,'Ytick',tickLoc);
            if ii<6
                set(gca,'Xticklabel',[])
            end
            if ii~=1 & ii~=6
                set(gca,'Yticklabel',[])
            end
            xlim([0.5 4.5])
            ylim([0.5 4.5])
            caxis([0 100]);
            title([units_str_hcr{ii},' (',num2str(size(lfUH,1)),')']);
        end

        if ii>5
            xlabel('PID')
            set(gca,'Xticklabel',tickLabel);
            xtickangle(0)
        end
        if ii==1 | ii==6
            ylabel('UWILD')
            set(gca,'Yticklabel',tickLabel);
        end
        text(0.75,4.15,panel{ii},'FontSize',11,'FontWeight','bold','Color',panelColor{ii});

        hitRate=0;

        for jj=1:length(Nnorm)
            if ~isnan(Nnorm(jj,jj))
                hitRate=hitRate+Nnorm(jj,jj);
            end
        end

        disp(['HitRate ',units_str_hcr{ii},' ',num2str(hitRate)]);
    end
end

liqFracHCRall=outTableAll.numLiqHCR./outTableAll.numAllHCR;
s{2}.Position=[0.36 0.56 0.15 0.39];
liqFracPallL=outTableAll.numLiqLargestP./outTableAll.numAllLargestP;

liqFracUW_HCR=cat(2,liqFracPallL,liqFracHCRall);
liqFracUW_HCR(any(isnan(liqFracUW_HCR),2),:)=[];

s{10}=subplot(2,5,10);

N2=hist3(liqFracUW_HCR,'Ctrs',centers);
N2(N2==0)=nan;
N2norm=N2'./size(liqFracUW_HCR,1).*100;
imagesc(N2norm,'AlphaData',~isnan(N2norm))

set(gca,'YDir','normal');
xlim([0.5 4.5])
ylim([0.5 4.5])
cb=colorbar;
caxis([0 100]);
title(['All (',num2str(size(liqFracUW_HCR,1)),')']);

xlabel('PID')
set(gca,'Yticklabel',[])

set(gca,'Xtick',tickLoc);
set(gca,'Ytick',tickLoc);
set(gca,'Xticklabel',tickLabel);
xtickangle(0)

text(0.75,4.15,'(j)','FontSize',11,'FontWeight','bold','Color','w');

s{1}.Position=[0.055 0.55 0.17 0.39];
s{2}.Position=[0.235 0.55 0.17 0.39];
s{3}.Position=[0.415 0.55 0.17 0.39];
s{4}.Position=[0.595 0.55 0.17 0.39];
s{5}.Position=[0.775 0.55 0.17 0.39];

s{6}.Position=[0.055 0.08 0.17 0.39];
s{7}.Position=[0.235 0.08 0.17 0.39];
s{8}.Position=[0.415 0.08 0.17 0.39];
s{9}.Position=[0.595 0.08 0.17 0.39];
s{10}.Position=[0.775 0.08 0.17 0.39];

cb.Position=[0.9515 0.2 0.017 0.6];
colorTitleHandle = get(cb,'Title');
set(colorTitleHandle ,'String','%');

set(gcf,'PaperPositionMode','auto')
print(fig1,[figdir,'heatMapOverlap.png'],'-dpng','-r0')

hitRate=0;

for ii=1:length(N2norm)
    hitRate=hitRate+N2norm(ii,ii);
end

disp(['HitRate All ',num2str(hitRate)]);