% Call cloud puzzle script

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project={'cset','socrates','otrec'}; %socrates, aristo, cset, otrec
freqData='10hz';
whichModel='era5';

figdir='/scr/snow2/rsfdata/projects/cset/hcr/qc3/cfradial/v3.0_full/cloudPropsProjects/paperFigs/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classTypes={'CloudLow','CloudMid','CloudHigh'};

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

for ii=1:length(project)
    if strcmp(project{ii},'cset')
        qcVersion='v3.0';
        quality='qc3';
    else
        qcVersion='v3.1';
        quality='qc3';
    end

    cfDir=HCRdir(project{ii},quality,qcVersion,freqData);

    indir=[cfDir(1:end-5),'cloudProps/'];

    in.(project{ii})=load([indir,project{ii},'_cloudProps.mat']);
end

plotVars=fields(in.cset);

for ii=1:length(plotVars)
    for jj=1:length(project)
        if ~(strcmp(project{jj},'spicule') & strcmp(plotVars{ii},'sstAll'))
            plotV.(plotVars{ii}).(project{jj})=in.(project{jj}).(plotVars{ii});
        end
    end
end

%% Plot

close all

disp('Plotting ...')


fig=figure('DefaultAxesFontSize',11,'position',[100,100,600,800],'visible','on','renderer','painters');

cols=[0,1,0; ...
    0,0,1; ...
    1,0,0];

titles={'(a) CloudLow','(b) CloudMid','(c) CloudHigh'};

projects=fields(plotV.cloudBaseAll);

edges=0.5:1:9.5;

s.s1=[];
s.s2=[];
s.s3=[];

sfields=fields(s);

for jj=1:3
    countsPerc=[];
    medians=[];
    nums=[];
    for ii=1:length(projects)
        thisVar=plotV.cloudLayersAll.(projects{ii}).(classTypes{jj});
        thisNum=sum(~isnan(thisVar));
        counts=histcounts(thisVar,edges);
        if thisNum<10
            counts(:)=0;
        end
        countsPerc=cat(1,countsPerc,counts./thisNum.*100);
        medians=cat(1,medians,median(thisVar,'omitnan'));
        nums=cat(1,nums,thisNum);
    end

    s.(sfields{jj})=subplot(3,1,jj);
    b=bar(edges(1:end-1)+(edges(2)-edges(1))/2,countsPerc,1,'FaceColor','flat');
    for kk=1:size(countsPerc,1)
        b(kk).CData=cols(kk,:);
    end

    grid on
    box on

    xlim([edges(1),edges(end)]);
    ylim([0,90]);
    
    ylabel('Percent (%)')
    set(gca, 'YGrid', 'on', 'XGrid', 'off')

    text(1.2,80,titles{jj},'FontSize',12,'FontWeight','bold');
end

xlabel('Number of cloud layers');
legend({'CSET','SOCRATES','OTREC'},'Location','southoutside','Orientation','horizontal');

s.s1.Position=[0.08,0.71,0.9,0.27];
s.s2.Position=[0.08,0.405,0.9,0.27];
s.s3.Position=[0.08,0.1,0.9,0.27];


set(gcf,'PaperPositionMode','auto')
print([figdir,'nonpLayers.png'],'-dpng','-r0');
