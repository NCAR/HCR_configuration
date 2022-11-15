% Plot HCR variables

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='otrec'; %socrates, aristo, cset, otrec
quality='qc3'; %field, qc1, or qc2
freqData='10hz';
qcVersion='v3.1';

% startTime=datetime(2019,8,16,14,52,0);
% endTime=datetime(2019,8,16,15,8,0);

startTime=datetime(2019,8,22,17,53,0);
endTime=datetime(2019,8,22,18,5,0);

indir=HCRdir(project,quality,qcVersion,freqData);

ylimUpper=14;

saveFig=1;
if saveFig
    outname=['meltingLayer_',datestr(startTime,'yyyymmdd_HHMMSS')];
    %figdir=['/scr/sci/romatsch/HCR/examplePlots/'];
    figdir=[indir(1:end-5),'examplePlots/'];
    if ~exist(figdir, 'dir')
        mkdir(figdir)
    end
end

%% Get data

disp('Reading data ...');

fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

data=[];

data.MELTING_LAYER=[];
data.DBZ_MASKED=[];

data=read_HCR(fileList,data,startTime,endTime);

disp('Plotting ...');

fig=figure('Position',[200 500 1200 600],'DefaultAxesFontSize',14);

elevenInds=find(data.MELTING_LAYER==11);
    twelveInds=find(data.MELTING_LAYER==12);
    thirteenInds=find(data.MELTING_LAYER==13);
    fourteenInds=find(data.MELTING_LAYER==14);
    
    twentyoneInds=find(data.MELTING_LAYER==21);
    twentytwoInds=find(data.MELTING_LAYER==22);
    twentythreeInds=find(data.MELTING_LAYER==23);
    twentyfourInds=find(data.MELTING_LAYER==24);

timeMat=repmat(data.time,size(data.DBZ_MASKED,1),1);

s1=surf(data.time,data.asl./1000,data.DBZ_MASKED,'EdgeColor','none');
view(2)
colorbar
s1=colMapDBZ(s1);

hold on

scatter(timeMat(elevenInds),data.asl(elevenInds)./1000,10,'k','filled');
    scatter(timeMat(fourteenInds),data.asl(fourteenInds)./1000,10,'g','filled');
    scatter(timeMat(thirteenInds),data.asl(thirteenInds)./1000,10,'c','filled');    
    scatter(timeMat(twelveInds),data.asl(twelveInds)./1000,10,'b','filled');
    
    scatter(timeMat(twentyoneInds),data.asl(twentyoneInds)./1000,10,'k','filled');
    scatter(timeMat(twentyfourInds),data.asl(twentyfourInds)./1000,10,'g','filled');
    scatter(timeMat(twentythreeInds),data.asl(twentythreeInds)./1000,10,'c','filled');    
    scatter(timeMat(twentytwoInds),data.asl(twentytwoInds)./1000,10,'b','filled');
    ax = gca;
    ax.SortMethod = 'childorder';

xlim([startTime,endTime]);
ylim([0 ylimUpper]);

ylabel('Altitude (km)')

grid on
box on

title('Reflectivity and melting layer')

if saveFig
    set(gcf,'PaperPositionMode','auto')
    print(fig,[figdir,outname,'.png'],'-dpng','-r0')
end