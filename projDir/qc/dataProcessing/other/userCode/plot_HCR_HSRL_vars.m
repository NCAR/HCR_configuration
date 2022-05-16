% Plot HCR variables

clear all
close all

%% Input variables

startTime=datetime(2021,6,21,1,53,5);
endTime=datetime(2021,6,21,1,53,35);

indir='/scr/sleet2/rsfdata/projects/spicule/hcr/qc1/cfradial/v1.1_full/10hz/';

ylimUpper=5;

% Variables to plot.
plotVars={'WIDTH'};

saveFig=0; % Set to 1 to save the figure. Adjust outname and figdir below. 
if saveFig
    outname='examplePlot' % Figure name
    figdir=['.']; % Output directory
end

%% Get data

disp('Reading data ...');

fileList=make_file_list(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

data=[];

for ii=1:length(plotVars)
    data.(plotVars{ii})=[];
end

data=read_HCR_HSRL(fileList,data,startTime,endTime);

disp('Plotting ...');

[fig,s]=do_plot_HCR_HSRL(data,ylimUpper);

% Plot can be adjusted with the figure handle s.<variable_name>
s.WIDTH.CLim=[0,2];

if saveFig
    set(gcf,'PaperPositionMode','auto')
    print(fig,[figdir,outname,'.png'],'-dpng','-r0')
end