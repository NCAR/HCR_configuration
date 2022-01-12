% Plot HCR variables

clear all
close all

%% Input variables

startTime=datetime(2018,1,29,1,10,0);
endTime=datetime(2018,1,29,1,20,0);

indir='/scr/snow2/rsfdata/projects/socrates/hcr/qc3/cfradial/hcr_hsrl_merge/v3.0/2hz/';

ylimUpper=18;

% Variables to plot.
plotVars={'HCR_DBZ','HSRL_Particle_Depolarization'};

saveFig=1; % Set to 1 to save the figure. Adjust outname and figdir below. 
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
s.HCR_DBZ.CLim=[-60,20];

if saveFig
    set(gcf,'PaperPositionMode','auto')
    print(fig,[figdir,outname,'.png'],'-dpng','-r0')
end