
% ECCO-V for NASA APR3 radar data.
% Author: Ulrike Romatschke, NCAR-EOL
% See https://doi.org/10.1175/JTECH-D-22-0019.1 for algorithm description

clear all; % Clear workspace
close all; % Close all figures

%% Input variables

% Directory path to RHI data. The data needs to be organized into
% subdirectories by date in the format yyyymmdd (20220906). The dataDir
% below needs to point to the parent directory of those subdirectories.
dataDir='/scr/cirrus3/rsfdata/projects/precip/grids/spol/radarPolar/qc2/rate/sband/v2.0/rhi/';

% Directory for output figures
figdir='/scr/cirrus3/rsfdata/projects/precip/grids/spol/radarPolar/qc2/rate/plots/ecco_v-RHIs/';

% Set showPlot below to either 'on' or 'off'. If 'on', the figures will pop up on
% the screen and also be saved. If 'off', they will be only saved but will
% not show on the screen while the code runs.
showPlot='on';

startTime=datetime(2022,6,8,0,0,0);
endTime=datetime(2022,6,8,0,15,0);

%% Tuning parameters

% These two tuning parameters affect the classification
pixRadDBZ=5; % Horizontal number of pixels over which reflectivity texture is calculated.
upperLimDBZ=30; % This affects how reflectivity texture translates into convectivity.

%% Loop through the files
fileList=makeFileList(dataDir,startTime,endTime,'xxxxxxxxxxxxxxxxxxxxxxxxxxxxx20YYMMDDxhhmmss',1);

for aa=1:length(fileList)

    disp(['File ',num2str(aa),' of ',num2str(length(fileList))]);

    %% Read data

    disp("Getting data ...");

    % Create a file list with files between the start and end time of the
    % case. Then load the data

    dataFile=[];
    dataFile.DBZ_F=[];
    dataFile.TEMP_FOR_PID=[];
    
    % Load data
    dataFile=read_spol(fileList{aa},dataFile);

    %% Loop through RHIs
    for ii=1:size(dataFile,2)
        
        dataPol=dataFile(ii);

        %% Interpolate data
        [phi,r]=meshgrid(deg2rad(dataPol.elevation),dataPol.range);
        [Xin,Yin]=pol2cart(phi,r);

        [X,Y]=meshgrid(dataPol.range(1):(dataPol.range(2)-dataPol.range(1))/4:dataPol.range(end),0:0.1:20);

        data.DBZ=griddata(double(Xin),double(Yin),dataPol.DBZ_F',double(X),double(Y),'nearest');
        data.TEMP=griddata(double(Xin),double(Yin),dataPol.TEMP_FOR_PID',double(X),double(Y),'nearest');

        %% Prepare data
               
        % % Remove specles
        % maskSub=~isnan(data.DBZ);
        % maskSub=bwareaopen(maskSub,10);
        % 
        % data.DBZ(maskSub==0)=nan;

        % Create a fake melting layer based on meltAlt
        data.MELTING_LAYER=nan(size(data.DBZ));
        data.MELTING_LAYER(data.TEMP<=0)=20;
        data.MELTING_LAYER(data.TEMP>0)=10;

        %% Texture from reflectivity

        disp('Calculating reflectivity texture ...');

        dbzBase=0; % Reflectivity base value which is subtracted from DBZ.
        dbzText=f_reflTexture(data.DBZ,pixRadDBZ,dbzBase);

        %% Convectivity

        convDBZ=1/upperLimDBZ.*dbzText;

        %% Basic classification

        disp('Basic classification ...');

        stratMixed=0.4; % Convectivity boundary between strat and mixed.
        mixedConv=0.5; % Convectivity boundary between mixed and conv.

        classBasic=f_classBasic(convDBZ,stratMixed,mixedConv,data.MELTING_LAYER);

        %% Sub classification

        disp('Sub classification ...');

        classSub=f_classSub(classBasic,data.asl,data.TOPO,data.MELTING_LAYER,data.TEMP);

        %% Plot

        disp('Plotting ...');

        % Change the default values of the subclassification to something that
        % is easier to plot
        classSubPlot=classSub;
        classSubPlot(classSub==14)=1;
        classSubPlot(classSub==16)=2;
        classSubPlot(classSub==18)=3;
        classSubPlot(classSub==25)=4;
        classSubPlot(classSub==30)=5;
        classSubPlot(classSub==32)=6;
        classSubPlot(classSub==34)=7;
        classSubPlot(classSub==36)=8;
        classSubPlot(classSub==38)=9;

        % Set up the 1D classification at the bottom of the plot
        stratConv1D=max(classSubPlot,[],1);
        time1D=data.time(~isnan(stratConv1D));
        stratConv1D=stratConv1D(~isnan(stratConv1D));

        % Set up color maps
        colmapSC=[0,0.1,0.6;
            0.38,0.42,0.96;
            0.65,0.74,0.86;
            0.32,0.78,0.59;
            1,0,0;
            1,0,1;
            1,1,0;
            0.99,0.77,0.22;
            0.7,0,0];

        col1D=colmapSC(stratConv1D,:);

        % Determine upper limit of y axis based on where the valid data ends
        ylimUpper=(max(data.asl(~isnan(data.DBZ)))./1000)+0.5;
        % Altitude of the labels within the subplots
        textAlt=ylimUpper-1;
        % Time of the labels within the subplots
        textDate=data.time(1)+seconds(5);

        close all

        f1 = figure('Position',[200 500 1600 1100],'DefaultAxesFontSize',12,'visible',showPlot);

        colormap('jet');

        % Plot reflectivity
        s1=subplot(4,1,1);

        hold on
        surf(data.time,data.asl./1000,data.DBZ,'edgecolor','none');
        view(2);
        ylabel('Altitude (km)');
        clim([-10 60]);
        ylim([0 ylimUpper]);
        xlim([data.time(1),data.time(end)]);
        set(gca,'XTickLabel',[]);
        cb1=colorbar;
        grid on
        box on

        text(textDate,textAlt,'(a) Reflectivity (dBZ)','FontSize',11,'FontWeight','bold');

        % Plot convectivity
        s2=subplot(4,1,2);

        hold on
        surf(data.time,data.asl./1000,convDBZ,'edgecolor','none');
        view(2);
        ylabel('Altitude (km)');
        clim([0 1]);
        ylim([0 ylimUpper]);
        xlim([data.time(1),data.time(end)]);
        cb2=colorbar;
        set(gca,'XTickLabel',[]);
        grid on
        box on
        text(textDate,textAlt,'(b) Convectivity','FontSize',11,'FontWeight','bold');

        % Plot the 1D classification at the very bottom (needs to be done
        % before the last plot for matlab specific reasons)
        s5=subplot(30,1,30);

        hold on
        scat1=scatter(time1D,ones(size(time1D)),10,col1D,'filled');
        set(gca,'clim',[0,1]);
        set(gca,'YTickLabel',[]);
        s5.Colormap=colmapSC;
        xlim([data.time(1),data.time(end)]);
        grid on
        box on

        % Plot classification
        s4=subplot(4,1,3);

        hold on
        surf(data.time,data.asl./1000,classSubPlot,'edgecolor','none');
        view(2);
        ylabel('Altitude (km)');
        clim([0 10]);
        ylim([0 ylimUpper]);
        xlim([data.time(1),data.time(end)]);
        s4.Colormap=colmapSC;
        clim([0.5 9.5]);
        cb4=colorbar;
        cb4.Ticks=1:9;
        cb4.TickLabels={'Strat Low','Strat Mid','Strat High','Mixed',...
            'Conv','Conv Elev','Conv Shallow','Conv Mid','Conv Deep'};
        set(gca,'XTickLabel',[]);
        grid on
        box on
        text(textDate,textAlt,'(c) Echo type','FontSize',11,'FontWeight','bold');

        % Matlab by default creates a lot of white space so we reposition the
        % panels to avoid that
        s1.Position=[0.049 0.69 0.82 0.29];
        s2.Position=[0.049 0.38 0.82 0.29];
        s4.Position=[0.049 0.07 0.82 0.29];
        s5.Position=[0.049 0.035 0.82 0.02];

        % The color bars also need to be repositioned
        cb1.Position=[0.875,0.69,0.02,0.29];
        cb2.Position=[0.875,0.38,0.02,0.29];
        cb4.Position=[0.875,0.07,0.02,0.29];

        % Save the figure based on the start and end time
        set(gcf,'PaperPositionMode','auto')
        print(f1,[figdir,'apr3_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS'),'.png'],'-dpng','-r0')
    end
end