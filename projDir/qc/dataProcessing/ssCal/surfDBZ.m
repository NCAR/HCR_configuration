% Ocean scan calibration for HCR data

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

% If 1, plots for individual calibration events will be made, if 0, only
% total plots will be made

project='socrates'; %socrates, aristo, cset, otrec
quality='qc2'; %raw, qc1, or qc2
addName=''; % Extra name part for output files. Default is ''.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/h/eol/romatsch/gitPriv/process_HCR/oceanScans/functions/');
addpath('/h/eol/romatsch/gitPriv/process_HCR/oceanScans/colormaps/');
addpath('/h/eol/romatsch/gitPriv/process_HCR/NSCAL/functions/');
addpath(genpath('/h/eol/romatsch/gitPriv/utils/'));

outName=[project,'_',quality,'_nadir_DBZ'];

directories.figdir=['/h/eol/romatsch/hcrCalib/oceanScans/figsCompare/',outName,'/'];

if ~exist(directories.figdir, 'dir')
    mkdir(directories.figdir)
end

if ~exist(directories.figdir, 'dir')
    mkdir(directories.figdir)
end

if strcmp(project,'otrec')
    oceanTemp=20; % Ocean temperature for sig0 model in C. Will be overwritten if model data is available.
    if strcmp(quality,'raw')
        directories.dataDir='/scr/snow2/rsfdata/projects/otrec/hcr/cfradial/moments/10hz/'; % raw data
    %elseif strcmp(quality,'qc1')
    %    directories.dataDir='/scr/snow2/rsfdata/projects/socrates/hcr/qc/cfradial/moments/10hz/'; %Final
    %elseif strcmp(quality,'qc2')
    %    directories.dataDir='/scr/snow2/rsfdata/projects/socrates/hcr/qc2/cfradial/moments/10hz/'; %Final
    end
end

if strcmp(project,'socrates')
    oceanTemp=20; % Ocean temperature for sig0 model in C. Will be overwritten if model data is available.
    if strcmp(quality,'raw')
        directories.dataDir='/scr/snow2/rsfdata/projects/socrates/hcr/cfradial/moments/10hz/'; % raw data
    elseif strcmp(quality,'qc1')
        directories.dataDir='/scr/snow2/rsfdata/projects/socrates/hcr/qc/cfradial/moments/10hz/'; %Final
    elseif strcmp(quality,'qc2')
        directories.dataDir='/scr/snow2/rsfdata/projects/socrates/hcr/qc2/cfradial/moments/10hz/'; %Final
    end
end

if strcmp(project,'aristo')
    oceanTemp=20; % Ocean temperature for sig0 model in C. Will be overwritten if model data is available.
    if strcmp(quality,'raw')
        directories.dataDir='/scr/snow2/rsfdata/projects/aristo-17/hcr/cfradial/moments/10hz/'; % raw data
    elseif strcmp(quality,'qc2')
        directories.dataDir=''; %Final
    end
end

if strcmp(project,'cset')
    oceanTemp=25; % Ocean temperature for sig0 model in C. Will be overwritten if model data is available.
    if strcmp(quality,'raw')
        directories.dataDir='/scr/snow2/rsfdata/projects/cset/hcr/cfradial/moments/10hz/'; %raw data
    elseif strcmp(quality,'qc1')
        directories.dataDir='/scr/snow2/rsfdata/projects/cset/hcr/qc/cfradial/moments/10hz/'; %qc data
    elseif strcmp(quality,'qc2')
        directories.dataDir='/scr/snow2/rsfdata/projects/cset/hcr/qc2/cfradial/moments/10hz/'; %qc data
    end
end

directories.modeldir=['/scr/sci/romatsch/data/reanalysis/ecmwf/era5/',project,'/'];

%File with start date, end date
infile=['/h/eol/romatsch/hcrCalib/oceanScans/biasInFiles/flights_',project,'.txt'];

caseList = table2array(readtable(infile));

%% Run processing

% Go through flights
for ii=1:size(caseList,1)
%for ii=1:1
    
    startTime=datetime(caseList(ii,1:6));
    endTimeAll=datetime(caseList(ii,7:12));
    
    % Go through hours
    while startTime<endTimeAll
        disp(startTime);
        
        endTime=startTime+hours(1);
        
        %% Get data
        [data frq]=f_load_sort_data_nadir(directories.dataDir,startTime,endTime);
        
        if ~isempty(data.refl)            
            %% Plot
            close all
            f1 = figure('Position',[200 500 2000 1200],'DefaultAxesFontSize',12);
            
            reflMeasGood=data.refl;
            reflMeasGood(data.reflMask==0)=nan;
            reflMeasBad=data.refl;
            reflMeasBad(data.reflMask==1)=nan;
            
            subplot(3,1,1)
            hold on
            l1=plot(data.time,reflMeasGood,'-b','linewidth',1);
            l2=plot(data.time,reflMeasBad,'color',[0.5 0.5 0.5],'linewidth',0.5);
            ylabel('Reflectivity (dBZ)');
            ylim([25 55]);
            grid on
            xlim([data.time(1),data.time(end)]);
            ax = gca;
            ax.YColor = 'k';
            
            yyaxis right
            l3=plot(data.time,reflMeasGood-(data.maxpwrv-data.maxpwrh),'-r','linewidth',1);
            l4=plot(data.time,reflMeasBad-(data.maxpwrv-data.maxpwrh),'color',[0.5 0.5 0.5],'linewidth',0.5);
            ylabel('Reflectivity (dBZ)');
            ylim([5 35]);
            ax = gca;
            ax.YColor = 'r';
            
            legend([l1 l3],{'Reflectivity','Reflectivity-(VC-HX)'},'location','southwest');
            
            title(['Ocean surface reflectivity: ',project,' ',datestr(data.time(1)),' to ',datestr(data.time(end))])
            
            subplot(3,1,2)
                      
            hold on
            plot(data.time,data.reflNum,'color',[0.5 0.5 0.5],'linewidth',0.5);
            plot(data.time,movmean(data.reflNum,50),'-k','linewidth',2);
            ylabel('Number of gates with DBZ');
            ylim([0 15]);
            xlim([data.time(1),data.time(end)]);
            yticks(0:3:15);
            ax = gca;
            ax.YColor = 'k';
            grid on
            
            legend({'NumReflGates','MeanNumReflGates'},'location','northeast');
            
            subplot(3,1,3)
            yyaxis right
            hold on
            plot(data.time,data.alt./1000,'-b','linewidth',2);
            plot(data.time,(data.alt-data.range)./1000,'-c','linewidth',2);
            ylabel('Altitude (km)');
            ylim([-0.2 15]);
            yticks(0:1.5:15);
            ax = gca;
            ax.YColor = 'k';
            grid on
            
            yyaxis left
            plot(data.time,data.elev,'-k','linewidth',2);
            ylabel('Elev (deg)');
            ylim([-0.4 30]);
            xlim([data.time(1),data.time(end)]);
            yticks(0:3:30);
            ax = gca;
            ax.YColor = 'k';
            
            legend({'Elevation','Altitude','Alt-sfcRange'},'location','northeast');
            
            set(gcf,'PaperPositionMode','auto')
            print([directories.figdir,project,'_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
        end
        
        startTime=endTime;
    end
end