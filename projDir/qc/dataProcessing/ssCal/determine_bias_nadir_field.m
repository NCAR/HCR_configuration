% Ocean scan calibration for HCR data

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

% If 1, plots for individual calibration events will be made, if 0, only
% total plots will be made

project='otrec'; %socrates, aristo, cset, otrec
quality='field'; %field, qc1, or qc2
dataFreq='10hz';
addName=''; % Extra name part for output files. Default is ''.

attenuationFrom='ecmwf'; %Which attenuation data should be used: ecmwf, era5

oceanTemp=20;
salinity=35; % Ocean salinity for sig0model in per mille (world wide default is 35) and sensitivity to that number is low
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('~/gitPriv/process_HCR/oceanScans/functions/');
addpath('~/gitPriv/process_HCR/oceanScans/colormaps/');
addpath('~/gitPriv/process_HCR/NSCAL/functions/');
addpath(genpath('~/gitPriv/utils/'));

outName=[project,'_',quality,'_nadir'];

directories.figdir=['/home/romatsch/plots/HCR/oceanCal/',outName,'/'];

if ~exist(directories.figdir, 'dir')
    mkdir(directories.figdir)
end

%directories.dataDir='/home/romatsch/data/HCR/otrec/10hz/';
directories.dataDir='/run/media/romatsch/RSF0043/rsf/qc0/cfradial/moments/10hz/';

if strcmp(attenuationFrom,'era5') | strcmp(attenuationFrom,'era5sonde')
    directories.modeldir=['/scr/sci/romatsch/data/reanalysis/ecmwf/era5interp/',project,'/',dataFreq,'/'];
elseif strcmp(attenuationFrom,'ecmwf') | strcmp(attenuationFrom,'ecmwfsonde')
    directories.modeldir=['/home/romatsch/data/reanalysis/ecmwfForeInterp/',project,'/'];
end

%File with start date, end date
infile=['/home/romatsch/infiles/flights_',project,'.txt'];

caseList = table2array(readtable(infile));

%% Run processing

attenuation.liebe={};
attenuation.itu={};
attenuation.windspeed={};
attenuation.u={};
attenuation.v={};
attenuation.sst={};

% Go through flights
for ii=1:size(caseList,1)
  
    startTime=datetime(caseList(ii,1:6));
    endTimeAll=datetime(caseList(ii,7:12));
    
    % Go through hours
    while startTime<endTimeAll
        disp(startTime);
        
        endTime=startTime+hours(1);
        
        %% Get data
        [data frq]=f_load_sort_data_nadir(directories.dataDir,startTime,endTime);
        
        if ~isempty(data.dbz)
            
            % Mask data over terrain
            modeldata.topo=[];
            modeldata=read_model(modeldata,directories.modeldir,data.time(1),data.time(end));
            
            data.reflMask(modeldata.topo>0)=0;
            
            if ~max(data.reflMask)==0
                %% Attenuation
                [attenuation.liebe{end+1},attenuation.itu{end+1},attenuation.windspeed{end+1},attenuation.u{end+1},attenuation.v{end+1},attenuation.sst{end+1}]= ...
                    f_atten_layers_modelInterp(directories.modeldir,frq/1e+9,data.time,1);
                
                %% Bias
                
                data=f_determine_bias_nadir(data,attenuation.liebe{1,end},attenuation.itu{1,end},frq);
                
                data=f_sigma0_model(data,attenuation.windspeed(1,end),frq,attenuation.sst{1,end},salinity);
                
                sst=attenuation.sst{1,end};
                sst(modeldata.topo>0)=nan;
                %% Plot
                close all
                f1 = figure('Position',[200 500 2000 1200],'DefaultAxesFontSize',12);
                
                sigMeasGood=data.sig0measured;
                sigMeasGood(data.reflMask==0)=nan;
                sigMeasBad=data.sig0measured;
                sigMeasBad(data.reflMask==1)=nan;
                
                subplot(3,1,1)
                hold on
                l1=plot(data.time,data.sig0model(:,2),'-c','linewidth',2);
                l2=plot(data.time,data.sig0model(:,5),'color',[0 0.5 0],'linewidth',2);
                l3=plot(data.time,data.sig0model(:,8),'-g','linewidth',2);
                l4=plot(data.time,sigMeasGood,'-b','linewidth',1);
                l5=plot(data.time,sigMeasBad,'color',[0.5 0.5 0.5],'linewidth',0.5);
                ylabel('Sig0 (dB)');
                if strcmp(project,'otrec')
                    ylim([-20 30]);
                else
                    ylim([-20 20]);
                end
                grid on
                
                yyaxis right
                l6=plot(data.time,sigMeasGood-(data.maxpwrv-data.maxpwrh),'-r','linewidth',1);
                l7=plot(data.time,sigMeasBad-(data.maxpwrv-data.maxpwrh),'color',[0.5 0.5 0.5],'linewidth',0.5);
                ylabel('Sig0 (dB)');
                if strcmp(project,'otrec')
                    ylim([-35 10]);
                else
                    ylim([-35 5]);
                end
                ax = gca;
                ax.YColor = 'r';
                xlim([data.time(1),data.time(end)]);
                
                legend([l1 l2 l3 l4 l6],{'FreilichVanhoff','Wu','CoxMunk','Measured','Measured-(VC-HX)'},'location','southwest');
                
                title(['Sigma0: ',project,' flight ',num2str(ii),' ',datestr(data.time(1)),' to ',datestr(data.time(end))])
                
                subplot(3,1,2)
                
                hold on
                l1=plot(data.time,data.reflNoOcean,'-k','linewidth',1);
                plot([data.time(1),data.time(end)],[0.8 0.8],...
                    '-r','linewidth',1);
                ylabel('Lin refl');
                ylim([0 100]);
                xlim([data.time(1),data.time(end)]);
                %yticks(0:3:15);
                ax = gca;
                ax.YColor = 'k';
                grid on
                set(gca, 'YScale', 'log')
                               
                legend(l1,{'ReflAboveOcean'},'location','northeast');
                
                subplot(3,1,3)
                yyaxis right
                hold on
                plot(data.time,data.alt./1000,'-b','linewidth',2);
                plot(data.time,(data.range-data.alt)./1000,'-c','linewidth',2);
                plot(data.time,attenuation.liebe{1,end}.*2,'color','m','linewidth',2);
                plot(data.time,data.elev,'-k','linewidth',2);
                ylabel('Alt (km), Att (dB), Elev (deg)');
                if strcmp(project,'cset') | strcmp(project,'otrec')
                    ylim([-0.2 16.5]);
                    yticks(0:1.5:15);
                else
                    ylim([-0.2 10]);
                end
                ax = gca;
                ax.YColor = 'k';
                grid on
                
                yyaxis left
                plot(data.time,attenuation.windspeed{1,end},'-g','linewidth',2);                
                plot(data.time,sst,'-r','linewidth',2);
                ylabel('Wdspd (m s^{-1}), SST (C)');
                if strcmp(project,'cset') | strcmp(project,'otrec')
                    ylim([-0.4 33]);
                    yticks(0:3:33);
                else
                    ylim([-0.6 30]);
                    yticks(0:3:30);
                end
                xlim([data.time(1),data.time(end)]);
                ax = gca;
                ax.YColor = 'k';
                
                legend({'WindSpd','SST','Altitude','SfcRange-Alt','2wayAtten','Elevation'},'location','northeast');
                
                set(gcf,'PaperPositionMode','auto')
                print(f1,[directories.figdir,project,'_Flight',num2str(ii),'_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
                
%                 f2 = figure('Position',[200 500 2000 600],'DefaultAxesFontSize',12);
%                 hold on
%                 plot(data.time,data.maxpwrv,'-b','linewidth',1);
%                 plot(data.time,data.maxpwrh,'-r','linewidth',1);
%                 plot(data.time,data.maxpwrv+attenuation.liebe{1,end}.*2,'-c','linewidth',1);
%                 plot(data.time,data.maxpwrh+attenuation.liebe{1,end}.*2,'-m','linewidth',1);
%                 
%                 title(['Power: ',project,' flight ',num2str(ii),' ',datestr(data.time(1)),' to ',datestr(data.time(end))])
%                 
%                 ylabel('Surf power (dB)');
%                 ylim([-90 -20]);
%                 
%                 xlim([data.time(1),data.time(end)]);
%                 ax = gca;
%                 ax.YColor = 'k';
%                 
%                 legend({'V','H','V+att','H+att'},'location','northeast');
%                 
%                 set(gcf,'PaperPositionMode','auto')
%                 print(f2,[directories.figdir,project,'_Flight',num2str(ii),'_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS'),'_power'],'-dpng','-r0')
            end
        end
        
        startTime=endTime;
    end
end